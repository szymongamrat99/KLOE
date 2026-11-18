import copy
import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "build/Subanalysis/CPFit/cpfit_standalone"
SCANNER = ROOT / "Subanalysis/CPFit/app/scan_cuts.py"


class CPFitStandaloneTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.path = Path(cls.temp.name)
        macro = ROOT / "tests/cpfit_test_data.C"
        subprocess.run(["root", "-n", "-l", "-b", "-q",
                        '{}("{}")'.format(macro, cls.path)], check=True)
        shutil.copyfile(str(cls.path / "events.root"), str(cls.path / "events2.root"))
        cls.config = json.loads((ROOT / "Subanalysis/CPFit/config/config.json").read_text())
        cls.config["deltaT"].update(rangeX=[-3.0, 3.0], xRangeDisplay=[-3.0, 3.0], resolution=1.0)
        cls.config["eventWeights"]["apply"] = dict.fromkeys(cls.config["eventWeights"]["apply"], False)
        cls.config["eventWeights"]["signalMaxRelativeError"] = 0.0
        cls.config["eventWeights"]["signalCorrectionFile"] = str(cls.path / "signal.root")
        cls.config["eventWeights"]["signalCorrectionHistogram"] = "weights"
        cls.config["regenerationShapeCorrection"]["weightFunctionsFile"] = str(cls.path / "regen.root")
        for cut in cls.config["cuts"].values():
            cut["enabled"] = False
        for name, parameter in cls.config["parameterSettings"].items():
            parameter["fixed"] = name != "A_signal"
        cls.config["fitSettings"]["printLevel"] = 0
        cls.sequence = 0

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def run_fit(self, config, files=1, expected=0):
        type(self).sequence += 1
        stem = self.path / ("run_{}".format(self.sequence))
        config_path = stem.with_suffix(".json")
        output_path = self.path / (stem.name + "_result.json")
        config_path.write_text(json.dumps(config))
        command = [str(RUNNER), str(config_path), str(output_path), "analysis"]
        command.append(str(self.path / "events.root"))
        if files == 2:
            command.append(str(self.path / "events2.root"))
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                universal_newlines=True)
        self.assertEqual(result.returncode, expected, result.stdout)
        return json.loads(output_path.read_text()) if expected == 0 else result.stdout

    def test_chain_and_ellipse_selection(self):
        nominal = self.run_fit(self.config, files=2)
        self.assertEqual(nominal["channels"]["Data"]["selected"], 1680)
        self.assertEqual(nominal["channels"]["Signal"]["selected"], 1120)
        self.assertEqual(nominal["channels"]["MC sum"]["selected"], 1680)
        self.assertGreater(nominal["physicsConstants"]["mK0"], 490)
        config = copy.deepcopy(self.config)
        config["cuts"]["ellipse"]["enabled"] = True
        selected = self.run_fit(config, files=2)
        self.assertEqual(selected["channels"]["Data"]["selected"], 840)
        self.assertEqual(selected["channels"]["Signal"]["selected"], 560)

    def test_independent_weights_and_regeneration_veto(self):
        baseline = self.run_fit(self.config)
        for key, channel, factor in [
            ("signal", "Signal", 0.96),
            ("semileptonic", "Semileptonic", 0.9272),
            ("threePi0", "3pi0", 0.9375),
            ("regeneration", "Regeneration", 0.8),
        ]:
            config = copy.deepcopy(self.config)
            config["eventWeights"]["apply"][key] = True
            result = self.run_fit(config)
            if key == "signal":
                self.assertAlmostEqual(result["signalWeights"]["sum"], 560 * factor)
            else:
                ratio = (result["channels"][channel]["fittedIntegral"] /
                         baseline["channels"][channel]["fittedIntegral"])
                self.assertAlmostEqual(ratio, factor)
        config = copy.deepcopy(self.config)
        veto = config["regenerationExclusion"]
        veto.update(enabled=True, beamPipeBoundary=0.0, sigmaMultiplier=1.5)
        veto["spherical"].update(meanCharged=4.8, meanNeutral=4.8, sigmaCharged=1.0, sigmaNeutral=1.0)
        result = self.run_fit(config)
        self.assertEqual(result["channels"]["Regeneration"]["selected"], 0)
        self.assertEqual(result["channels"]["Data"]["selected"], 840)

    def test_signal_weight_relative_error_cut(self):
        config = copy.deepcopy(self.config)
        config["eventWeights"]["apply"]["signal"] = True
        config["eventWeights"]["signalMaxRelativeError"] = 0.1
        weights = self.run_fit(config)["signalWeights"]
        replaced = weights["largeRelativeErrorReplaced"]
        self.assertGreater(replaced, 0)
        self.assertLess(replaced, 560)
        self.assertAlmostEqual(weights["sum"], (560 - replaced) * 0.96 + replaced)

    def test_diagnostics_profile(self):
        config = copy.deepcopy(self.config)
        config["diagnostics"] = {"profile": {"enabled": True, "points": 5, "nSigma": 2.0}}
        config["parameterSettings"]["A_bcg_group"].update(enabled=True, fixed=False)
        for name in ("A_three", "A_semileptonic", "A_other"):
            config["parameterSettings"][name]["enabled"] = False
        config["parameterSettings"]["A_signal"]["fixed"] = False
        diagnostics = self.run_fit(config)["diagnostics"]
        self.assertNotIn("error", diagnostics)
        self.assertIn("conditionNumber", diagnostics["degeneracy"])
        profile = diagnostics["profile"]
        self.assertTrue(profile)
        for scan in profile.values():
            self.assertEqual(len(scan), 5)
            self.assertGreater(max(point["deltaChi2"] for point in scan), -1e-3)

    def test_missing_regeneration_file(self):
        config = copy.deepcopy(self.config)
        config["eventWeights"]["signalCorrectionFile"] = str(self.path / "missing_signal.root")
        config["regenerationShapeCorrection"]["weightFunctionsFile"] = str(self.path / "missing.root")
        self.run_fit(config)
        config["eventWeights"]["apply"]["signal"] = True
        self.assertIn("Could not open signal correction file", self.run_fit(config, expected=1))
        config["eventWeights"]["apply"]["signal"] = False
        config["eventWeights"]["apply"]["regeneration"] = True
        self.assertIn("Could not open regeneration weight file", self.run_fit(config, expected=1))

    def test_existing_output_is_not_overwritten(self):
        config = self.path / "existing_config.json"
        output = self.path / "existing_result.json"
        config.write_text(json.dumps(self.config))
        output.write_text("keep")
        result = subprocess.run(
            [str(RUNNER), str(config), str(output), "analysis", str(self.path / "events.root")],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        self.assertEqual(result.returncode, 1)
        self.assertIn("Output file already exists", result.stdout)
        self.assertEqual(output.read_text(), "keep")

    def test_scan_records_failure(self):
        config_path = self.path / "scan_base.json"
        config_path.write_text(json.dumps(self.config))
        output = self.path / "scan"
        result = subprocess.run([
            "python3", str(SCANNER), "--runner", str(RUNNER), "--config", str(config_path),
            "--output-dir", str(output), "--tree", "analysis", "--key", "cuts.ellipse.nSigma",
            "--values", "[-1, 2.5]", str(self.path / "events.root"),
        ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        self.assertEqual(result.returncode, 1, result.stdout)
        self.assertIn("Invalid ellipse", (output / "variant_000/run.log").read_text())
        self.assertTrue((output / "variant_001/result.json").is_file())
        self.assertTrue((output / "summary.csv").is_file())


if __name__ == "__main__":
    unittest.main()
