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
        cls.config["eventWeights"]["crossChannelMeanWeights"] = {
            "threePi0FromKsPm": False,
            "semileptonicFromKs00": False,
        }
        cls.config["eventWeights"]["signalMaxRelativeError"] = 0.0
        cls.config["signalNormalizationInRangeOnly"] = False
        for name, component in cls.config["eventWeights"]["signalComponents"].items():
            component["enabled"] = name in ("ks_pm", "kl_pm")
            component["correctionFile"] = str(cls.path / "signal.root")
            component["correctionHistogram"] = "{}_weights".format(name)
            component["ratioFile"] = str(cls.path / "signal_ratios.root")
            component["ratioSpline"] = "spline_1" if name.startswith("ks_") else "spline_2"
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
            ("semileptonic", "Semileptonic", self.config["eventWeights"]["semileptonic"]),
            ("threePi0", "3pi0", self.config["eventWeights"]["threePi0"]),
            ("regeneration", "Regeneration", 0.8),
        ]:
            config = copy.deepcopy(self.config)
            config["eventWeights"]["apply"][key] = True
            if key == "signal":
                config["eventWeights"]["signalComponents"]["kl_pm"]["enabled"] = False
            result = self.run_fit(config)
            if key == "signal":
                self.assertAlmostEqual(result["signalWeights"]["sum"], 560 * factor * 1.1)
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
        config["eventWeights"]["signalComponents"]["kl_pm"]["enabled"] = False
        config["eventWeights"]["signalMaxRelativeError"] = 0.1
        weights = self.run_fit(config)["signalWeights"]
        replaced = weights["largeRelativeErrorReplaced"]
        self.assertGreater(replaced, 0)
        self.assertLess(replaced, 560)
        self.assertEqual(weights["largeRelativeErrorReplacedByComponent"]["ks_pm"], replaced)
        self.assertAlmostEqual(weights["sum"], (560 - replaced) * 0.96 * 1.1 + replaced * 1.1)

    def test_signal_weight_error_replaces_only_the_affected_component(self):
        config = copy.deepcopy(self.config)
        config["eventWeights"]["apply"]["signal"] = True
        config["eventWeights"]["signalMaxRelativeError"] = 0.1
        weights = self.run_fit(config)["signalWeights"]
        replaced = weights["largeRelativeErrorReplaced"]
        self.assertGreater(replaced, 0)
        self.assertLess(replaced, 560)
        self.assertEqual(weights["largeRelativeErrorReplacedByComponent"]["ks_pm"], replaced)
        self.assertNotIn("kl_pm", weights["largeRelativeErrorReplacedByComponent"])
        nominal_component_sum = 0.96 * 1.1 + 0.8 * 0.6
        substituted_component_sum = 1.1 + 0.8 * 0.6
        self.assertAlmostEqual(
            weights["sum"], (560 - replaced) * nominal_component_sum + replaced * substituted_component_sum)

    def test_cross_channel_weights_use_each_channel_kinematics(self):
        baseline_config = copy.deepcopy(self.config)
        baseline_config["eventWeights"]["apply"]["threePi0"] = True
        baseline_config["eventWeights"]["apply"]["semileptonic"] = True
        baseline = self.run_fit(baseline_config)

        ks_pm_config = copy.deepcopy(baseline_config)
        ks_pm_config["eventWeights"]["crossChannelMeanWeights"]["threePi0FromKsPm"] = True
        ks_pm = self.run_fit(ks_pm_config)
        ks_pm_mean = 0.96 * 1.05
        self.assertAlmostEqual(ks_pm["crossChannelMeanWeights"]["values"]["threePi0FromKsPm"], ks_pm_mean)
        self.assertAlmostEqual(
            ks_pm["channels"]["3pi0"]["fittedIntegral"] / baseline["channels"]["3pi0"]["fittedIntegral"],
            ks_pm_mean)
        self.assertAlmostEqual(
            ks_pm["channels"]["Semileptonic"]["fittedIntegral"],
            baseline["channels"]["Semileptonic"]["fittedIntegral"])

        ks_00_config = copy.deepcopy(baseline_config)
        ks_00_config["eventWeights"]["crossChannelMeanWeights"]["semileptonicFromKs00"] = True
        ks_00 = self.run_fit(ks_00_config)
        ks_00_mean = 0.7 * 1.05
        self.assertAlmostEqual(
            ks_00["crossChannelMeanWeights"]["values"]["semileptonicFromKs00"], ks_00_mean)
        self.assertAlmostEqual(
            ks_00["channels"]["Semileptonic"]["fittedIntegral"] /
            baseline["channels"]["Semileptonic"]["fittedIntegral"],
            ks_00_mean)
        self.assertAlmostEqual(
            ks_00["channels"]["3pi0"]["fittedIntegral"],
            baseline["channels"]["3pi0"]["fittedIntegral"])

    def test_signal_component_weight_composition(self):
        config = copy.deepcopy(self.config)
        config["eventWeights"]["apply"]["signal"] = True
        config["eventWeights"]["signalMaxRelativeError"] = 0.0
        combined = self.run_fit(config)["signalWeights"]["sum"]
        self.assertAlmostEqual(combined, 560 * (0.96 * 1.1 + 0.8 * 0.6))

        config["eventWeights"]["signalComponents"]["kl_pm"]["enabled"] = False
        ks_only = self.run_fit(config)["signalWeights"]["sum"]
        self.assertAlmostEqual(ks_only, 560 * 0.96 * 1.1)

        config["eventWeights"]["signalComponents"]["ks_pm"]["enabled"] = False
        config["eventWeights"]["signalComponents"]["kl_pm"]["enabled"] = False
        no_components = self.run_fit(config)["signalWeights"]["sum"]
        self.assertAlmostEqual(no_components, 560.0)

        config["eventWeights"]["signalComponents"]["ks_00"].update(
            enabled=True, correctionFile="", correctionHistogram="", ratioFile="", ratioSpline="")
        self.assertIn("incomplete ROOT object configuration: ks_00", self.run_fit(config, expected=1))

    def test_luminosity_factor(self):
        base = copy.deepcopy(self.config)
        base["parameterSettings"]["A_signal"]["fixed"] = True
        base["parameterSettings"]["Re"]["fixed"] = False
        baseline = self.run_fit(base)
        config = copy.deepcopy(base)
        config["luminosityFactor"] = 2.0 * baseline["luminosityFactor"]
        result = self.run_fit(config)
        self.assertAlmostEqual(result["luminosityFactor"], 2.0 * baseline["luminosityFactor"])
        ratio = result["channels"]["Signal"]["fittedIntegral"] / baseline["channels"]["Signal"]["fittedIntegral"]
        self.assertAlmostEqual(ratio, 2.0)

    def test_signal_normalization_in_range_only(self):
        config = copy.deepcopy(self.config)
        config["eventWeights"]["apply"]["signal"] = True
        config["parameterSettings"]["A_signal"]["fixed"] = True
        config["parameterSettings"]["Re"]["fixed"] = False
        config["deltaT"]["rangeX"] = [-1.0, 1.0]
        config["deltaT"]["xRangeDisplay"] = [-1.0, 1.0]
        full = self.run_fit(config)
        config["signalNormalizationInRangeOnly"] = True
        in_range = self.run_fit(config)
        self.assertTrue(in_range["signalNormalizationInRangeOnly"])
        self.assertLess(in_range["channels"]["Signal"]["fittedIntegral"],
                        full["channels"]["Signal"]["fittedIntegral"])

    def test_branching_ratio_correction_switch(self):
        base = copy.deepcopy(self.config)
        base["parameterSettings"]["A_signal"]["fixed"] = True
        base["parameterSettings"]["Re"]["fixed"] = False
        on = self.run_fit(base)
        self.assertTrue(all(on["branchingRatioCorrectionApplied"].values()))
        base["branchingRatioCorrection"] = {"Signal": False}
        partial = self.run_fit(base)
        self.assertFalse(partial["branchingRatioCorrectionApplied"]["Signal"])
        self.assertTrue(partial["branchingRatioCorrectionApplied"]["3pi0"])
        self.assertEqual(partial["branchingRatioCorrections"]["Signal"], 1.0)
        self.assertEqual(partial["branchingRatioCorrections"]["3pi0"], on["branchingRatioCorrections"]["3pi0"])
        base["branchingRatioCorrection"] = False
        off = self.run_fit(base)
        self.assertEqual(set(off["branchingRatioCorrections"].values()), {1.0})
        base["branchingRatioCorrection"] = {"Bogus": False}
        self.assertIn("Unknown channel", self.run_fit(base, expected=1))

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
        config["eventWeights"]["signalComponents"]["ks_pm"]["correctionFile"] = str(self.path / "missing_signal.root")
        config["regenerationShapeCorrection"]["weightFunctionsFile"] = str(self.path / "missing.root")
        self.run_fit(config)
        config["eventWeights"]["apply"]["signal"] = True
        self.assertIn("Could not open signal correction file for ks_pm", self.run_fit(config, expected=1))
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
