#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Run a one-parameter CP-fit systematic scan")
    parser.add_argument("--runner", required=True, help="Path to cpfit_standalone")
    parser.add_argument("--config", required=True, help="Base CP-fit JSON configuration")
    parser.add_argument("--output-dir", required=True, help="New directory for all scan outputs")
    parser.add_argument("--tree", required=True, help="Input ROOT tree name")
    parser.add_argument("--key", required=True, help="Dot-separated configuration key")
    parser.add_argument("--values", required=True, help='JSON array, e.g. "[2.0, 2.5, 3.0]"')
    parser.add_argument("files", nargs="+", help="Input ROOT files")
    args = parser.parse_args()

    with open(args.config, encoding="utf-8") as source:
        base = json.load(source)
    values = json.loads(args.values)
    if not isinstance(values, list) or not values:
        parser.error("--values must be a non-empty JSON array")
    try:
        json.dumps(values, allow_nan=False)
    except ValueError:
        parser.error("--values cannot contain non-finite numbers")
    keys = args.key.split(".")
    if not all(keys):
        parser.error("--key must contain non-empty key names")
    node = base
    for key in keys:
        if not isinstance(node, dict) or key not in node:
            parser.error("Configuration key does not exist: " + args.key)
        node = node[key]
    for value in values:
        numeric = isinstance(node, (int, float)) and not isinstance(node, bool)
        if type(value) is not type(node) and not (
            numeric and isinstance(value, (int, float)) and not isinstance(value, bool)
        ):
            parser.error("Scan values must have the same JSON type as the configuration value")

    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=False)
    manifest = {
        "baseConfig": str(Path(args.config).resolve()),
        "key": args.key,
        "values": values,
        "tree": args.tree,
        "files": args.files,
        "runner": str(Path(args.runner).resolve()),
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    rows = []
    for index, value in enumerate(values):
        config = json.loads(json.dumps(base))
        target = config
        for key in keys[:-1]:
            target = target[key]
        target[keys[-1]] = value
        run = output / ("variant_{:03d}".format(index))
        run.mkdir()
        config_text = json.dumps(config, indent=2, allow_nan=False) + "\n"
        config_path = run / "config.json"
        config_path.write_text(config_text, encoding="utf-8")
        result_path = run / "result.json"
        with (run / "run.log").open("w", encoding="utf-8") as log:
            completed = subprocess.run(
                [args.runner, str(config_path), str(result_path), args.tree] + args.files,
                stdout=log, stderr=subprocess.STDOUT, check=False,
            )
        row = {
            "variant": run.name,
            "value": json.dumps(value),
            "config_sha256": hashlib.sha256(config_text.encode("utf-8")).hexdigest(),
            "returncode": completed.returncode,
            "chi2": "",
            "Re": "",
            "Re_error": "",
            "Im": "",
            "Im_error": "",
        }
        if completed.returncode == 0:
            with result_path.open(encoding="utf-8") as result_file:
                fit = json.load(result_file)["fit"]
            row["chi2"] = fit["chi2"]
            for name in ("Re", "Im"):
                if name in fit["parameters"]:
                    row[name] = fit["parameters"][name]["value"]
                    row[name + "_error"] = fit["parameters"][name]["error"]
        rows.append(row)

    with (output / "summary.csv").open("w", newline="", encoding="utf-8") as summary:
        writer = csv.DictWriter(summary, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    failures = sum(row["returncode"] != 0 for row in rows)
    if failures:
        print("{} scan variant(s) failed; see run.log files".format(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
