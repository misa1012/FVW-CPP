#!/usr/bin/env python3
"""Generate cutoffParam convergence study configs and a submission script.

Usage:
  python tools/python/studies/cutoff/generate_cases.py \
      --base tutorials/NTNU/config.json \
      --out-dir tutorials/NTNU/cutoff_study \
      --values 0.01 0.05 0.1 0.25 0.5 0.75
"""

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Generate cutoffParam study configs.")
    parser.add_argument("--base", required=True, help="Base config.json")
    parser.add_argument("--out-dir", required=True, help="Output directory for configs")
    parser.add_argument("--values", nargs="+", type=float, required=True, help="cutoffParam values")
    args = parser.parse_args()

    base_path = Path(args.base)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(base_path, "r") as f:
        base = json.load(f)

    base_model = base.get("turbine", {}).get("model", "NTNU")
    output_root = str(out_dir)

    configs = []
    for val in args.values:
        cfg = json.loads(json.dumps(base))
        cfg["simulation"]["cutoffParam"] = val
        cfg["turbine"]["model"] = base_model
        cfg["caseName"] = f"{base_model}_cutoff_{val:.2f}"
        cfg["outputRoot"] = output_root
        out_path = out_dir / f"config_cutoff_{val:.2f}.json"
        with open(out_path, "w") as f:
            json.dump(cfg, f, indent=4)
        configs.append(out_path)

    # Write a run list for convenience
    run_list = out_dir / "cutoff_cases.txt"
    with open(run_list, "w") as f:
        for p in configs:
            f.write(str(p) + "\n")

    # Write a submit script using sbatch and CONFIG env var
    submit_script = out_dir / "submit_all.sh"
    with open(submit_script, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("set -euo pipefail\n")
        f.write("ROOT=\"$(cd \"$(dirname \"$0\")/../..\" && pwd)\"\n")
        f.write("while read -r cfg; do\n")
        f.write("  [ -z \"$cfg\" ] && continue\n")
        f.write("  echo \"Submitting $cfg\"\n")
        f.write("  (cd \"$ROOT\" && CONFIG=\"$cfg\" sbatch submit_job.slurm)\n")
        f.write("done < \"$ROOT/$(python - <<'PY'\nfrom pathlib import Path\nprint(Path('tutorials/NTNU/cutoff_study/cutoff_cases.txt'))\nPY)\"\n")
    submit_script.chmod(0o755)

    print("Wrote configs to:", out_dir)
    print("Run list:", run_list)
    print("Submit script:", submit_script)


if __name__ == "__main__":
    main()
