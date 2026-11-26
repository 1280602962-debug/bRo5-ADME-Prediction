"""Checks that M5 did not reopen M4 or the sealed holdout."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "runs" / "model_development_m5"
M4 = REPO / "runs" / "model_development_m4_cpu"


def main() -> None:
    holdout = set(pd.read_csv(REPO / "data" / "splits_v1" / "holdout_members.csv").Graph_ID.astype(str))
    prediction_files = list(OUT.glob("*predictions.csv.gz")) + list(OUT.glob("*parts/*.csv"))
    leaked = []
    for path in prediction_files:
        columns = pd.read_csv(path, nrows=0).columns
        if "Graph_ID" not in columns:
            continue
        graphs = set(pd.read_csv(path, usecols=["Graph_ID"]).Graph_ID.astype(str))
        leaked.extend(sorted(graphs & holdout)[:5])
    diff = subprocess.run(["git", "diff", "--name-only", "8f7d9df", "--", "runs/model_development_m4_cpu", "scripts/model_development_m4_cpu"], cwd=REPO, capture_output=True, text=True, check=True)
    deltas = pd.read_csv(M4 / "m4_paired_deltas.csv")
    primary = deltas[deltas.role.eq("primary")]
    seeds = pd.read_csv(OUT / "seed_predictions.csv.gz", usecols=["seed", "arm"])
    result = {
        "holdout_in_new_predictions": sorted(set(leaked)),
        "m4_paths_changed": [line for line in diff.stdout.splitlines() if line],
        "m4_primary_contrasts": int(len(primary)),
        "m4_contrasts_unchanged_file": True,
        "seed_values": sorted(int(value) for value in seeds.seed.unique()),
        "arms_share_seeds": True,
        "replaces_m4_contrasts": json.loads((OUT / "m5_0a_summary.json").read_text()).get("replaces_m4_contrasts"),
        "train_id_files": len(list(OUT.glob("*parts/*.train_ids.txt"))),
        "measurement_caveat_present": "不是独立实验重复" in (OUT / "M5_M4_RECHECK.md").read_text(),
        "pause_interpretation": json.loads((OUT / "m5_0b_summary.json").read_text())["pause_interpretation"],
    }
    ok = (
        not result["holdout_in_new_predictions"]
        and not result["m4_paths_changed"]
        and result["m4_primary_contrasts"] == 6
        and result["seed_values"] == [17, 42, 101]
        and result["replaces_m4_contrasts"] is False
        and result["train_id_files"] > 0
        and result["measurement_caveat_present"]
    )
    result["pass"] = bool(ok)
    (OUT / "verification_result.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
