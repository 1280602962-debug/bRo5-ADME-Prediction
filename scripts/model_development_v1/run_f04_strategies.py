"""Fixed RF-F04 inner fits for G2 and G3, on the same validation labels as G1."""

from __future__ import annotations

import json

import pandas as pd

from model_lib import OUT, load_store
from run_fixed import inner_frames, run_job

PARTS_LOG = OUT / "fixed_timing_f04_strategies.jsonl"


def main() -> None:
    store = load_store()
    splits = pd.read_csv(OUT / "evaluation_splits.csv.gz")
    for strategy in ["G2", "G3"]:
        for outer in range(5):
            for inner in range(3):
                train, valid = inner_frames(splits, outer, inner, strategy)
                record = run_job(store, train, valid, "RF", "F04", strategy, outer, inner, "fixed")
                record.update({"outer_fold": outer, "inner_fold": inner, "model": "RF", "feature": "F04", "strategy": strategy})
                with PARTS_LOG.open("a") as handle:
                    handle.write(json.dumps(record) + "\n")
                print(json.dumps(record), flush=True)


if __name__ == "__main__":
    main()
