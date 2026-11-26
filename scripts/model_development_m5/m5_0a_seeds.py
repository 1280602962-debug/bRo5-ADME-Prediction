"""Sixty fixed-forest fits on the M4 ensemble cohort. Seeds are shared inside a fold."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from m5_lib import OUT, SEEDS, arm_matrix, fit_fold, folds, load_tables, paired_interval, regression
from model_lib import load_store

PARTS = OUT / "seed_parts"
ARMS = ("R0", "R1", "R2", "R3")
CONTRASTS = (("R3-R0", "R3", "R0"), ("R3-R1", "R3", "R1"), ("R3-R2", "R3", "R2"))


def with_error(frame: pd.DataFrame) -> pd.DataFrame:
    work = frame.copy()
    work["abs_error"] = (work.y_pred - work.y).abs()
    return work


def average_seeds(frames: list[pd.DataFrame]) -> pd.DataFrame:
    stacked = pd.concat(frames, ignore_index=True)
    grouped = stacked.groupby("label_id", as_index=False).agg(
        y_pred=("y_pred", "mean"),
        abs_error=("abs_error", "mean"),
        y=("y", "first"),
        Graph_ID=("Graph_ID", "first"),
        Study_key=("Study_key", "first"),
        Scaffold_connectivity_group=("Scaffold_connectivity_group", "first"),
        cohort=("cohort", "first"),
        MW=("MW", "first"),
    )
    return grouped


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    PARTS.mkdir(parents=True, exist_ok=True)
    store = load_store()
    tables = load_tables(store)
    cohort = folds(store, tables["ensemble_ids"])
    predictions = []
    for arm in ARMS:
        matrix, names = arm_matrix(store, tables, arm)
        for outer, train, test in cohort:
            for seed in SEEDS:
                path = PARTS / f"o{outer}_{arm}_s{seed}.csv"
                pred = fit_fold(store, train, test, matrix, names, seed, path)
                pred["arm"] = arm
                pred["outer_fold"] = outer
                predictions.append(pred)
                print(json.dumps({"arm": arm, "outer": outer, "seed": seed, "train": int(len(train)), "test": int(len(test))}), flush=True)
    pred = pd.concat(predictions, ignore_index=True)
    pred.to_csv(OUT / "seed_predictions.csv.gz", index=False)
    rows = []
    for (arm, seed), part in pred.groupby(["arm", "seed"]):
        row = {"arm": arm, "seed": int(seed), "summary": "one_seed"}
        row.update(regression(part))
        rows.append(row)
    averaged = {}
    for arm in ARMS:
        averaged[arm] = average_seeds([with_error(pred[(pred.arm == arm) & (pred.seed == seed)]) for seed in SEEDS])
        row = {"arm": arm, "seed": "mean", "summary": "seed_mean_prediction"}
        row.update(regression(averaged[arm]))
        rows.append(row)
        # identical-data repeat range
        maes = [regression(pred[(pred.arm == arm) & (pred.seed == seed)])["mae_source"] for seed in SEEDS]
        rows.append({"arm": arm, "seed": "range", "summary": "repeat_range", "mae_source_min": min(maes), "mae_source_max": max(maes), "mae_source_span": max(maes) - min(maes)})
    metrics = pd.DataFrame(rows)
    contrast_rows = []
    for name, left_name, right_name in CONTRASTS:
        for seed in SEEDS:
            left = with_error(pred[(pred.arm == left_name) & (pred.seed == seed)])
            right = with_error(pred[(pred.arm == right_name) & (pred.seed == seed)])
            row = {"contrast": name, "seed": int(seed), "summary": "one_seed"}
            row.update(paired_interval(left, right))
            contrast_rows.append(row)
        left_mean = averaged[left_name]
        right_mean = averaged[right_name]
        row = {"contrast": name, "seed": "mean", "summary": "seed_mean_absolute_error"}
        row.update(paired_interval(left_mean, right_mean))
        contrast_rows.append(row)
    contrasts = pd.DataFrame(contrast_rows)
    metrics.to_csv(OUT / "seed_sensitivity_metrics.csv", index=False)
    contrasts.to_csv(OUT / "seed_sensitivity.csv", index=False)
    summary = {
        "fits": 60,
        "seeds": list(SEEDS),
        "same_seed_within_fold": True,
        "seed_not_derived_from_arm_name": True,
        "replaces_m4_contrasts": False,
        "contrasts": contrasts[contrasts.seed.astype(str).eq("mean")][["contrast", "estimate", "ci95_low", "ci95_high"]].to_dict(orient="records"),
    }
    (OUT / "m5_0a_summary.json").write_text(json.dumps(summary, indent=2))
    print(contrasts[contrasts.seed.astype(str).eq("mean")].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
