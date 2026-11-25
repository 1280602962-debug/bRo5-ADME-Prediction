"""Fixed-parameter inner-validation fits. Resumes from completed part files."""

from __future__ import annotations

import argparse
import json
import resource
import time
from pathlib import Path

import numpy as np
import pandas as pd

from model_lib import FIXED, OUT, feature_matrix, fit_predict, load_store, model_seed, preprocess, rows_for

PARTS = OUT / "parts"
SCREEN = [
    (model, feature, "screen12")
    for model in ["RF", "XGBoost", "LightGBM", "SVR"]
    for feature in ["F02", "F03", "F05"]
]
BASELINES = [(model, "F01", "baseline_f01") for model in ["RF", "XGBoost", "LightGBM", "SVR"]]
BASELINES += [("RF", "F04", "baseline_f04"), ("Median", "F01", "baseline_simple"), ("Ridge", "F01", "baseline_simple")]


def rss_mb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def inner_frames(splits: pd.DataFrame, outer: int, inner: int, strategy: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    part = splits[splits.evaluation.eq("C_inner") & splits.outer_fold.eq(outer) & splits.inner_fold.eq(inner)]
    valid = part[part.role.eq("validation")].copy()
    train = part[part.role.eq("train_candidate")].copy()
    if strategy == "G2":
        train = train[train.g2]
    elif strategy == "G3":
        train = train[train.g3]
    elif strategy != "G1":
        raise KeyError(strategy)
    return train, valid


def run_job(store: dict, train: pd.DataFrame, valid: pd.DataFrame, model: str, feature: str, strategy: str, outer: int, inner: int, phase: str) -> dict:
    path = PARTS / f"{phase}_{strategy}_o{outer}_i{inner}_{model}_{feature}.csv"
    if path.exists() and len(pd.read_csv(path)) == len(valid):
        return {"file": path.name, "resumed": True}
    matrix, names = feature_matrix(store, feature)
    train_x, y, weights = rows_for(train, store, matrix)
    valid_x, y_valid, _ = rows_for(valid, store, matrix)
    scale = model in {"Ridge", "SVR"}
    train_ready, valid_ready, kept = preprocess(train_x, valid_x, names, scale)
    started = time.perf_counter()
    seed = model_seed(phase, strategy, outer, inner, model, feature)
    valid_pred, train_pred = fit_predict(model, train_ready, y, weights, [valid_ready, train_ready], FIXED.get(model, {}), seed)
    elapsed = time.perf_counter() - started
    output = valid[["label_id", "Study_key", "Graph_ID", "y"]].copy()
    output["y_pred"] = valid_pred
    output["model"] = model
    output["feature"] = feature
    output["strategy"] = strategy
    output["outer_fold"] = outer
    output["inner_fold"] = inner
    output["phase"] = phase
    output.to_csv(path, index=False)
    return {
        "file": path.name,
        "resumed": False,
        "seconds": round(elapsed, 3),
        "rss_mb": round(rss_mb(), 1),
        "n_train": int(len(train)),
        "n_valid": int(len(valid)),
        "n_features": int(train_ready.shape[1]),
        "train_mae_label": float(np.mean(np.abs(train_pred - y))),
        "valid_mae_label": float(np.mean(np.abs(valid_pred - y_valid))),
        "kept_note": f"{len(kept)} columns",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", action="store_true")
    parser.add_argument("--selected", action="store_true")
    parser.add_argument("--strategy", default="G1")
    parser.add_argument("--phase", default="fixed")
    args = parser.parse_args()
    PARTS.mkdir(parents=True, exist_ok=True)
    store = load_store()
    splits = pd.read_csv(OUT / "evaluation_splits.csv.gz")
    if args.benchmark:
        jobs = [("RF", "F05", "benchmark"), ("XGBoost", "F05", "benchmark"), ("LightGBM", "F05", "benchmark"), ("SVR", "F05", "benchmark")]
        outers = [0]
        inners = [0]
    elif args.selected:
        choice = pd.read_csv(OUT / "selected_pairs.csv")
        jobs = [(row.model, row.feature, "selected") for row in choice.drop_duplicates(["model", "feature"]).itertuples(index=False)]
        outers = range(5)
        inners = range(3)
    else:
        jobs = SCREEN + BASELINES
        outers = range(5)
        inners = range(3)
    log_path = OUT / ("benchmark_timing.jsonl" if args.benchmark else f"fixed_timing_{args.strategy}.jsonl")
    for outer in outers:
        for inner in inners:
            train, valid = inner_frames(splits, outer, inner, args.strategy)
            if train.empty or valid.empty:
                raise RuntimeError(f"empty split {args.strategy} {outer} {inner}")
            for model, feature, _tag in jobs:
                record = run_job(store, train, valid, model, feature, args.strategy, outer, inner, "benchmark" if args.benchmark else args.phase)
                record.update({"outer_fold": outer, "inner_fold": inner, "model": model, "feature": feature, "strategy": args.strategy})
                with log_path.open("a") as handle:
                    handle.write(json.dumps(record) + "\n")
                print(json.dumps(record), flush=True)


if __name__ == "__main__":
    main()
