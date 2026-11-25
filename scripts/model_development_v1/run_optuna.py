"""Inner-validation Optuna search for the two fixed-screen winners in each outer fold."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import optuna
import pandas as pd

from model_lib import OUT, feature_matrix, fit_predict, load_store, model_seed, preprocess, rows_for
from run_fixed import inner_frames

import sys

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "stage_b"))
from score_stage_b import source_table  # noqa: E402

optuna.logging.set_verbosity(optuna.logging.WARNING)


def suggest(trial: optuna.Trial, model: str) -> dict:
    if model == "RF":
        depth = trial.suggest_int("max_depth", 0, 24)
        return {
            "n_estimators": trial.suggest_int("n_estimators", 300, 1000),
            "max_depth": None if depth == 0 else depth,
            "min_samples_leaf": trial.suggest_int("min_samples_leaf", 1, 10),
            "min_samples_split": trial.suggest_int("min_samples_split", 2, 20),
            "max_features": trial.suggest_float("max_features", 0.2, 1.0),
        }
    if model == "XGBoost":
        return {
            "n_estimators": trial.suggest_int("n_estimators", 200, 1500),
            "max_depth": trial.suggest_int("max_depth", 2, 8),
            "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.15),
            "min_child_weight": trial.suggest_float("min_child_weight", 1.0, 15.0),
            "subsample": trial.suggest_float("subsample", 0.6, 1.0),
            "colsample_bytree": trial.suggest_float("colsample_bytree", 0.3, 1.0),
            "reg_alpha": trial.suggest_float("reg_alpha", 1e-5, 10.0, log=True),
            "reg_lambda": trial.suggest_float("reg_lambda", 1e-3, 30.0, log=True),
        }
    if model == "LightGBM":
        depth = trial.suggest_categorical("max_depth", [-1, 3, 4, 5, 6, 7, 8, 9, 10])
        leaves = trial.suggest_int("num_leaves", 8, 63)
        if depth > 0:
            leaves = min(leaves, 2**depth - 1)
        return {
            "n_estimators": trial.suggest_int("n_estimators", 200, 1500),
            "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.15),
            "num_leaves": leaves,
            "max_depth": depth,
            "min_child_samples": trial.suggest_int("min_child_samples", 5, 50),
            "subsample": trial.suggest_float("subsample", 0.6, 1.0),
            "colsample_bytree": trial.suggest_float("colsample_bytree", 0.3, 1.0),
            "reg_alpha": trial.suggest_float("reg_alpha", 1e-5, 10.0, log=True),
            "reg_lambda": trial.suggest_float("reg_lambda", 1e-3, 30.0, log=True),
        }
    if model == "SVR":
        return {
            "C": trial.suggest_float("C", 0.1, 1000.0, log=True),
            "gamma": trial.suggest_float("gamma", 1e-5, 1.0, log=True),
            "epsilon": trial.suggest_float("epsilon", 0.01, 0.3),
        }
    raise KeyError(model)


def pooled_mae(store, splits, outer: int, model: str, feature: str, params: dict, strategy: str) -> float:
    matrix, names = feature_matrix(store, feature)
    frames = []
    scale = model in {"Ridge", "SVR"}
    for inner in range(3):
        train, valid = inner_frames(splits, outer, inner, strategy)
        train_x, y, weights = rows_for(train, store, matrix)
        valid_x, y_valid, _ = rows_for(valid, store, matrix)
        train_ready, valid_ready, _ = preprocess(train_x, valid_x, names, scale)
        seed = model_seed("optuna", strategy, outer, model, feature, inner, json.dumps(params, sort_keys=True, default=str))
        prediction = fit_predict(model, train_ready, y, weights, [valid_ready], params, seed)[0]
        part = valid[["label_id", "Study_key", "Graph_ID", "y"]].copy()
        part["y_pred"] = prediction
        frames.append(part)
    pred = pd.concat(frames, ignore_index=True).rename(columns={"Study_key": "study_id", "Graph_ID": "graph_id"})
    pred["abs_error"] = (pred.y - pred.y_pred).abs()
    return float(source_table(pred, "abs_error").mean())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--strategy", default="G1")
    parser.add_argument("--trials", type=int, default=30)
    args = parser.parse_args()
    store = load_store()
    splits = pd.read_csv(OUT / "evaluation_splits.csv.gz")
    choice = pd.read_csv(OUT / "selected_pairs.csv")
    records = []
    path = OUT / "hyperparameter_trials.csv"
    if path.exists():
        records = pd.read_csv(path).to_dict("records")
    for row in choice.itertuples(index=False):
        name = f"{args.strategy}_o{int(row.outer_fold)}_{row.model}_{row.feature}"
        study = optuna.create_study(
            study_name=name,
            direction="minimize",
            sampler=optuna.samplers.TPESampler(seed=20261125),
            storage="sqlite:///runs/model_development_v1/optuna.db",
            load_if_exists=True,
        )
        done = len([trial for trial in study.trials if trial.state == optuna.trial.TrialState.COMPLETE])
        remaining = max(0, args.trials - done)

        def objective(trial: optuna.Trial, outer=int(row.outer_fold), model=row.model, feature=row.feature) -> float:
            params = suggest(trial, model)
            value = pooled_mae(store, splits, outer, model, feature, params, args.strategy)
            trial.set_user_attr("params_resolved", json.dumps(params, default=str))
            return value

        if remaining:
            study.optimize(objective, n_trials=remaining, n_jobs=1)
        best = study.best_trial
        for trial in study.trials:
            if trial.state != optuna.trial.TrialState.COMPLETE:
                continue
            records.append(
                {
                    "strategy": args.strategy,
                    "outer_fold": int(row.outer_fold),
                    "model": row.model,
                    "feature": row.feature,
                    "trial": int(trial.number),
                    "mae_source": float(trial.value),
                    "params": trial.user_attrs.get("params_resolved", json.dumps(trial.params, default=str)),
                    "best": trial.number == best.number,
                }
            )
        pd.DataFrame(records).drop_duplicates(["strategy", "outer_fold", "model", "feature", "trial"]).to_csv(path, index=False)
        print(name, "best", best.value, "trials", len(study.trials), flush=True)
    frame = pd.DataFrame(records).drop_duplicates(["strategy", "outer_fold", "model", "feature", "trial"])
    frame["best"] = frame.groupby(["strategy", "outer_fold", "model", "feature"]).mae_source.transform(lambda values: values.eq(values.min()))
    frame.to_csv(path, index=False)


if __name__ == "__main__":
    main()
