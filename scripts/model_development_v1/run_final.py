"""Refit selected and fixed models on outer training pools and on evaluations A/B."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from model_lib import FIXED, OUT, feature_matrix, fit_predict, load_store, model_seed, preprocess, rows_for

PARTS = OUT / "final_parts"


def best_params() -> dict:
    trials = pd.read_csv(OUT / "hyperparameter_trials.csv")
    trials = trials[trials.strategy.eq("G1") & trials.best.astype(str).str.lower().eq("true")]
    trials = trials.sort_values("trial")
    chosen = {}
    for row in trials.itertuples(index=False):
        chosen.setdefault((int(row.outer_fold), row.model, row.feature), json.loads(row.params))
    return chosen


def predict_split(store, train: pd.DataFrame, test: pd.DataFrame, model: str, feature: str, params: dict, seed_key: str) -> pd.DataFrame:
    matrix, names = feature_matrix(store, feature)
    train_x, y, weights = rows_for(train, store, matrix)
    test_x, _, _ = rows_for(test, store, matrix)
    train_ready, test_ready, _ = preprocess(train_x, test_x, names, model in {"Ridge", "SVR"})
    seed = model_seed("final", seed_key, model, feature)
    test_pred, train_pred = fit_predict(model, train_ready, y, weights, [test_ready, train_ready], params, seed)
    output = test[["label_id", "Graph_ID", "Study_key", "y", "cohort", "MW"]].copy()
    output["y_pred"] = test_pred
    output["model"] = model
    output["feature"] = feature
    output["train_mae_label"] = float(np.mean(np.abs(train_pred - y)))
    output["n_train_labels"] = int(len(train))
    output["n_train_graphs"] = int(train.Graph_ID.nunique())
    output["n_train_studies"] = int(train.Study_key.nunique())
    return output


def main() -> None:
    PARTS.mkdir(parents=True, exist_ok=True)
    store = load_store()
    dev = store["dev"]
    splits = pd.read_csv(OUT / "evaluation_splits.csv.gz")
    test = pd.read_csv(Path(__file__).resolve().parents[2] / "data" / "splits_v1" / "test_labels.csv.gz")
    chosen = best_params()
    frames = []
    g1 = dev[dev.g1].copy()
    g1["assigned_fold"] = g1.assigned_fold.astype(int)
    for outer in range(5):
        train = g1[g1.assigned_fold.ne(outer)].copy()
        evaluation = test[test.assigned_fold.eq(outer)].copy()
        jobs = [
            ("RF", "F04", FIXED["RF"], "fixed_g1_rf"),
            ("Median", "F01", {}, "fixed_median"),
            ("Ridge", "F01", FIXED["Ridge"], "fixed_ridge"),
        ]
        for (fold, model, feature), params in chosen.items():
            if fold == outer:
                jobs.append((model, feature, params, "optuna"))
                jobs.append((model, feature, FIXED[model], "fixed_selected"))
        for model, feature, params, source in jobs:
            key = f"C_o{outer}_{source}_{model}_{feature}"
            path = PARTS / f"{key}.csv"
            if path.exists():
                frames.append(pd.read_csv(path))
                continue
            pred = predict_split(store, train, evaluation, model, feature, params, key)
            pred["evaluation"] = "C"
            pred["outer_fold"] = outer
            pred["parameter_source"] = source
            pred["exploratory"] = False
            pred.to_csv(path, index=False)
            frames.append(pred)
            print(key, "train_mae", round(float(pred.train_mae_label.iloc[0]), 4), flush=True)
    for evaluation_name in ["A", "B"]:
        part = splits[splits.evaluation.eq(evaluation_name)].copy()
        for fold in range(5):
            train = part[part.outer_fold.ne(fold)].copy()
            evaluation = part[part.outer_fold.eq(fold)].copy()
            for model in ["RF", "XGBoost", "LightGBM"]:
                key = f"{evaluation_name}_o{fold}_fixed_{model}_F05"
                path = PARTS / f"{key}.csv"
                if path.exists():
                    frames.append(pd.read_csv(path))
                    continue
                pred = predict_split(store, train, evaluation, model, "F05", FIXED[model], key)
                pred["evaluation"] = evaluation_name
                pred["outer_fold"] = fold
                pred["parameter_source"] = "fixed_ab"
                pred["exploratory"] = False
                pred.to_csv(path, index=False)
                frames.append(pred)
                print(key, flush=True)
    all_pred = pd.concat(frames, ignore_index=True)
    all_pred.to_csv(OUT / "oof_predictions.csv.gz", index=False)
    summary = (
        all_pred.groupby(["evaluation", "outer_fold", "model", "feature", "parameter_source"], as_index=False)
        .agg(train_mae_label=("train_mae_label", "first"), n_train_labels=("n_train_labels", "first"), n_train_graphs=("n_train_graphs", "first"), n_train_studies=("n_train_studies", "first"), n_test=("label_id", "size"))
    )
    summary.to_csv(OUT / "training_summary.csv", index=False)
    params_out = {f"{fold}|{model}|{feature}": params for (fold, model, feature), params in chosen.items()}
    (OUT / "selected_parameters.json").write_text(json.dumps(params_out, indent=2))
    print("FINAL_DONE", len(all_pred))


if __name__ == "__main__":
    main()
