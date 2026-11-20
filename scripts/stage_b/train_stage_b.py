"""Fit the locked Stage B models. Reuses identical training sets and resumes finished parts."""

from __future__ import annotations

import hashlib
import json
import os
import resource
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stageb_lib import DESC_NAMES, derive_seed  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
SPLIT = REPO / "data" / "splits_v1"
RUN = REPO / "runs" / "stage_b"
PARTS = RUN / "pred_parts"
N_JOBS = int(os.environ.get("STAGEB_NJOBS", "8"))


def cpu_time() -> float:
    self_use = resource.getrusage(resource.RUSAGE_SELF)
    child = resource.getrusage(resource.RUSAGE_CHILDREN)
    return (self_use.ru_utime + self_use.ru_stime) + (child.ru_utime + child.ru_stime)


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values, kind="mergesort")
    values = values[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    return float(values[min(np.searchsorted(cumulative, 0.5 * weights.sum(), side="left"), len(values) - 1)])


def signature(model: str, rf_seed: int, labels: pd.DataFrame, test_ids: list[str]) -> str:
    digest = hashlib.sha256()
    digest.update(f"{model}|{rf_seed}|{N_JOBS}".encode())
    for row in labels.sort_values("label_id").itertuples(index=False):
        digest.update(f"{row.label_id}:{row.sample_weight:.10f}|".encode())
    for label_id in sorted(test_ids):
        digest.update(label_id.encode())
    return digest.hexdigest()


def training_frame(usable: pd.DataFrame, selected: pd.DataFrame, fold: int, block: str, dropped: dict, pilot_studies: set[str]) -> pd.DataFrame:
    frames = []
    for component, part in selected.groupby("component"):
        graphs = set(part.Graph_ID)
        labels = usable[usable.Graph_ID.isin(graphs) & usable.assigned_fold.ne(fold)].copy()
        if block == "B04" and component == "M":
            labels = labels[labels.Study_key.ne(dropped[str(fold)])]
        if block == "PILOT":
            labels = labels[~labels.Study_key.isin(pilot_studies)]
        frames.append(labels)
    train = pd.concat(frames, ignore_index=True).drop_duplicates("label_id").sort_values("label_id")
    counts = train.groupby("Graph_ID").label_id.transform("count")
    train["sample_weight"] = 1.0 / counts.to_numpy()
    totals = train.groupby("Graph_ID").sample_weight.sum()
    if not np.allclose(totals.to_numpy(), 1.0):
        raise RuntimeError("graph weights do not sum to 1")
    if train.Graph_ID.nunique() != selected.Graph_ID.nunique():
        missing = set(selected.Graph_ID) - set(train.Graph_ID)
        raise RuntimeError(f"{block} fold {fold} lost graphs {list(missing)[:5]}")
    return train


def impute(train_x: np.ndarray, test_x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    train_x = train_x.astype(np.float64)
    test_x = test_x.astype(np.float64)
    train_x[~np.isfinite(train_x)] = np.nan
    test_x[~np.isfinite(test_x)] = np.nan
    median = np.nanmedian(train_x, axis=0)
    median = np.where(np.isfinite(median), median, 0.0)
    return np.where(np.isfinite(train_x), train_x, median), np.where(np.isfinite(test_x), test_x, median)


def predict_model(model: str, train: pd.DataFrame, test: pd.DataFrame, desc: np.ndarray, fp: np.ndarray, index: dict[str, int], rf_seed: int) -> np.ndarray:
    y = train.y.to_numpy(np.float64)
    weights = train.sample_weight.to_numpy(np.float64)
    if model == "Median":
        baseline = weighted_median(y, weights)
        return np.full(len(test), baseline)
    train_idx = train.Graph_ID.map(index).to_numpy()
    test_idx = test.Graph_ID.map(index).to_numpy()
    if model == "kNN":
        graph_y = train.groupby("Graph_ID").y.median()
        graph_ids = graph_y.index.tolist()
        train_fp = fp[[index[gid] for gid in graph_ids]].astype(np.float32)
        query_fp = fp[test_idx].astype(np.float32)
        intersection = query_fp @ train_fp.T
        denom = query_fp.sum(1, keepdims=True) + train_fp.sum(1) - intersection
        similarity = intersection / np.maximum(denom, 1e-12)
        medians = graph_y.to_numpy(np.float64)
        fallback = weighted_median(y, weights)
        output = np.empty(len(test), dtype=np.float64)
        order_ids = np.array(graph_ids)
        for row_i in range(len(test)):
            sims = similarity[row_i]
            order = sorted(range(len(graph_ids)), key=lambda j: (-float(sims[j]), order_ids[j]))
            chosen = [j for j in order if sims[j] > 0][:5]
            if not chosen:
                output[row_i] = fallback
            else:
                chosen_w = sims[chosen]
                output[row_i] = float(np.dot(chosen_w, medians[chosen]) / chosen_w.sum())
        return output
    raw_train = desc[train_idx]
    raw_test = desc[test_idx]
    x_train, x_test = impute(raw_train, raw_test)
    if model == "Ridge":
        scaler = StandardScaler()
        fitted = Ridge(alpha=1.0, fit_intercept=True)
        fitted.fit(scaler.fit_transform(x_train), y, sample_weight=weights)
        return fitted.predict(scaler.transform(x_test))
    if model == "RF":
        bits_train = fp[train_idx].astype(np.float64)
        bits_test = fp[test_idx].astype(np.float64)
        fitted = RandomForestRegressor(
            n_estimators=500,
            max_features=0.5,
            min_samples_leaf=2,
            max_depth=None,
            bootstrap=True,
            random_state=rf_seed,
            n_jobs=N_JOBS,
        )
        fitted.fit(np.hstack([x_train, bits_train]), y, sample_weight=weights)
        return fitted.predict(np.hstack([x_test, bits_test]))
    raise RuntimeError(model)


def jobs_for(decision: dict) -> list[tuple]:
    folds = range(decision["n_folds"])
    seeds = decision["seeds"]
    planned = [("PILOT", [0, 1], [17], ["RF"], ["B", "BS", "BM", "BN"])]
    planned.append(("B01", folds, seeds, ["Median", "kNN", "Ridge", "RF"], ["B", "BS", "BM", "BN", "ALL_DEV"]))
    if decision["b02"]["run"]:
        planned.append(("B02", folds, seeds, ["RF"], ["B", "BS", "BM", "BN"]))
    if decision["b03"]["run"]:
        planned.append(("B03", folds, seeds, ["RF"], ["B", "BS", "BM", "BN"]))
    if decision["b04"]["run"]:
        planned.append(("B04", folds, seeds, ["RF"], ["B", "BM", "BN"]))
    jobs = []
    for block, block_folds, block_seeds, models, arms in planned:
        for fold in block_folds:
            if block == "PILOT" and fold >= decision["n_folds"]:
                continue
            if block == "PILOT" and not any(item.get("outer_fold") == fold and item.get("run") for item in decision["pilot"]):
                continue
            for seed in block_seeds:
                for model in models:
                    for arm in arms:
                        if arm == "ALL_DEV":
                            continue
                        jobs.append((block, fold, seed, model, arm))
            if "ALL_DEV" in arms:
                for model in models:
                    jobs.append((block, fold, -1, model, "ALL_DEV"))
    return jobs


def main() -> None:
    PARTS.mkdir(parents=True, exist_ok=True)
    decision = json.loads((SPLIT / "budget_decision.json").read_text())
    features = np.load(RUN / "features_dev.npz", allow_pickle=True)
    index = {str(gid): i for i, gid in enumerate(features["graph_id"])}
    desc = features["desc"].astype(np.float64)
    fp = features["fp"]
    usable = pd.read_csv(SPLIT / "outer_members.csv.gz")
    usable = usable[usable.role.isin(["test", "train_eligible"])].copy()
    usable["assigned_fold"] = usable.assigned_fold.astype(int)
    test = pd.read_csv(SPLIT / "test_labels.csv.gz")
    test["assigned_fold"] = test.assigned_fold.astype(int)
    members = pd.read_csv(SPLIT / "sampling_members.csv.gz")
    dropped = decision.get("b04", {}).get("dropped_study_by_fold", {})
    pilot_studies = {}
    for fold in (0, 1):
        path = SPLIT / f"pilot_test_fold{fold}.csv.gz"
        pilot_studies[fold] = set(pd.read_csv(path, usecols=["Study_key"]).Study_key) if path.exists() else set()
        budgets = {
        "B01": (decision["b"], decision["k"]),
        "B02": (decision["b"], 25),
        "B03": (decision["b03"].get("budget_b"), decision["b03"].get("budget_k")),
        "B04": (decision["b04"].get("budget_b"), decision["b04"].get("budget_k")),
    }
    cache: dict[str, Path] = {}
    timing = []
    job_list = jobs_for(decision)
    print(f"JOBS {len(job_list)}", flush=True)
    started = time.perf_counter()
    cpu_started = cpu_time()
    for block, fold, seed, model, arm in job_list:
        path = PARTS / f"{block}_f{fold}_s{seed}_{model}_{arm}.csv"
        member_seed = -1 if arm == "ALL_DEV" else seed
        selected = members[
            members.block.eq(block) & members.fold.eq(fold) & members.seed.eq(member_seed) & members.arm.eq(arm)
        ]
        if selected.empty:
            raise RuntimeError(f"missing members for {path.name}")
        if block == "PILOT":
            evaluation = pd.read_csv(SPLIT / f"pilot_test_fold{fold}.csv.gz")
        else:
            evaluation = test[test.assigned_fold.eq(fold)].copy()
        train = training_frame(usable, selected, fold, block, dropped, pilot_studies.get(fold, set()))
        if block == "PILOT":
            pilot_budget = next(item for item in decision["pilot"] if item.get("outer_fold") == fold and item.get("run"))
            budget = (pilot_budget["b"], pilot_budget["k"])
        else:
            budget = budgets[block]
        rf_seed = derive_seed("Caco-2:AB", fold, seed, arm, model)
        sig = signature(model, rf_seed, train, evaluation.label_id.astype(str).tolist())
        if path.exists() and len(pd.read_csv(path)) == len(evaluation):
            cache.setdefault(sig, path)
            continue
        wall = time.perf_counter()
        cpu = cpu_time()
        if sig in cache:
            previous = pd.read_csv(cache[sig])
            prediction = evaluation[["label_id"]].merge(previous[["label_id", "y_pred"]], on="label_id", how="left")
            y_pred = prediction.y_pred.to_numpy()
            if not np.isfinite(y_pred).all():
                raise RuntimeError(f"cache miss alignment for {path.name}")
        else:
            y_pred = predict_model(model, train, evaluation, desc, fp, index, rf_seed)
        output = pd.DataFrame(
            {
                "endpoint": "Caco-2:AB",
                "block": block,
                "fold": fold,
                "seed": seed,
                "model": model,
                "arm": arm,
                "budget_b": budget[0],
                "budget_k": budget[1],
                "label_id": evaluation.label_id.to_numpy(),
                "graph_id": evaluation.Graph_ID.to_numpy(),
                "study_id": evaluation.Study_key.to_numpy(),
                "condition_id": evaluation.condition_id.to_numpy(),
                "y": evaluation.y.to_numpy(),
                "y_pred": y_pred,
            }
        )
        output.to_csv(path, index=False)
        cache[sig] = path
        record = {
            "file": path.name,
            "wall_s": round(time.perf_counter() - wall, 3),
            "cpu_s": round(cpu_time() - cpu, 3),
            "n_train_labels": int(len(train)),
            "n_train_graphs": int(train.Graph_ID.nunique()),
            "n_test": int(len(evaluation)),
        }
        timing.append(record)
        print(f"FIT {path.name} wall={record['wall_s']} cpu={record['cpu_s']} graphs={record['n_train_graphs']}", flush=True)
    summary = {
        "fits_recorded": len(timing),
        "wall_s": round(time.perf_counter() - started, 3),
        "cpu_s": round(cpu_time() - cpu_started, 3),
        "n_jobs": N_JOBS,
        "descriptor_names": DESC_NAMES,
        "parts": timing,
    }
    (RUN / "timing.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"TRAIN_DONE fits={len(timing)} wall={summary['wall_s']} cpu={summary['cpu_s']}", flush=True)


if __name__ == "__main__":
    main()
