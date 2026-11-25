"""Shared matrices, preprocessing, and model fitting for bRo5-Perm."""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVR
from xgboost import XGBRegressor

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "stage_b"))
from stageb_lib import DESC_NAMES, derive_seed  # noqa: E402

OUT = REPO / "runs" / "model_development_v1"
N_JOBS = 4
FIXED = {
    "RF": {"n_estimators": 500, "max_depth": None, "min_samples_leaf": 2, "max_features": 0.5},
    "XGBoost": {"n_estimators": 500, "max_depth": 6, "learning_rate": 0.05, "min_child_weight": 1, "subsample": 0.8, "colsample_bytree": 0.8, "reg_alpha": 0.0, "reg_lambda": 1.0},
    "LightGBM": {"n_estimators": 500, "learning_rate": 0.05, "num_leaves": 31, "max_depth": -1, "min_child_samples": 20, "subsample": 0.8, "colsample_bytree": 0.8, "reg_alpha": 0.0, "reg_lambda": 0.0},
    "SVR": {"C": 10.0, "gamma": "scale", "epsilon": 0.1},
    "Ridge": {"alpha": 1.0},
}


def load_store() -> dict:
    dev = pd.read_csv(OUT / "development_pool.csv.gz")
    rdkit = pd.read_csv(OUT / "descriptors_rdkit2d.csv.gz")
    stage = np.load(REPO / "runs" / "stage_b" / "features_dev.npz", allow_pickle=True)
    stage_index = {str(gid): i for i, gid in enumerate(stage["graph_id"])}
    order = rdkit.Graph_ID.astype(str).tolist()
    if any(gid not in stage_index for gid in order):
        missing = [gid for gid in order if gid not in stage_index]
        raise RuntimeError(f"{len(missing)} graphs lack Stage B features")
    stage_rows = [stage_index[gid] for gid in order]
    desc_names = [column for column in rdkit.columns if column != "Graph_ID"]
    return {
        "dev": dev,
        "graph_id": order,
        "row": {gid: i for i, gid in enumerate(order)},
        "f01": stage["desc"][stage_rows].astype(np.float64),
        "f01_names": list(DESC_NAMES),
        "f02": rdkit[desc_names].to_numpy(np.float64),
        "f02_names": desc_names,
        "fp": stage["fp"][stage_rows].astype(np.float64),
    }


def feature_matrix(store: dict, name: str) -> tuple[np.ndarray, list[str]]:
    fp_names = [f"ECFP4_{i}" for i in range(store["fp"].shape[1])]
    if name == "F01":
        return store["f01"], store["f01_names"]
    if name == "F02":
        return store["f02"], store["f02_names"]
    if name == "F03":
        return store["fp"], fp_names
    if name == "F04":
        return np.hstack([store["f01"], store["fp"]]), store["f01_names"] + fp_names
    if name == "F05":
        return np.hstack([store["f02"], store["fp"]]), store["f02_names"] + fp_names
    raise KeyError(name)


def preprocess(train_x: np.ndarray, test_x: np.ndarray, names: list[str], scale: bool) -> tuple[np.ndarray, np.ndarray, list[str]]:
    train_x = train_x.astype(np.float64, copy=True)
    test_x = test_x.astype(np.float64, copy=True)
    train_x[~np.isfinite(train_x)] = np.nan
    test_x[~np.isfinite(test_x)] = np.nan
    missing = np.mean(~np.isfinite(train_x), axis=0)
    keep = missing <= 0.20
    train_x = train_x[:, keep]
    test_x = test_x[:, keep]
    names = [name for name, flag in zip(names, keep) if flag]
    median = np.nanmedian(train_x, axis=0)
    median = np.where(np.isfinite(median), median, 0.0)
    train_filled = np.where(np.isfinite(train_x), train_x, median)
    test_filled = np.where(np.isfinite(test_x), test_x, median)
    std = train_filled.std(axis=0)
    vary = std >= 1e-6
    train_filled = train_filled[:, vary]
    test_filled = test_filled[:, vary]
    names = [name for name, flag in zip(names, vary) if flag]
    fingerprint = np.array([name.startswith("ECFP4_") for name in names])
    if scale:
        desc = ~fingerprint
        if desc.any():
            center = train_filled[:, desc].mean(axis=0)
            spread = train_filled[:, desc].std(axis=0)
            spread = np.where(spread >= 1e-6, spread, 1.0)
            train_filled[:, desc] = (train_filled[:, desc] - center) / spread
            test_filled[:, desc] = (test_filled[:, desc] - center) / spread
    if not np.isfinite(train_filled).all() or not np.isfinite(test_filled).all():
        raise RuntimeError("non-finite values remained after preprocessing")
    return train_filled, test_filled, names


def rows_for(frame: pd.DataFrame, store: dict, matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    counts = frame.groupby("Graph_ID").label_id.transform("count").to_numpy(np.float64)
    weights = 1.0 / counts
    idx = frame.Graph_ID.map(store["row"]).to_numpy()
    if np.isnan(idx.astype(float)).any():
        raise RuntimeError("missing feature row")
    return matrix[idx.astype(int)], frame.y.to_numpy(np.float64), weights


def fit_predict(model: str, train_x: np.ndarray, y: np.ndarray, weights: np.ndarray, test_matrices: list[np.ndarray], params: dict, seed: int) -> list[np.ndarray]:
    if model == "Median":
        order = np.argsort(y, kind="mergesort")
        cumulative = np.cumsum(weights[order])
        value = float(y[order][min(np.searchsorted(cumulative, 0.5 * weights.sum(), side="left"), len(y) - 1)])
        return [np.full(len(matrix), value) for matrix in test_matrices]
    if model == "Ridge":
        from sklearn.linear_model import Ridge

        fitted = Ridge(alpha=float(params.get("alpha", 1.0)))
        fitted.fit(train_x, y, sample_weight=weights)
        return [fitted.predict(matrix) for matrix in test_matrices]
    if model == "RF":
        depth = params.get("max_depth", None)
        fitted = RandomForestRegressor(
            n_estimators=int(params.get("n_estimators", 500)),
            max_depth=None if depth in (None, "None") else int(depth),
            min_samples_leaf=int(params.get("min_samples_leaf", 2)),
            min_samples_split=int(params.get("min_samples_split", 2)),
            max_features=float(params.get("max_features", 0.5)),
            bootstrap=True,
            random_state=seed,
            n_jobs=N_JOBS,
        )
        fitted.fit(train_x, y, sample_weight=weights)
        return [fitted.predict(matrix) for matrix in test_matrices]
    if model == "XGBoost":
        fitted = XGBRegressor(
            n_estimators=int(params["n_estimators"]),
            max_depth=int(params["max_depth"]),
            learning_rate=float(params["learning_rate"]),
            min_child_weight=float(params["min_child_weight"]),
            subsample=float(params["subsample"]),
            colsample_bytree=float(params["colsample_bytree"]),
            reg_alpha=float(params["reg_alpha"]),
            reg_lambda=float(params["reg_lambda"]),
            objective="reg:squarederror",
            tree_method="hist",
            random_state=seed,
            n_jobs=N_JOBS,
            verbosity=0,
        )
        fitted.fit(train_x, y, sample_weight=weights)
        return [fitted.predict(matrix) for matrix in test_matrices]
    if model == "LightGBM":
        depth = int(params["max_depth"])
        leaves = int(params["num_leaves"])
        if depth > 0:
            leaves = min(leaves, 2**depth - 1)
        fitted = LGBMRegressor(
            n_estimators=int(params["n_estimators"]),
            learning_rate=float(params["learning_rate"]),
            num_leaves=leaves,
            max_depth=depth,
            min_child_samples=int(params["min_child_samples"]),
            subsample=float(params["subsample"]),
            subsample_freq=1,
            colsample_bytree=float(params["colsample_bytree"]),
            reg_alpha=float(params["reg_alpha"]),
            reg_lambda=float(params["reg_lambda"]),
            random_state=seed,
            n_jobs=N_JOBS,
            verbose=-1,
        )
        fitted.fit(train_x, y, sample_weight=weights)
        return [fitted.predict(matrix) for matrix in test_matrices]
    if model == "SVR":
        gamma = params.get("gamma", "scale")
        fitted = SVR(kernel="rbf", C=float(params["C"]), gamma=gamma if gamma == "scale" else float(gamma), epsilon=float(params["epsilon"]))
        fitted.fit(train_x, y, sample_weight=weights)
        return [fitted.predict(matrix) for matrix in test_matrices]
    raise KeyError(model)


def model_seed(*parts: object) -> int:
    return derive_seed("bRo5-Perm", *parts)
