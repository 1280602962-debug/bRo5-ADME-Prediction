"""Shared matrices and fixed random forests for M5. Seeds are explicit."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "model_development_v1"))
sys.path.insert(0, str(REPO / "scripts" / "stage_b"))
from model_lib import feature_matrix, load_store, preprocess, rows_for  # noqa: E402
from score_stage_b import bootstrap_weights, components, graph_mean, source_table  # noqa: E402
from scipy.stats import pearsonr, spearmanr

M4 = REPO / "runs" / "model_development_m4_cpu"
OUT = REPO / "runs" / "model_development_m5"
SEEDS = (17, 42, 101)
RF_PARAMS = {
    "n_estimators": 500,
    "max_features": 0.5,
    "min_samples_leaf": 2,
    "max_depth": None,
    "bootstrap": True,
    "n_jobs": 4,
}
STATIC_COLUMNS = [
    "static_Rgyr",
    "static_NPR1",
    "static_NPR2",
    "static_MolecularVolume",
    "static_MaximumDimension",
    "static_TotalSASA",
    "static_PolarSASA_NOH",
    "static_PolarFraction",
    "static_IMHB_count",
]
ENSEMBLE_COLUMNS = [
    "ens_Rgyr_median",
    "ens_Rgyr_p10",
    "ens_Rgyr_iqr",
    "ens_PolarSASA_median",
    "ens_PolarSASA_p10",
    "ens_PolarSASA_iqr",
    "ens_PolarFraction_median",
    "ens_PolarFraction_p10",
    "ens_IMHB_mean",
    "ens_IMHB_max",
    "ens_IMHB_fraction",
    "ens_NPR1_median",
    "ens_NPR2_median",
    "ens_Volume_median",
    "ens_MaximumDimension_median",
]
GEOMETRY_COLUMNS = [
    "ens_Rgyr_median",
    "ens_Rgyr_p10",
    "ens_Rgyr_iqr",
    "ens_NPR1_median",
    "ens_NPR2_median",
    "ens_Volume_median",
    "ens_MaximumDimension_median",
]
POLAR_COLUMNS = [
    "ens_PolarSASA_median",
    "ens_PolarSASA_p10",
    "ens_PolarSASA_iqr",
    "ens_PolarFraction_median",
    "ens_PolarFraction_p10",
]
IMHB_COLUMNS = ["ens_IMHB_mean", "ens_IMHB_max", "ens_IMHB_fraction"]


def load_tables(store: dict) -> dict:
    static = pd.read_csv(M4 / "l1_static_features.csv.gz")
    ensemble = pd.read_csv(M4 / "l1_ensemble_features.csv.gz")
    for frame in (static, ensemble):
        frame["Graph_ID"] = frame.Graph_ID.astype(str)
    holdout = set(pd.read_csv(REPO / "data" / "splits_v1" / "holdout_members.csv").Graph_ID.astype(str))
    if holdout & set(static.Graph_ID) or holdout & set(ensemble.Graph_ID):
        raise RuntimeError("holdout graph is in an M4 feature table")
    order = store["graph_id"]

    def matrix(frame: pd.DataFrame, columns: list[str]) -> np.ndarray:
        kept = frame.drop_duplicates("Graph_ID").set_index("Graph_ID")
        return kept.reindex(order)[columns].to_numpy(np.float64)

    return {
        "ensemble_ids": set(ensemble.Graph_ID),
        "static_ids": set(static.Graph_ID),
        "static": matrix(static, STATIC_COLUMNS),
        "ensemble": matrix(ensemble, ENSEMBLE_COLUMNS),
    }


def arm_matrix(store: dict, tables: dict, arm: str) -> tuple[np.ndarray, list[str]]:
    groups = {
        "R0": ("F04", []),
        "R1": ("F05", []),
        "R2": ("F04", STATIC_COLUMNS),
        "R3": ("F04", ENSEMBLE_COLUMNS),
        "geometry": ("F04", GEOMETRY_COLUMNS),
        "polar": ("F04", POLAR_COLUMNS),
        "imhb": ("F04", IMHB_COLUMNS),
    }
    base, columns = groups[arm]
    matrix, names = feature_matrix(store, base)
    if not columns:
        return matrix, names
    extra = tables["static"] if columns == STATIC_COLUMNS else tables["ensemble"]
    positions = [ENSEMBLE_COLUMNS.index(name) for name in columns] if columns != STATIC_COLUMNS else list(range(len(STATIC_COLUMNS)))
    if columns == STATIC_COLUMNS:
        block = extra
    else:
        block = extra[:, positions]
    return np.hstack([matrix, block]), names + columns


def folds(store: dict, graphs: set[str]) -> list[tuple[int, pd.DataFrame, pd.DataFrame]]:
    dev = store["dev"].copy()
    dev["Graph_ID"] = dev.Graph_ID.astype(str)
    dev["assigned_fold"] = dev.assigned_fold.astype(int)
    g1 = dev[dev.g1].copy()
    test = pd.read_csv(REPO / "data" / "splits_v1" / "test_labels.csv.gz")
    test["Graph_ID"] = test.Graph_ID.astype(str)
    test["assigned_fold"] = test.assigned_fold.astype(int)
    rows = []
    for outer in range(5):
        train = g1[g1.assigned_fold.ne(outer) & g1.Graph_ID.isin(graphs)].sort_values("label_id")
        evaluation = test[test.assigned_fold.eq(outer) & test.Graph_ID.isin(graphs)].sort_values("label_id")
        rows.append((outer, train, evaluation))
    return rows


def fit_fold(store, train, test, matrix, names, seed: int, path: Path) -> pd.DataFrame:
    if path.exists():
        return pd.read_csv(path)
    train = train.sort_values("label_id")
    test = test.sort_values("label_id")
    train_x, y, weights = rows_for(train, store, matrix)
    test_x, _, _ = rows_for(test, store, matrix)
    train_ready, test_ready, _ = preprocess(train_x, test_x, names, False)
    model = RandomForestRegressor(
        n_estimators=RF_PARAMS["n_estimators"],
        max_features=RF_PARAMS["max_features"],
        min_samples_leaf=RF_PARAMS["min_samples_leaf"],
        max_depth=RF_PARAMS["max_depth"],
        bootstrap=RF_PARAMS["bootstrap"],
        random_state=int(seed),
        n_jobs=RF_PARAMS["n_jobs"],
    )
    model.fit(train_ready, y, sample_weight=weights)
    output = test[["label_id", "Graph_ID", "Study_key", "Scaffold_connectivity_group", "y", "cohort", "MW"]].copy()
    output["y_pred"] = model.predict(test_ready)
    output["seed"] = int(seed)
    output["n_train_labels"] = int(len(train))
    path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(path, index=False)
    train_ids = train.label_id.astype(str).tolist()
    (path.with_suffix(".train_ids.txt")).write_text("\n".join(train_ids))
    return output


def regression(frame: pd.DataFrame) -> dict:
    work = frame.rename(columns={"Graph_ID": "graph_id", "Study_key": "study_id"})
    work["abs_error"] = (work.y_pred - work.y).abs()
    work["signed_error"] = work.y_pred - work.y
    error = work.signed_error.to_numpy()
    total = float(np.sum((work.y - work.y.mean()) ** 2))
    return {
        "mae_source": float(source_table(work, "abs_error").mean()),
        "mae_graph": graph_mean(work, "abs_error"),
        "mae_label": float(work.abs_error.mean()),
        "rmse": float(np.sqrt(np.mean(error ** 2))),
        "r2": float(1 - np.sum(error ** 2) / total) if total else np.nan,
        "spearman": float(spearmanr(work.y, work.y_pred).statistic) if len(work) > 2 else np.nan,
        "pearson": float(pearsonr(work.y, work.y_pred).statistic) if len(work) > 2 else np.nan,
        "n_labels": int(work.label_id.nunique()),
        "n_graphs": int(work.graph_id.nunique()),
        "n_studies": int(work.study_id.nunique()),
    }


def paired_interval(left: pd.DataFrame, right: pd.DataFrame, value: str = "abs_error") -> dict:
    merged = left[["label_id", "y", "Study_key", "Graph_ID", "Scaffold_connectivity_group", value]].merge(
        right[["label_id", value]], on="label_id", suffixes=("_left", "_right")
    )
    if len(merged) != len(left) or len(merged) != len(right):
        raise RuntimeError("paired frames do not share labels")
    named = merged.rename(columns={"Graph_ID": "graph_id", "Study_key": "study_id"})
    left_source = source_table(named.assign(abs_error=named[f"{value}_left"]), "abs_error")
    right_source = source_table(named.assign(abs_error=named[f"{value}_right"]), "abs_error")
    studies = left_source.index.intersection(right_source.index)
    point = float(left_source.loc[studies].mean() - right_source.loc[studies].mean())
    study_component = components(merged)
    rng = np.random.default_rng(271828)
    study_list, weights = bootstrap_weights(study_component, rng)
    vector = (left_source.loc[study_list] - right_source.loc[study_list]).to_numpy(np.float64)
    draw = (weights * vector).sum(1) / weights.sum(1)
    low, high = np.quantile(draw, [0.025, 0.975])
    return {
        "estimate": point,
        "ci95_low": float(low),
        "ci95_high": float(high),
        "ci95_excludes_zero": bool(low > 0 or high < 0),
        "n_labels": int(len(merged)),
        "n_studies": int(len(studies)),
        "n_components": int(study_component.max()) + 1,
    }
