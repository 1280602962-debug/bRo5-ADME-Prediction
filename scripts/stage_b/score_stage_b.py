"""Score saved Stage B predictions. Does not refit models or change members."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[2]
SPLIT = REPO / "data" / "splits_v1"
RUN = REPO / "runs" / "stage_b"
SEEDS = [17, 29, 42, 71, 101, 137, 271, 314, 577, 911]
CONTRASTS = [
    ("P01", "BN", "BM"),
    ("BS_minus_B", "BS", "B"),
    ("BM_minus_B", "BM", "B"),
    ("BN_minus_B", "BN", "B"),
    ("BN_minus_BS", "BN", "BS"),
    ("BM_minus_BS", "BM", "BS"),
]


def load_predictions() -> pd.DataFrame:
    frames = [pd.read_csv(path) for path in sorted((RUN / "pred_parts").glob("*.csv"))]
    pred = pd.concat(frames, ignore_index=True)
    alldev = pred[pred.arm.eq("ALL_DEV")].copy()
    if len(alldev):
        expanded = []
        for seed in SEEDS:
            part = alldev.copy()
            part["seed"] = seed
            expanded.append(part)
        pred = pd.concat([pred[pred.arm.ne("ALL_DEV")], *expanded], ignore_index=True)
    pred["abs_error"] = (pred.y - pred.y_pred).abs()
    pred["signed_error"] = pred.y_pred - pred.y
    return pred


def source_table(frame: pd.DataFrame, value: str) -> pd.Series:
    step = frame.groupby(["study_id", "graph_id"], as_index=False)[value].mean()
    return step.groupby("study_id")[value].mean()


def graph_mean(frame: pd.DataFrame, value: str) -> float:
    step = frame.groupby(["graph_id", "study_id"], as_index=False)[value].mean()
    return float(step.groupby("graph_id")[value].mean().mean())


def components(test: pd.DataFrame) -> pd.Series:
    parent: dict[tuple, tuple] = {}

    def find(node: tuple) -> tuple:
        parent.setdefault(node, node)
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(left: tuple, right: tuple) -> None:
        root_left, root_right = find(left), find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    for row in test.itertuples(index=False):
        union(("study", row.Study_key), ("sc", row.Scaffold_connectivity_group))
    studies = sorted(set(test.Study_key))
    roots = {study: find(("study", study)) for study in studies}
    codes = {root: i for i, root in enumerate(sorted(set(roots.values()), key=str))}
    return pd.Series({study: codes[root] for study, root in roots.items()}, name="component")


def bootstrap_weights(study_component: pd.Series, rng: np.random.Generator) -> tuple[list[str], np.ndarray]:
    studies = list(study_component.index)
    comp = study_component.loc[studies].to_numpy()
    n_comp = int(comp.max()) + 1
    draws = rng.integers(0, n_comp, size=(10000, n_comp))
    counts = np.vstack([np.bincount(row, minlength=n_comp) for row in draws])
    return studies, counts[:, comp]


def regression_stats(part: pd.DataFrame) -> dict:
    error = part.y_pred.to_numpy() - part.y.to_numpy()
    total = float(np.sum((part.y - part.y.mean()) ** 2))
    return {
        "rmse": float(np.sqrt(np.mean(error ** 2))),
        "r2": float(1 - np.sum(error ** 2) / total) if total > 0 else np.nan,
        "spearman": float(spearmanr(part.y, part.y_pred).statistic),
        "signed_label": float(np.mean(error)),
    }


def chemical_space(desc: pd.DataFrame, fp: np.ndarray, ids: list[str], members: pd.DataFrame, test: pd.DataFrame, outer: pd.DataFrame) -> pd.DataFrame:
    index = {gid: i for i, gid in enumerate(ids)}
    rows = []
    selected = members[members.block.eq("B01") & members.component.isin(["M", "N"]) & members.arm.isin(["BM", "BN"])]
    selected = selected.merge(desc, left_on="Graph_ID", right_on="Graph_ID", how="left")
    for component, part in selected.groupby("component"):
        rows.append(
            {
                "record": "donor_distribution",
                "component": component,
                "n_draws": int(len(part)),
                "n_unique_graphs": int(part.Graph_ID.nunique()),
                "MW_median": float(part.MW_x.median()) if "MW_x" in part.columns else float(part.MW.median()),
                "TPSA_median": float(part.TPSA.median()),
                "MolLogP_median": float(part.MolLogP.median()),
                "HBD_median": float(part.HBD.median()),
                "HBA_median": float(part.HBA.median()),
                "RotB_median": float(part.NumRotatableBonds.median()),
            }
        )
    usable = outer[outer.role.isin(["test", "train_eligible"])]
    nn_rows = []
    primary = members[members.block.eq("B01") & members.arm.isin(["B", "BS", "BM", "BN"])]
    for (fold, seed, arm), part in primary.groupby(["fold", "seed", "arm"]):
        train_ids = [gid for gid in part.Graph_ID.unique() if gid in index]
        test_ids = [gid for gid in test.loc[test.assigned_fold.eq(fold), "Graph_ID"].unique() if gid in index]
        if not train_ids or not test_ids:
            continue
        train_fp = fp[[index[gid] for gid in train_ids]].astype(np.float32)
        query_fp = fp[[index[gid] for gid in test_ids]].astype(np.float32)
        intersection = query_fp @ train_fp.T
        denom = query_fp.sum(1, keepdims=True) + train_fp.sum(1) - intersection
        nearest = (intersection / np.maximum(denom, 1e-12)).max(1)
        nn_rows.append({"fold": fold, "seed": seed, "arm": arm, "nearest_mean": float(nearest.mean()), "nearest_median": float(np.median(nearest))})
    nn = pd.DataFrame(nn_rows)
    for arm, part in nn.groupby("arm"):
        rows.append(
            {
                "record": "train_test_nearest_tanimoto",
                "arm": arm,
                "nearest_mean": float(part.nearest_mean.mean()),
                "nearest_mean_seed_std": float(part.groupby("seed").nearest_mean.mean().std(ddof=1)),
                "nearest_median": float(part.nearest_median.mean()),
            }
        )
    label_y = usable.groupby("Graph_ID").y.mean()
    for component, part in selected.groupby("component"):
        rows.append({"record": "donor_label_y", "component": component, "y_median": float(label_y.loc[part.Graph_ID.unique()].median()), "y_mean": float(label_y.loc[part.Graph_ID.unique()].mean())})
    return pd.DataFrame(rows)


def main() -> None:
    pred = load_predictions()
    formal = pred[pred.block.ne("PILOT")].copy()
    formal.to_csv(RUN / "oof_predictions.csv.gz", index=False)
    test = pd.read_csv(SPLIT / "test_labels.csv.gz")
    study_component = components(test)
    studies, weights = bootstrap_weights(study_component, np.random.default_rng(271828))
    pooled_rows = []
    fold_rows = []
    delta_rows = []
    component_rows = []
    n_comp = int(study_component.max()) + 1
    largest = study_component.value_counts().iloc[0]
    largest_id = int(study_component.value_counts().index[0])
    for (block, model, arm), part in formal.groupby(["block", "model", "arm"]):
        seed_mean = part.groupby("label_id", as_index=False).agg(
            abs_error=("abs_error", "mean"),
            signed_error=("signed_error", "mean"),
            y=("y", "first"),
            y_pred=("y_pred", "mean"),
            study_id=("study_id", "first"),
            graph_id=("graph_id", "first"),
            fold=("fold", "first"),
        )
        per_seed = [regression_stats(chunk) for _, chunk in part.groupby("seed")]
        stat = pd.DataFrame(per_seed)
        ensemble = seed_mean.assign(abs_error=(seed_mean.y - seed_mean.y_pred).abs(), signed_error=seed_mean.y_pred - seed_mean.y)
        ensemble_stat = regression_stats(seed_mean.rename(columns={"y_pred": "y_pred"}))
        pooled_rows.append(
            {
                "block": block,
                "model": model,
                "arm": arm,
                "mae_source": float(source_table(seed_mean, "abs_error").mean()),
                "mae_graph": graph_mean(seed_mean, "abs_error"),
                "signed_source": float(source_table(seed_mean, "signed_error").mean()),
                "rmse_seed_mean": float(stat.rmse.mean()),
                "rmse_seed_std": float(stat.rmse.std(ddof=1)),
                "r2_seed_mean": float(stat.r2.mean()),
                "r2_seed_std": float(stat.r2.std(ddof=1)),
                "spearman_seed_mean": float(stat.spearman.mean()),
                "spearman_seed_std": float(stat.spearman.std(ddof=1)),
                "signed_label_seed_mean": float(stat.signed_label.mean()),
                "ensemble_mae_source": float(source_table(ensemble, "abs_error").mean()),
                "ensemble_rmse": ensemble_stat["rmse"],
                "ensemble_r2": ensemble_stat["r2"],
                "ensemble_spearman": ensemble_stat["spearman"],
                "n_labels": int(seed_mean.label_id.nunique()),
                "n_graphs": int(seed_mean.graph_id.nunique()),
                "n_studies": int(seed_mean.study_id.nunique()),
            }
        )
        for (fold, seed), chunk in part.groupby(["fold", "seed"]):
            fold_rows.append(
                {
                    "block": block,
                    "model": model,
                    "arm": arm,
                    "fold": int(fold),
                    "seed": int(seed),
                    "mae_source": float(source_table(chunk.assign(abs_error=chunk.abs_error), "abs_error").mean()),
                    "rmse_label": float(np.sqrt(np.mean((chunk.y_pred - chunk.y) ** 2))),
                    "n_labels": int(len(chunk)),
                    "n_graphs": int(chunk.graph_id.nunique()),
                    "n_studies": int(chunk.study_id.nunique()),
                }
            )
    pooled = pd.DataFrame(pooled_rows)
    by_fold = pd.DataFrame(fold_rows)
    for (block, model), part in formal.groupby(["block", "model"]):
        seed_mean = part.groupby(["arm", "label_id"], as_index=False).agg(
            abs_error=("abs_error", "mean"), study_id=("study_id", "first"), graph_id=("graph_id", "first")
        )
        study_mae = {arm: source_table(chunk, "abs_error") for arm, chunk in seed_mean.groupby("arm")}
        arm_draws = {}
        for arm, values in study_mae.items():
            vector = values.loc[studies].to_numpy(np.float64)
            arm_draws[arm] = (weights * vector).sum(1) / weights.sum(1)
        point = {arm: float(values.mean()) for arm, values in study_mae.items()}
        for name, left, right in CONTRASTS:
            if left not in arm_draws or right not in arm_draws:
                continue
            draw = arm_draws[left] - arm_draws[right]
            estimate = point[left] - point[right]
            low, high = np.quantile(draw, [0.025, 0.975])
            delta_rows.append(
                {
                    "block": block,
                    "model": model,
                    "comparison": name,
                    "weighting": "source_equal_seed_mean_abs_error",
                    "estimate": estimate,
                    "ci_low": float(low),
                    "ci_high": float(high),
                    "ci_width": float(high - low),
                    "n_studies": int(len(study_component)),
                    "n_components": n_comp,
                    "max_component_studies": int(largest),
                    "max_component_study_fraction": float(largest / len(study_component)),
                    "interval_status": "exploratory_few_components" if n_comp < 10 else "reported_interval",
                    "primary": bool(block == "B01" and model == "RF" and name == "P01"),
                }
            )
        if "BN" in point and "BM" in point:
            kept = study_component[study_component.ne(largest_id)]
            reduced = {arm: values.loc[values.index.isin(kept.index)].mean() for arm, values in study_mae.items()}
            component_rows.append(
                {
                    "block": block,
                    "model": model,
                    "n_studies": int(len(study_component)),
                    "n_components": n_comp,
                    "max_component_studies": int(largest),
                    "max_component_study_fraction": float(largest / len(study_component)),
                    "p01_full": point["BN"] - point["BM"],
                    "p01_without_largest_component": reduced["BN"] - reduced["BM"],
                }
            )
    deltas = pd.DataFrame(delta_rows)
    components_out = pd.DataFrame(component_rows)
    pooled.to_csv(RUN / "metrics_pooled.csv", index=False)
    by_fold.to_csv(RUN / "metrics_by_fold.csv", index=False)
    deltas.to_csv(RUN / "paired_deltas.csv", index=False)
    components_out.to_csv(RUN / "source_component_summary.csv", index=False)
    features = np.load(RUN / "features_dev.npz", allow_pickle=True)
    desc = pd.read_csv(RUN / "descriptors_dev.csv.gz")
    members = pd.read_csv(SPLIT / "sampling_members.csv.gz")
    outer = pd.read_csv(SPLIT / "outer_members.csv.gz")
    chemical_space(desc, features["fp"], [str(gid) for gid in features["graph_id"]], members, test, outer).to_csv(
        RUN / "chemical_space.csv", index=False
    )
    primary = deltas[(deltas.block.eq("B01")) & (deltas.model.eq("RF")) & (deltas.comparison.eq("P01"))]
    print(primary.to_string(index=False))
    print(pooled[(pooled.block.eq("B01")) & (pooled.model.eq("RF"))][["arm", "mae_source", "mae_graph", "rmse_seed_mean", "r2_seed_mean", "spearman_seed_mean"]].to_string(index=False))


if __name__ == "__main__":
    main()
