"""Full-coverage reliability. Thresholds come from training-side scores only."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge

from m5_lib import OUT, RF_PARAMS, regression
from model_lib import feature_matrix, load_store, preprocess, rows_for

PARTS = OUT / "reliability_parts"
SEED = 17
COVERAGES = (1.0, 0.8, 0.6, 0.4)
FEATURES = ["risk_neighbor", "risk_tree", "risk_gap", "prediction"]


def fit_feature(store, train, test, feature: str, path):
    if path.exists():
        saved = pd.read_csv(path)
        if "tree_std" in saved.columns:
            return saved
    matrix, names = feature_matrix(store, feature)
    train = train.sort_values("label_id")
    test = test.sort_values("label_id")
    train_x, y, weights = rows_for(train, store, matrix)
    test_x, _, _ = rows_for(test, store, matrix)
    train_ready, test_ready, _ = preprocess(train_x, test_x, names, False)
    model = RandomForestRegressor(random_state=SEED, **RF_PARAMS)
    model.fit(train_ready, y, sample_weight=weights)
    stacked = np.vstack([tree.predict(test_ready) for tree in model.estimators_])
    output = test[["label_id", "Graph_ID", "Study_key", "Scaffold_connectivity_group", "y", "cohort", "MW"]].copy()
    output["y_pred"] = model.predict(test_ready)
    output["tree_std"] = stacked.std(axis=0)
    output["seed"] = SEED
    output["n_train_labels"] = int(len(train))
    path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(path, index=False)
    path.with_suffix(".train_ids.txt").write_text("\n".join(train.label_id.astype(str)))
    return output


def neighbor_risk(store, train_graphs: set[str], query_graphs: list[str]) -> dict[str, float]:
    fp = store["fp"] > 0
    row = store["row"]
    train_idx = [int(row[graph]) for graph in train_graphs]
    train_fp = fp[train_idx]
    train_count = train_fp.sum(axis=1)
    scores = {}
    for graph in query_graphs:
        query = fp[int(row[graph])]
        intersection = train_fp[:, query].sum(axis=1)
        union = train_count + query.sum() - intersection
        tanimoto = np.divide(intersection, union, out=np.zeros(len(union)), where=union > 0)
        scores[graph] = float(1 - tanimoto.max()) if len(tanimoto) else 1.0
    return scores


def standardize(train: np.ndarray, test: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    center = np.nanmean(train, axis=0)
    spread = np.nanstd(train, axis=0)
    spread[spread < 1e-8] = 1.0
    return (train - center) / spread, (test - center) / spread


def graph_table(frame: pd.DataFrame, risk: str) -> pd.DataFrame:
    return frame.groupby("Graph_ID", as_index=False).agg(risk=(risk, "mean"), Study_key=("Study_key", "first"))


def source_metrics(frame: pd.DataFrame) -> dict:
    if frame.empty:
        return {"mae_source": np.nan, "spearman": np.nan, "n_labels": 0, "n_graphs": 0, "n_studies": 0}
    row = regression(frame)
    return {key: row[key] for key in ["mae_source", "spearman", "n_labels", "n_graphs", "n_studies"]}


def main() -> None:
    store = load_store()
    graphs = set(store["graph_id"])
    dev = store["dev"].copy()
    dev["Graph_ID"] = dev.Graph_ID.astype(str)
    dev["assigned_fold"] = dev.assigned_fold.astype(int)
    g1 = dev[dev.g1].copy()
    test = pd.read_csv(OUT.parents[1] / "data" / "splits_v1" / "test_labels.csv.gz")
    test["Graph_ID"] = test.Graph_ID.astype(str)
    test["assigned_fold"] = test.assigned_fold.astype(int)
    splits = pd.read_csv(OUT.parents[1] / "runs" / "model_development_v1" / "evaluation_splits.csv.gz")
    splits = splits[(splits.evaluation.eq("C_inner")) & splits.Graph_ID.astype(str).isin(graphs)].copy()
    splits["Graph_ID"] = splits.Graph_ID.astype(str)
    outer_rows = []
    inner_rows = []
    for outer in range(5):
        train = g1[g1.assigned_fold.ne(outer) & g1.Graph_ID.isin(graphs)].sort_values("label_id")
        evaluation = test[test.assigned_fold.eq(outer) & test.Graph_ID.isin(graphs)].sort_values("label_id")
        fitted = {}
        for feature in ("F04", "F05"):
            fitted[feature] = fit_feature(store, train, evaluation, feature, PARTS / f"outer_o{outer}_{feature}.csv")
            print(json.dumps({"split": "outer", "outer": outer, "feature": feature, "test": int(len(evaluation))}), flush=True)
        neighbors = neighbor_risk(store, set(train.Graph_ID), evaluation.Graph_ID.drop_duplicates().tolist())
        merged = fitted["F05"].merge(fitted["F04"][["label_id", "y_pred", "tree_std"]], on="label_id", suffixes=("_f05", "_f04"))
        merged["outer_fold"] = outer
        merged["risk_neighbor"] = merged.Graph_ID.map(neighbors)
        merged["risk_tree"] = merged.tree_std_f05
        merged["risk_gap"] = (merged.y_pred_f04 - merged.y_pred_f05).abs()
        merged["prediction"] = merged.y_pred_f05
        merged["y_pred"] = merged.y_pred_f05
        merged["abs_error"] = (merged.y_pred - merged.y).abs()
        outer_rows.append(merged)
        for inner in range(3):
            part = splits[(splits.outer_fold.eq(outer)) & (splits.inner_fold.eq(inner))]
            inner_train = part[part.role.eq("train_candidate")].sort_values("label_id")
            validation = part[part.role.eq("validation")].sort_values("label_id")
            inner_fit = {}
            for feature in ("F04", "F05"):
                inner_fit[feature] = fit_feature(store, inner_train, validation, feature, PARTS / f"inner_o{outer}_i{inner}_{feature}.csv")
            inner_neighbors = neighbor_risk(store, set(inner_train.Graph_ID), validation.Graph_ID.drop_duplicates().tolist())
            one = inner_fit["F05"].merge(inner_fit["F04"][["label_id", "y_pred", "tree_std"]], on="label_id", suffixes=("_f05", "_f04"))
            one["outer_fold"] = outer
            one["inner_fold"] = inner
            one["risk_neighbor"] = one.Graph_ID.map(inner_neighbors)
            one["risk_tree"] = one.tree_std_f05
            one["risk_gap"] = (one.y_pred_f04 - one.y_pred_f05).abs()
            one["prediction"] = one.y_pred_f05
            one["abs_error"] = (one.y_pred_f05 - one.y).abs()
            inner_rows.append(one)
            print(json.dumps({"split": "inner", "outer": outer, "inner": inner, "validation": int(len(validation))}), flush=True)
    outer = pd.concat(outer_rows, ignore_index=True)
    inner = pd.concat(inner_rows, ignore_index=True)
    outer.to_csv(OUT / "reliability_outer_predictions.csv.gz", index=False)
    scored = []
    for outer_id, part in inner.groupby("outer_fold"):
        risks = []
        for held in range(3):
            train_part = part[part.inner_fold.ne(held)]
            test_part = part[part.inner_fold.eq(held)].copy()
            train_x, test_x = standardize(train_part[FEATURES].to_numpy(np.float64), test_part[FEATURES].to_numpy(np.float64))
            model = Ridge(alpha=1.0)
            model.fit(train_x, train_part.abs_error.to_numpy(np.float64))
            test_part["crossfit_risk"] = model.predict(test_x)
            risks.append(test_part)
        cross = pd.concat(risks, ignore_index=True)
        train_x, test_x = standardize(cross[FEATURES].to_numpy(np.float64), outer.loc[outer.outer_fold.eq(outer_id), FEATURES].to_numpy(np.float64))
        model = Ridge(alpha=1.0)
        model.fit(train_x, cross.abs_error.to_numpy(np.float64))
        block = outer.loc[outer.outer_fold.eq(outer_id)].copy()
        block["model_risk"] = model.predict(test_x)
        graphs = graph_table(cross, "crossfit_risk")
        for coverage in COVERAGES:
            test_graphs = graph_table(block, "model_risk")
            test_threshold = float(graphs.risk.quantile(coverage))
            kept_graphs = set(test_graphs.loc[test_graphs.risk.le(test_threshold), "Graph_ID"])
            block[f"keep_{coverage}"] = block.Graph_ID.isin(kept_graphs)
        for feature in ("risk_neighbor", "risk_tree", "risk_gap"):
            feature_graphs = graph_table(cross, feature)
            for coverage in COVERAGES:
                cutoff = float(feature_graphs.risk.quantile(coverage))
                kept = set(graph_table(block, feature).loc[lambda frame: frame.risk.le(cutoff), "Graph_ID"])
                block[f"keep_{feature}_{coverage}"] = block.Graph_ID.isin(kept)
        association = float(pd.Series(cross.crossfit_risk).corr(cross.abs_error, method="spearman"))
        block["training_risk_error_spearman"] = association
        scored.append(block)
        print(json.dumps({"outer": int(outer_id), "training_risk_error_spearman": association}), flush=True)
    scored = pd.concat(scored, ignore_index=True)
    rows = []
    full = source_metrics(scored)
    rows.append({"method": "F05_full_coverage", "coverage_target": 1.0, "realized_molecule_coverage": 1.0, **full, "interpretation": "full_coverage_baseline"})
    rng = np.random.default_rng(271828)
    for method, column in [("ridge_crossfit", "keep_"), ("nearest_neighbor", "keep_risk_neighbor_"), ("tree_dispersion", "keep_risk_tree_"), ("f04_f05_gap", "keep_risk_gap_")]:
        for coverage in COVERAGES:
            key = f"{column}{coverage}"
            kept = scored[scored[key]]
            metrics = source_metrics(kept)
            realized = kept.Graph_ID.nunique() / scored.Graph_ID.nunique()
            rows.append({"method": method, "coverage_target": coverage, "realized_molecule_coverage": realized, **metrics, "interpretation": "rejected_subset_not_full_coverage_accuracy"})
            random_maes = []
            random_ranks = []
            random_labels = []
            random_graphs = []
            random_studies = []
            for _ in range(150):
                chosen = []
                for _, fold in scored.groupby("outer_fold"):
                    molecules = fold.Graph_ID.drop_duplicates().to_numpy()
                    n_keep = int(fold.loc[fold[key], "Graph_ID"].nunique())
                    pick = rng.choice(molecules, size=min(n_keep, len(molecules)), replace=False)
                    chosen.extend(pick.tolist())
                sample = scored[scored.Graph_ID.isin(chosen)]
                sample_metrics = source_metrics(sample)
                random_maes.append(sample_metrics["mae_source"])
                random_ranks.append(sample_metrics["spearman"])
                random_labels.append(sample_metrics["n_labels"])
                random_graphs.append(sample_metrics["n_graphs"])
                random_studies.append(sample_metrics["n_studies"])
            rows.append(
                {
                    "method": f"random_same_coverage::{method}",
                    "coverage_target": coverage,
                    "realized_molecule_coverage": realized,
                    "mae_source": float(np.mean(random_maes)),
                    "mae_source_low": float(np.quantile(random_maes, 0.025)),
                    "mae_source_high": float(np.quantile(random_maes, 0.975)),
                    "spearman": float(np.nanmean(random_ranks)),
                    "n_labels": float(np.mean(random_labels)),
                    "n_graphs": float(np.mean(random_graphs)),
                    "n_studies": float(np.mean(random_studies)),
                    "random_repeats": 150,
                    "interpretation": "random_rejection_reference_counts_are_means",
                }
            )
    p4_path = OUT / "seed_predictions.csv.gz"
    if p4_path.exists():
        p4 = pd.read_csv(p4_path)
        p4 = p4[(p4.arm.eq("R3")) & (p4.seed.eq(SEED))]
        both = scored.merge(p4[["label_id", "y_pred"]], on="label_id", how="left", suffixes=("", "_p4"))
        covered = both[both.y_pred_p4.notna()]
        rows.append(
            {
                "method": "P4_minus_F05_coverage_only",
                "coverage_target": np.nan,
                "realized_molecule_coverage": covered.Graph_ID.nunique() / scored.Graph_ID.nunique(),
                "n_labels": int(covered.label_id.nunique()),
                "n_graphs": int(covered.Graph_ID.nunique()),
                "n_studies": int(covered.Study_key.nunique()),
                "mean_abs_prediction_difference": float((covered.y_pred_p4 - covered.y_pred).abs().mean()) if len(covered) else np.nan,
                "interpretation": "feature_coverage_reported_separately_not_imputed",
            }
        )
    metrics = pd.DataFrame(rows)
    metrics.to_csv(OUT / "reliability_metrics.csv", index=False)
    by_fold = {str(int(key)): float(value) for key, value in scored.groupby("outer_fold").training_risk_error_spearman.first().items()}
    (OUT / "m5_5_summary.json").write_text(json.dumps({"training_risk_error_spearman_by_fold": by_fold, "full_coverage_mae": full["mae_source"]}, indent=2))
    print(metrics.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
