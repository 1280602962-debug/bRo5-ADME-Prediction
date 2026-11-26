"""Matched evaluation-C ablation for ionization and conformer features.

P4 minus P0 is the primary scientific contrast. The six-contrast family stays
fixed even when a contrast cannot be estimated. P4 minus P1 is exploratory.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "model_development_v1"))
sys.path.insert(0, str(REPO / "scripts" / "stage_b"))
import model_lib
from model_lib import FIXED, feature_matrix, fit_predict, load_store, model_seed, preprocess, rows_for  # noqa: E402
from score_stage_b import bootstrap_weights, components, graph_mean, source_table  # noqa: E402
from scipy.stats import pearsonr, spearmanr

OUT = REPO / "runs" / "model_development_m4_cpu"
PARTS = OUT / "model_parts"
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
ION_COLUMNS = [
    "molgpka_n_acid",
    "molgpka_n_base",
    "molgpka_has_acid",
    "molgpka_has_base",
    "molgpka_permanent",
    "molgpka_min_acid_pka_x_has",
    "molgpka_max_base_pka_x_has",
]
PRIMARY = [
    ("P2-P0", "P2", "P0"),
    ("P3-P0", "P3", "P0"),
    ("P4-P0", "P4", "P0"),
    ("P4-P3", "P4", "P3"),
    ("P5-P4", "P5", "P4"),
    ("P5-P2", "P5", "P2"),
]


def load_feature_tables(store: dict) -> dict:
    static = pd.read_csv(OUT / "l1_static_features.csv.gz")
    ensemble = pd.read_csv(OUT / "l1_ensemble_features.csv.gz")
    ion = pd.read_csv(OUT / "ionization_features.csv.gz")
    for frame in (static, ensemble, ion):
        frame["Graph_ID"] = frame.Graph_ID.astype(str)
    holdout = set(pd.read_csv(REPO / "data" / "splits_v1" / "holdout_members.csv").Graph_ID.astype(str))
    if holdout & set(static.Graph_ID) or holdout & set(ensemble.Graph_ID) or holdout & set(ion.Graph_ID):
        raise RuntimeError("holdout graph is present in an M4 feature table")
    order = store["graph_id"]

    def matrix(frame: pd.DataFrame, columns: list[str], mask: pd.Series) -> np.ndarray:
        kept = frame.loc[mask, ["Graph_ID", *columns]].drop_duplicates("Graph_ID").set_index("Graph_ID")
        return kept.reindex(order)[columns].to_numpy(np.float64)

    static_ok = static.static_available.astype(str).str.lower().eq("true") if "static_available" in static.columns else pd.Series(True, index=static.index)
    ensemble_ok = ensemble.ensemble_reliable.astype(str).str.lower().eq("true") if "ensemble_reliable" in ensemble.columns else pd.Series(True, index=ensemble.index)
    ion_ok = ion.molgpka_status.eq("ok")
    return {
        "static_ids": set(static.loc[static_ok, "Graph_ID"]),
        "ensemble_ids": set(ensemble.loc[ensemble_ok, "Graph_ID"]),
        "ion_ids": set(ion.loc[ion_ok, "Graph_ID"]),
        "static": matrix(static, STATIC_COLUMNS, static_ok),
        "ensemble": matrix(ensemble, ENSEMBLE_COLUMNS, ensemble_ok),
        "ion": matrix(ion, ION_COLUMNS, ion_ok),
        "sampling_status": ensemble.set_index("Graph_ID")["sampling_status"].to_dict() if "sampling_status" in ensemble else {},
    }


def spec_for(name: str, tables: dict) -> tuple[str, list[np.ndarray], set[str]] | None:
    ion = tables["ion_ids"]
    static_ids = tables["static_ids"]
    ensemble_ids = tables["ensemble_ids"]
    if name == "P0":
        return "F04", [], set(tables["order"])
    if name == "P1":
        return "F05", [], set(tables["order"])
    if name == "P2":
        if not ion:
            return None
        return "F04", [tables["ion"]], ion
    if name == "P3":
        return "F04", [tables["static"]], static_ids
    if name == "P4":
        return "F04", [tables["ensemble"]], ensemble_ids
    if name == "P5":
        both = ion & ensemble_ids
        if not both:
            return None
        return "F04", [tables["ion"], tables["ensemble"]], both
    raise KeyError(name)


def design(base_name: str, extras: list[np.ndarray], names: list[str]) -> tuple[np.ndarray, list[str]]:
    pieces = [feature_matrix_cache[base_name][0]]
    labels = list(feature_matrix_cache[base_name][1])
    for extra in extras:
        pieces.append(extra)
    labels.extend(names)
    return np.hstack(pieces), labels


def fit_frame(store, train, test, base_name, extras, extra_names, tag, model, params) -> pd.DataFrame:
    path = PARTS / f"{tag}.csv"
    if path.exists():
        return pd.read_csv(path)
    matrix, names = feature_matrix(store, base_name)
    if extras:
        matrix = np.hstack([matrix, *extras])
        names = names + extra_names
    train_x, y, weights = rows_for(train, store, matrix)
    test_x, _, _ = rows_for(test, store, matrix)
    if not np.isfinite(train_x).all() or not np.isfinite(test_x).all():
        missing = ~np.isfinite(train_x).all(axis=0)
        if missing.any() and model == "RF":
            pass
    train_ready, test_ready, _ = preprocess(train_x, test_x, names, False)
    seed = model_seed("m4", tag)
    predicted = fit_predict(model, train_ready, y, weights, [test_ready], params, seed)[0]
    output = test[["label_id", "Graph_ID", "Study_key", "Scaffold_connectivity_group", "y", "cohort", "MW"]].copy()
    output["y_pred"] = predicted
    output["n_train_labels"] = int(len(train))
    output["n_test_labels"] = int(len(test))
    output.to_csv(path, index=False)
    print(json.dumps({"fit": tag, "train": int(len(train)), "test": int(len(test))}), flush=True)
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
        "pearson": float(pearsonr(work.y, work.y_pred).statistic) if len(work) > 2 else np.nan,
        "spearman": float(spearmanr(work.y, work.y_pred).statistic) if len(work) > 2 else np.nan,
        "signed_source": float(source_table(work, "signed_error").mean()),
        "n_labels": int(work.label_id.nunique()),
        "n_graphs": int(work.graph_id.nunique()),
        "n_studies": int(work.study_id.nunique()),
    }


def subset_labels(frame: pd.DataFrame, graphs: set[str]) -> pd.DataFrame:
    return frame[frame.Graph_ID.astype(str).isin(graphs)].copy()


def contrast_graphs(left: str, right: str, tables: dict) -> set[str] | None:
    specs = []
    for name in (left, right):
        spec = spec_for(name, tables)
        if spec is None:
            return None
        specs.append(spec[2])
    graphs = set.intersection(*specs)
    if not graphs:
        return None
    return graphs


def paired_interval(left: pd.DataFrame, right: pd.DataFrame) -> dict:
    merged = left[["label_id", "y", "y_pred", "Study_key", "Graph_ID", "Scaffold_connectivity_group"]].merge(
        right[["label_id", "y_pred"]], on="label_id", suffixes=("_left", "_right")
    )
    if len(merged) != len(left) or len(merged) != len(right):
        raise RuntimeError("paired predictions do not share test labels")
    merged["abs_left"] = (merged.y_pred_left - merged.y).abs()
    merged["abs_right"] = (merged.y_pred_right - merged.y).abs()
    named = merged.rename(columns={"Graph_ID": "graph_id", "Study_key": "study_id"})
    left_source = source_table(named.assign(abs_error=named.abs_left), "abs_error")
    right_source = source_table(named.assign(abs_error=named.abs_right), "abs_error")
    studies = left_source.index.intersection(right_source.index)
    point = float(left_source.loc[studies].mean() - right_source.loc[studies].mean())
    study_component = components(merged)
    rng = np.random.default_rng(271828)
    study_list, weights = bootstrap_weights(study_component, rng)
    vector = (left_source.loc[study_list] - right_source.loc[study_list]).to_numpy(np.float64)
    draw = (weights * vector).sum(1) / weights.sum(1)
    low95, high95 = np.quantile(draw, [0.025, 0.975])
    low_b, high_b = np.quantile(draw, [0.05 / 12, 1 - 0.05 / 12])
    return {
        "estimate": point,
        "ci95_low": float(low95),
        "ci95_high": float(high95),
        "ci95_excludes_zero": bool(low95 > 0 or high95 < 0),
        "ci_bonferroni_low": float(low_b),
        "ci_bonferroni_high": float(high_b),
        "ci_bonferroni_excludes_zero": bool(low_b > 0 or high_b < 0),
        "n_labels": int(len(merged)),
        "n_studies": int(len(studies)),
        "n_components": int(study_component.max()) + 1,
        "family_size": 6,
        "bonferroni_nominal": 1 - 0.05 / 6,
    }


def extra_names(columns_groups: list[list[str]]) -> list[str]:
    names = []
    for group in columns_groups:
        names.extend(group)
    return names


def arm_inputs(name: str, tables: dict) -> tuple[str, list[np.ndarray], list[str]]:
    spec = spec_for(name, tables)
    if spec is None:
        raise KeyError(name)
    base, extras, _graphs = spec
    groups = []
    if name in {"P2", "P5"}:
        groups.append(ION_COLUMNS)
    if name in {"P3"}:
        groups.append(STATIC_COLUMNS)
    if name in {"P4", "P5"}:
        groups.append(ENSEMBLE_COLUMNS)
    return base, extras, extra_names(groups)


def fit_arm(store, train, test, name, tables, tag, model="RF") -> pd.DataFrame:
    base, extras, names = arm_inputs(name, tables)
    params = FIXED[model]
    return fit_frame(store, train, test, base, extras, names, tag, model, params)


def inner_gate(store, splits, tables) -> pd.DataFrame:
    rows = []
    graphs = tables["ensemble_ids"]
    for outer in range(5):
        chunks = {"P0": [], "P4": []}
        for inner in range(3):
            part = splits[splits.evaluation.eq("C_inner") & splits.outer_fold.eq(outer) & splits.inner_fold.eq(inner)]
            train = subset_labels(part[part.role.eq("train_candidate")], graphs)
            valid = subset_labels(part[part.role.eq("validation")], graphs)
            if train.empty or valid.empty:
                rows.append({"outer_fold": outer, "inner_fold": inner, "status": "empty"})
                continue
            for name in ("P0", "P4"):
                pred = fit_arm(store, train, valid, name, tables, f"inner_o{outer}_i{inner}_{name}")
                chunks[name].append(pred)
        if not chunks["P0"] or not chunks["P4"]:
            continue
        left = pd.concat(chunks["P4"], ignore_index=True)
        right = pd.concat(chunks["P0"], ignore_index=True)
        stats_left = regression(left)
        stats_right = regression(right)
        rows.append(
            {
                "outer_fold": outer,
                "inner_fold": "pooled",
                "status": "ok",
                "mae_p4": stats_left["mae_source"],
                "mae_p0": stats_right["mae_source"],
                "delta_p4_minus_p0": stats_left["mae_source"] - stats_right["mae_source"],
                "improved": stats_left["mae_source"] < stats_right["mae_source"],
            }
        )
    return pd.DataFrame(rows)


def permutation_note(store, train, test, tables) -> pd.DataFrame:
    base, extras, names = arm_inputs("P4", tables)
    matrix, base_names = feature_matrix(store, base)
    matrix = np.hstack([matrix, *extras])
    all_names = base_names + names
    train_x, y, weights = rows_for(train, store, matrix)
    test_x, y_test, _ = rows_for(test, store, matrix)
    train_ready, test_ready, kept = preprocess(train_x, test_x, all_names, False)
    seed = model_seed("m4", "C_o0_P4_ensemble")
    model = RandomForestRegressor(
        n_estimators=500,
        max_depth=None,
        min_samples_leaf=2,
        max_features=0.5,
        bootstrap=True,
        random_state=seed,
        n_jobs=4,
    )
    model.fit(train_ready, y, sample_weight=weights)
    baseline = float(np.mean(np.abs(model.predict(test_ready) - y_test)))
    focus = [name for name in kept if name in set(ENSEMBLE_COLUMNS + ["MW", "TPSA", "MolLogP", "NumRotatableBonds"])]
    rng = np.random.default_rng(271828)
    rows = []
    for name in focus:
        column = kept.index(name)
        increases = []
        for _ in range(5):
            shuffled = test_ready.copy()
            shuffled[:, column] = rng.permutation(shuffled[:, column])
            increases.append(float(np.mean(np.abs(model.predict(shuffled) - y_test)) - baseline))
        rows.append({"feature": name, "label_mae_increase_mean": float(np.mean(increases)), "baseline_label_mae": baseline, "n_repeats": 5})
    return pd.DataFrame(rows)


def main() -> None:
    import os
    model_lib.N_JOBS = int(os.environ.get("M4_MODEL_JOBS", "4"))
    PARTS.mkdir(parents=True, exist_ok=True)
    store = load_store()
    tables = load_feature_tables(store)
    tables["order"] = set(store["graph_id"])
    dev = store["dev"]
    dev["Graph_ID"] = dev.Graph_ID.astype(str)
    g1 = dev[dev.g1].copy()
    g1["assigned_fold"] = g1.assigned_fold.astype(int)
    test = pd.read_csv(REPO / "data" / "splits_v1" / "test_labels.csv.gz")
    test["Graph_ID"] = test.Graph_ID.astype(str)
    splits = pd.read_csv(REPO / "runs" / "model_development_v1" / "evaluation_splits.csv.gz")
    gate = inner_gate(store, splits, tables)
    gate.to_csv(OUT / "m4_inner_gate.csv", index=False)
    pooled = gate[gate.inner_fold.astype(str).eq("pooled")]
    use_xgb = int(pooled.improved.fillna(False).sum()) >= 4 if len(pooled) else False
    predictions = []
    cohort_rows = []
    coverage_rows = []
    arms = {
        "P0": None,
        "P1": None,
        "P2": "ion",
        "P3": "static",
        "P4": "ensemble",
        "P5": "both",
    }
    fitted = {}
    for outer in range(5):
        train_all = g1[g1.assigned_fold.ne(outer)].copy()
        test_all = test[test.assigned_fold.eq(outer)].copy()
        universes = {
            "all": set(store["graph_id"]),
            "ion": tables["ion_ids"],
            "static": tables["static_ids"],
            "ensemble": tables["ensemble_ids"],
            "both": tables["ion_ids"] & tables["ensemble_ids"],
        }
        needed = {
            "P0_all": ("P0", "all"),
            "P0_ion": ("P0", "ion"),
            "P0_static": ("P0", "static"),
            "P0_ensemble": ("P0", "ensemble"),
            "P0_both": ("P0", "both"),
            "P1_ensemble": ("P1", "ensemble"),
            "P2_ion": ("P2", "ion"),
            "P2_both": ("P2", "both"),
            "P3_static": ("P3", "static"),
            "P3_ensemble": ("P3", "ensemble"),
            "P4_ensemble": ("P4", "ensemble"),
            "P4_both": ("P4", "both"),
            "P5_both": ("P5", "both"),
        }
        local = {}
        for key, (arm, universe_name) in needed.items():
            graphs = universes[universe_name]
            if arm in {"P2", "P5"} and not tables["ion_ids"]:
                continue
            if universe_name == "both" and not universes["both"]:
                continue
            train = subset_labels(train_all, graphs)
            evaluation = subset_labels(test_all, graphs)
            if train.empty or evaluation.empty:
                continue
            pred = fit_arm(store, train, evaluation, arm, tables, f"C_o{outer}_{key}")
            pred["evaluation"] = "C"
            pred["outer_fold"] = outer
            pred["arm"] = key
            pred["model"] = "RF"
            pred["parameter_source"] = "fixed_rf"
            local[key] = pred
            predictions.append(pred)
        if "P4_ensemble" in local and "P0_all" in local:
            enhanced = local["P4_ensemble"].set_index("label_id")
            fallback = local["P0_all"].set_index("label_id")
            stitched = fallback.copy()
            overlap = stitched.index.intersection(enhanced.index)
            stitched.loc[overlap, "y_pred"] = enhanced.loc[overlap, "y_pred"]
            stitched = stitched.reset_index()
            stitched["arm"] = "P4_fallback_full"
            stitched["model"] = "RF"
            stitched["parameter_source"] = "fixed_rf_fallback"
            stitched["used_enhanced"] = stitched.label_id.isin(overlap)
            predictions.append(stitched)
            coverage_rows.append(
                {
                    "outer_fold": outer,
                    "test_labels": int(len(test_all)),
                    "ensemble_test_labels": int(len(local["P4_ensemble"])),
                    "static_test_labels": int(len(local.get("P3_static", []))),
                    "ion_test_labels": int(len(local["P2_ion"])) if "P2_ion" in local else 0,
                    "fallback_uses_enhanced": int(stitched.used_enhanced.sum()),
                }
            )
        fitted[outer] = local
        for key, pred in local.items():
            cohort_rows.append({"evaluation": "C", "outer_fold": outer, "arm": key, **{k: pred.iloc[0][k] for k in ("n_train_labels", "n_test_labels")}})
    if use_xgb:
        for outer in range(5):
            train_all = g1[g1.assigned_fold.ne(outer)].copy()
            test_all = test[test.assigned_fold.eq(outer)].copy()
            graphs = tables["ensemble_ids"]
            train = subset_labels(train_all, graphs)
            evaluation = subset_labels(test_all, graphs)
            for arm in ("P0", "P4"):
                pred = fit_arm(store, train, evaluation, arm, tables, f"C_o{outer}_xgb_{arm}", model="XGBoost")
                pred["evaluation"] = "C"
                pred["outer_fold"] = outer
                pred["arm"] = f"XGB_{arm}_ensemble"
                pred["model"] = "XGBoost"
                pred["parameter_source"] = "fixed_xgb_inner_gate"
                predictions.append(pred)
    for evaluation_name in ("A", "B"):
        part = splits[splits.evaluation.eq(evaluation_name)].copy()
        for fold in range(5):
            train = subset_labels(part[part.outer_fold.ne(fold)], tables["ensemble_ids"])
            evaluation = subset_labels(part[part.outer_fold.eq(fold)], tables["ensemble_ids"])
            if train.empty or evaluation.empty:
                continue
            for arm in ("P0", "P4"):
                pred = fit_arm(store, train, evaluation, arm, tables, f"{evaluation_name}_o{fold}_{arm}")
                pred["evaluation"] = evaluation_name
                pred["outer_fold"] = fold
                pred["arm"] = arm
                pred["model"] = "RF"
                pred["parameter_source"] = "fixed_ab"
                pred["exploratory"] = False
                predictions.append(pred)
    pred = pd.concat(predictions, ignore_index=True)
    pred.to_csv(OUT / "m4_oof_predictions.csv.gz", index=False)
    holdout = set(pd.read_csv(REPO / "data" / "splits_v1" / "holdout_members.csv").Graph_ID.astype(str))
    if holdout & set(pred.Graph_ID.astype(str)):
        raise RuntimeError("holdout graph entered M4 predictions")
    metric_rows = []
    for keys, part in pred.groupby(["evaluation", "model", "arm", "parameter_source"], dropna=False):
        evaluation, model, arm, source = keys
        row = {"evaluation": evaluation, "model": model, "arm": arm, "parameter_source": source}
        row.update(regression(part))
        metric_rows.append(row)
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(OUT / "m4_metrics.csv", index=False)
    pd.DataFrame(coverage_rows).to_csv(OUT / "m4_coverage_metrics.csv", index=False)
    pd.DataFrame(cohort_rows).to_csv(OUT / "matched_cohort_summary.csv", index=False)
    c = pred[pred.evaluation.eq("C") & pred.model.eq("RF") & pred.parameter_source.eq("fixed_rf")]
    contrast_map = {
        "P2-P0": ("P2_ion", "P0_ion"),
        "P3-P0": ("P3_static", "P0_static"),
        "P4-P0": ("P4_ensemble", "P0_ensemble"),
        "P4-P3": ("P4_ensemble", "P3_ensemble"),
        "P5-P4": ("P5_both", "P4_both"),
        "P5-P2": ("P5_both", "P2_both"),
    }
    delta_rows = []
    for name, _left, _right in PRIMARY:
        left_name, right_name = contrast_map[name]
        left = c[c.arm.eq(left_name)]
        right = c[c.arm.eq(right_name)]
        row = {"contrast": name, "role": "primary", "candidate": left_name, "reference": right_name}
        if left.empty or right.empty or set(left.label_id) != set(right.label_id):
            row.update({"status": "not_executable", "family_size": 6, "bonferroni_nominal": 1 - 0.05 / 6})
        else:
            row.update({"status": "estimated"})
            row.update(paired_interval(left, right))
        delta_rows.append(row)
    left = c[c.arm.eq("P4_ensemble")]
    right = c[c.arm.eq("P1_ensemble")]
    exploratory = {"contrast": "P4-P1", "role": "exploratory", "candidate": "P4_ensemble", "reference": "P1_ensemble", "family_size": 6}
    if left.empty or right.empty or set(left.label_id) != set(right.label_id):
        exploratory["status"] = "not_executable"
    else:
        interval = paired_interval(left, right)
        interval.pop("ci_bonferroni_low")
        interval.pop("ci_bonferroni_high")
        interval.pop("ci_bonferroni_excludes_zero")
        exploratory.update({"status": "estimated", **interval})
        exploratory["bonferroni_applied"] = False
    delta_rows.append(exploratory)
    deltas = pd.DataFrame(delta_rows)
    deltas.to_csv(OUT / "m4_paired_deltas.csv", index=False)
    train0 = subset_labels(g1[g1.assigned_fold.ne(0)], tables["ensemble_ids"])
    test0 = subset_labels(test[test.assigned_fold.eq(0)], tables["ensemble_ids"])
    importance = permutation_note(store, train0, test0, tables)
    importance.to_csv(OUT / "m4_permutation_importance.csv", index=False)
    summary = {
        "primary_contrast": "P4-P0",
        "family_size": 6,
        "xgb_added": bool(use_xgb),
        "inner_improved_folds": int(pooled.improved.fillna(False).sum()) if len(pooled) else 0,
        "ensemble_graphs": len(tables["ensemble_ids"]),
        "static_graphs": len(tables["static_ids"]),
        "ion_graphs": len(tables["ion_ids"]),
        "ab_parameter_source": "fixed_ab",
        "ab_inherits_evaluation_c_selection": False,
    }
    (OUT / "m4_3_summary.json").write_text(json.dumps(summary, indent=2))
    print(deltas.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
