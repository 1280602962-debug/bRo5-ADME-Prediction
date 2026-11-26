"""Exploratory group ablation. Base arms are reused from the shared-seed fits."""

from __future__ import annotations

import json

import pandas as pd

from m5_0a_seeds import with_error
from m5_lib import OUT, SEEDS, arm_matrix, fit_fold, folds, load_tables, paired_interval, regression
from model_lib import load_store

PARTS = OUT / "group_parts"
NEW_ARMS = ("geometry", "polar", "imhb")


def main() -> None:
    store = load_store()
    tables = load_tables(store)
    cohort = folds(store, tables["ensemble_ids"])
    rows = []
    for arm in NEW_ARMS:
        matrix, names = arm_matrix(store, tables, arm)
        for outer, train, test in cohort:
            for seed in SEEDS:
                pred = fit_fold(store, train, test, matrix, names, seed, PARTS / f"o{outer}_{arm}_s{seed}.csv")
                pred["arm"] = arm
                pred["outer_fold"] = outer
                pred["seed"] = seed
                rows.append(pred)
                print(json.dumps({"arm": arm, "outer": outer, "seed": seed}), flush=True)
    pred = pd.concat(rows, ignore_index=True)
    base = pd.read_csv(OUT / "seed_predictions.csv.gz")
    base = base[base.arm.isin(["R0", "R1", "R3"])]
    all_pred = pd.concat([base, pred], ignore_index=True)
    all_pred.to_csv(OUT / "group_predictions.csv.gz", index=False)
    metric_rows = []
    contrasts = []
    references = {"F04": "R0", "F05": "R1", "all15": "R3"}
    for arm in list(NEW_ARMS) + ["R0", "R1", "R3"]:
        parts = [with_error(all_pred[(all_pred.arm == arm) & (all_pred.seed == seed)]) for seed in SEEDS]
        stacked = pd.concat(parts, ignore_index=True)
        averaged = stacked.groupby("label_id", as_index=False).agg(
            y_pred=("y_pred", "mean"),
            abs_error=("abs_error", "mean"),
            y=("y", "first"),
            Graph_ID=("Graph_ID", "first"),
            Study_key=("Study_key", "first"),
            Scaffold_connectivity_group=("Scaffold_connectivity_group", "first"),
        )
        row = {"arm": arm, "exploratory": arm in NEW_ARMS}
        row.update(regression(averaged.assign(y_pred=averaged.y_pred)))
        metric_rows.append(row)
        if arm in NEW_ARMS:
            reference = with_error(all_pred[(all_pred.arm == "R0") & (all_pred.seed == SEEDS[0])])
            # seed-mean absolute error versus F04
            ref_parts = [with_error(all_pred[(all_pred.arm == "R0") & (all_pred.seed == seed)]) for seed in SEEDS]
            ref_avg = pd.concat(ref_parts, ignore_index=True).groupby("label_id", as_index=False).agg(
                abs_error=("abs_error", "mean"),
                y=("y", "first"),
                Graph_ID=("Graph_ID", "first"),
                Study_key=("Study_key", "first"),
                Scaffold_connectivity_group=("Scaffold_connectivity_group", "first"),
            )
            contrast = {"contrast": f"{arm}-F04", "exploratory": True}
            contrast.update(paired_interval(averaged, ref_avg))
            contrasts.append(contrast)
            del reference
    ensemble = pd.read_csv(REPO_ENSEMBLE())
    groups = pd.read_csv(OUT.parents[1] / "data" / "modeling_v1" / "molecule_groups.csv.gz")
    eligibility = pd.read_csv(OUT.parents[1] / "runs" / "model_development_m4_cpu" / "conformer_eligibility.csv")
    corr = correlation_table(ensemble, groups, eligibility)
    corr.to_csv(OUT / "feature_group_correlations.csv", index=False)
    pd.DataFrame(metric_rows).to_csv(OUT / "feature_group_metrics.csv", index=False)
    pd.DataFrame(contrasts).to_csv(OUT / "feature_group_ablation.csv", index=False)
    print(pd.DataFrame(contrasts).to_string(index=False), flush=True)


def REPO_ENSEMBLE() -> str:
    return str(OUT.parents[1] / "runs" / "model_development_m4_cpu" / "l1_ensemble_features.csv.gz")


def correlation_table(ensemble: pd.DataFrame, groups: pd.DataFrame, eligibility: pd.DataFrame) -> pd.DataFrame:
    from scipy.stats import spearmanr

    frame = ensemble.merge(groups[["Graph_ID", "MW", "HBD", "HBA"]], on="Graph_ID", how="left")
    frame = frame.merge(eligibility[["Graph_ID", "TPSA"]], on="Graph_ID", how="left")
    pairs = [
        ("MW", "ens_Volume_median"),
        ("MW", "ens_Rgyr_median"),
        ("TPSA", "ens_PolarSASA_median"),
        ("HBD", "ens_IMHB_mean"),
        ("HBA", "ens_IMHB_mean"),
    ]
    rows = []
    for left, right in pairs:
        clean = frame[[left, right]].dropna()
        coef = float(spearmanr(clean[left], clean[right]).statistic) if len(clean) > 2 else float("nan")
        rows.append({"x": left, "y": right, "spearman": coef, "n": int(len(clean)), "causal": False})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    main()
