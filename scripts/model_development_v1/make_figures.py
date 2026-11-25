"""Figures for optimization, evaluation C, A/B/C, MW regions, and G1/G2/G3."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "stage_b"))
from score_stage_b import source_table  # noqa: E402

OUT = REPO / "runs" / "model_development_v1"
FIG = OUT / "figures"


def source_mae(frame: pd.DataFrame) -> float:
    work = frame.rename(columns={"Study_key": "study_id", "Graph_ID": "graph_id"})
    work["abs_error"] = (work.y - work.y_pred).abs()
    return float(source_table(work, "abs_error").mean())


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    screen = pd.read_csv(OUT / "model_screening.csv")
    trials = pd.read_csv(OUT / "hyperparameter_trials.csv")
    best = trials[trials.best.astype(str).str.lower().eq("true")].copy()
    metrics = pd.read_csv(OUT / "metrics_pooled.csv")
    pred = pd.read_csv(OUT / "oof_predictions.csv.gz")
    inner = pd.read_csv(OUT / "strategy_inner_summary.csv")

    labels = []
    fixed_vals = []
    opt_vals = []
    for row in best.sort_values(["outer_fold", "model", "feature"]).itertuples(index=False):
        match = screen[
            screen.strategy.eq("G1")
            & screen.outer_fold.eq(row.outer_fold)
            & screen.model.eq(row.model)
            & screen.feature.eq(row.feature)
        ]
        labels.append(f"o{int(row.outer_fold)}\n{row.model}-{row.feature}")
        fixed_vals.append(float(match.mae_source.iloc[0]))
        opt_vals.append(float(row.mae_source))
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(10, 4.2))
    ax.bar(x - 0.18, fixed_vals, width=0.36, label="fixed")
    ax.bar(x + 0.18, opt_vals, width=0.36, label="Optuna best")
    ax.set_xticks(x, labels, fontsize=7)
    ax.set_ylabel("Inner source-equal MAE")
    ax.set_title("G1 inner validation before and after Optuna")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "optimization_inner.png", dpi=160)
    plt.close(fig)

    outer = metrics[
        metrics.evaluation.eq("C")
        & metrics.subset.eq("nonmacro_500_700")
        & metrics.model.eq("RF")
        & metrics.feature.eq("F05")
        & metrics.parameter_source.isin(["fixed_selected", "optuna"])
    ]
    fig, ax = plt.subplots(figsize=(4.2, 3.6))
    order = ["fixed_selected", "optuna"]
    values = [float(outer.loc[outer.parameter_source.eq(name), "mae_source"].iloc[0]) for name in order]
    ax.bar(["fixed RF-F05", "Optuna RF-F05"], values, color=["#4C78A8", "#F58518"])
    ax.set_ylabel("Outer source-equal MAE")
    ax.set_title("Evaluation C, 530 labels")
    ax.set_ylim(0.5, 0.65)
    fig.tight_layout()
    fig.savefig(FIG / "optimization_outer.png", dpi=160)
    plt.close(fig)

    scatter = pred[
        pred.evaluation.eq("C") & pred.model.eq("RF") & pred.feature.eq("F05") & pred.parameter_source.eq("fixed_selected")
    ]
    fig, ax = plt.subplots(figsize=(4.4, 4.2))
    ax.scatter(scatter.y, scatter.y_pred, s=12, alpha=0.55, linewidths=0)
    lo = min(scatter.y.min(), scatter.y_pred.min())
    hi = max(scatter.y.max(), scatter.y_pred.max())
    ax.plot([lo, hi], [lo, hi], color="0.3", linewidth=1)
    ax.set_xlabel("Experimental log10(Papp [cm/s])")
    ax.set_ylabel("Predicted log10(Papp [cm/s])")
    ax.set_title("Evaluation C, fixed RF-F05")
    fig.tight_layout()
    fig.savefig(FIG / "pred_vs_obs_c.png", dpi=160)
    plt.close(fig)

    target = metrics[metrics.subset.eq("nonmacro_500_700")].copy()
    keep = target.parameter_source.isin(["fixed_ab", "fixed_selected"]) & target.model.isin(["RF", "XGBoost", "LightGBM"])
    keep &= ~((target.evaluation.eq("C")) & ~((target.model.eq("RF")) & (target.feature.eq("F05"))))
    shown = target[keep].copy()
    shown["label"] = shown.evaluation + " " + shown.model
    shown = shown.sort_values(["evaluation", "mae_source"])
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.8))
    axes[0].barh(shown.label, shown.mae_source, color="#4C78A8")
    axes[0].set_xlabel("Source-equal MAE")
    axes[1].barh(shown.label, shown.r2, color="#54A24B")
    axes[1].set_xlabel("R2")
    fig.suptitle("500–700 Da nonmacro subset")
    fig.tight_layout()
    fig.savefig(FIG / "evaluation_abc.png", dpi=160)
    plt.close(fig)

    regional = pred[pred.evaluation.eq("A") & pred.parameter_source.eq("fixed_ab") & pred.model.eq("RF") & pred.feature.eq("F05")]
    region_rows = []
    for cohort, part in regional.groupby("cohort"):
        region_rows.append({"cohort": cohort, "mae_source": source_mae(part), "n_labels": int(len(part))})
    c_part = pred[pred.evaluation.eq("C") & pred.model.eq("RF") & pred.feature.eq("F05") & pred.parameter_source.eq("fixed_selected")]
    region_rows.append({"cohort": "C nonmacro_500_700", "mae_source": source_mae(c_part), "n_labels": int(len(c_part))})
    regions = pd.DataFrame(region_rows)
    regions.to_csv(OUT / "error_by_region.csv", index=False)
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.bar(regions.cohort, regions.mae_source, color="#4C78A8")
    ax.set_ylabel("Source-equal MAE")
    ax.set_title("RF-F05 error by cohort")
    ax.tick_params(axis="x", labelrotation=20)
    fig.tight_layout()
    fig.savefig(FIG / "error_by_mw_region.png", dpi=160)
    plt.close(fig)

    focus = inner[inner.model.eq("RF") & inner.feature.isin(["F04", "F05"])].copy()
    outer_focus = metrics[
        metrics.evaluation.eq("C")
        & metrics.subset.eq("nonmacro_500_700")
        & metrics.model.eq("RF")
        & metrics.feature.isin(["F04", "F05"])
        & metrics.parameter_source.isin(["fixed_g1_rf", "fixed_selected", "fixed_G2", "fixed_G3"])
    ].copy()
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6))
    for ax, feature in zip(axes, ["F04", "F05"]):
        part = focus[focus.feature.eq(feature)].set_index("strategy").reindex(["G1", "G2", "G3"])
        ax.bar(part.index, part.mae_source_mean_of_outer, color="#4C78A8")
        ax.set_title(f"Inner RF-{feature}")
        ax.set_ylabel("Mean of five outer-fold MAEs")
    fig.tight_layout()
    fig.savefig(FIG / "strategy_inner.png", dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    names = []
    vals = []
    for feature, source, name in [
        ("F04", "fixed_g1_rf", "G1 F04"),
        ("F04", "fixed_G2", "G2 F04"),
        ("F04", "fixed_G3", "G3 F04"),
        ("F05", "fixed_selected", "G1 F05"),
        ("F05", "fixed_G2", "G2 F05"),
        ("F05", "fixed_G3", "G3 F05"),
    ]:
        row = outer_focus[outer_focus.feature.eq(feature) & outer_focus.parameter_source.eq(source)]
        names.append(name)
        vals.append(float(row.mae_source.iloc[0]))
    ax.bar(names, vals, color="#F58518")
    ax.set_ylabel("Outer source-equal MAE")
    ax.set_title("Evaluation C, fixed parameters")
    ax.tick_params(axis="x", labelrotation=20)
    fig.tight_layout()
    fig.savefig(FIG / "strategy_outer.png", dpi=160)
    plt.close(fig)
    print("figures", len(list(FIG.glob('*.png'))))


if __name__ == "__main__":
    main()
