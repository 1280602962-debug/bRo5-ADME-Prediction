"""Score final predictions and compare evaluation C with the frozen Stage B RF baselines."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "stage_b"))
from score_stage_b import bootstrap_weights, components, graph_mean, source_table  # noqa: E402

OUT = REPO / "runs" / "model_development_v1"
FIG = OUT / "figures"
SEEDS = [17, 29, 42, 71, 101, 137, 271, 314, 577, 911]


def regression(frame: pd.DataFrame) -> dict:
    error = frame.y_pred.to_numpy() - frame.y.to_numpy()
    total = float(np.sum((frame.y - frame.y.mean()) ** 2))
    return {
        "mae_source": float(source_table(frame, "abs_error").mean()),
        "mae_graph": graph_mean(frame, "abs_error"),
        "mae_label": float(frame.abs_error.mean()),
        "rmse": float(np.sqrt(np.mean(error ** 2))),
        "r2": float(1 - np.sum(error ** 2) / total) if total else np.nan,
        "pearson": float(pearsonr(frame.y, frame.y_pred).statistic),
        "spearman": float(spearmanr(frame.y, frame.y_pred).statistic),
        "signed_source": float(source_table(frame, "signed_error").mean()),
        "n_labels": int(frame.label_id.nunique()),
        "n_graphs": int(frame.graph_id.nunique()),
        "n_studies": int(frame.study_id.nunique()),
    }


def stage_b() -> dict[str, pd.DataFrame]:
    parts = REPO / "runs" / "stage_b" / "pred_parts"
    frames = {"B": [], "BM": [], "ALL_DEV": []}
    for fold in range(5):
        for seed in SEEDS:
            for arm in ["B", "BM"]:
                path = parts / f"B01_f{fold}_s{seed}_RF_{arm}.csv"
                frames[arm].append(pd.read_csv(path))
        frames["ALL_DEV"].append(pd.read_csv(parts / f"B01_f{fold}_s-1_RF_ALL_DEV.csv"))
    out = {}
    for arm, chunks in frames.items():
        pred = pd.concat(chunks, ignore_index=True)
        pred["abs_error"] = (pred.y - pred.y_pred).abs()
        if arm == "ALL_DEV":
            out[arm] = pred
        else:
            out[arm] = pred.groupby("label_id", as_index=False).agg(
                abs_error=("abs_error", "mean"), y=("y", "first"), study_id=("study_id", "first"), graph_id=("graph_id", "first")
            )
    return out


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    pred = pd.read_csv(OUT / "oof_predictions.csv.gz")
    pred = pred.rename(columns={"Graph_ID": "graph_id", "Study_key": "study_id"})
    pred["abs_error"] = (pred.y - pred.y_pred).abs()
    pred["signed_error"] = pred.y_pred - pred.y
    rows = []
    fold_rows = []
    for keys, part in pred.groupby(["evaluation", "model", "feature", "parameter_source"]):
        evaluation, model, feature, source = keys
        for subset, chunk in [("all", part), ("nonmacro_500_700", part[part.cohort.eq("nonmacro_500_700")])]:
            if chunk.empty:
                continue
            row = {"evaluation": evaluation, "model": model, "feature": feature, "parameter_source": source, "subset": subset}
            row.update(regression(chunk))
            row["mw_median"] = float(chunk.MW.median())
            rows.append(row)
        for fold, chunk in part.groupby("outer_fold"):
            fold_rows.append(
                {
                    "evaluation": evaluation,
                    "model": model,
                    "feature": feature,
                    "parameter_source": source,
                    "outer_fold": int(fold),
                    "mae_source": float(source_table(chunk, "abs_error").mean()),
                    "n_labels": int(len(chunk)),
                }
            )
    metrics = pd.DataFrame(rows)
    metrics.to_csv(OUT / "metrics_pooled.csv", index=False)
    pd.DataFrame(fold_rows).to_csv(OUT / "metrics_by_fold.csv", index=False)
    test = pd.read_csv(REPO / "data" / "splits_v1" / "test_labels.csv.gz")
    study_component = components(test)
    studies, weights = bootstrap_weights(study_component, np.random.default_rng(271828))
    baselines = stage_b()
    base_draw = {}
    base_point = {}
    for arm, frame in baselines.items():
        vector = source_table(frame.assign(abs_error=frame.abs_error), "abs_error").loc[studies].to_numpy(np.float64)
        base_draw[arm] = (weights * vector).sum(1) / weights.sum(1)
        base_point[arm] = float(vector.mean())
    c_pred = pred[pred.evaluation.eq("C")]
    delta_rows = []
    for keys, part in c_pred.groupby(["model", "feature", "parameter_source"]):
        model, feature, source = keys
        if part.label_id.nunique() != 530 or part.study_id.nunique() != len(studies):
            continue
        vector = source_table(part, "abs_error").loc[studies].to_numpy(np.float64)
        draw = (weights * vector).sum(1) / weights.sum(1)
        point = float(vector.mean())
        references = {"stageB_RF_B": "B", "stageB_RF_BM": "BM", "stageB_RF_ALL_DEV": "ALL_DEV"}
        if not (model == "RF" and feature == "F04" and source == "fixed_g1_rf"):
            references["g1_fixed_RF_F04"] = None
        for name, arm in references.items():
            if arm is None:
                fixed = c_pred[c_pred.model.eq("RF") & c_pred.feature.eq("F04") & c_pred.parameter_source.eq("fixed_g1_rf")]
                other_vector = source_table(fixed, "abs_error").loc[studies].to_numpy(np.float64)
                other_draw = (weights * other_vector).sum(1) / weights.sum(1)
                other_point = float(other_vector.mean())
            else:
                other_draw = base_draw[arm]
                other_point = base_point[arm]
            diff = draw - other_draw
            low, high = np.quantile(diff, [0.025, 0.975])
            delta_rows.append(
                {
                    "model": model,
                    "feature": feature,
                    "parameter_source": source,
                    "reference": name,
                    "estimate": point - other_point,
                    "ci_low": float(low),
                    "ci_high": float(high),
                    "interval_excludes_zero": bool(low > 0 or high < 0),
                    "n_studies": int(len(studies)),
                    "n_components": int(study_component.max()) + 1,
                }
            )
    deltas = pd.DataFrame(delta_rows)
    deltas.to_csv(OUT / "paired_deltas.csv", index=False)
    screen = pd.read_csv(OUT / "model_screening.csv")
    heat = screen[screen.feature.isin(["F01", "F02", "F03", "F04", "F05"])].groupby(["model", "feature"]).mae_source.mean().unstack()
    fig, ax = plt.subplots(figsize=(7, 3.8))
    image = ax.imshow(heat.to_numpy(), cmap="viridis_r", aspect="auto")
    ax.set_xticks(range(heat.shape[1]), heat.columns)
    ax.set_yticks(range(heat.shape[0]), heat.index)
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            value = heat.to_numpy()[i, j]
            if np.isfinite(value):
                ax.text(j, i, f"{value:.3f}", ha="center", va="center", color="white", fontsize=8)
    ax.set_title("Inner source-equal MAE")
    fig.colorbar(image, ax=ax, fraction=0.046)
    fig.tight_layout()
    fig.savefig(FIG / "feature_model_heatmap.png", dpi=160)
    plt.close(fig)
    print(metrics[metrics.subset.eq("nonmacro_500_700")].sort_values("mae_source").head(20).to_string(index=False))
    print(deltas.head(12).to_string(index=False))


if __name__ == "__main__":
    main()
