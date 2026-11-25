"""Pool inner-fold predictions and rank fixed G1 model-feature pairs within each outer fold."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "stage_b"))
from score_stage_b import source_table  # noqa: E402

OUT = REPO / "runs" / "model_development_v1"
PARTS = OUT / "parts"


def metrics(frame: pd.DataFrame) -> dict:
    error = frame.y_pred.to_numpy() - frame.y.to_numpy()
    total = float(np.sum((frame.y - frame.y.mean()) ** 2))
    return {
        "mae_source": float(source_table(frame, "abs_error").mean()),
        "mae_label": float(frame.abs_error.mean()),
        "rmse": float(np.sqrt(np.mean(error ** 2))),
        "r2": float(1 - np.sum(error ** 2) / total) if total else np.nan,
        "pearson": float(pearsonr(frame.y, frame.y_pred).statistic),
        "spearman": float(spearmanr(frame.y, frame.y_pred).statistic),
        "signed_source": float(source_table(frame, "signed_error").mean()),
        "n_labels": int(frame.label_id.nunique()),
        "n_studies": int(frame.study_id.nunique()),
    }


def main() -> None:
    files = sorted(PARTS.glob("fixed_G1_*.csv"))
    if not files:
        raise RuntimeError("no fixed G1 parts")
    frame = pd.concat([pd.read_csv(path) for path in files], ignore_index=True)
    frame = frame.rename(columns={"Study_key": "study_id", "Graph_ID": "graph_id"})
    frame["abs_error"] = (frame.y - frame.y_pred).abs()
    frame["signed_error"] = frame.y_pred - frame.y
    rows = []
    for keys, part in frame.groupby(["model", "feature", "outer_fold"]):
        model, feature, outer = keys
        row = {"model": model, "feature": feature, "outer_fold": int(outer), "phase": "fixed", "strategy": "G1"}
        row.update(metrics(part))
        fold_maes = []
        for _, chunk in part.groupby("inner_fold"):
            fold_maes.append(float(source_table(chunk, "abs_error").mean()))
        row["mean_of_inner_fold_mae"] = float(np.mean(fold_maes))
        rows.append(row)
    pooled = pd.DataFrame(rows)
    pooled.to_csv(OUT / "model_screening.csv", index=False)
    screen = pooled[pooled.feature.isin(["F02", "F03", "F05"]) & pooled.model.isin(["RF", "XGBoost", "LightGBM", "SVR"])]
    selected = []
    for outer, part in screen.groupby("outer_fold"):
        best = part.sort_values(["mae_source", "model", "feature"]).head(2)
        for rank, row in enumerate(best.itertuples(index=False), start=1):
            selected.append({"outer_fold": int(outer), "rank": rank, "model": row.model, "feature": row.feature, "mae_source": row.mae_source})
    choice = pd.DataFrame(selected)
    choice.to_csv(OUT / "selected_pairs.csv", index=False)
    print(choice.to_string(index=False))
    print("screen rows", len(pooled))


if __name__ == "__main__":
    main()
