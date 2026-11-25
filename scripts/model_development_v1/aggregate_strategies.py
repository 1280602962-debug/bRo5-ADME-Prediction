"""Pool G1/G2/G3 inner predictions on the shared validation labels.

Each outer fold is pooled across its three inner folds. Outer folds are not
concatenated, because the same label can appear in more than one outer pool.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "stage_b"))
from score_stage_b import source_table  # noqa: E402

OUT = REPO / "runs" / "model_development_v1"
PARTS = OUT / "parts"


def main() -> None:
    files = sorted(PARTS.glob("fixed_G*_*.csv"))
    frame = pd.concat([pd.read_csv(path) for path in files], ignore_index=True)
    frame = frame.rename(columns={"Study_key": "study_id", "Graph_ID": "graph_id"})
    frame["abs_error"] = (frame.y - frame.y_pred).abs()
    rows = []
    for keys, part in frame.groupby(["strategy", "model", "feature", "outer_fold"]):
        strategy, model, feature, outer = keys
        rows.append(
            {
                "strategy": strategy,
                "model": model,
                "feature": feature,
                "outer_fold": int(outer),
                "mae_source": float(source_table(part, "abs_error").mean()),
                "n_labels": int(part.label_id.nunique()),
                "n_studies": int(part.study_id.nunique()),
            }
        )
    pooled = pd.DataFrame(rows).sort_values(["model", "feature", "strategy", "outer_fold"])
    pooled.to_csv(OUT / "strategy_inner.csv", index=False)
    summary = (
        pooled.groupby(["strategy", "model", "feature"], as_index=False)
        .agg(mae_source_mean_of_outer=("mae_source", "mean"), n_outer=("outer_fold", "nunique"), n_labels_min=("n_labels", "min"), n_labels_max=("n_labels", "max"))
    )
    summary.to_csv(OUT / "strategy_inner_summary.csv", index=False)
    focus = summary[summary.model.eq("RF") & summary.feature.isin(["F04", "F05"])]
    print(focus.sort_values(["feature", "strategy"]).to_string(index=False))


if __name__ == "__main__":
    main()
