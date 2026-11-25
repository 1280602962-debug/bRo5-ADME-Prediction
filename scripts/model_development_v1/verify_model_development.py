"""Independent checks for the completed bRo5-Perm M0–M3 round."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "stage_b"))
from score_stage_b import source_table  # noqa: E402

OUT = REPO / "runs" / "model_development_v1"
errors = []


def check(condition: bool, message: str) -> None:
    if not condition:
        errors.append(message)


def source_mae(frame: pd.DataFrame) -> float:
    work = frame.rename(columns={"Study_key": "study_id", "Graph_ID": "graph_id"}).copy()
    work["abs_error"] = (work.y - work.y_pred).abs()
    return float(source_table(work, "abs_error").mean())


def main() -> None:
    identity = pd.read_csv(OUT / "g1_alldev_identity.csv")
    check(bool(identity.identical.all()), "G1 graphs differ from ALL_DEV")
    check(int(identity.only_g1.sum()) == 0 and int(identity.only_alldev.sum()) == 0, "G1/ALL_DEV set difference is not empty")
    pred = pd.read_csv(OUT / "oof_predictions.csv.gz")
    c_model = pred[pred.evaluation.eq("C") & pred.model.eq("RF") & pred.feature.eq("F05") & pred.parameter_source.eq("fixed_selected")]
    check(c_model.label_id.nunique() == 530, f"fixed RF-F05 covers {c_model.label_id.nunique()} labels")
    check(c_model.Study_key.nunique() == 153, "fixed RF-F05 does not cover 153 studies")
    metrics = pd.read_csv(OUT / "metrics_pooled.csv")
    published = metrics[
        metrics.evaluation.eq("C")
        & metrics.model.eq("RF")
        & metrics.feature.eq("F05")
        & metrics.parameter_source.eq("fixed_selected")
        & metrics.subset.eq("nonmacro_500_700")
    ]
    check(abs(source_mae(c_model) - float(published.mae_source.iloc[0])) < 1e-8, "recomputed source MAE disagrees with metrics_pooled")
    holdout = pd.read_csv(REPO / "data" / "splits_v1" / "holdout_members.csv")
    hold_graphs = set(holdout.Graph_ID.astype(str))
    check(len(hold_graphs) == 106, f"holdout graph count is {len(hold_graphs)}")
    check(hold_graphs.isdisjoint(set(pred.Graph_ID.astype(str))), "holdout graph appears in predictions")
    dev = pd.read_csv(OUT / "development_pool.csv.gz")
    check(hold_graphs.isdisjoint(set(dev.Graph_ID.astype(str))), "holdout graph appears in the development pool")
    trials = pd.read_csv(OUT / "hyperparameter_trials.csv")
    check(len(trials) == 300, f"Optuna rows are {len(trials)}")
    check(set(trials.strategy) == {"G1"}, "Optuna contains a strategy other than G1")
    ab = pred[pred.evaluation.isin(["A", "B"])]
    check(set(ab.parameter_source) == {"fixed_ab"}, "A/B uses a parameter source other than fixed_ab")
    check(not bool(ab.exploratory.astype(bool).any()), "A/B contains an exploratory prediction")
    inner = pd.read_csv(OUT / "evaluation_splits.csv.gz")
    inner = inner[inner.evaluation.eq("C_inner")]
    for outer in range(5):
        counts = []
        for inner_fold in range(3):
            part = inner[inner.outer_fold.eq(outer) & inner.inner_fold.eq(inner_fold) & inner.role.eq("validation")]
            counts.append(set(part.label_id))
        check(len(counts) == 3 and all(counts), f"outer {outer} is missing an inner validation fold")
    g_parts = list((OUT / "parts").glob("fixed_G*_o0_i0_RF_F05.csv"))
    labels = []
    for path in sorted(g_parts):
        labels.append(set(pd.read_csv(path).label_id))
    check(len(labels) == 3 and labels[0] == labels[1] == labels[2], "G1/G2/G3 inner validation labels differ")
    report = (OUT / "MODEL_DEVELOPMENT_REPORT.md").read_text()
    for phrase in ["log10(Papp [cm/s])", "没有观察到 TL-N0", "0 次试验", "42400"]:
        check(phrase in report, f"report missing {phrase}")
    result = {"status": "PASS" if not errors else "FAIL", "errors": errors}
    (OUT / "verification_result.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
