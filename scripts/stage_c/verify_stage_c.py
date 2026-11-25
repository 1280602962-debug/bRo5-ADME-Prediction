"""Independent check of the archived C02 predictions and the recomputed archive."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "stage_b"))
from score_stage_b import source_table  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
RUN = REPO / "runs" / "stage_c"
PRED = RUN / "c02" / "pred"
ARMS = ["TL-B", "TL-N0", "TL-JS", "TL-JM", "TL-JN", "TL-FS", "TL-FM", "TL-FN"]


def main() -> None:
    errors = []
    files = sorted(PRED.glob("*.csv"))
    if len(files) != 400:
        errors.append(f"prediction files {len(files)}")
    frames = [pd.read_csv(path) for path in files]
    pred = pd.concat(frames, ignore_index=True)
    pred["abs_error"] = (pred.y - pred.y_pred).abs()
    if len(pred) != 42400:
        errors.append(f"rows {len(pred)}")
    test = pd.read_csv(REPO / "data" / "splits_v1" / "test_labels.csv.gz")
    expected = set(test.label_id)
    if set(test.corrected_task) != {"Caco-2:AB"} or set(test.analysis_set) != {"caco2_ab_main"}:
        errors.append("endpoint is not Caco-2 A to B")
    holdout = set(pd.read_csv(REPO / "data" / "splits_v1" / "holdout_members.csv").Graph_ID)
    if set(pred.graph_id) & holdout:
        errors.append("holdout graph in predictions")
    if sorted(pred.arm.unique()) != sorted(ARMS):
        errors.append("arm set mismatch")
    for arm, part in pred.groupby("arm"):
        if len(part) != 5300:
            errors.append(f"{arm} rows {len(part)}")
        for _, chunk in part.groupby("seed"):
            if set(chunk.label_id) != expected or chunk.label_id.duplicated().any():
                errors.append(f"{arm} seed coverage")
                break
    merged = pred.drop_duplicates("label_id")[["label_id", "y", "study_id", "graph_id"]].merge(
        test[["label_id", "y", "Study_key", "Graph_ID"]], on="label_id"
    )
    if len(merged) != 530 or not np.allclose(merged.y_x, merged.y_y) or not merged.study_id.eq(merged.Study_key).all():
        errors.append("key or y mismatch")
    recomputed = []
    for arm, part in pred.groupby("arm"):
        seed_mean = part.groupby("label_id", as_index=False).agg(
            abs_error=("abs_error", "mean"), study_id=("study_id", "first"), graph_id=("graph_id", "first")
        )
        recomputed.append({"arm": arm, "mae_source": float(source_table(seed_mean, "abs_error").mean())})
    fresh = pd.DataFrame(recomputed)
    saved = pd.read_csv(RUN / "c02" / "metrics_recomputed.csv")
    published = pd.read_csv(RUN / "c02" / "metrics_pooled.csv")
    check = fresh.merge(saved[["arm", "mae_source"]], on="arm", suffixes=("_fresh", "_saved")).merge(
        published[["arm", "mae_source"]].rename(columns={"mae_source": "mae_published"}), on="arm"
    )
    if not np.allclose(check.mae_source_fresh, check.mae_source_saved, atol=1e-8):
        errors.append("recomputed file does not match an independent aggregation")
    if not np.allclose(check.mae_source_fresh, check.mae_published, atol=1e-8):
        errors.append("published source MAE differs from the independent aggregation")
    records = [json.loads(line) for line in (RUN / "c02" / "training_records.jsonl").read_text().splitlines() if line.strip()]
    if len(records) != 400:
        errors.append("training record count")
    n0 = [row for row in records if row["arm"] == "TL-N0"]
    if len(n0) != 50:
        errors.append("TL-N0 count")
    if any(row.get("encoder_frozen") or row.get("reused_from") or row.get("steps") != 60 or row.get("encoder_trainable") != 4 for row in n0):
        errors.append("TL-N0 record flags")
    if len({row["encoder_hash_before"] for row in n0}) != 50 or any(row["encoder_hash_before"] == row["encoder_hash_after"] for row in n0):
        errors.append("TL-N0 hash evidence")
    for arm in ["TL-FS", "TL-FM", "TL-FN"]:
        rows = [row for row in records if row["arm"] == arm]
        if any(row["encoder_hash_before"] != row["encoder_hash_after"] or not row["encoder_frozen"] for row in rows):
            errors.append(f"{arm} encoder was not held fixed")
    reused = [row for row in records if row["arm"] == "TL-B"]
    if len(reused) != 50 or any("c01/dmpnn" not in str(row.get("reused_from")) for row in reused):
        errors.append("TL-B is not recorded as reused C01 B predictions")
    deltas = pd.read_csv(RUN / "c02" / "paired_deltas_recomputed.csv")
    tl01 = deltas[deltas.comparison.eq("TL01_FN_minus_JN")]
    if len(tl01) != 1 or tl01.iloc[0].family != "FTRANSFER":
        errors.append("TL01 family")
    report = (RUN / "STAGE_C_REPORT.md").read_text()
    for phrase in ["log10(Papp [cm/s])", "不表示这些网络训练已经被独立复现", "C01 XGBoost", "FTRANSFER"]:
        if phrase not in report:
            errors.append(f"report missing {phrase}")
    status = "PASS" if not errors else "FAIL"
    payload = {"status": status, "errors": errors, "n_rows": int(len(pred)), "mae_source": fresh.set_index("arm").mae_source.to_dict()}
    (RUN / "verification_result.json").write_text(json.dumps(payload, indent=2))
    print(status, errors)
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
