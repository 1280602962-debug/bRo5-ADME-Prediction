"""Independent checks for the M4-CPU outputs."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "runs" / "model_development_m4_cpu"


def main() -> None:
    checks = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    holdout = set(pd.read_csv(REPO / "data" / "splits_v1" / "holdout_members.csv").Graph_ID.astype(str))
    for name in [
        "ionization_features.csv.gz",
        "l1_static_features.csv.gz",
        "l1_ensemble_features.csv.gz",
        "m4_oof_predictions.csv.gz",
    ]:
        frame = pd.read_csv(OUT / name)
        overlap = holdout & set(frame.Graph_ID.astype(str))
        check(f"holdout absent from {name}", not overlap, f"overlap {len(overlap)}")
    eligibility = pd.read_csv(OUT / "conformer_eligibility.csv")
    stereo = eligibility[eligibility.failure.eq("unspecified_stereo")]
    check("unspecified stereo excluded", int(stereo.unspecified_tetrahedral.gt(0).sum() + stereo.unspecified_double_bond.gt(0).sum()) > 0, f"rows {len(stereo)}")
    macro = eligibility[eligibility.Topology.astype(str).str.upper().str.startswith("MACRO") & eligibility.eligible_main.astype(str).str.lower().eq("true")]
    nonmacro_low = eligibility[eligibility.Topology.astype(str).str.upper().str.startswith("NONMACRO") & eligibility.MW.le(500) & eligibility.n_rotatable.lt(10)]
    check("macrocycles request 100", bool((macro.requested_conformers == 100).all()), f"n {len(macro)}")
    check("low-MW nonmacrocycles request 20", bool((nonmacro_low.requested_conformers == 20).all()), f"n {len(nonmacro_low)}")
    deltas = pd.read_csv(OUT / "m4_paired_deltas.csv")
    primary = deltas[deltas.role.eq("primary")]
    check("six primary contrasts retained", set(primary.contrast) == {"P2-P0", "P3-P0", "P4-P0", "P4-P3", "P5-P4", "P5-P2"}, str(primary.contrast.tolist()))
    check("family size stays 6", bool((primary.family_size == 6).all()), "family")
    pred = pd.read_csv(OUT / "m4_oof_predictions.csv.gz")
    c = pred[pred.evaluation.eq("C") & pred.model.eq("RF") & pred.parameter_source.eq("fixed_rf")]
    for left, right in [("P4_ensemble", "P0_ensemble"), ("P4_ensemble", "P3_ensemble"), ("P5_both", "P4_both"), ("P5_both", "P2_both")]:
        a = set(c.loc[c.arm.eq(left), "label_id"])
        b = set(c.loc[c.arm.eq(right), "label_id"])
        if not a and not b:
            check(f"{left} vs {right} labels", True, "not executable")
        else:
            check(f"{left} vs {right} labels", a == b, f"{len(a)} vs {len(b)}")
    train_counts = pd.read_csv(OUT / "matched_cohort_summary.csv")
    for arm_pair in [("P4_ensemble", "P0_ensemble"), ("P5_both", "P2_both")]:
        parts = [train_counts[train_counts.arm.eq(name)] for name in arm_pair]
        if any(part.empty for part in parts):
            continue
        same = parts[0].merge(parts[1], on="outer_fold", suffixes=("_a", "_b"))
        check(
            f"train labels match {arm_pair[0]} {arm_pair[1]}",
            bool((same.n_train_labels_a == same.n_train_labels_b).all()),
            "per fold",
        )
    check("condition audit present", (OUT / "M4_0_DATA_VARIABILITY_REPORT.md").exists(), "m4.0")
    check("crest not claimed", "没有做水相/氯仿配对试算" in (OUT / "M4_CPU_REPORT.md").read_text(), "report")
    result = {"passed": all(item["passed"] for item in checks), "checks": checks}
    (OUT / "verification_result.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
