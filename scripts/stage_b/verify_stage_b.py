"""Independent recount of isolation, budgets, prediction coverage, MAE, and P01."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
SPLIT = REPO / "data" / "splits_v1"
RUN = REPO / "runs" / "stage_b"


def source_mae(frame: pd.DataFrame) -> float:
    step = frame.groupby(["study_id", "graph_id"], as_index=False).abs_error.mean()
    return float(step.groupby("study_id").abs_error.mean().mean())


def main() -> None:
    errors = []
    decision = json.loads((SPLIT / "budget_decision.json").read_text())
    outer = pd.read_csv(SPLIT / "outer_members.csv.gz")
    members = pd.read_csv(SPLIT / "sampling_members.csv.gz")
    test = pd.read_csv(SPLIT / "test_labels.csv.gz")
    pred = pd.read_csv(RUN / "oof_predictions.csv.gz")
    if outer.label_id.duplicated().any():
        errors.append("outer members repeat a label")
    if set(outer.loc[outer.role.eq("holdout"), "Graph_ID"]) & set(test.Graph_ID):
        errors.append("holdout graph entered the test set")
    for fold, part in test.groupby("assigned_fold"):
        studies = set(part.Study_key)
        graphs = set(part.Graph_ID)
        sc = set(part.Scaffold_connectivity_group)
        conn = set(part.Connectivity_SMILES)
        selected = members[members.block.eq("B01") & members.fold.eq(fold) & members.seed.eq(17) & members.arm.eq("BM")]
        train_graphs = set(selected.Graph_ID)
        usable = outer[outer.role.isin(["test", "train_eligible"]) & outer.Graph_ID.isin(train_graphs) & outer.assigned_fold.ne(fold)]
        if train_graphs & graphs:
            errors.append(f"fold {fold} shares a graph")
        if set(usable.Study_key) & studies:
            errors.append(f"fold {fold} shares a study")
        if set(usable.Scaffold_connectivity_group) & sc:
            errors.append(f"fold {fold} shares an SC")
        if set(usable.Connectivity_SMILES) & conn:
            errors.append(f"fold {fold} shares a connectivity")
        if part.Study_key.nunique() < 2:
            errors.append(f"fold {fold} has fewer than two test studies")
    b01 = members[members.block.eq("B01") & members.seed.eq(17)]
    for fold, part in b01.groupby("fold"):
        m = part[(part.arm.eq("BM")) & (part.component.eq("M"))]
        n = part[(part.arm.eq("BN")) & (part.component.eq("N"))]
        b = part[(part.arm.eq("B")) & (part.component.eq("B"))]
        s = part[(part.arm.eq("BS")) & (part.component.eq("S"))]
        if b.Graph_ID.nunique() != decision["b"] or s.Graph_ID.nunique() != decision["k"]:
            errors.append(f"fold {fold} base budget mismatch")
        if m.Graph_ID.nunique() != decision["k"] or n.Graph_ID.nunique() != decision["k"]:
            errors.append(f"fold {fold} donor budget mismatch")
        if m.bin50.value_counts().sort_index().to_dict() != n.bin50.value_counts().sort_index().to_dict():
            errors.append(f"fold {fold} MW quotas differ between M and N")
        groups = [set(part.loc[part.component.eq(name), "Connectivity_SMILES"]) for name in ("B", "S", "M", "N") if (part.component.eq(name)).any()]
        # component column repeats across arms; use one arm each
        groups = []
        for arm, component in (("B", "B"), ("BS", "S"), ("BM", "M"), ("BN", "N")):
            groups.append(set(part.loc[part.arm.eq(arm) & part.component.eq(component), "Connectivity_SMILES"]))
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                if groups[i] & groups[j]:
                    errors.append(f"fold {fold} connectivity reused across budget slots")
    if decision["b02"]["run"]:
        for fold in members.fold.unique():
            large = set(members.loc[members.block.eq("B01") & members.fold.eq(fold) & members.seed.eq(17) & members.arm.eq("BM") & members.component.eq("M"), "Graph_ID"])
            small = set(members.loc[members.block.eq("B02") & members.fold.eq(fold) & members.seed.eq(17) & members.arm.eq("BM") & members.component.eq("M"), "Graph_ID"])
            if not small <= large:
                errors.append(f"fold {fold} nested M is not inside the frozen 50")
    formal = pred[pred.block.ne("PILOT")]
    expected_labels = set(test.label_id.astype(str))
    for (block, model, arm, seed), part in formal.groupby(["block", "model", "arm", "seed"]):
        if arm == "ALL_DEV":
            continue
        if set(part.label_id.astype(str)) != expected_labels:
            errors.append(f"incomplete predictions {block} {model} {arm} seed {seed}")
            break
    rf = formal[formal.block.eq("B01") & formal.model.eq("RF")].copy()
    rf["abs_error"] = (rf.y - rf.y_pred).abs()
    averaged = rf.groupby(["arm", "label_id"], as_index=False).agg(abs_error=("abs_error", "mean"), study_id=("study_id", "first"), graph_id=("graph_id", "first"))
    maes = {arm: source_mae(chunk) for arm, chunk in averaged.groupby("arm")}
    delta = maes["BN"] - maes["BM"]
    saved = pd.read_csv(RUN / "paired_deltas.csv")
    row = saved[(saved.block.eq("B01")) & (saved.model.eq("RF")) & (saved.comparison.eq("P01"))].iloc[0]
    if abs(delta - float(row.estimate)) > 1e-8:
        errors.append(f"recomputed P01 {delta} != saved {row.estimate}")
    result = {"status": "PASS" if not errors else "FAIL", "errors": errors, "p01": delta, "mae": maes}
    (RUN / "verification_result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
