"""Independent checks on exported decisions, labels and split memberships."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

root = Path(__file__).resolve().parent
folder = root / "current_design_audit_20261118"
decisions = pd.read_csv(folder / "occurrence_decisions.csv", keep_default_na=False)
labels = pd.read_csv(folder / "working_graph_study_labels.csv", keep_default_na=False)
membership = pd.read_csv(folder / "split_membership.csv", keep_default_na=False)
groups = pd.read_csv(next(root.glob("bro5_*/**/molecule_groups.csv")), keep_default_na=False)
raw = pd.read_csv(next(root.glob("bro5_*/**/point_candidates_after_primary_overlay.csv")), keep_default_na=False, low_memory=False)
expected = raw[raw.Task.isin(["Caco-2:AB", "PAMPA:DIRECTION_NOT_APPLICABLE"])]
assert len(decisions) == len(expected)
assert set(decisions.Occurrence_ID) == set(expected.Occurrence_ID)
chosen = decisions[decisions.Working_inclusion]
assert not chosen.Graph_ID.isna().any()
assert not (chosen.Log10_cm_s > -3).any()
assert not chosen.Study_key.str.startswith("DOI:10.6019/").any()
assert not chosen.Hold_reason.ne("").any()
keys = ["Task", "Study_key", "Graph_ID"]
assert not labels.duplicated(keys).any()
numeric = chosen.assign(Value_signature=chosen.Log10_cm_s.round(6)).drop_duplicates(keys + ["Value_signature"])
numeric["Value"] = 10.0 ** numeric.Log10_cm_s
recomputed = numeric.groupby(keys).Value.mean().pipe(np.log10).rename("Recomputed")
check = labels.set_index(keys).join(recomputed)
assert np.allclose(check.Log10_cm_s, check.Recomputed, atol=1e-12)
membership = membership.merge(groups[["Graph_ID", "Scaffold_connectivity_group", "Scaffold_FP06_component", "MW"]], on="Graph_ID", validate="many_to_one")
checks = 0
for (task, experiment, fold, split), rows in membership.groupby(["Task", "Experiment", "Fold", "Split"]):
    train = rows[rows.Role == "TRAIN_AVAILABLE"]
    test = rows[rows.Role == "TEST"]
    fields = ["Graph_ID", "Study_key"]
    if split == "SOURCE_SC":
        fields.append("Scaffold_connectivity_group")
    if split == "SOURCE_FP06":
        fields.append("Scaffold_FP06_component")
    for field in fields:
        assert not set(train[field]) & set(test[field]), (task, experiment, fold, split, field)
        checks += 1
    if experiment.startswith("STRICT"):
        assert train.MW.le(700).all()
        assert test.MW.gt(700).all()
unique_test = membership[membership.Role == "TEST"].groupby(["Task", "Experiment", "Split", "Graph_ID"]).Fold.nunique()
assert unique_test.max() == 1
manifest = []
for name in ["point_candidates_after_primary_overlay.csv", "molecule_groups.csv", "source_sample_identity_conflicts.csv"]:
    path = next(root.glob(f"bro5_*/**/{name}"))
    manifest.append(dict(Path=str(path.relative_to(root)), SHA256=hashlib.sha256(path.read_bytes()).hexdigest()))
result = dict(Status="PASS", Retained_input_main_rows=len(decisions), Working_graph_study_labels=len(labels),
              Explicit_source_identity_chemical_group_checks=checks, Recomputed_label_values="PASS",
              Test_graph_assigned_to_one_fold="PASS", Models_trained=False, Inputs=manifest)
(folder / "independent_verification.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, ensure_ascii=True))
