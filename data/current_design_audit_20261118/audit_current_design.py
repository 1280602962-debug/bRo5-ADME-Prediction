"""Recount and construct reviewable working cohorts from the existing stage-3 overlay.

This audit does not certify experimental identity or infer independent replicates.
"""
from pathlib import Path
import hashlib
import json
import platform
import zipfile

import numpy as np
import pandas as pd
from rdkit import Chem, rdBase

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "current_design_audit_20261118"
OUT.mkdir(exist_ok=True)


def read(name):
    return pd.read_csv(next(ROOT.glob(f"bro5_*/**/{name}")), keep_default_na=False, low_memory=False)


def save(frame, name):
    frame.to_csv(OUT / name, index=False, encoding="utf-8-sig")


points = read("point_candidates_after_primary_overlay.csv")
molecules = read("molecule_groups.csv")
conflicts = read("source_sample_identity_conflicts.csv")
conflict_graphs = set("|".join(conflicts.Graph_IDs).split("|"))
data = points.merge(molecules, on="Graph_ID", how="left", validate="many_to_one")
assert data.MW.notna().all()
data["MW_bin"] = pd.cut(data.MW, [0, 500, 700, np.inf], labels=["LE500", "GT500_LE700", "GT700"]).astype(str)
data["Traceable_task"] = data.Split_eligible.astype(str).str.lower().eq("true")
data["Identity_conflict"] = data.Graph_ID.isin(conflict_graphs)
data["Outside_upper"] = data.Log10_cm_s > -3
data["Below_floor"] = data.Log10_cm_s < -9
data["Nonfinite"] = ~np.isfinite(data.Log10_cm_s)
data["Invalid_structure"] = False
invalid_graphs = []
for row in molecules.itertuples():
    mol = Chem.MolFromSmiles(row.SMILES)
    if mol is None or any(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()):
        invalid_graphs.append(row.Graph_ID)
data["Invalid_structure"] = data.Graph_ID.isin(invalid_graphs)
main = data[data.Task.isin(["PAMPA:DIRECTION_NOT_APPLICABLE", "Caco-2:AB"])].copy()
main["Deposited_dataset_reference"] = main.Study_key.str.startswith("DOI:10.6019/")

# Large within-source disagreement is a review trigger, not evidence of bad experiments.
# Keep every conflicting row in the review file. The conservative cohort temporarily
# holds the entire graph-study-task group instead of silently averaging it.
keys = ["Task", "Study_key", "Graph_ID"]
spans = main[main.Traceable_task].groupby(keys, as_index=False).agg(
    Min_log=("Log10_cm_s", "min"), Max_log=("Log10_cm_s", "max"),
    Occurrences=("Occurrence_ID", "size"), Protocol_count=("Protocol", "nunique"),
    Source_count=("Source", "nunique"))
spans["Span_log10"] = spans.Max_log - spans.Min_log
main = main.merge(spans[keys + ["Span_log10"]], on=keys, how="left", validate="many_to_one")
main["Disagreement_gt05"] = main.Span_log10 > 0.5
main["Eligible_before_disagreement"] = (
    main.Traceable_task & ~main.Identity_conflict & ~main.Outside_upper
    & ~main.Nonfinite & ~main.Invalid_structure & ~main.Deposited_dataset_reference)
main["Working_inclusion"] = main.Eligible_before_disagreement & ~main.Disagreement_gt05
main["Hold_reason"] = ""
for column, reason in [("Identity_conflict", "SOURCE_SAMPLE_IDENTITY_CONFLICT"),
                       ("Outside_upper", "LOG10_CM_S_GT_MINUS3"),
                       ("Nonfinite", "NONFINITE"), ("Invalid_structure", "INVALID_OR_RADICAL"),
                       ("Deposited_dataset_reference", "DATASET_DEPOSIT_ORIGINAL_STUDY_UNRESOLVED"),
                       ("Disagreement_gt05", "WITHIN_STUDY_SPAN_GT_0.5_REVIEW")]:
    main.loc[main[column], "Hold_reason"] += reason + "|"
main.loc[~main.Traceable_task, "Hold_reason"] += "SOURCE_OR_TASK_UNRESOLVED|"
save(main[["Occurrence_ID", "Graph_ID", "Study_key", "Task", "Source", "Source_row_ID",
           "Log10_cm_s", "Protocol", "MW", "Topology", "MW_bin", "Working_inclusion",
           "Hold_reason", "Below_floor", "Span_log10"]], "occurrence_decisions.csv")
save(spans[spans.Span_log10 > 0.5], "within_study_disagreement.csv")
save(main[(main.Hold_reason != "") & main.Traceable_task], "targeted_review_rows.csv")

summary = []
layers = {"RAW_POINT": main,
          "TRACEABLE_TASK": main[main.Traceable_task],
          "BASIC_QC": main[main.Eligible_before_disagreement],
          "WORKING_CONSERVATIVE": main[main.Working_inclusion]}
for layer, frame in layers.items():
    for (task, top, weight), subset in frame.groupby(["Task", "Topology", "MW_bin"], observed=True):
        summary.append(dict(Layer=layer, Task=task, Topology=top, MW_bin=weight,
                            Structures=subset.Graph_ID.nunique(), Occurrences=len(subset),
                            Studies=subset.loc[subset.Traceable_task, "Study_key"].nunique(),
                            SC_groups=subset.Scaffold_connectivity_group.nunique(),
                            FP06_components=subset.Scaffold_FP06_component.nunique(),
                            Unassigned_stereo_structures=subset.loc[subset.Unassigned_chiral_centers > 0, "Graph_ID"].nunique()))
summary = pd.DataFrame(summary)
save(summary, "cohort_counts.csv")

work = main[main.Working_inclusion].copy()
work["Value_rounded_for_copy_collapse"] = work.Log10_cm_s.round(6)
# A shared rounded value is a numerical copy signature only, never an event ID.
unique = work.drop_duplicates(keys + ["Value_rounded_for_copy_collapse"]).copy()
unique["Linear_cm_s"] = 10.0 ** unique.Log10_cm_s
cohort = unique.groupby(keys, as_index=False).agg(
    Linear_mean_cm_s=("Linear_cm_s", "mean"), Unique_numeric_values=("Linear_cm_s", "size"))
cohort["Log10_cm_s"] = np.log10(cohort.Linear_mean_cm_s)
provenance = work.groupby(keys, as_index=False).agg(
    Source=("Source", lambda x: "|".join(sorted(set(x)))),
    Occurrence_IDs=("Occurrence_ID", lambda x: "|".join(sorted(set(x)))),
    Protocols=("Protocol", lambda x: "|".join(sorted(set(x)))),
    Description=("Assay_description", lambda x: " | ".join(sorted(set(x)))),
    Inhibitor_status=("Inhibitor_status", lambda x: "|".join(sorted(set(x)))),
    pH_reported=("pH_reported", lambda x: "|".join(sorted(set(x)))))
cohort = cohort.merge(provenance, on=keys, validate="one_to_one").merge(molecules, on="Graph_ID", validate="many_to_one")
cohort["MW_bin"] = pd.cut(cohort.MW, [0, 500, 700, np.inf], labels=["LE500", "GT500_LE700", "GT700"]).astype(str)
cohort["Below_minus9"] = cohort.Log10_cm_s < -9
cohort["Has_inhibitor_annotation"] = cohort.Description.str.contains(r"inhibitor|inhibition|\+inh|with inh|presence of.*(?:verapamil|elacridar|zosuquidar|tariquidar)", case=False, regex=True) | cohort.Protocols.str.contains(r"\+inh|_inh|inhibitor", case=False, regex=True)
cohort["Aggregation_note"] = "Provisional linear mean of distinct six-decimal log-value signatures within graph-study-task; protocols may differ; no independent-event claim"
save(cohort, "working_graph_study_labels.csv")

concentration = []
for (task, top, weight, study), part in cohort.groupby(["Task", "Topology", "MW_bin", "Study_key"]):
    concentration.append(dict(Task=task, Topology=top, MW_bin=weight, Study_key=study,
                              Structures=part.Graph_ID.nunique(), SC_groups=part.Scaffold_connectivity_group.nunique(),
                              FP06_components=part.Scaffold_FP06_component.nunique()))
save(pd.DataFrame(concentration), "study_composition.csv")
save(cohort.groupby(["Task", "Source"], as_index=False).agg(Structures=("Graph_ID", "nunique"), Labels=("Graph_ID", "size")), "source_composition.csv")
save(data.groupby("Task", as_index=False).agg(Occurrences=("Graph_ID", "size"), Structures=("Graph_ID", "nunique"), Traceable_rows=("Traceable_task", "sum")), "all_endpoint_inventory.csv")

# Reversible, performance-blind feasibility splits. Fold assignment depends only on
# the study identifier. Each variant purges training rows overlapping test identity
# or the selected chemical group, for every training arm together.
cohort["Study_fold"] = cohort.Study_key.map(lambda x: int(hashlib.sha256(("20261118|" + x).encode()).hexdigest()[:16], 16) % 5)
save(cohort[keys + ["Study_fold"]], "provisional_source_folds.csv")
fold_counts = cohort.groupby(["Task", "Graph_ID"], as_index=False).Study_fold.nunique()
cross_fold = fold_counts[fold_counts.Study_fold > 1][["Task", "Graph_ID"]]
save(cross_fold, "cross_source_fold_graphs.csv")
split_cohort = cohort.merge(cross_fold.assign(Cross_fold=True), on=["Task", "Graph_ID"], how="left")
split_cohort = split_cohort[split_cohort.Cross_fold.isna()].drop(columns="Cross_fold")
save(split_cohort, "split_eligible_working_labels.csv")
fold_results = []
assignments = []
largest_donor = {}
for task, frame in split_cohort.groupby("Task"):
    donor = frame[(frame.Topology == "MACRO_SCREEN") & frame.MW.between(500, 700, inclusive="right")]
    largest_donor[task] = donor.groupby("Study_key").Graph_ID.nunique().idxmax()
    for experiment in ["GT500_NONMACRO", "STRICT_GT700_NONMACRO"]:
        scope = (frame.Topology == "NONMACRO_SCREEN") & (frame.MW > (500 if experiment == "GT500_NONMACRO" else 700))
        global_test_studies = set(frame.loc[scope, "Study_key"])
        for fold in range(5):
            test = frame[scope & (frame.Study_fold == fold)]
            train = frame[frame.Study_fold != fold]
            if experiment.startswith("STRICT"):
                train = train[(train.MW <= 700) & ~train.Study_key.isin(global_test_studies)]
            for mode, field in [("SOURCE_GRAPH", "Graph_ID"), ("SOURCE_SC", "Scaffold_connectivity_group"), ("SOURCE_FP06", "Scaffold_FP06_component")]:
                candidate = train[~train[field].isin(test[field]) & ~train.Graph_ID.isin(test.Graph_ID)]
                assert not set(candidate.Graph_ID) & set(test.Graph_ID)
                assert not set(candidate.Study_key) & set(test.Study_key)
                assert not set(candidate[field]) & set(test[field])
                for variant in ["ALL_DONORS", "DROP_LARGEST_MACRO_STUDY"]:
                    current = candidate
                    if variant.startswith("DROP"):
                        current = current[current.Study_key != largest_donor[task]]
                    # Select one graph per training arm regardless of its source count.
                    graph = current.drop_duplicates("Graph_ID")
                    small = graph[(graph.Topology == "NONMACRO_SCREEN") & (graph.MW <= 500)]
                    macro = graph[(graph.Topology == "MACRO_SCREEN") & graph.MW.between(500, 700, inclusive="right")]
                    nonmacro = graph[(graph.Topology == "NONMACRO_SCREEN") & graph.MW.between(500, 700, inclusive="right")]
                    matched = 0
                    for low, high in [(500, 550), (550, 600), (600, 650), (650, 700)]:
                        matched += min(macro.MW.between(low, high, inclusive="right").sum(), nonmacro.MW.between(low, high, inclusive="right").sum())
                    fold_results.append(dict(Task=task, Experiment=experiment, Fold=fold, Split=mode, Variant=variant,
                        Test_graphs=test.Graph_ID.nunique(), Test_studies=test.Study_key.nunique(),
                        Test_SC=test.Scaffold_connectivity_group.nunique(), Test_FP06=test.Scaffold_FP06_component.nunique(),
                        Small_nonmacro_pool=len(small), Potential_common_base=len(small)//2,
                        Additional_small_pool=len(small)-len(small)//2,
                        Macro_500_700=len(macro), Nonmacro_500_700=len(nonmacro), MW_matched_capacity=int(matched),
                        Macro_studies=current.loc[current.Graph_ID.isin(macro.Graph_ID), "Study_key"].nunique(),
                        Nonmacro_studies=current.loc[current.Graph_ID.isin(nonmacro.Graph_ID), "Study_key"].nunique()))
                for role, part in [("TRAIN_AVAILABLE", candidate), ("TEST", test)]:
                    assignments.extend(dict(Task=task, Experiment=experiment, Fold=fold, Split=mode, Role=role,
                                            Graph_ID=r.Graph_ID, Study_key=r.Study_key) for r in part.itertuples())
folds = pd.DataFrame(fold_results)
save(folds, "split_feasibility.csv")
save(pd.DataFrame(assignments), "split_membership.csv")

checks = dict(
    Original_point_rows=len(points), Original_point_graphs=points.Graph_ID.nunique(),
    Working_rows=len(work), Working_graphs=work.Graph_ID.nunique(), Working_graph_study_labels=len(cohort),
    Cross_source_fold_task_graphs_temporarily_held=len(cross_fold),
    Identity_conflict_graphs_total=len(conflict_graphs), Invalid_or_radical_graphs_in_point_overlay=sorted(set(invalid_graphs) & set(points.Graph_ID)),
    Largest_macro_500_700_source=largest_donor,
    All_split_identity_source_and_group_checks="PASS", Models_trained=False,
    Python=platform.python_version(), RDKit=rdBase.rdkitVersion,
    Decisions={"identity_conflicts":"Temporarily hold all graph representations in the existing ten named-sample conflicts",
               "upper_bound":"Hold log10(cm/s) > -3 using reference-paper upper threshold as a provisional rule",
               "lower_bound":"Retain original values below -9, flagged; no clipping in this audit",
               "disagreement":"Hold graph-study-task groups with range >0.5 log units; conservative sensitivity cohort, not mandatory universal cutoff",
               "stereo":"Retain reported stereochemistry and flag unassigned centers; use connectivity grouping for split sensitivity",
               "protocol":"Do not require identical protocol; mixing coefficients, pH, inhibitors and ambiguous units requires endpoint-specific resolution",
               "deposit_reference":"Hold ChEMBL deposited-dataset DOI 10.6019 from the original-study benchmark until its underlying studies are resolved",
               "topology":"Existing RDKit ring-size >=12 screen retained; chemical classification remains computational"})
(OUT / "audit_checks.json").write_text(json.dumps(checks, ensure_ascii=False, indent=2), encoding="utf-8")
print(summary[summary.Layer.isin(["TRACEABLE_TASK", "WORKING_CONSERVATIVE"])].to_string(index=False))
print("CHECKS", json.dumps(checks, ensure_ascii=True))
print("SPLIT RANGES")
print(folds.groupby(["Task", "Experiment", "Split", "Variant"])[["Test_graphs", "Test_studies", "Small_nonmacro_pool", "MW_matched_capacity", "Macro_studies"]].agg(["min", "max"]).to_string())
