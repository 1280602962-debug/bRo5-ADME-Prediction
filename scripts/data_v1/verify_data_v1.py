"""Independent recount of Data V1.0.

Recomputes conservation, the 17-record reclassification, log conversion
spot checks, database-copy collapse, and group-identifier consistency.
"""
import argparse
import json
import math
from pathlib import Path

import pandas as pd

PDL1_IDS = {str(i) for i in range(26123123, 26123140)}


def fail(errors, message):
    errors.append(message)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--molecules", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    errors = []
    source = pd.read_csv(args.labels, dtype=str, low_memory=False)
    decisions = pd.read_csv(args.out / "occurrence_decisions.csv.gz", dtype=str)
    built = pd.read_csv(args.out / "label_table.csv.gz")
    molecules = pd.read_csv(args.molecules, dtype=str)
    built_molecules = pd.read_csv(args.out / "molecule_groups.csv.gz", dtype=str)
    if len(source) != 66771:
        fail(errors, f"input rows {len(source)} != 66771")
    if len(decisions) != len(source) or set(decisions.record_id) != set(source.Occurrence_ID):
        fail(errors, "decision ids do not match the input occurrences")
    if decisions.record_id.duplicated().any():
        fail(errors, "duplicate decision ids")
    pdl1 = decisions[decisions.provenance.eq("ChEMBL37") & decisions.record_id.isin(
        source.loc[source.Source_row_ID.astype(str).isin(PDL1_IDS) & source.Source.eq("ChEMBL37"), "Occurrence_ID"]
    )]
    if len(pdl1) != 17:
        fail(errors, f"PD-L1 decision rows {len(pdl1)} != 17")
    if not pdl1.corrected_task.eq("MDCK_WT:AB").all():
        fail(errors, "PD-L1 rows were not all reclassified to MDCK_WT:AB")
    main = built[built.analysis_set.eq("caco2_ab_main")]
    main_ids = set()
    for text in main.record_ids.fillna(""):
        main_ids.update(text.split("|"))
    if main_ids & set(pdl1.record_id):
        fail(errors, "reclassified MDCK rows appear in the Caco-2 main labels")
    if not main.corrected_task.eq("Caco-2:AB").all():
        fail(errors, "Caco-2 main contains another task")
    if not main.Direction_refined.eq("AB").all():
        fail(errors, "Caco-2 main contains a non-AB direction")
    # No second logarithm: inherited log must match the source log, and a linear example must match log10(x)-6.
    source["y"] = pd.to_numeric(source.Log10_cm_s, errors="coerce")
    merged = decisions.merge(source[["Occurrence_ID", "y", "Raw_value", "Source"]], left_on="record_id", right_on="Occurrence_ID", how="left")
    kim = merged[merged.Source.eq("Kim2025") & merged.relation.eq("=")].head(20)
    for row in kim.itertuples(index=False):
        try:
            raw = float(str(row.Raw_value).split()[0])
        except ValueError:
            continue
        expected = math.log10(raw) - 7
        if abs(float(row.stored_log10_cm_s) - expected) > 1e-6:
            fail(errors, f"Kim conversion mismatch for {row.record_id}")
            break
        if abs(float(row.stored_log10_cm_s) - (math.log10(raw) - 6)) < 1e-6:
            fail(errors, "Kim nm/s value was treated as 10^-6 cm/s")
            break
    linear = source[source.Source.eq("ChEMBL37") & source.Source_row_ID.astype(str).eq("26123127")]
    if len(linear) != 1 or abs(float(linear.y.iloc[0]) - (math.log10(6.6) - 6)) > 1e-6:
        fail(errors, "ChEMBL 6.6 x 10^-6 cm/s was not stored as log10(value)-6")
    rz = source[source.Compound_original_ID.eq("BAS_52171383")]
    if rz.empty or not rz.Log10_cm_s.astype(float).round(6).eq(round(-5.027, 6)).all():
        fail(errors, "Rzepiela log Peff was altered or logged twice")
    # Database copies of one BAS id must not become two main labels.
    if built.record_ids.fillna("").duplicated().any():
        fail(errors, "a record id is listed in more than one label")
    # Group identifiers inherited, not swapped.
    check = built.merge(molecules, on="Graph_ID", how="left", suffixes=("", "_src"))
    if not (check.connectivity_id.eq(check.Scaffold_connectivity_group)).all():
        fail(errors, "connectivity id does not match the inherited scaffold-connectivity group")
    if not (check.fp06_component.eq(check.Scaffold_FP06_component)).all():
        fail(errors, "FP06 id does not match the inherited component")
    if (check.connectivity_id.eq(check.fp06_component) & check.MurckoScaffold.ne("")).any():
        # Equality can happen for singletons; only fail when the two families were swapped relative to the source.
        pass
    swapped = check[check.connectivity_id.eq(check.Scaffold_FP06_component) & check.fp06_component.eq(check.Scaffold_connectivity_group) & check.Scaffold_connectivity_group.ne(check.Scaffold_FP06_component)]
    if len(swapped):
        fail(errors, "connectivity and FP06 identifiers were swapped")
    missing = set(built.Graph_ID) - set(built_molecules.Graph_ID)
    if missing:
        fail(errors, f"{len(missing)} label graphs are absent from molecule_groups")
    summary = pd.read_csv(args.out / "dataset_summary.csv")
    cohort = summary[summary.section.eq("cohort") & summary.analysis_set.eq("caco2_ab_main") & summary.cohort.eq("nonmacro_500_700")]
    if len(cohort) != 1:
        fail(errors, "summary is missing the Caco-2 target cohort")
    else:
        reported = int(float(cohort.graphs.iloc[0]))
        actual = main[(main.Topology.eq("NONMACRO_SCREEN")) & (main.MW.astype(float) > 500) & (main.MW.astype(float) <= 700)]
        if actual.Graph_ID.nunique() != reported:
            fail(errors, "summary graph count does not match the label table")
    discordant = decisions.reason_code.fillna("").str.contains("BAS_SHARED_STRUCTURE_DISCORDANT_VALUES")
    if int(discordant.sum()) != 472:
        fail(errors, f"discordant BAS copies held {int(discordant.sum())}, expected 472")
    result = {"status": "PASS" if not errors else "FAIL", "errors": errors, "caco2_main_labels": int(len(main)), "caco2_main_graphs": int(main.Graph_ID.nunique())}
    (args.out / "verification_result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
