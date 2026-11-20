"""Build Data V1.0 labels from the audited occurrence table.

This script does not retrain models, generate conformers, or recompute
Graph_ID, connectivity, Murcko, or FP06 groups. Those identifiers are
inherited from the frozen molecule table.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

PDL1_ACTIVITY_IDS = {str(i) for i in range(26123123, 26123140)}
PDL1_DOI = "10.1016/j.ejmech.2024.116730"
RZEPIELA_DOI = "10.1021/acs.jmedchem.1c02090"
KIM_DOI = "10.1021/acs.jcim.5c01600"
STAGE3_GEOMETRIC = {
    "O_10f5b2bdce46f5db05f45b2f",
    "O_b458993ccf65ba03c426fd7b",
    "O_cf64fbc4f23edb7013e325ac",
}
STAGE3_ALIAS = {"O_a0d9efaa8694e5481000c297": "G_e9503599b2a73e145750fbf8"}
STAGE3_BOUND = {"O_b42b5ff0dce0790fd0d9fd98", "O_1bdf3a153125dae9ca92d1d1"}
STAGE3_PRECISION = {"O_02873b2d7247db0140df18bc": -6.56}
NAMED_INHIBITOR = r"elacridar|tariquidar|zosuquidar|verapamil|gf120918|cyclosporin"
POSITIVE_INHIBITOR = re.compile(
    rf"(?:in the presence of|treated with|co-?incubat\w* with|addition of|\bwith\b|\bplus\b|\+).{{0,40}}(?:{NAMED_INHIBITOR})",
    re.I,
)


def as_bool(value):
    return str(value).strip().lower() == "true"


def clean_text(value):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "null"} else text


def hid(prefix, value):
    return prefix + "_" + hashlib.sha256(str(value).encode()).hexdigest()[:16]


def mw_cohort(mw, topology):
    nonmacro = topology == "NONMACRO_SCREEN"
    macro = topology == "MACRO_SCREEN"
    if nonmacro and mw <= 500:
        return "nonmacro_le500"
    if nonmacro and 500 < mw <= 700:
        return "nonmacro_500_700"
    if macro and 500 < mw <= 700:
        return "macro_500_700"
    if nonmacro and mw > 700:
        return "nonmacro_gt700"
    return "other"


def mw_bin_50(mw):
    edges = [(500, 550), (550, 600), (600, 650), (650, 700)]
    for low, high in edges:
        if low < mw <= high:
            return f"{low}_{high}"
    return ""


def pampa_coefficient(endpoint, description, protocol):
    ep = clean_text(endpoint).upper()
    if ep == "PAMPA_INTRINSIC":
        return "PAMPA_INTRINSIC"
    if ep == "PAMPA_APPARENT":
        return "PAMPA_APPARENT"
    text = f"{endpoint} {description} {protocol}".lower()
    if "log peff" in text or re.search(r"\bpeff\b|p_eff", text):
        return "PAMPA_PEFF"
    if "intrinsic" in text:
        return "PAMPA_INTRINSIC"
    if re.search(r"\bp0\b", text):
        return "PAMPA_P0"
    if "papp" in text or "apparent" in text:
        return "PAMPA_PAPP"
    return "PAMPA_COEFFICIENT_UNSPECIFIED"


def intervention_class(protocol, description, existing):
    text = f"{protocol} {description}"
    if "+inh" in text.lower() or POSITIVE_INHIBITOR.search(text):
        return "KNOWN_INTERVENTION"
    if clean_text(existing) == "MENTION_REQUIRES_REVIEW" or re.search(NAMED_INHIBITOR, text, re.I):
        return "MENTION_NOT_CONFIRMED"
    return "NOT_REPORTED"


def specific_protocol(protocol, description):
    protocol = clean_text(protocol)
    description = clean_text(description)
    if protocol.startswith("CHEMBL_ASSAY:"):
        return protocol
    if len(description) >= 40 and "log peff" in description.lower():
        return "STATED_METHOD:" + hashlib.sha256(description.encode()).hexdigest()[:12]
    if protocol and not protocol.startswith("DATABASE_ASSAY") and len(protocol) >= 24:
        return "STATED_PROTOCOL:" + hashlib.sha256(protocol.encode()).hexdigest()[:12]
    return ""


def load_rzepiela_pairs(path):
    raw = pd.read_csv(path, dtype=str)
    raw = raw[raw.Source.fillna("").str.contains(RZEPIELA_DOI, regex=False)].copy()
    if len(raw) != 3348:
        raise RuntimeError(f"Rzepiela matched rows are {len(raw)}, expected 3348")
    counts = raw.groupby("SMILES").Original_ID.nunique()
    shared = counts[counts > 1]
    if len(shared) != 120:
        raise RuntimeError(f"Shared-SMILES BAS groups are {len(shared)}, expected 120")
    pair = raw[raw.SMILES.isin(shared.index)].copy()
    pair["value"] = pd.to_numeric(pair.Standardized_Value, errors="coerce")
    discordant = []
    concordant = []
    for smiles, part in pair.groupby("SMILES"):
        ids = sorted(set(part.Original_ID))
        nval = part["value"].round(6).nunique(dropna=True)
        (discordant if nval > 1 else concordant).extend(ids)
    sentinel = set(raw.loc[pd.to_numeric(raw.Standardized_Value, errors="coerce").eq(-8), "Original_ID"])
    full = pd.read_csv(path, dtype=str)
    archive_ids = {str(i + 2): clean_text(original) for i, original in enumerate(full.Original_ID)}
    return {
        "discordant_ids": set(discordant),
        "concordant_ids": set(concordant),
        "sentinel_ids": sentinel,
        "n_rows": len(raw),
        "n_shared_smiles": int(len(shared)),
        "n_discordant_ids": len(set(discordant)),
        "n_concordant_ids": len(set(concordant)),
        "n_sentinel_ids": len(sentinel),
        "archive_ids": archive_ids,
    }


def load_conflicts(path):
    frame = pd.read_csv(path, dtype=str)
    pairs = []
    for row in frame.itertuples(index=False):
        for graph in clean_text(row.Graph_IDs).split("|"):
            if graph:
                pairs.append((row.Study_key, graph, row.Resolution))
    return pd.DataFrame(pairs, columns=["Study_key", "Graph_ID", "Conflict_resolution"])


def apply_stage3(frame):
    frame = frame.copy()
    frame["Stage3_action"] = ""
    geometric = frame.Occurrence_ID.isin(STAGE3_GEOMETRIC)
    frame.loc[geometric, "Task"] = "Caco-2:GEOMETRIC_MEAN_AB_BA"
    frame.loc[geometric, "Direction_refined"] = "GEOMETRIC_MEAN_AB_BA"
    frame.loc[geometric, "Stage3_action"] = "RESOLVE_GEOMETRIC_MEAN_TASK"
    for occurrence, graph in STAGE3_ALIAS.items():
        hit = frame.Occurrence_ID.eq(occurrence)
        frame.loc[hit, "Graph_ID"] = graph
        frame.loc[hit, "Stage3_action"] = "CORRECT_GRAPH_ALIAS"
    bound = frame.Occurrence_ID.isin(STAGE3_BOUND)
    frame.loc[bound, "Relation"] = "<"
    frame.loc[bound, "Stage3_action"] = "POINT_TO_PRIMARY_BOUND"
    for occurrence, value in STAGE3_PRECISION.items():
        hit = frame.Occurrence_ID.eq(occurrence)
        frame.loc[hit, "Log10_cm_s"] = value
        frame.loc[hit, "Stage3_action"] = "NORMALIZE_PRIMARY_REPORTED_PRECISION"
    expected = {"RESOLVE_GEOMETRIC_MEAN_TASK": 3, "CORRECT_GRAPH_ALIAS": 1, "POINT_TO_PRIMARY_BOUND": 2, "NORMALIZE_PRIMARY_REPORTED_PRECISION": 1}
    found = frame.Stage3_action.value_counts().to_dict()
    for key, count in expected.items():
        if found.get(key, 0) != count:
            raise RuntimeError(f"Stage3 action {key} count {found.get(key, 0)} != {count}")
    return frame


def assign_tasks(frame, rzepiela):
    frame = frame.copy()
    frame["y"] = pd.to_numeric(frame.Log10_cm_s, errors="coerce")
    frame["Compound_original_ID"] = frame.Compound_original_ID.map(clean_text)
    archive = frame.Source.eq("NPMMP_archive") & frame.Compound_original_ID.eq("")
    frame.loc[archive, "Compound_original_ID"] = frame.loc[archive, "Source_row_ID"].map(rzepiela["archive_ids"]).fillna("")
    frame["Protocol"] = frame.Protocol.map(clean_text)
    frame["Assay_description"] = frame.Assay_description.map(clean_text)
    frame["Study_key"] = frame.Study_key.map(clean_text)
    pdl1 = frame.Source.eq("ChEMBL37") & frame.Source_row_ID.astype(str).isin(PDL1_ACTIVITY_IDS)
    if int(pdl1.sum()) != 17:
        raise RuntimeError(f"PD-L1 block has {int(pdl1.sum())} rows, expected 17")
    if not frame.loc[pdl1, "Study_key"].eq(f"DOI:{PDL1_DOI}").all():
        raise RuntimeError("PD-L1 block is not confined to the expected DOI")
    frame["corrected_task"] = frame.Task
    frame.loc[pdl1, "corrected_task"] = "MDCK_WT:AB"
    frame.loc[pdl1, "Cell_model"] = "MDCKII-WT"
    frame["correction_note"] = ""
    frame.loc[pdl1, "correction_note"] = "PDL1_17_RECLASSIFIED_MDCKII_WT"
    frame["coefficient"] = ""
    pampa = frame.Endpoint.fillna("").str.startswith("PAMPA") | frame.corrected_task.str.startswith("PAMPA")
    frame.loc[pampa, "coefficient"] = [
        pampa_coefficient(ep, desc, protocol)
        for ep, desc, protocol in zip(frame.loc[pampa, "Endpoint"], frame.loc[pampa, "Assay_description"], frame.loc[pampa, "Protocol"])
    ]
    frame.loc[pampa, "corrected_task"] = "PAMPA:" + frame.loc[pampa, "coefficient"]
    frame["intervention"] = [
        intervention_class(protocol, description, existing)
        for protocol, description, existing in zip(frame.Protocol, frame.Assay_description, frame.Inhibitor_status)
    ]
    frame["rzepiela"] = frame.Study_key.eq(f"DOI:{RZEPIELA_DOI}")
    frame["identity_issue"] = ""
    discordant = frame.rzepiela & frame.Compound_original_ID.isin(rzepiela["discordant_ids"])
    concordant = frame.rzepiela & frame.Compound_original_ID.isin(rzepiela["concordant_ids"])
    frame.loc[discordant, "identity_issue"] = "BAS_SHARED_STRUCTURE_DISCORDANT_VALUES"
    frame.loc[concordant, "identity_issue"] = "BAS_SHARED_STRUCTURE_SAME_VALUE_COLLAPSE"
    conflict = load_conflicts  # set by caller via column if merged
    return frame, pdl1


def decide(frame, conflicts, sentinel_ids):
    frame = frame.merge(conflicts, on=["Study_key", "Graph_ID"], how="left", validate="many_to_one")
    named = frame.Conflict_resolution.fillna("")
    frame.loc[named.eq("CONNECTIVITY_CONFLICT_REVIEW"), "identity_issue"] = "NAMED_SAMPLE_CONNECTIVITY_CONFLICT"
    stereo_name = named.eq("STEREO_REPRESENTATION_REVIEW") & frame.identity_issue.eq("")
    frame.loc[stereo_name, "identity_issue"] = "NAMED_SAMPLE_STEREO_ASSIGNMENT_AMBIGUOUS"
    reasons = []
    for row in frame.itertuples(index=False):
        code = []
        if row.correction_note:
            code.append(row.correction_note)
        if row.Stage3_action:
            code.append(row.Stage3_action)
        if not as_bool(row.Stage2_point_eligible):
            code.append("STAGE2_HOLD:" + (clean_text(row.Stage2_hold_reason) or "UNSPECIFIED"))
        if clean_text(row.Relation) != "=":
            code.append("QUALIFIER_NOT_EXACT")
        if not np.isfinite(row.y):
            code.append("NONFINITE")
        elif row.y > -3:
            code.append("LOG10_CM_S_GT_MINUS3")
        if row.identity_issue in {
            "BAS_SHARED_STRUCTURE_DISCORDANT_VALUES",
            "NAMED_SAMPLE_CONNECTIVITY_CONFLICT",
            "NAMED_SAMPLE_STEREO_ASSIGNMENT_AMBIGUOUS",
        }:
            code.append(row.identity_issue)
        if row.rzepiela and row.Compound_original_ID in sentinel_ids and np.isfinite(row.y) and abs(row.y + 8) < 1e-9:
            code.append("RZEPIELA_MINUS8_CENSORING_UNRESOLVED")
        if row.Study_key.startswith("DOI:10.6019/"):
            code.append("DATASET_DEPOSIT_DOI_NOT_ORIGINAL_STUDY")
        reasons.append("|".join(code))
    frame["reason_code"] = reasons
    blocking = {
        "STAGE2_HOLD",
        "QUALIFIER_NOT_EXACT",
        "NONFINITE",
        "LOG10_CM_S_GT_MINUS3",
        "BAS_SHARED_STRUCTURE_DISCORDANT_VALUES",
        "NAMED_SAMPLE_CONNECTIVITY_CONFLICT",
        "NAMED_SAMPLE_STEREO_ASSIGNMENT_AMBIGUOUS",
        "RZEPIELA_MINUS8_CENSORING_UNRESOLVED",
        "ORIGINAL_CENSORED_NOT_EXACT_POINT",
        "ORIGINAL_CENSORED_POINT_OR_POTENTIAL_CROSS_SOURCE_COPY_REVIEW",
        "PRIMARY_LOW_RECOVERY_AND_CENSORED_RESULT_REVIEW",
    }
    # Stage2 hold reasons are included whole. Also block any reason containing these tokens.
    def blocked(text):
        parts = [part for part in text.split("|") if part]
        if any(part in blocking or part.startswith("STAGE2_HOLD") for part in parts):
            return True
        if any(token in text for token in [
            "CENSORED", "LOW_RECOVERY", "NOT_EXACT", "NONFINITE", "GT_MINUS3",
            "DISCORDANT", "CONNECTIVITY_CONFLICT", "STEREO_ASSIGNMENT", "MINUS8",
        ]):
            return True
        return False
    frame["blocked"] = frame.reason_code.map(blocked)
    # A stage2-eligible exact row with only a reclassification note is not blocked.
    frame.loc[frame.reason_code.eq("PDL1_17_RECLASSIFIED_MDCKII_WT"), "blocked"] = False
    frame.loc[frame.reason_code.eq("NORMALIZE_PRIMARY_REPORTED_PRECISION"), "blocked"] = False
    frame.loc[frame.reason_code.eq("CORRECT_GRAPH_ALIAS"), "blocked"] = False
    frame.loc[frame.reason_code.eq(""), "blocked"] = False
    return frame


def measurement_id(row):
    compound = clean_text(row.Compound_original_ID)
    rounded = "" if not np.isfinite(row.y) else f"{row.y:.6f}"
    if clean_text(row.identity_issue) == "BAS_SHARED_STRUCTURE_SAME_VALUE_COLLAPSE":
        return hid("M", ("BAS_SAME_VALUE", row.Graph_ID, row.Study_key, row.corrected_task, row.Direction_refined, row.Relation, rounded))
    if compound:
        return hid("M", (compound, row.Graph_ID, row.Study_key, row.corrected_task, row.Direction_refined, row.Relation, rounded))
    snapshot = clean_text(row.Snapshot_record_key)
    if snapshot:
        return hid("M", ("SNAP", snapshot, row.corrected_task))
    return hid("M", ("OCC", row.Occurrence_ID))


def condition_of(row, measurement):
    stated = specific_protocol(row.Protocol, row.Assay_description)
    body = stated + "|" + row.intervention if stated else "UNSTATED|" + measurement
    deposit = clean_text(row.Study_key).startswith("DOI:10.6019/")
    resolution = "NONSTRICT" if deposit or not as_bool(row.Split_eligible) else "STRICT"
    return resolution + "|" + body


def analysis_set(row, auxiliary):
    if clean_text(row.Study_key).startswith("DOI:10.6019/"):
        return "audit_deposit_doi"
    strict = as_bool(row.Split_eligible)
    if row.intervention == "KNOWN_INTERVENTION" and row.corrected_task == "Caco-2:AB" and strict:
        return "caco2_known_intervention"
    if row.corrected_task == "Caco-2:AB" and strict:
        return "caco2_ab_main"
    if row.corrected_task.startswith("PAMPA:"):
        coeff = row.corrected_task.split(":", 1)[1]
        role = "auxiliary" if coeff == auxiliary else "exploratory"
        if not strict:
            role = "exploratory_nonstrict"
        return f"pampa_{coeff}_{role}"
    if not strict:
        return "separated_nonstrict_source"
    return "separated_other_endpoint"


def aggregate(exact):
    exact = exact.copy()
    exact["_stated"] = ~exact.condition_id.astype(str).str.contains("UNSTATED|", regex=False)
    measurements = (
        exact.sort_values(["_stated", "Occurrence_ID"], ascending=[False, True])
        .groupby("measurement_id", as_index=False)
        .agg(
            Graph_ID=("Graph_ID", "first"),
            Study_key=("Study_key", "first"),
            corrected_task=("corrected_task", "first"),
            Direction_refined=("Direction_refined", "first"),
            condition_id=("condition_id", "first"),
            y=("y", "first"),
            n_database_occurrences=("Occurrence_ID", "size"),
            Occurrence_IDs=("Occurrence_ID", lambda s: "|".join(sorted(set(s)))),
            Sources=("Source", lambda s: "|".join(sorted(set(s)))),
            intervention=("intervention", "first"),
            Protocol=("Protocol", "first"),
            below_floor=("below_floor", "max"),
            protocol_unknown=("protocol_unknown", "max"),
            analysis_set=("analysis_set", "first"),
        )
    )
    grouped = []
    for key, part in measurements.groupby(
        ["Graph_ID", "Study_key", "corrected_task", "Direction_refined", "condition_id"], sort=False
    ):
        values = part.y.to_numpy(dtype=float)
        span = float(values.max() - values.min()) if len(values) else np.nan
        shift = len(values) == 2 and abs(span - round(span)) < 0.02 and round(span) in {1, 2, 3}
        stated = "UNSTATED|" not in str(key[4])
        if span > 0.5 and shift:
            status = "HOLD_POSSIBLE_UNIT_SHIFT"
        elif span > 0.5 and not stated:
            status = "HOLD_UNEXPLAINED_RANGE"
        else:
            status = "LABEL"
        grouped.append({
            "Graph_ID": key[0], "Study_key": key[1], "corrected_task": key[2],
            "Direction_refined": key[3], "condition_id": key[4],
            "y": float(np.median(values)),
            "y_mean": float(values.mean()),
            "y_sd": float(values.std(ddof=1)) if len(values) > 1 else np.nan,
            "y_range": span,
            "n_distinct_measurements": int(len(part)),
            "n_database_occurrences": int(part.n_database_occurrences.sum()),
            "record_ids": "|".join(part.Occurrence_IDs),
            "sources": "|".join(sorted(set("|".join(part.Sources).split("|")))),
            "intervention": part.intervention.iloc[0],
            "below_floor": bool(part.below_floor.max()),
            "protocol_unknown": bool(part.protocol_unknown.max()),
            "analysis_set": part.analysis_set.iloc[0],
            "variance_status": "RETAINED_COMPATIBLE_RANGE_GT_0.5" if status == "LABEL" and span > 0.5 else status,
        })
    labels = pd.DataFrame(grouped)
    if labels.empty:
        return labels, measurements
    held = labels.variance_status.isin(["HOLD_POSSIBLE_UNIT_SHIFT", "HOLD_UNEXPLAINED_RANGE"])
    labels["label_status"] = np.where(held, "HELD_AFTER_AGGREGATION", "EXACT_LABEL")
    labels["label_id"] = [
        hid("L", (r.Graph_ID, r.Study_key, r.corrected_task, r.Direction_refined, r.condition_id))
        for r in labels.itertuples(index=False)
    ]
    return labels, measurements


def choose_auxiliary(labels, molecules):
    use = labels[labels.label_status.eq("EXACT_LABEL") & labels.corrected_task.str.startswith("PAMPA:")].copy()
    if use.empty:
        return "", []
    use = use.merge(molecules[["Graph_ID", "MW", "Topology"]], on="Graph_ID", how="left", validate="many_to_one")
    rows = []
    for task, part in use.groupby("corrected_task"):
        coeff = task.split(":", 1)[1]
        target = part[(part.Topology == "NONMACRO_SCREEN") & (part.MW > 500) & (part.MW <= 700)]
        studies = target.Study_key.nunique()
        graphs = target.Graph_ID.nunique()
        largest = target.groupby("Study_key").Graph_ID.nunique().max() if len(target) else 0
        fraction = (largest / graphs) if graphs else np.nan
        defined = coeff != "PAMPA_COEFFICIENT_UNSPECIFIED"
        eligible = bool(defined and graphs >= 30 and studies >= 5 and fraction <= 0.80)
        rows.append({
            "coefficient": coeff, "defined": defined, "labels": len(part),
            "graphs": part.Graph_ID.nunique(), "studies": part.Study_key.nunique(),
            "nonmacro_500_700_graphs": graphs, "nonmacro_500_700_studies": studies,
            "largest_target_study_fraction": fraction, "eligible_auxiliary": eligible,
            "macro_500_700_graphs": int(((part.Topology == "MACRO_SCREEN") & part.MW.between(500, 700, inclusive="right")).sum() and part.loc[(part.Topology == "MACRO_SCREEN") & part.MW.between(500, 700, inclusive="right"), "Graph_ID"].nunique()),
            "nonmacro_gt700_graphs": int(part.loc[(part.Topology == "NONMACRO_SCREEN") & (part.MW > 700), "Graph_ID"].nunique()),
        })
    table = pd.DataFrame(rows).sort_values(["eligible_auxiliary", "nonmacro_500_700_studies", "nonmacro_500_700_graphs"], ascending=False)
    eligible = table[table.eligible_auxiliary]
    auxiliary = "" if eligible.empty else "PAMPA:" + eligible.iloc[0].coefficient
    return auxiliary, table


def summarize(labels, molecules, decisions):
    exact = labels[labels.label_status.eq("EXACT_LABEL")].merge(
        molecules[["Graph_ID", "MW", "Topology", "Scaffold_connectivity_group", "Scaffold_FP06_component", "MurckoScaffold", "Unassigned_chiral_centers"]],
        on="Graph_ID", how="left", validate="many_to_one",
    )
    exact["cohort"] = [mw_cohort(mw, top) for mw, top in zip(exact.MW, exact.Topology)]
    exact["bin50"] = exact.MW.map(mw_bin_50)
    rows = []
    for key, part in exact.groupby(["analysis_set", "cohort"], dropna=False):
        rows.append({
            "section": "cohort", "analysis_set": key[0], "cohort": key[1],
            "labels": len(part), "graphs": part.Graph_ID.nunique(),
            "connectivity_groups": part.Scaffold_connectivity_group.nunique(),
            "studies": part.Study_key.nunique(),
            "fp06_groups": part.Scaffold_FP06_component.nunique(),
            "murcko": part.MurckoScaffold.nunique(),
            "measurements": int(part.n_distinct_measurements.sum()),
            "database_occurrences": int(part.n_database_occurrences.sum()),
            "stereo_defined_graphs": int(part.loc[part.Unassigned_chiral_centers.eq(0)].Graph_ID.nunique()),
            "stereo_unassigned_graphs": int(part.loc[part.Unassigned_chiral_centers.gt(0)].Graph_ID.nunique()),
            "preliminary_3d_graphs": int(part.loc[part.Unassigned_chiral_centers.eq(0)].Graph_ID.nunique()),
        })
    bins = []
    target = exact[exact.cohort.isin(["nonmacro_500_700", "macro_500_700"])]
    for key, part in target.groupby(["analysis_set", "cohort", "bin50"]):
        bins.append({"section": "bin50", "analysis_set": key[0], "cohort": key[1], "bin50": key[2], "graphs": part.Graph_ID.nunique(), "labels": len(part)})
    sources = []
    for key, part in exact.groupby(["analysis_set", "cohort", "Study_key"]):
        sources.append({"section": "study", "analysis_set": key[0], "cohort": key[1], "study_id": key[2], "graphs": part.Graph_ID.nunique(), "labels": len(part)})
    held = decisions.decision.value_counts().rename_axis("decision").reset_index(name="occurrences")
    held["section"] = "decision"
    return exact, pd.DataFrame(rows), pd.DataFrame(bins), pd.DataFrame(sources), held


def write_report(path, info, cohort, bins, sources, decisions, descriptor_summary):
    main = cohort[(cohort.analysis_set == "caco2_ab_main")] if len(cohort) else cohort
    def cell(frame, cohort_name, column):
        hit = frame[frame.cohort == cohort_name]
        return int(hit[column].sum()) if len(hit) else 0
    lines = []
    lines.append("# Data V1.0 数据整理报告")
    lines.append("")
    lines.append("日期：2026-11-20。本文件是阶段A的正式数据结果，不是训练、划分或三维计算结果。")
    lines.append("")
    lines.append("## 纳入标准")
    lines.append("")
    lines.append("Caco-2 主标签只保留能明确识别的 A→B Papp，数值为精确点值，单位已是 log10(cm/s)。进入主集还要求 stage2 已把来源解析为单一 DOI 或专利号。MolMeDB 中由多条参考文献拼出的 DOI 不进入主集。数据库存储 DOI（10.6019/）不作为原始研究。已知外加抑制剂或诱导的记录单独成集。温度或浓度缺失本身不删除终点、方向、单位、结构和来源已经明确的记录，但缺失条件的记录不会自动并入同一实验组。")
    lines.append("")
    lines.append("聚合只发生在同一分子图、同一原始研究、同一终点、同一方向和同一已确认实验协议内。主标签是组内 log10 值的中位数。跨研究的数值不合并。同一数据库记录被两个库收录时只计一次测量。极差超过 0.5 log 单位时先检查单位错位；已确认属于同一协议的测量保留并报告极差，无法解释的冲突暂缓。0.5 是本项目的操作性核查规则，不是文献里的标准差阈值。")
    lines.append("")
    lines.append("## 已知问题")
    lines.append("")
    lines.append(f"- PD-L1：DOI {PDL1_DOI}、ChEMBL assay 2489294 的 17 条活性记录（activity 26123123–26123139，化合物 CB13–CB31）从 Caco-2:AB 改归 MDCKII-WT A→B。17 是出现记录数。CB30 与 CB31 的 ChEMBL 原值分别是 6.9 和 3.9 × 10^-6 cm/s，属于这 17 条之中，已随整组改归，没有重复补录。")
    lines.append("- CB30/CB31 的 25.8 与 21.7 × 10^-6 cm/s：底表中没有这两条数值。本次没有重新读到原文表格，结构、条件和独立 Caco-2 测量不能同时复核，因此不新增记录。")
    lines.append(f"- Rzepiela（DOI {RZEPIELA_DOI}）：匹配 CSV 3,348 条，官方库与归档库各收录一次，共 {info['rzepiela_occurrences']} 条出现。数值字段已是 log10(Peff)，单位 cm/s，没有再次取对数。其中 {info['sentinel_ids']} 个编号的记录值为 −8；检测限语义仍未由原文闭合，不进入精确回归。归档说明中的 54 条与这次匹配到的 {info['sentinel_ids']} 个编号不一致，没有补造差额记录。")
    lines.append(f"- BAS：120 个未写明立体化学的结构各对应 2 个 BAS 编号。其中 {info['discordant_graphs']} 个结构的两个编号数值不同，相关出现记录暂缓，不把编号当成两个分子，也不把两个数值合成一个标签。{info['concordant_graphs']} 个结构的两个编号数值相同，合并为一次测量，手性按未知处理，可以进入二维标签。")
    lines.append(f"- Kim（DOI {KIM_DOI}）：来源列单位为 nm/s，换算是 log10(数值)−7，与列单位一致；原文实验条件仍未复核，协议标记为未知。方向为 A→B 的精确记录保留在主集，条件缺失的记录不互相合并。B→A 与主集分开。")
    lines.append("- 同研究同名样品冲突沿用 stage3 的 10 组。连接关系冲突和无法确定立体归属的记录暂缓。仅仅是手性中心未指定、且记录能对应到唯一结构的分子，保留在二维数据中，并标为三维资格受限。")
    lines.append("- log10(cm/s) > −3 的记录暂缓。log10(cm/s) < −9 的记录保留原值并标记，不截断。限定符不是等号的记录不进入精确回归。")
    lines.append("")
    lines.append("## 终点")
    lines.append("")
    lines.append(f"PAMPA 按系数分开。辅助任务选择看的是系数是否有明确实验定义、500<MW≤700 非大环是否至少有 30 个图和 5 个原始研究，以及最大来源占比是否不超过 80%。未指定系数的混合 PAMPA 不因为行数多而被选中。本次辅助任务：{info['auxiliary'] or '未指定，各系数保留为探索性任务'}。Caco-2 主集不依赖这个选择。")
    lines.append("")
    lines.append("## 正式规模")
    lines.append("")
    lines.append("下表是精确标签。数据库出现次数、测量记录数、聚合标签数、分子图数和原始研究数分开计数。测量记录数不是独立实验重复次数。")
    lines.append("")
    lines.append("| 集合 | 分组 | 标签 | 分子图 | 连接组 | 原始研究 | 数据库出现 |")
    lines.append("| --- | --- | ---: | ---: | ---: | ---: | ---: |")
    show = cohort.sort_values(["analysis_set", "cohort"]) if len(cohort) else cohort
    for row in show.itertuples(index=False):
        if row.analysis_set in {"caco2_ab_main", "caco2_known_intervention"} or str(row.analysis_set).startswith("pampa_"):
            lines.append(f"| {row.analysis_set} | {row.cohort} | {row.labels} | {row.graphs} | {row.connectivity_groups} | {row.studies} | {row.database_occurrences} |")
    lines.append("")
    lines.append("500<MW≤700 的 50 Da 箱按分子图计数：")
    lines.append("")
    if len(bins):
        lines.append("| 集合 | 分组 | 箱 | 分子图 | 标签 |")
        lines.append("| --- | --- | --- | ---: | ---: |")
        for row in bins.sort_values(["analysis_set", "cohort", "bin50"]).itertuples(index=False):
            if row.analysis_set == "caco2_ab_main" or "auxiliary" in str(row.analysis_set):
                lines.append(f"| {row.analysis_set} | {row.cohort} | {row.bin50} | {row.graphs} | {row.labels} |")
    lines.append("")
    lines.append("## 来源与旧候选数")
    lines.append("")
    lines.append("2026-11-19 的候选数（Caco-2 非大环 500–700 为 600，大环 207，>700 为 88；PAMPA 为 764、1,481、144）是修订前的工作集汇总。当时把组内极差超过 0.5 的记录一律暂缓，并把不同 PAMPA 系数放在一起计数。本版只在无法确认同一实验协议时暂缓冲突，并把 Peff、intrinsic、apparent 和未指明系数分开。因此正式分子图数不要求与那些数字相同。")
    lines.append("")
    focus = sources[(sources.analysis_set == "caco2_ab_main") & (sources.cohort == "nonmacro_500_700")] if len(sources) else sources
    if len(focus):
        total = focus.graphs.sum()
        ranked = focus.sort_values("graphs", ascending=False)
        top = ranked.iloc[0]
        lines.append(f"Caco-2 主集、500<MW≤700 非大环最大原始研究是 {top.study_id}，{int(top.graphs)} 个分子图，占该组按研究计数的分子图合计 {top.graphs / total:.1%}。同一分子图可以出现在多个研究中，因此这个合计可以大于唯一分子图数。")
        lines.append("")
        lines.append("该组按分子图数排列的前五个原始研究：")
        lines.append("")
        for row in ranked.head(5).itertuples(index=False):
            lines.append(f"- {row.study_id}：{int(row.graphs)} 个分子图")
    target = main[main.cohort == "nonmacro_500_700"] if len(main) else main
    if len(target):
        row = target.iloc[0]
        lines.append("")
        lines.append(
            f"该组精确标签 {int(row.labels)} 个，不同测量记录 {int(row.measurements)} 条，数据库出现 {int(row.database_occurrences)} 次，"
            f"分子图 {int(row.graphs)} 个，连接组 {int(row.connectivity_groups)} 个，Murcko 骨架 {int(row.murcko)} 个，"
            f"FP06 分量 {int(row.fp06_groups)} 个，原始研究 {int(row.studies)} 个。连接组、Murcko 骨架和 FP06 分量分别保存。"
        )
    lines.append("")
    lines.append("MW>700 非大环精确标签只标记为外推候选，本阶段不封存、不划分、不抽样。")
    lines.append("")
    lines.append("## 结构与三维资格")
    lines.append("")
    lines.append("未改结构继承 2026-11-18 的 Graph_ID、Connectivity_ID（Scaffold_connectivity_group）、Murcko 和 FP06。FP06 是 Morgan 半径 2、2048 位、不含手性、Tanimoto≥0.6 的单链接分量，并包含骨架/连接约束；当时 RDKit 为 2026.03.3。本次环境的 RDKit 是 2025.09.6，只用于确认缓存与 34,150 个 SMILES 一致，没有重算分组。stage3 将含 Pb 的错误图改指到已有 MZ3 图，没有新增立体异构体。")
    lines.append("")
    lines.append("初步三维资格只表示：结构对应关系明确，且继承的未指定手性中心数为 0。它不是构象生成或力场检查。手性未知的记录仍可留在二维主集。")
    if len(main):
        hit = main[main.cohort == "nonmacro_500_700"]
        if len(hit):
            row = hit.iloc[0]
            lines.append(f"Caco-2 主集 500<MW≤700 非大环中，手性标记完整的分子图 {int(row.stereo_defined_graphs)} 个，手性未指定的 {int(row.stereo_unassigned_graphs)} 个。初步三维资格按手性标记完整计数，为 {int(row.preliminary_3d_graphs)} 个图。")
    lines.append("")
    lines.append("## 描述符")
    lines.append("")
    if descriptor_summary:
        lines.append("下表使用冻结分子表中已有的 MW、cLogP、HBD、HBA，以及标签 y。TPSA 和可旋转键不在冻结表中，留到阶段B。")
        lines.append("")
        lines.append("| 集合 | 分组 | n | MW中位 | cLogP中位 | TPSA | HBD中位 | HBA中位 | RotB | y中位 |")
        lines.append("| --- | --- | ---: | ---: | ---: | --- | ---: | ---: | --- | ---: |")
        for row in descriptor_summary:
            lines.append(f"| {row['analysis_set']} | {row['cohort']} | {row['n']} | {row['MW']:.1f} | {row['cLogP']:.2f} | 未在冻结表 | {row['HBD']:.1f} | {row['HBA']:.1f} | 未在冻结表 | {row['y']:.3f} |")
    else:
        lines.append("冻结分子表没有附带完整二维描述符，TPSA、可旋转键和分布汇总留到阶段B。")
    lines.append("")
    lines.append("## 暂缓与排除")
    lines.append("")
    lines.append("| 决定 | 出现记录 |")
    lines.append("| --- | ---: |")
    counts = decisions.decision.value_counts()
    for key, value in counts.items():
        lines.append(f"| {key} | {int(value)} |")
    lines.append("")
    lines.append(f"输入出现记录 {info['input_rows']} 条，决定表 {info['decision_rows']} 条，两者相同。没有为 CB30/CB31 新增行。")
    lines.append("")
    lines.append("暂缓原因按出现记录计数。一条记录可以同时有多个原因，所以各项之和可以大于暂缓记录数。")
    lines.append("")
    lines.append("| 原因 | 出现记录 |")
    lines.append("| --- | ---: |")
    held = decisions[decisions.decision.ne("IN_EXACT_LABEL")]
    tokens = held.reason_code.fillna("").str.split("|").explode()
    tokens = tokens[tokens.astype(str).str.len().gt(0)]
    for key, value in tokens.value_counts().items():
        lines.append(f"| {key} | {int(value)} |")
    lines.append("")
    lines.append("## 进入阶段B需要的信息")
    lines.append("")
    lines.append("阶段B可以从 `label_table.csv.gz` 的 `analysis_set=caco2_ab_main` 取开发标签，用 `cohort=nonmacro_gt700` 标记外推候选，用 `connectivity_id` 和 `fp06_component` 做隔离。PAMPA 若辅助任务为空，则各 `pampa_*_exploratory` 集合都不要并入 Caco-2。本阶段没有生成折、预算或封存成员。")
    lines.append("")
    lines.append("## 参考文献与本项目规则的差别")
    lines.append("")
    lines.append("Esaki 等，Journal of Cheminformatics，2024，DOI 10.1186/s13321-024-00826-z：本项目沿用其按终点、单位、方向和细胞体系保留已审核记录的做法。Feng 等，Scientific Data，2025，DOI 10.1038/s41597-024-04302-z：本项目参考其对大环、终点类别、单位和限定符的区分；大环阈值使用本项目已有的最大计算环≥12。Wang 等，Journal of Cheminformatics，2025，DOI 10.1186/s13321-025-00947-z：参考重复测量和 RDKit 标准化的必要性；组内中位数是本项目预先规定的标签，不是该文的原句。Pharmaceutics，2022，DOI 10.3390/pharmaceutics14101998：该文讨论测量变异和标准差，本项目的 0.5 是 log 极差核查规则，二者不是同一个阈值。Fourches 等，JCIM，2010，DOI 10.1021/ci100176x：结构和重复记录需要可追溯整理。CycPeptMPDB，JCIM，2023，DOI 10.1021/acs.jcim.2c01573：立体化学、实验条件和检测限要分开保存。Steinbauer 等，JCIM，2024，DOI 10.1021/acs.jcim.4c01083：相同实验条件才适合放在一起比较；本项目据此不把缺失条件自动视为兼容。")
    lines.append("")
    lines.append("## 复核")
    lines.append("")
    lines.append("独立复核脚本是 `scripts/data_v1/verify_data_v1.py`。它重新检查出现记录守恒、17 条不再进入 Caco-2、单位没有二次取对数、不一致的 BAS 副本没有进入精确标签、测量没有按数据库副本重复计数，以及分组标识与分子表一致。结果写在 `verification_result.json`。")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--molecules", type=Path, required=True)
    parser.add_argument("--rzepiela-raw", type=Path, required=True)
    parser.add_argument("--conflicts", type=Path, required=True)
    parser.add_argument("--descriptors", type=Path, default=None)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    labels = pd.read_csv(args.labels, dtype=str, low_memory=False)
    input_rows = len(labels)
    if labels.Occurrence_ID.nunique() != input_rows:
        raise RuntimeError("Occurrence_ID is not unique in the input")
    labels = apply_stage3(labels)
    molecules = pd.read_csv(args.molecules, dtype=str)
    molecules["MW"] = pd.to_numeric(molecules.MW, errors="coerce")
    molecules["Unassigned_chiral_centers"] = pd.to_numeric(molecules.Unassigned_chiral_centers, errors="coerce").fillna(-1).astype(int)
    if molecules.Graph_ID.duplicated().any():
        raise RuntimeError("molecule Graph_ID is not unique")
    rzepiela = load_rzepiela_pairs(args.rzepiela_raw)
    frame, _ = assign_tasks(labels, rzepiela)
    before = len(frame)
    frame = frame.merge(molecules, on="Graph_ID", how="left", validate="many_to_one")
    if len(frame) != before:
        raise RuntimeError("Molecule join changed the row count")
    if frame.MW.isna().any():
        missing = frame.loc[frame.MW.isna(), "Graph_ID"].unique()[:5]
        raise RuntimeError(f"Missing molecule metadata for {missing}")
    conflicts = load_conflicts(args.conflicts)
    conflicts = conflicts.drop_duplicates(["Study_key", "Graph_ID"])
    frame = decide(frame, conflicts, rzepiela["sentinel_ids"])
    if len(frame) != input_rows:
        raise RuntimeError("Decision table row count drifted")
    frame["below_floor"] = frame.y < -9
    frame["protocol_unknown"] = (frame.Source.eq("Kim2025")) | (frame.Protocol.eq("") & frame.Assay_description.eq(""))
    frame["measurement_id"] = [measurement_id(row) for row in frame.itertuples(index=False)]
    frame["condition_id"] = [condition_of(row, measurement) for row, measurement in zip(frame.itertuples(index=False), frame.measurement_id)]
    # First pass uses no auxiliary so PAMPA sets can be counted, then labels are tagged.
    frame["analysis_set"] = "pending"
    exact_mask = ~frame.blocked
    exact = frame.loc[exact_mask].copy()
    labels_out, measurements = aggregate(exact)
    auxiliary, pampa_table = choose_auxiliary(labels_out, molecules)
    auxiliary_coeff = auxiliary.split(":", 1)[1] if auxiliary else ""
    frame["analysis_set"] = [analysis_set(row, auxiliary_coeff) for row in frame.itertuples(index=False)]
    exact = frame.loc[exact_mask].copy()
    labels_out, measurements = aggregate(exact)
    # Propagation of aggregation holds back to occurrences.
    hold_ids = set()
    if len(labels_out):
        held_labels = labels_out[labels_out.label_status.ne("EXACT_LABEL")]
        for record_ids in held_labels.record_ids:
            hold_ids.update(record_ids.split("|"))
    frame["decision"] = np.where(frame.blocked, "HELD", "IN_EXACT_LABEL")
    frame.loc[frame.Occurrence_ID.isin(hold_ids), "decision"] = "HELD_AFTER_AGGREGATION"
    needs_reason = frame.decision.eq("HELD_AFTER_AGGREGATION") & ~frame.reason_code.str.contains("AGGREGATION_CONFLICT", na=False)
    frame.loc[needs_reason, "reason_code"] = frame.loc[needs_reason, "reason_code"].map(lambda text: (text + "|" if text else "") + "AGGREGATION_CONFLICT")
    usable = labels_out[labels_out.label_status.eq("EXACT_LABEL")].copy() if len(labels_out) else labels_out
    usable = usable.merge(
        molecules[["Graph_ID", "SMILES", "Connectivity_SMILES", "MurckoScaffold", "Scaffold_connectivity_group", "Scaffold_FP06_component", "MW", "Topology", "Unassigned_chiral_centers"]],
        on="Graph_ID", how="left", validate="many_to_one",
    )
    usable["connectivity_id"] = usable.Scaffold_connectivity_group
    usable["fp06_component"] = usable.Scaffold_FP06_component
    usable["cohort"] = [mw_cohort(mw, top) for mw, top in zip(usable.MW, usable.Topology)]
    usable["bin50"] = usable.MW.map(mw_bin_50)
    usable["stereo_flag"] = np.where(usable.Unassigned_chiral_centers.eq(0), "stereo_assigned_or_absent", "stereo_unassigned")
    usable["preliminary_3d"] = usable.Unassigned_chiral_centers.eq(0)
    usable["structure_class"] = "unknown"
    usable["extrapolation_candidate"] = usable.cohort.eq("nonmacro_gt700")
    label_columns = [
        "label_id", "Graph_ID", "connectivity_id", "fp06_component", "Study_key", "corrected_task",
        "Direction_refined", "condition_id", "analysis_set", "y", "y_mean", "y_sd", "y_range",
        "n_distinct_measurements", "n_database_occurrences", "record_ids", "sources", "MW", "Topology",
        "cohort", "bin50", "stereo_flag", "preliminary_3d", "structure_class", "intervention",
        "below_floor", "protocol_unknown", "variance_status", "extrapolation_candidate", "SMILES",
        "MurckoScaffold",
    ]
    # Structure annotation uses source libraries already attached to the label, not a new chemical judgement.
    source_text = usable.sources.fillna("")
    usable.loc[source_text.str.contains("PROTAC", case=False), "structure_class"] = "annotated_protac_source"
    peptide = source_text.str.contains("CycPept", case=False)
    usable.loc[peptide & usable.structure_class.eq("unknown"), "structure_class"] = "annotated_cyclic_peptide_database"
    usable.loc[peptide & usable.structure_class.eq("annotated_protac_source"), "structure_class"] = "annotated_protac_and_cyclic_peptide_source"
    usable[label_columns].to_csv(args.out / "label_table.csv.gz", index=False, compression="gzip")
    used_graphs = set(frame.Graph_ID)
    molecule_out = molecules[molecules.Graph_ID.isin(used_graphs)].copy()
    molecule_out["group_inheritance"] = "inherited_20261118_molecule_groups_no_recompute"
    if args.descriptors and args.descriptors.exists():
        desc = pd.read_csv(args.descriptors, dtype=str, usecols=lambda c: c in {"Graph_ID", "cLogP", "HBD", "HBA", "Max_ring"})
        desc = desc.drop_duplicates("Graph_ID")
        for column in ["cLogP", "HBD", "HBA", "Max_ring"]:
            if column in desc.columns:
                desc[column] = pd.to_numeric(desc[column], errors="coerce")
        before_m = len(molecule_out)
        molecule_out = molecule_out.merge(desc, on="Graph_ID", how="left", validate="one_to_one")
        if len(molecule_out) != before_m:
            raise RuntimeError("Descriptor join changed molecule rows")
    molecule_out.to_csv(args.out / "molecule_groups.csv.gz", index=False, compression="gzip")
    decisions = pd.DataFrame({
        "record_id": frame.Occurrence_ID,
        "provenance": frame.Source,
        "original_value": frame.Raw_value,
        "stored_log10_cm_s": frame.y,
        "relation": frame.Relation,
        "original_task": frame.Task,
        "corrected_task": frame.corrected_task,
        "graph_id": frame.Graph_ID,
        "study_id": frame.Study_key,
        "condition_id": frame.condition_id,
        "decision": frame.decision,
        "reason_code": frame.reason_code,
        "evidence_locator": frame.Evidence_locator,
        "transformation": np.where(frame.Source.eq("Kim2025"), "column_unit_nm_s_log10(value)-7_no_second_log", "inherited_log10_cm_s_no_second_log"),
        "measurement_id": frame.measurement_id,
        "analysis_set": frame.analysis_set,
    })
    decisions.to_csv(args.out / "occurrence_decisions.csv.gz", index=False, compression="gzip")
    registry = (
        frame.groupby("Study_key", as_index=False)
        .agg(databases=("Source", lambda s: "|".join(sorted(set(s)))), occurrences=("Occurrence_ID", "size"), graphs=("Graph_ID", "nunique"))
    )
    registry["identifier_type"] = np.where(registry.Study_key.str.startswith("DOI:10.6019/"), "dataset_deposit_doi", np.where(registry.Study_key.str.startswith("DOI:"), "doi", np.where(registry.Study_key.str.startswith("PATENT:"), "patent", "unresolved_or_alias")))
    registry.to_csv(args.out / "source_registry.csv", index=False)
    exact_labels, cohort, bins, study_counts, decision_counts = summarize(labels_out, molecules, decisions)
    summary_parts = [cohort, bins, study_counts, decision_counts]
    if len(pampa_table):
        pampa_table = pampa_table.copy()
        pampa_table["section"] = "pampa_coefficient"
        summary_parts.append(pampa_table)
    summary = pd.concat(summary_parts, ignore_index=True, sort=False)
    summary.to_csv(args.out / "dataset_summary.csv", index=False)
    descriptor_summary = []
    if "cLogP" in molecule_out.columns:
        joined = exact_labels.merge(molecule_out[["Graph_ID", "cLogP", "HBD", "HBA"]], on="Graph_ID", how="left", validate="many_to_one")
        focus = joined[joined.analysis_set.eq("caco2_ab_main") & joined.cohort.isin(["nonmacro_le500", "nonmacro_500_700", "macro_500_700", "nonmacro_gt700"])]
        for (analysis, group), part in focus.groupby(["analysis_set", "cohort"]):
            descriptor_summary.append({
                "analysis_set": analysis, "cohort": group, "n": int(part.Graph_ID.nunique()),
                "MW": float(part.drop_duplicates("Graph_ID").MW.median()),
                "cLogP": float(part.drop_duplicates("Graph_ID").cLogP.median()),
                "HBD": float(part.drop_duplicates("Graph_ID").HBD.median()),
                "HBA": float(part.drop_duplicates("Graph_ID").HBA.median()),
                "y": float(part.y.median()),
            })
    info = {
        "input_rows": input_rows,
        "decision_rows": len(decisions),
        "rzepiela_occurrences": int(frame.rzepiela.sum()),
        "sentinel_ids": rzepiela["n_sentinel_ids"],
        "discordant_graphs": rzepiela["n_discordant_ids"] // 2,
        "concordant_graphs": rzepiela["n_concordant_ids"] // 2,
        "auxiliary": auxiliary or "",
    }
    write_report(args.out / "DATA_PREPARATION_REPORT.md", info, cohort, bins, study_counts, decisions, descriptor_summary)
    manifest = {
        "input_rows": input_rows,
        "decision_rows": len(decisions),
        "exact_labels": int(labels_out.label_status.eq("EXACT_LABEL").sum()) if len(labels_out) else 0,
        "held_after_aggregation": int(labels_out.label_status.ne("EXACT_LABEL").sum()) if len(labels_out) else 0,
        "auxiliary": auxiliary,
        "pdl1_reclassified": 17,
        "rzepiela": {key: rzepiela[key] for key in ["n_rows", "n_shared_smiles", "n_discordant_ids", "n_concordant_ids", "n_sentinel_ids"]},
        "prospective_not_added": [
            {"compound": "CB30", "reported_papp_1e-6_cm_s": 25.8, "reason": "not present in base table; primary table not re-read; existing ChEMBL value 6.9 is inside the 17 MDCKII-WT rows"},
            {"compound": "CB31", "reported_papp_1e-6_cm_s": 21.7, "reason": "not present in base table; primary table not re-read; existing ChEMBL value 3.9 is inside the 17 MDCKII-WT rows"},
        ],
    }
    (args.out / "build_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
