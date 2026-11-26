"""M4.0 measurement heterogeneity for the frozen Caco-2 A→B development pool.

condition_id is STRICT|source_token|protocol_token. The source token is often a
ChEMBL assay id, not the literature study. Cross-study counts are therefore
reported both with the full condition_id and with the protocol token alone.
Identical missing tokens are not treated as the same experimental condition.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "runs" / "model_development_m4_cpu"
MISSING_PROTOCOL = {"NOT_REPORTED", "UNKNOWN", "MENTION_NOT_CONFIRMED", "MENTION_REQUIRES_REVIEW", ""}


def pairwise_abs(values: np.ndarray) -> float:
    if len(values) < 2:
        return np.nan
    diffs = [abs(float(values[i] - values[j])) for i in range(len(values)) for j in range(i + 1, len(values))]
    return float(np.mean(diffs))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    dev = pd.read_csv(REPO / "runs" / "model_development_v1" / "development_pool.csv.gz")
    holdout = set(pd.read_csv(REPO / "data" / "splits_v1" / "holdout_members.csv").Graph_ID.astype(str))
    if holdout & set(dev.Graph_ID.astype(str)):
        raise RuntimeError("holdout graph entered the development pool")
    graphs = set(dev.Graph_ID.astype(str))
    meta = dev.drop_duplicates("Graph_ID")[["Graph_ID", "cohort", "MW"]].copy()
    meta["Graph_ID"] = meta.Graph_ID.astype(str)
    occ = pd.read_csv(REPO / "data" / "modeling_v1" / "occurrence_decisions.csv.gz", low_memory=False)
    main = occ[occ.analysis_set.eq("caco2_ab_main") & occ.corrected_task.eq("Caco-2:AB")].copy()
    main["graph_id"] = main.graph_id.astype(str)
    main = main[main.graph_id.isin(graphs)].copy()
    parts = main.condition_id.fillna("").str.split("|", expand=True)
    main["strict_token"] = parts[0]
    main["source_token"] = parts[1].fillna("")
    main["protocol_token"] = parts[2].fillna("")
    main["source_kind"] = np.where(main.source_token.str.startswith("CHEMBL_ASSAY"), "chembl_assay_id", np.where(main.source_token.eq("UNSTATED"), "unstated", np.where(main.source_token.eq("STATED_PROTOCOL"), "stated_protocol", "other")))
    main["protocol_missing"] = main.protocol_token.isin(MISSING_PROTOCOL)
    main["y"] = main.stored_log10_cm_s.astype(float)
    group_cols = ["graph_id", "study_id", "condition_id"]
    rows = []
    for keys, part in main.groupby(group_cols, dropna=False):
        graph_id, study_id, condition_id = keys
        values = part.y.to_numpy(np.float64)
        n_meas = int(part.measurement_id.nunique())
        same_value = int(part.y.nunique()) == 1
        if n_meas == 1 and len(part) > 1:
            kind = "database_reentry"
        elif n_meas > 1 and same_value:
            kind = "independent_ids_same_value"
        elif n_meas > 1:
            kind = "independent_repeat"
        else:
            kind = "single_record"
        rows.append(
            {
                "graph_id": graph_id,
                "study_id": study_id,
                "condition_id": condition_id,
                "source_token": part.source_token.iloc[0],
                "source_kind": part.source_kind.iloc[0],
                "protocol_token": part.protocol_token.iloc[0],
                "protocol_missing": bool(part.protocol_missing.iloc[0]),
                "n_records": int(len(part)),
                "n_measurement_ids": n_meas,
                "n_provenances": int(part.provenance.nunique()),
                "n_distinct_values": int(part.y.nunique()),
                "record_kind": kind,
                "y_mean": float(np.mean(values)),
                "y_sd": float(np.std(values, ddof=1)) if len(values) > 1 and part.y.nunique() > 1 else (0.0 if len(values) > 1 else np.nan),
                "y_range": float(np.max(values) - np.min(values)),
                "pairwise_abs": pairwise_abs(np.unique(values) if kind == "database_reentry" else values),
            }
        )
    groups = pd.DataFrame(rows).merge(meta, left_on="graph_id", right_on="Graph_ID", how="left")
    groups.to_csv(OUT / "measurement_reproducibility.csv", index=False)

    informative = main[~main.protocol_missing]
    coverage_rows = [
        {"item": "development_graphs", "n": int(dev.Graph_ID.nunique())},
        {"item": "development_labels", "n": int(dev.label_id.nunique())},
        {"item": "caco2_ab_main_records_on_development_graphs", "n": int(len(main))},
        {"item": "graphs_with_a_main_record", "n": int(main.graph_id.nunique())},
        {"item": "studies", "n": int(main.study_id.nunique())},
        {"item": "full_condition_id_values", "n": int(main.condition_id.nunique())},
        {"item": "source_token_is_chembl_assay", "n": int(main.source_kind.eq("chembl_assay_id").sum())},
        {"item": "source_token_unstated", "n": int(main.source_kind.eq("unstated").sum())},
        {"item": "protocol_token_missing", "n": int(main.protocol_missing.sum())},
        {"item": "protocol_token_informative", "n": int((~main.protocol_missing).sum())},
        {"item": "database_reentry_groups", "n": int(groups.record_kind.eq("database_reentry").sum())},
        {"item": "independent_repeat_groups", "n": int(groups.record_kind.eq("independent_repeat").sum())},
        {"item": "independent_ids_same_value_groups", "n": int(groups.record_kind.eq("independent_ids_same_value").sum())},
    ]
    full_cross = main.groupby(["graph_id", "condition_id"]).study_id.nunique()
    protocol_cross = informative.groupby(["graph_id", "protocol_token"]).study_id.nunique() if len(informative) else pd.Series(dtype=int)
    missing_cross = main[main.protocol_missing].groupby(["graph_id", "protocol_token"]).study_id.nunique()
    coverage_rows.extend(
        [
            {"item": "cross_study_groups_full_condition_id", "n": int((full_cross > 1).sum())},
            {"item": "cross_study_groups_informative_protocol_token", "n": int((protocol_cross > 1).sum()) if len(protocol_cross) else 0},
            {"item": "cross_study_groups_identical_missing_token_not_used", "n": int((missing_cross > 1).sum())},
            {"item": "graphs_in_more_than_one_study", "n": int((main.groupby("graph_id").study_id.nunique() > 1).sum())},
            {"item": "graphs_with_more_than_one_condition_id", "n": int((main.groupby("graph_id").condition_id.nunique() > 1).sum())},
        ]
    )
    work = pd.read_csv(REPO / "data" / "curated" / "working_graph_study_labels.csv.gz")
    work = work[work.Task.eq("Caco-2:AB") & work.Graph_ID.astype(str).isin(graphs)]
    coverage_rows.extend(
        [
            {"item": "curated_rows_ph_unknown", "n": int(work.pH_reported.astype(str).eq("UNKNOWN").sum())},
            {"item": "curated_rows_ph_numeric", "n": int((~work.pH_reported.astype(str).eq("UNKNOWN")).sum())},
            {"item": "curated_rows_inhibitor_not_reported", "n": int(work.Inhibitor_status.astype(str).str.contains("NOT_REPORTED").sum())},
            {"item": "curated_rows", "n": int(len(work))},
        ]
    )
    pd.DataFrame(coverage_rows).to_csv(OUT / "assay_condition_coverage.csv", index=False)

    summary = []
    independent = groups[groups.record_kind.eq("independent_repeat")]
    for cohort, part in groups.groupby(groups.cohort.fillna("unassigned")):
        ind = part[part.record_kind.eq("independent_repeat")]
        summary.append(
            {
                "cohort": cohort,
                "graphs": int(part.graph_id.nunique()),
                "records": int(part.n_records.sum()),
                "independent_repeat_groups": int(len(ind)),
                "database_reentry_groups": int(part.record_kind.eq("database_reentry").sum()),
                "independent_sd_median": float(ind.y_sd.median()) if len(ind) else np.nan,
                "independent_pairwise_abs_median": float(ind.pairwise_abs.median()) if len(ind) else np.nan,
                "independent_range_median": float(ind.y_range.median()) if len(ind) else np.nan,
            }
        )
    summary.append(
        {
            "cohort": "ALL_DEVELOPMENT",
            "graphs": int(groups.graph_id.nunique()),
            "records": int(groups.n_records.sum()),
            "independent_repeat_groups": int(len(independent)),
            "database_reentry_groups": int(groups.record_kind.eq("database_reentry").sum()),
            "independent_sd_median": float(independent.y_sd.median()) if len(independent) else np.nan,
            "independent_pairwise_abs_median": float(independent.pairwise_abs.median()) if len(independent) else np.nan,
            "independent_range_median": float(independent.y_range.median()) if len(independent) else np.nan,
        }
    )
    pd.DataFrame(summary).to_csv(OUT / "measurement_variability_summary.csv", index=False)

    oof = pd.read_csv(REPO / "runs" / "model_development_v1" / "oof_predictions.csv.gz")
    oof = oof[oof.evaluation.eq("C") & oof.parameter_source.isin(["fixed_g1_rf", "fixed_selected"])]
    oof["abs_error"] = (oof.y - oof.y_pred).abs()
    link_rows = []
    for row in independent.itertuples(index=False):
        for source, feature in [("fixed_g1_rf", "F04"), ("fixed_selected", "F05")]:
            pred = oof[oof.Graph_ID.astype(str).eq(row.graph_id) & oof.parameter_source.eq(source) & oof.feature.eq(feature)]
            link_rows.append(
                {
                    "graph_id": row.graph_id,
                    "cohort": row.cohort,
                    "repeat_sd": row.y_sd,
                    "repeat_range": row.y_range,
                    "model": f"RF-{feature}",
                    "n_test_labels": int(len(pred)),
                    "mean_abs_error": float(pred.abs_error.mean()) if len(pred) else np.nan,
                }
            )
    link = pd.DataFrame(link_rows)
    link.to_csv(OUT / "repeat_error_link.csv", index=False)
    note = {
        "independent_repeat_groups": int(len(independent)),
        "independent_groups_with_eval_c_prediction": int(link.dropna(subset=["mean_abs_error"]).graph_id.nunique()) if len(link) else 0,
        "variance_components": "not_fit",
        "reason": "Cross-study compatible repeats are absent after removing missing-condition tokens, and independent within-study repeats are only a handful of groups.",
    }
    (OUT / "m4_0_summary.json").write_text(json.dumps(note, indent=2))
    med_sd = float(independent.y_sd.median()) if len(independent) else float("nan")
    text = f"""# M4.0 实验测量异质性

终点保持冻结的 Caco-2 A→B log10(Papp [cm/s])。统计只使用开发池中的 2511 个图，不包含 106 个封存图。没有改写冻结标签。

`condition_id` 的结构是 `STRICT|source_token|protocol_token`。在开发池的主终点记录中，source token 多数是 `CHEMBL_ASSAY:` 实验编号，不是文献 DOI。研究来源在单独的 `study_id` 字段。若只用完整 `condition_id` 寻找跨研究重复，assay 编号会把不同研究拆开。因此另外按去掉 assay 编号后的 protocol token 计数。`NOT_REPORTED`、`UNKNOWN` 和其他缺失标记单独统计；相同缺失标记不视为实验条件相同。

开发池主终点原始记录 {len(main)} 条，覆盖 {main.graph_id.nunique()} 个图和 {main.study_id.nunique()} 个研究。同一分子、同一研究、同一 condition_id 下，有 {int(groups.record_kind.eq('database_reentry').sum())} 组是同一 `measurement_id` 的多条数据库记录，记为重复收录，不作为独立测量。测量编号不同且数值不同的组有 {len(independent)} 个。测量编号不同但数值相同的组有 {int(groups.record_kind.eq('independent_ids_same_value').sum())} 个。

按完整 condition_id，跨研究组有 {int((full_cross > 1).sum())} 个。去掉缺失 protocol 后，按 protocol token 的跨研究组有 {int((protocol_cross > 1).sum()) if len(protocol_cross) else 0} 个。缺失标记本身如果被错误当成相同条件，会得到 {int((missing_cross > 1).sum())} 个跨研究组，这些组没有使用。有 {int((main.groupby('graph_id').study_id.nunique() > 1).sum())} 个分子出现在一个以上的研究中，但没有足够的已知兼容条件把它们配成重复测量。

因此不能充分量化研究间变异，没有做方差分解或混合效应模型。也不能把实验异质性写成已经证实的模型失败原因。独立重复组的标准差中位数是 {med_sd:.3f}。这些组里能连到评价 C 的 RF-F04/F05 预测的分子数记录在 `repeat_error_link.csv`。数量很少，只作描述，不解释为因果，也不写成全体分子的预测误差下限。

整理表中的 pH 与抑制剂注释只说明注释覆盖率。数值 pH 很少，绝大多数为 UNKNOWN。pH 6.5 和 7.4 不是这些原始实验已被证实的统一 pH。
"""
    (OUT / "M4_0_DATA_VARIABILITY_REPORT.md").write_text(text)
    print(json.dumps(note, indent=2))
    print("independent", len(independent), "reentry", int(groups.record_kind.eq("database_reentry").sum()))


if __name__ == "__main__":
    main()
