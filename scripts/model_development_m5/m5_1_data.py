"""Capacity of five candidate sources. Does not change the frozen training data."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import urllib.request
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.inchi import MolToInchiKey

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "runs" / "model_development_m5"
EXTERNAL = OUT / "external"


def connectivity_key(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(str(smiles)) if isinstance(smiles, str) else None
    if mol is None:
        return None
    return MolToInchiKey(mol).split("-")[0]


def max_ring(mol) -> int:
    rings = mol.GetRingInfo().AtomRings()
    return max((len(ring) for ring in rings), default=0)


def describe(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"parse_ok": False}
    ring = max_ring(mol)
    mw = float(Descriptors.MolWt(mol))
    return {
        "parse_ok": True,
        "connectivity_key": MolToInchiKey(mol).split("-")[0],
        "MW": mw,
        "max_ring": int(ring),
        "topology": "MACRO_SCREEN" if ring >= 12 else "NONMACRO_SCREEN",
        "sealed_nonmacro_gt700": bool(ring < 12 and mw > 700),
    }


def index_graphs() -> pd.DataFrame:
    groups = pd.read_csv(REPO / "data" / "modeling_v1" / "molecule_groups.csv.gz", usecols=["Graph_ID", "SMILES", "MW", "Topology"])
    keys = groups.SMILES.map(connectivity_key)
    groups = groups.assign(connectivity_key=keys).dropna(subset=["connectivity_key"])
    return groups.drop_duplicates("connectivity_key")


def pubchem_smiles(name: str) -> str:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{urllib.request.quote(name)}/property/SMILES,MolecularWeight/JSON"
    with urllib.request.urlopen(url, timeout=40) as response:
        payload = json.load(response)
    return str(payload["PropertyTable"]["Properties"][0]["SMILES"])


def inganas_records() -> pd.DataFrame:
    overview = pd.read_excel(EXTERNAL / "inganias_2025_supporting.xlsx", sheet_name="Overview")
    caco = pd.read_excel(EXTERNAL / "inganias_2025_supporting.xlsx", sheet_name="Caco-2", header=1)
    caco = caco.rename(columns={caco.columns[0]: "ID", caco.columns[1]: "raw_papp", caco.columns[2]: "recovery"})
    overview["ID"] = overview.ID.astype(int)
    caco["ID"] = caco.ID.astype(int)
    frame = caco.merge(overview[["ID", "SMILES", "DOI"]], on="ID", how="left")
    rows = []
    for row in frame.itertuples(index=False):
        text = str(row.raw_papp).strip()
        censored = text.startswith("<") or text.startswith(">")
        limit = float(text[1:]) if censored else float(text)
        rows.append(
            {
                "source_id": "10.1021/acs.jmedchem.5c01499",
                "record_id": f"inganias_{int(row.ID)}",
                "endpoint": "Caco-2:AB",
                "inhibitor_status": "efflux_inhibitors_present",
                "raw_value": text,
                "raw_unit": "1e-6 cm/s",
                "relation": "<" if text.startswith("<") else (">" if text.startswith(">") else "="),
                "exact_value": np.nan if censored else limit,
                "log10_cm_s": np.nan if censored else float(np.log10(limit) - 6),
                "recovery": float(row.recovery),
                "SMILES": row.SMILES,
                "literature_locator": "Zenodo 10.5281/zenodo.17708199 Caco-2 sheet",
                "structure_source_doi": row.DOI,
            }
        )
    return pd.DataFrame(rows)


def hudec_records() -> pd.DataFrame:
    frame = pd.read_csv(EXTERNAL / "hudec_combined_protacs.csv").rename(columns={"Article DOI": "article_doi"})
    rows = []
    for row in frame.itertuples(index=False):
        text = str(row.original_pampa).strip()
        censored = "<" in text or ">" in text
        rows.append(
            {
                "source_id": "arxiv:2610.07998",
                "record_id": f"hudec_{int(row.protac_id)}",
                "endpoint": "PAMPA",
                "inhibitor_status": "not_applicable_cell_free",
                "raw_value": text,
                "raw_unit": "as released; logPAMPA column was clipped and is not used",
                "relation": "censored" if censored else "=",
                "exact_value": np.nan if censored else row.PAMPA,
                "log10_cm_s": np.nan,
                "recovery": np.nan,
                "SMILES": row.SMILES,
                "literature_locator": str(row.article_doi),
                "structure_source_doi": str(row.article_doi),
            }
        )
    return pd.DataFrame(rows)


def murakami_records() -> pd.DataFrame:
    rows = []
    for name in ["KT-474", "NX-2127", "MS4078", "dBET6", "dBET57"]:
        rows.append(
            {
                "source_id": "10.1021/jacsau.6c00033",
                "record_id": f"murakami_{name}",
                "endpoint": "Caco-2:AB",
                "inhibitor_status": "not_stated_in_retrieved_main_text",
                "raw_value": np.nan,
                "raw_unit": "not_retrieved",
                "relation": "value_not_retrieved",
                "exact_value": np.nan,
                "log10_cm_s": np.nan,
                "recovery": np.nan,
                "SMILES": pubchem_smiles(name),
                "literature_locator": "JACS Au 2026 main text names five compounds; Table S1 values were not retrieved",
                "structure_source_doi": "PubChem name lookup, not the assay table",
                "compound_name": name,
            }
        )
    return pd.DataFrame(rows)


def kim_row(holdout: set[str]) -> dict:
    occurrences = pd.read_csv(REPO / "data" / "modeling_v1" / "occurrence_decisions.csv.gz", usecols=["study_id", "graph_id", "analysis_set"])
    hit = occurrences[occurrences.study_id.astype(str).str.contains("10.1021/acs.jcim.5c01600", na=False)]
    main = hit[hit.analysis_set.eq("caco2_ab_main")]
    graphs = set(main.graph_id.astype(str))
    return {
        "source_id": "10.1021/acs.jcim.5c01600",
        "endpoint": "Caco-2:AB",
        "raw_records": int(len(main)),
        "qualified_exact_labels": 0,
        "new_graphs": 0,
        "existing_graphs": int(len(graphs)),
        "sealed_graphs": int(len(graphs & holdout)),
        "usable_training_graphs_after_isolation": 0,
        "decision": "already_in_frozen_holdout_not_reextracted",
    }


def attach(records: pd.DataFrame, graphs: pd.DataFrame, holdout: set[str], development: set[str]) -> pd.DataFrame:
    described = records.SMILES.map(describe)
    extra = pd.DataFrame(list(described))
    frame = pd.concat([records.reset_index(drop=True), extra], axis=1)
    lookup = graphs.set_index("connectivity_key")
    frame["existing_graph_id"] = frame.connectivity_key.map(lambda key: lookup.Graph_ID.get(key) if isinstance(key, str) and key in lookup.index else pd.NA)
    frame["already_in_v1"] = frame.existing_graph_id.notna()
    frame["in_development_pool"] = frame.existing_graph_id.astype(str).isin(development)
    frame["in_holdout"] = frame.existing_graph_id.astype(str).isin(holdout)
    frame["new_graph"] = ~frame.already_in_v1 & frame.parse_ok.fillna(False)
    return frame


def capacity_row(frame: pd.DataFrame, source: str, endpoint: str, decision: str) -> dict:
    exact = frame.relation.eq("=") & frame.exact_value.notna()
    new_exact = frame[exact & frame.new_graph]
    developable = new_exact[~new_exact.sealed_nonmacro_gt700.fillna(True)]
    return {
        "source_id": source,
        "endpoint": endpoint,
        "raw_records": int(len(frame)),
        "qualified_exact_labels": int(exact.sum()),
        "censored_or_unretrieved": int((~exact).sum()),
        "new_graphs": int(frame.new_graph.sum()),
        "existing_graphs": int(frame.already_in_v1.sum()),
        "new_exact_graphs": int(new_exact.connectivity_key.nunique()) if len(new_exact) else 0,
        "new_exact_not_sealed": int(developable.connectivity_key.nunique()) if len(developable) else 0,
        "usable_training_graphs_after_isolation": 0,
        "decision": decision,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    graphs = index_graphs()
    holdout = set(pd.read_csv(REPO / "data" / "splits_v1" / "holdout_members.csv").Graph_ID.astype(str))
    development = set(pd.read_csv(REPO / "runs" / "model_development_v1" / "development_pool.csv.gz", usecols=["Graph_ID"]).Graph_ID.astype(str))
    inganas = attach(inganas_records(), graphs, holdout, development)
    hudec = attach(hudec_records(), graphs, holdout, development)
    murakami = attach(murakami_records(), graphs, holdout, development)
    records = pd.concat([inganas, hudec, murakami], ignore_index=True)
    records.to_csv(OUT / "new_experimental_records.csv.gz", index=False)
    overlap = records[
        [
            "source_id",
            "record_id",
            "endpoint",
            "relation",
            "existing_graph_id",
            "already_in_v1",
            "in_development_pool",
            "in_holdout",
            "new_graph",
            "MW",
            "topology",
            "sealed_nonmacro_gt700",
        ]
    ]
    overlap.to_csv(OUT / "overlap_audit.csv", index=False)
    rows = [
        kim_row(holdout),
        capacity_row(inganas, "10.1021/acs.jmedchem.5c01499", "Caco-2:AB_efflux_inhibited", "not_added_condition_differs_and_freeze_unchanged"),
        capacity_row(murakami, "10.1021/jacsau.6c00033", "Caco-2:AB", "values_not_retrieved_freeze_unchanged"),
        capacity_row(hudec, "arxiv:2610.07998", "PAMPA", "different_endpoint_censored_values_not_used_as_exact"),
        {
            "source_id": "10.1021/acs.jcim.6c00597",
            "endpoint": "multiple_public_assays",
            "raw_records": 0,
            "qualified_exact_labels": 0,
            "new_graphs": 0,
            "existing_graphs": 0,
            "usable_training_graphs_after_isolation": 0,
            "decision": "row_level_compilation_not_retrieved_no_multitask",
        },
    ]
    capacity = pd.DataFrame(rows)
    capacity.to_csv(OUT / "endpoint_capacity.csv", index=False)
    domain = (
        records.groupby(["source_id", "endpoint", "topology"], dropna=False)
        .agg(records=("record_id", "size"), median_mw=("MW", "median"), sealed=("sealed_nonmacro_gt700", "sum"), new_graphs=("new_graph", "sum"))
        .reset_index()
    )
    domain.to_csv(OUT / "target_domain_capacity.csv", index=False)
    inventory = pd.DataFrame(
        [
            {"source_id": "PMC12750995", "resolved": "10.1021/acs.jmedchem.5c01499", "local_registry_before": 0, "retrieved": "Zenodo Caco-2 sheet, 10 rows"},
            {"source_id": "PMC12933330", "resolved": "10.1021/jacsau.6c00033", "local_registry_before": 0, "retrieved": "five compound names; assay numbers not in the retrieved text"},
            {"source_id": "10.1021/acs.jcim.5c01600", "resolved": "10.1021/acs.jcim.5c01600", "local_registry_before": 1, "retrieved": "local frozen records only"},
            {"source_id": "10.1021/acs.jcim.6c00597", "resolved": "PubMed 42283410", "local_registry_before": 0, "retrieved": "abstract only; no row file"},
            {"source_id": "arxiv:2610.07998", "resolved": "github xhudec2/protac-permeability", "local_registry_before": 0, "retrieved": "combined_protacs.csv, 87 PAMPA rows"},
        ]
    )
    inventory.to_csv(OUT / "candidate_source_inventory.csv", index=False)
    text = f"""# M5.1 数据增量

冻结训练数据没有修改，也没有做多任务训练。Kim 2025 的 Caco-2 A→B 记录已经在封存集中，没有重新提取标签，也没有计算新特征。

Inganäs 2025 的机器可读表有 10 条 Caco-2 A→B 记录，其中精确值 {int((inganas.relation.eq('=')).sum())} 条，删失值 {int(inganas.relation.eq('<').sum())} 条。测定使用了外排抑制剂，供体 pH 6.5。这不是当前冻结终点的同一实验条件。新图 {int(inganas.new_graph.sum())} 个，已有图 {int(inganas.already_in_v1.sum())} 个。隔离后实际用于训练的图是 0。

Murakami 2026 正文点名 5 个 Caco-2 化合物。补充表中的数值没有取到，因此合格标签按 0 计。结构只用于重叠核对。

Hudec 2026 释放了 87 条 PAMPA 记录，其中原始字段含不等号的有 {int(hudec.relation.eq('censored').sum())} 条。这些不等号没有被改成精确值。PAMPA 不并入 Caco-2。可用 Caco-2 训练图是 0。

Ohlsson 2026 是已有公共测定的多任务模型论文。这次没有取得逐行数据，没有把它当作新的实验监督，也没有训练多任务模型。

结论：没有找到可以在不改冻结规则的前提下加入开发集的新 Caco-2 A→B 监督。
"""
    (OUT / "M5_1_DATA_OPPORTUNITY_REPORT.md").write_text(text)
    print(capacity.to_string(index=False))
    print(text)


if __name__ == "__main__":
    main()
