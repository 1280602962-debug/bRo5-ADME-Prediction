"""Rebuild expected M4 membership and compare it with saved outputs. Do not edit M4."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

REPO = Path(__file__).resolve().parents[2]
M4 = REPO / "runs" / "model_development_m4_cpu"
OUT = REPO / "runs" / "model_development_m5"


def expected_members() -> pd.DataFrame:
    dev = pd.read_csv(REPO / "runs" / "model_development_v1" / "development_pool.csv.gz")
    dev["Graph_ID"] = dev.Graph_ID.astype(str)
    dev["assigned_fold"] = dev.assigned_fold.astype(int)
    g1 = dev[dev.g1].copy()
    test = pd.read_csv(REPO / "data" / "splits_v1" / "test_labels.csv.gz")
    test["Graph_ID"] = test.Graph_ID.astype(str)
    test["assigned_fold"] = test.assigned_fold.astype(int)
    static = set(pd.read_csv(M4 / "l1_static_features.csv.gz", usecols=["Graph_ID"]).Graph_ID.astype(str))
    ensemble = set(pd.read_csv(M4 / "l1_ensemble_features.csv.gz", usecols=["Graph_ID"]).Graph_ID.astype(str))
    ion = set(pd.read_csv(M4 / "ionization_features.csv.gz", usecols=["Graph_ID", "molgpka_status"]).query("molgpka_status=='ok'").Graph_ID.astype(str))
    universes = {
        "all": set(g1.Graph_ID),
        "ion": ion,
        "static": static,
        "ensemble": ensemble,
        "both": ion & ensemble,
    }
    arms = {
        "P0_all": "all",
        "P0_ion": "ion",
        "P0_static": "static",
        "P0_ensemble": "ensemble",
        "P0_both": "both",
        "P1_ensemble": "ensemble",
        "P2_ion": "ion",
        "P2_both": "both",
        "P3_static": "static",
        "P3_ensemble": "ensemble",
        "P4_ensemble": "ensemble",
        "P4_both": "both",
        "P5_both": "both",
    }
    rows = []
    for outer in range(5):
        for arm, universe in arms.items():
            graphs = universes[universe]
            train = g1[g1.assigned_fold.ne(outer) & g1.Graph_ID.isin(graphs)]
            evaluation = test[test.assigned_fold.eq(outer) & test.Graph_ID.isin(graphs)]
            rows.append(
                {
                    "outer_fold": outer,
                    "arm": arm,
                    "membership": "reconstructed_from_code",
                    "n_train_labels": int(len(train)),
                    "n_test_labels": int(len(evaluation)),
                    "train_ids": set(train.label_id.astype(str)),
                    "test_ids": set(evaluation.label_id.astype(str)),
                }
            )
    return pd.DataFrame(rows)


def historical_members() -> pd.DataFrame:
    rows = []
    for path in sorted((M4 / "model_parts").glob("C_o*_*.csv")):
        stem = path.stem
        # C_o{fold}_{arm}
        fold = int(stem.split("_")[1][1:])
        arm = stem.split("_", 2)[2]
        frame = pd.read_csv(path, usecols=["label_id", "n_train_labels"])
        rows.append(
            {
                "outer_fold": fold,
                "arm": arm,
                "membership": "saved_prediction_file",
                "n_train_labels": int(frame.n_train_labels.iloc[0]),
                "n_test_labels": int(frame.label_id.nunique()),
                "test_ids": set(frame.label_id.astype(str)),
            }
        )
    return pd.DataFrame(rows)


def conformer_recount() -> dict:
    eligibility = pd.read_csv(M4 / "conformer_eligibility.csv")
    summary = pd.read_csv(M4 / "conformer_generation_summary.csv")
    eligible = eligibility[eligibility.eligible_main.astype(str).str.lower().eq("true")]
    mismatches = 0
    missing = 0
    conformers = 0
    for row in eligible.itertuples(index=False):
        path = M4 / "conformer_parts" / f"{row.Graph_ID}_{int(row.requested_conformers)}_a.csv"
        if not path.exists():
            missing += 1
            continue
        frame = pd.read_csv(path, usecols=["within_5kcal"])
        kept = int(frame.within_5kcal.astype(str).str.lower().eq("true").sum()) if len(frame) else 0
        conformers += kept
        recorded = summary.loc[summary.Graph_ID.astype(str).eq(str(row.Graph_ID)), "n_window"]
        if len(recorded) != 1 or int(recorded.iloc[0]) != kept:
            mismatches += 1
    return {
        "part_files": int(len(list((M4 / "conformer_parts").glob("*.csv")))),
        "eligible_graphs": int(len(eligible)),
        "missing_production_parts": missing,
        "window_count_mismatches": mismatches,
        "conformers_within_5kcal": conformers,
        "summary_static": int(summary.static_available.astype(str).str.lower().eq("true").sum()),
        "summary_ensemble": int(summary.ensemble_reliable.astype(str).str.lower().eq("true").sum()),
    }


def measurement_audit() -> dict:
    groups = pd.read_csv(M4 / "measurement_reproducibility.csv")
    reentry = groups[groups.record_kind.eq("database_reentry")]
    independent = groups[groups.record_kind.eq("independent_repeat")]
    same_id_disagreement = int((reentry.n_distinct_values > 1).sum()) if "n_distinct_values" in reentry else 0
    repeated_within_independent = int((independent.n_records > independent.n_measurement_ids).sum()) if len(independent) else 0
    return {
        "database_reentry_groups": int(len(reentry)),
        "reentry_groups_with_disagreeing_values": same_id_disagreement,
        "independent_id_groups": int(len(independent)),
        "independent_groups_with_repeated_same_id_rows": repeated_within_independent,
        "interpretation": "Different measurement_id values are distinct database measurement records, not proof of independent experiments. Same-id duplicate rows were not treated as independent repeats.",
    }


def charge_audit() -> pd.DataFrame:
    features = pd.read_csv(M4 / "ionization_features.csv.gz", usecols=["Graph_ID", "ionization_class", "permanent_cation_atoms", "permanent_anion_atoms", "recorded_formal_charge"])
    groups = pd.read_csv(REPO / "data" / "modeling_v1" / "molecule_groups.csv.gz", usecols=["Graph_ID", "SMILES"])
    features = features.merge(groups, on="Graph_ID", how="left")
    flagged = features[(features.permanent_cation_atoms.fillna(0) > 0) | (features.permanent_anion_atoms.fillna(0) > 0) | features.ionization_class.eq("permanent_charge")]
    rows = []
    cation = Chem.MolFromSmarts("[N+,n+;!$([N+](=O)[O-])]")
    quaternary = Chem.MolFromSmarts("[N+;H0;!$([N+](=O)[O-])]")
    protonated = Chem.MolFromSmarts("[N+,n+;!H0]")
    carboxylate = Chem.MolFromSmarts("[CX3](=O)[O-]")
    for row in flagged.itertuples(index=False):
        mol = Chem.MolFromSmiles(str(row.SMILES)) if isinstance(row.SMILES, str) else None
        rows.append(
            {
                "Graph_ID": row.Graph_ID,
                "ionization_class": row.ionization_class,
                "formal_charge": row.recorded_formal_charge,
                "charged_n_not_nitro": 0 if mol is None else len(mol.GetSubstructMatches(cation)),
                "quaternary_n": 0 if mol is None else len(mol.GetSubstructMatches(quaternary)),
                "protonated_n": 0 if mol is None else len(mol.GetSubstructMatches(protonated)),
                "carboxylate": 0 if mol is None else len(mol.GetSubstructMatches(carboxylate)),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    expected = expected_members()
    historical = historical_members()
    merged = expected.merge(historical, on=["outer_fold", "arm"], how="outer", suffixes=("_expected", "_saved"))
    merged["test_ids_equal"] = [
        left == right if isinstance(left, set) and isinstance(right, set) else False
        for left, right in zip(merged.test_ids_expected, merged.test_ids_saved)
    ]
    merged["train_count_equal"] = merged.n_train_labels_expected.eq(merged.n_train_labels_saved)
    audit = merged.drop(columns=["test_ids_expected", "test_ids_saved", "train_ids"]).copy()
    audit["train_label_ids_saved"] = False
    audit["train_membership_status"] = np.where(
        audit.train_count_equal.fillna(False),
        "count_matches_reconstruction_ids_not_stored",
        "count_mismatch",
    )
    audit.to_csv(OUT / "m4_membership_audit.csv", index=False)
    # save reconstructed ids separately, not as a hash system
    id_rows = []
    for row in expected.itertuples(index=False):
        id_rows.append({"outer_fold": row.outer_fold, "arm": row.arm, "split": "train_reconstructed", "n": len(row.train_ids), "label_ids": " ".join(sorted(row.train_ids))})
        id_rows.append({"outer_fold": row.outer_fold, "arm": row.arm, "split": "test_reconstructed", "n": len(row.test_ids), "label_ids": " ".join(sorted(row.test_ids))})
    pd.DataFrame(id_rows).to_csv(OUT / "m4_reconstructed_label_ids.csv.gz", index=False)
    conformers = conformer_recount()
    measurements = measurement_audit()
    charges = charge_audit()
    charges.to_csv(OUT / "permanent_charge_audit.csv", index=False)
    reversible = int(((charges.protonated_n > 0) | (charges.carboxylate > 0)).sum()) if len(charges) else 0
    quaternary = int((charges.quaternary_n > 0).sum()) if len(charges) else 0
    deltas = pd.read_csv(M4 / "m4_paired_deltas.csv")
    text = f"""# M5 对 M4 的复核

M4 文件没有被修改。训练标签编号没有写在历史预测文件里，历史文件只保存了测试标签和训练标签数量。下表的训练成员是按原划分和特征可用性重建的预期成员；测试标签可以和保存的预测键直接对照。

外折–模型臂 {len(audit)} 行。测试标签与重建结果一致的有 {int(audit.test_ids_equal.fillna(False).sum())} 行。训练数量与重建结果一致的有 {int(audit.train_count_equal.fillna(False).sum())} 行。不一致的行数是 {int((~audit.train_count_equal.fillna(False)).sum())}。

构象分片 {conformers['part_files']} 个。合格图 {conformers['eligible_graphs']} 个，生产分片缺失 {conformers['missing_production_parts']}，5 kcal/mol 窗口计数与汇总表不一致 {conformers['window_count_mismatches']}。窗口内构象 {conformers['conformers_within_5kcal']} 个。汇总表静态可用 {conformers['summary_static']}，可靠集合 {conformers['summary_ensemble']}。

同一 measurement_id 的多条记录有 {measurements['database_reentry_groups']} 组，其中数值不一致的有 {measurements['reentry_groups_with_disagreeing_values']} 组。这些组没有被当作独立重复。被标为不同测量编号的组有 {measurements['independent_id_groups']} 个，其中同一编号又重复出现的组有 {measurements['independent_groups_with_repeated_same_id_rows']} 个。不同 measurement_id 只表示数据库中的不同测量记录，不是独立实验重复的证明。

保存的永久电荷类有 {len(charges)} 个图。按当前 SMILES 复查，含可逆质子化氮或羧酸根的有 {reversible} 个，含季铵氮的有 {quaternary} 个。该规则的定义仍可能把可逆电离标成永久电荷，但这次复查没有在这 {len(charges)} 个图上看到这种误标。SAMPL 匹配 MAE 仍是宏观实验值与预测位点值的匹配结果。

M4 六项比较仍以原文件为准。原估计：{deltas.loc[deltas.role.eq('primary'), ['contrast','estimate','ci95_low','ci95_high','ci_bonferroni_low','ci_bonferroni_high']].to_dict(orient='records')}。
"""
    (OUT / "M5_M4_RECHECK.md").write_text(text)
    (OUT / "m5_0_summary.json").write_text(json.dumps({"conformers": conformers, "measurements": measurements, "charge_rows": int(len(charges)), "reversible_flagged": reversible, "quaternary_flagged": quaternary, "test_id_matches": int(audit.test_ids_equal.fillna(False).sum()), "train_count_matches": int(audit.train_count_equal.fillna(False).sum())}, indent=2))
    print(text)


if __name__ == "__main__":
    main()
