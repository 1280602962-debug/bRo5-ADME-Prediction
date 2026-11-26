"""Same-fold series and matched-pair direction. Pair rules are fixed before predictions are scored."""

from __future__ import annotations

import json
from collections import defaultdict

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdMMPA

REPO_TEST = None
OUT = None
MISSING = {"NOT_REPORTED", "UNKNOWN", "MENTION_NOT_CONFIRMED", "MENTION_REQUIRES_REVIEW", ""}


def protocol(condition_id: str) -> str:
    parts = str(condition_id).split("|")
    return parts[2] if len(parts) > 2 else ""


def fragments(mol) -> list[str]:
    found = []
    try:
        pieces = rdMMPA.FragmentMol(mol, maxCuts=1, resultsAsMols=False)
    except Exception:
        return found
    for core, side in pieces:
        if core and side:
            found.append(str(core))
    return found


def build_pairs(test: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    work = test.copy()
    work["protocol"] = work.condition_id.map(protocol)
    work = work[~work.protocol.isin(MISSING)].copy()
    counts = work.groupby(["assigned_fold", "Graph_ID"]).label_id.transform("count")
    work = work[counts.eq(1)].copy()
    mols = {}
    fps = {}
    cores = defaultdict(list)
    for row in work.itertuples(index=False):
        mol = Chem.MolFromSmiles(str(row.SMILES))
        if mol is None:
            continue
        mols[row.Graph_ID] = mol
        fps[row.Graph_ID] = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        for core in fragments(mol):
            cores[(int(row.assigned_fold), row.Study_key, row.protocol, core)].append(row.Graph_ID)
    mmp_keys = set()
    mmp_rows = []
    for key, graphs in cores.items():
        unique = sorted(set(graphs))
        if len(unique) < 2:
            continue
        fold, study, proto, _core = key
        for i, left in enumerate(unique):
            for right in unique[i + 1 :]:
                pair = (fold, study, proto, left, right)
                if pair in mmp_keys:
                    continue
                mmp_keys.add(pair)
                mmp_rows.append({"assigned_fold": fold, "Study_key": study, "protocol": proto, "reference_id": left, "modified_id": right, "pair_type": "mmp"})
    sim_rows = []
    for (fold, study, proto), part in work.groupby(["assigned_fold", "Study_key", "protocol"]):
        ids = [gid for gid in part.Graph_ID.unique() if gid in fps]
        scaffolds = part.drop_duplicates("Graph_ID").set_index("Graph_ID").Scaffold_connectivity_group
        masses = part.drop_duplicates("Graph_ID").set_index("Graph_ID").MW
        for i, left in enumerate(ids):
            for right in ids[i + 1 :]:
                if scaffolds.get(left) != scaffolds.get(right):
                    continue
                if abs(float(masses.get(left, np.nan)) - float(masses.get(right, np.nan))) > 80:
                    continue
                if DataStructs.TanimotoSimilarity(fps[left], fps[right]) < 0.6:
                    continue
                ordered = tuple(sorted((left, right)))
                if (fold, study, proto, ordered[0], ordered[1]) in mmp_keys:
                    continue
                sim_rows.append({"assigned_fold": int(fold), "Study_key": study, "protocol": proto, "reference_id": ordered[0], "modified_id": ordered[1], "pair_type": "similarity"})
    return pd.DataFrame(mmp_rows), pd.DataFrame(sim_rows)


def attach_labels(pairs: pd.DataFrame, test: pd.DataFrame) -> pd.DataFrame:
    lookup = test.drop_duplicates(["assigned_fold", "Graph_ID"]).set_index(["assigned_fold", "Graph_ID"])
    rows = []
    for row in pairs.itertuples(index=False):
        try:
            left = lookup.loc[(row.assigned_fold, row.reference_id)]
            right = lookup.loc[(row.assigned_fold, row.modified_id)]
        except KeyError:
            continue
        rows.append({**row._asdict(), "y_reference": float(left.y), "y_modified": float(right.y), "label_reference": left.label_id, "label_modified": right.label_id})
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame["delta_y"] = frame.y_modified - frame.y_reference
    return frame


def score(frame: pd.DataFrame, predictions: pd.DataFrame, name: str) -> dict:
    if frame.empty:
        return {"model": name, "pair_type": None, "status": "insufficient", "n_pairs": 0, "n_studies": 0}
    merged = frame.merge(predictions.rename(columns={"label_id": "label_reference", "y_pred": "pred_reference"}), on="label_reference", how="inner")
    merged = merged.merge(predictions.rename(columns={"label_id": "label_modified", "y_pred": "pred_modified"}), on="label_modified", how="inner")
    if merged.empty:
        return {"model": name, "pair_type": frame.pair_type.iloc[0], "status": "insufficient", "n_pairs": 0, "n_studies": 0}
    merged["delta_pred"] = merged.pred_modified - merged.pred_reference
    merged["abs_delta_error"] = (merged.delta_pred - merged.delta_y).abs()
    merged["pred_zero"] = merged.delta_pred.abs() <= 1e-8
    study_rows = []
    for study, part in merged.groupby("Study_key"):
        usable = part[~part.pred_zero]
        direction = usable[usable.delta_y.abs() >= 0.3]
        large = usable[usable.delta_y.abs() >= 1.0]
        study_rows.append(
            {
                "Study_key": study,
                "n_pairs": int(len(part)),
                "mae": float(part.abs_delta_error.mean()),
                "zero_baseline_mae": float(part.delta_y.abs().mean()),
                "direction_n": int(len(direction)),
                "direction_correct": int((np.sign(direction.delta_pred) == np.sign(direction.delta_y)).sum()) if len(direction) else 0,
                "large_n": int(len(large)),
                "large_correct": int((np.sign(large.delta_pred) == np.sign(large.delta_y)).sum()) if len(large) else 0,
                "n_pred_zero": int(part.pred_zero.sum()),
            }
        )
    studies = pd.DataFrame(study_rows)
    n_pairs = int(len(merged))
    n_studies = int(studies.Study_key.nunique())
    status = "estimated" if n_pairs >= 20 and n_studies >= 5 else "insufficient"
    direction_n = int(studies.direction_n.sum())
    direction_studies = studies[studies.direction_n > 0]
    large_studies = studies[studies.large_n > 0]
    return {
        "model": name,
        "pair_type": frame.pair_type.iloc[0],
        "status": status,
        "n_pairs": n_pairs,
        "n_studies": n_studies,
        "delta_mae_study_mean": float(studies.mae.mean()),
        "zero_change_baseline_mae": float(studies.zero_baseline_mae.mean()),
        "direction_accuracy_0.3": float((direction_studies.direction_correct / direction_studies.direction_n).mean()) if len(direction_studies) else np.nan,
        "direction_pairs_0.3": direction_n,
        "direction_studies_0.3": int(len(direction_studies)),
        "direction_accuracy_1.0": float((large_studies.large_correct / large_studies.large_n).mean()) if len(large_studies) else np.nan,
        "direction_pairs_1.0": int(studies.large_n.sum()),
        "predicted_zero_pairs": int(studies.n_pred_zero.sum()),
        "random_direction_reference": 0.5,
        "application_thresholds_are_not_experimental_error_bounds": True,
    }


def main() -> None:
    from pathlib import Path

    repo = Path(__file__).resolve().parents[2]
    out = repo / "runs" / "model_development_m5"
    out.mkdir(parents=True, exist_ok=True)
    test = pd.read_csv(repo / "data" / "splits_v1" / "test_labels.csv.gz")
    test["Graph_ID"] = test.Graph_ID.astype(str)
    test["assigned_fold"] = test.assigned_fold.astype(int)
    mmp, sim = build_pairs(test)
    mmp = attach_labels(mmp, test)
    sim = attach_labels(sim, test)
    pairs = pd.concat([mmp, sim], ignore_index=True) if len(mmp) or len(sim) else pd.DataFrame()
    pairs.to_csv(out / "series_pairs.csv.gz", index=False)
    m3 = pd.read_csv(repo / "runs" / "model_development_v1" / "oof_predictions.csv.gz")
    m4 = pd.read_csv(repo / "runs" / "model_development_m4_cpu" / "m4_oof_predictions.csv.gz")
    sources = {
        "M3_RF_F05": m3[(m3.evaluation.eq("C")) & (m3.model.eq("RF")) & (m3.feature.eq("F05")) & (m3.parameter_source.eq("fixed_selected"))][["label_id", "y_pred"]],
        "M3_RF_F04": m3[(m3.evaluation.eq("C")) & (m3.model.eq("RF")) & (m3.feature.eq("F04")) & (m3.parameter_source.eq("fixed_g1_rf"))][["label_id", "y_pred"]],
        "M4_P4": m4[(m4.evaluation.eq("C")) & (m4.arm.eq("P4_ensemble")) & (m4.model.eq("RF"))][["label_id", "y_pred"]],
    }
    rows = []
    for pair_frame, pair_type in ((mmp, "mmp"), (sim, "similarity")):
        if pair_frame is None or pair_frame.empty:
            for name in sources:
                rows.append({"model": name, "pair_type": pair_type, "status": "insufficient", "n_pairs": 0, "n_studies": 0, "reason": "no same-fold pair had a compatible non-missing protocol"})
            continue
        for name, pred in sources.items():
            rows.append(score(pair_frame, pred, name))
    metrics = pd.DataFrame(rows)
    metrics.to_csv(out / "series_ranking_metrics.csv", index=False)
    excluded_studies = 0
    excluded_pairs = 0
    for (_, _), part in test.groupby(["assigned_fold", "Study_key"]):
        n_graphs = int(part.Graph_ID.nunique())
        if n_graphs >= 2:
            excluded_studies += 1
            excluded_pairs += n_graphs * (n_graphs - 1) // 2
    (out / "m5_2b_summary.json").write_text(json.dumps({
        "mmp_pairs": int(len(mmp)),
        "similarity_pairs": int(len(sim)),
        "informative_protocol_test_labels": int((~test.condition_id.map(protocol).isin(MISSING)).sum()),
        "same_fold_same_study_graph_pairs_excluded_by_protocol": excluded_pairs,
        "same_fold_same_study_studies_excluded_by_protocol": excluded_studies,
        "conditions_not_relaxed": True,
    }, indent=2))
    print(metrics.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
