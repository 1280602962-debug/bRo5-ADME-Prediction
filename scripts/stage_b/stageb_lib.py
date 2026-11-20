"""Stage B split, feature, and sampling library.

Group semantics, fixed before any model is fit:
- Graph_ID: standardized molecular graph.
- Connectivity_SMILES: full molecular connectivity.
- connectivity_id / Scaffold_connectivity_group: SC group from connectivity and Murcko.
- Scaffold_FP06_component: inherited fingerprint component. Not recomputed.
Budget counts unique Graph_ID values. One Connectivity_SMILES occupies one budget
slot, represented by the lexicographically smallest Graph_ID. Distinct graphs in
the same SC may appear together in one training set.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

SEEDS = [17, 29, 42, 71, 101, 137, 271, 314, 577, 911]
BINS = ["500_550", "550_600", "600_650", "650_700"]
BUDGETS = [(500, 50), (500, 25), (250, 50), (250, 25)]
DESC_NAMES = [
    "MW",
    "HeavyAtomCount",
    "MolLogP",
    "TPSA",
    "HBD",
    "HBA",
    "NumRotatableBonds",
    "FractionCSP3",
    "RingCount",
    "NumAromaticRings",
    "FormalCharge",
    "MaxRingSize",
    "TPSA_over_MW",
]
ARMS = {
    "B": ["B"],
    "BS": ["B", "S"],
    "BM": ["B", "M"],
    "BN": ["B", "N"],
}


def derive_seed(*parts: object) -> int:
    digest = hashlib.sha256("||".join(map(str, parts)).encode()).digest()
    return int.from_bytes(digest[:4], "little")


def tie_key(study: str) -> str:
    return hashlib.sha256(f"20261119{study}".encode()).hexdigest()


def largest_remainder(total: int, caps: list[int]) -> list[int] | None:
    caps = [int(c) for c in caps]
    if total < 0 or sum(caps) < total:
        return None
    if total == 0:
        return [0 for _ in caps]
    weight = sum(caps)
    base = [(total * c) // weight for c in caps]
    left = total - sum(base)
    rems = [(total * c) % weight for c in caps]
    order = sorted(range(len(caps)), key=lambda i: (-rems[i], i))
    for i in order:
        if left == 0:
            break
        if base[i] < caps[i]:
            base[i] += 1
            left -= 1
    if left != 0:
        return None
    return base


def assert_quota_examples() -> None:
    if largest_remainder(50, [10, 10, 10, 30]) != [9, 8, 8, 25]:
        raise RuntimeError("largest remainder tie example failed")
    if largest_remainder(5, [2, 2, 2, 2]) != [2, 1, 1, 1]:
        raise RuntimeError("largest remainder unit example failed")


def sample_ids(candidates: list[str], n: int, *seed_parts: object) -> list[str]:
    ordered = sorted(set(candidates))
    if n > len(ordered):
        raise RuntimeError(f"cannot sample {n} from {len(ordered)}")
    if n == 0:
        return []
    rng = np.random.default_rng(derive_seed(*seed_parts))
    picked = rng.choice(len(ordered), size=n, replace=False)
    return [ordered[int(i)] for i in picked]


def load_main(repo: Path) -> pd.DataFrame:
    labels = pd.read_csv(repo / "data/modeling_v1/label_table.csv.gz")
    labels = labels[labels.analysis_set.eq("caco2_ab_main")].copy()
    if len(labels) != 3224 or labels.Graph_ID.nunique() != 2916:
        raise RuntimeError("main-set census does not match Data V1.0")
    molecules = pd.read_csv(repo / "data/modeling_v1/molecule_groups.csv.gz", dtype=str)
    keep = [
        "Graph_ID",
        "SMILES",
        "Connectivity_SMILES",
        "Scaffold_connectivity_group",
        "Scaffold_FP06_component",
    ]
    frame = labels.merge(molecules[keep], on="Graph_ID", how="left", validate="many_to_one", suffixes=("", "_mol"))
    if frame.Connectivity_SMILES.isna().any():
        raise RuntimeError("labels are missing molecule connectivity")
    if not frame.connectivity_id.eq(frame.Scaffold_connectivity_group).all():
        raise RuntimeError("connectivity_id is not the inherited SC group")
    if not frame.fp06_component.eq(frame.Scaffold_FP06_component).all():
        raise RuntimeError("fp06_component is not the inherited FP06 group")
    frame["MW"] = pd.to_numeric(frame.MW, errors="raise")
    frame["y"] = pd.to_numeric(frame.y, errors="raise")
    frame["smiles_std"] = frame.SMILES_mol.fillna(frame.SMILES)
    return frame


def bin_of(mw: float) -> str:
    if 500 < mw <= 550:
        return "500_550"
    if 550 < mw <= 600:
        return "550_600"
    if 600 < mw <= 650:
        return "600_650"
    if 650 < mw <= 700:
        return "650_700"
    return ""


def compute_descriptor_row(mol) -> list[float]:
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

    values = {
        "MW": Descriptors.MolWt(mol),
        "HeavyAtomCount": float(Lipinski.HeavyAtomCount(mol)),
        "MolLogP": Crippen.MolLogP(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "HBD": float(rdMolDescriptors.CalcNumHBD(mol)),
        "HBA": float(rdMolDescriptors.CalcNumHBA(mol)),
        "NumRotatableBonds": float(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        "FractionCSP3": rdMolDescriptors.CalcFractionCSP3(mol),
        "RingCount": float(rdMolDescriptors.CalcNumRings(mol)),
        "NumAromaticRings": float(rdMolDescriptors.CalcNumAromaticRings(mol)),
        "FormalCharge": float(Chem.GetFormalCharge(mol)),
    }
    rings = mol.GetRingInfo().AtomRings()
    values["MaxRingSize"] = float(max((len(ring) for ring in rings), default=0))
    values["TPSA_over_MW"] = values["TPSA"] / values["MW"] if values["MW"] else np.nan
    return [float(values[name]) for name in DESC_NAMES]


def featurize(graph_smiles: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem

    RDLogger.DisableLog("rdApp.warning")

    ids: list[str] = []
    desc_rows: list[list[float]] = []
    fps: list[np.ndarray] = []
    failed: list[str] = []
    for row in graph_smiles.itertuples(index=False):
        mol = Chem.MolFromSmiles(str(row.smiles_std))
        if mol is None:
            failed.append(row.Graph_ID)
            continue
        try:
            values = compute_descriptor_row(mol)
            bitvect = AllChem.GetMorganFingerprintAsBitVect(
                mol, 2, nBits=2048, useChirality=False, useFeatures=False
            )
            array = np.zeros((2048,), dtype=np.uint8)
            AllChem.DataStructs.ConvertToNumpyArray(bitvect, array)
        except Exception:
            failed.append(row.Graph_ID)
            continue
        if not np.isfinite(values).any():
            failed.append(row.Graph_ID)
            continue
        ids.append(row.Graph_ID)
        desc_rows.append(values)
        fps.append(array)
    desc = pd.DataFrame(desc_rows, columns=DESC_NAMES)
    desc.insert(0, "Graph_ID", ids)
    fp = np.vstack(fps) if fps else np.zeros((0, 2048), dtype=np.uint8)
    return desc, fp, failed


def assign_studies(study_rows: list[dict], n_folds: int) -> list[list[str]]:
    folds = [[] for _ in range(n_folds)]
    target_sets = [set() for _ in range(n_folds)]
    graph_sets = [set() for _ in range(n_folds)]
    targeted = [row for row in study_rows if row["n_target"] > 0]
    targeted.sort(key=lambda row: (-row["n_target"], tie_key(row["study"])))
    for row in targeted:
        fold = min(range(n_folds), key=lambda i: (len(target_sets[i]), i))
        folds[fold].append(row["study"])
        target_sets[fold].update(row["target_graphs"])
        graph_sets[fold].update(row["graphs"])
    plain = [row for row in study_rows if row["n_target"] == 0]
    plain.sort(key=lambda row: (-len(row["graphs"]), tie_key(row["study"])))
    for row in plain:
        fold = min(range(n_folds), key=lambda i: (len(graph_sets[i]), i))
        folds[fold].append(row["study"])
        graph_sets[fold].update(row["graphs"])
    return folds


def representatives(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame
    ordered = frame.sort_values(["Connectivity_SMILES", "Graph_ID"])
    return ordered.groupby("Connectivity_SMILES", as_index=False).first()


def pool_caps(train: pd.DataFrame) -> dict:
    small = representatives(train[train.cohort.eq("nonmacro_le500")])
    macro = representatives(train[train.cohort.eq("macro_500_700")])
    donor = representatives(train[train.cohort.eq("nonmacro_500_700")])
    macro_bins = [int((macro["bin50"] == name).sum()) if len(macro) else 0 for name in BINS]
    donor_bins = [int((donor["bin50"] == name).sum()) if len(donor) else 0 for name in BINS]
    common = [min(a, b) for a, b in zip(macro_bins, donor_bins)]
    return {
        "n_small": int(len(small)),
        "macro_bins": macro_bins,
        "donor_bins": donor_bins,
        "common_bins": common,
        "common_total": int(sum(common)),
    }


def budget_fits(caps: list[dict], b: int, k: int) -> bool:
    for cap in caps:
        if cap["n_small"] < b + k or cap["common_total"] < k:
            return False
        if largest_remainder(k, cap["common_bins"]) is None:
            return False
    return True


def build_fold_state(dev: pd.DataFrame, n_folds: int) -> dict | None:
    study_rows = []
    for study, part in dev.groupby("Study_key"):
        target = set(part.loc[part.cohort.eq("nonmacro_500_700"), "Graph_ID"])
        study_rows.append(
            {
                "study": study,
                "n_target": len(target),
                "target_graphs": target,
                "graphs": set(part.Graph_ID),
            }
        )
    assignment = assign_studies(study_rows, n_folds)
    study_fold = {study: fold for fold, studies in enumerate(assignment) for study in studies}
    dev = dev.copy()
    dev["assigned_fold"] = dev.Study_key.map(study_fold)
    held = set()
    for sc, part in dev.groupby("Scaffold_connectivity_group"):
        present = set(int(v) for v in part.assigned_fold)
        if len(present) > 1:
            held.add(sc)
    dev["held_sc"] = dev.Scaffold_connectivity_group.isin(held)
    test_mask = dev.cohort.eq("nonmacro_500_700") & ~dev.held_sc
    test = dev.loc[test_mask].copy()
    test_studies = {
        fold: set(test.loc[test.assigned_fold.eq(fold), "Study_key"]) for fold in range(n_folds)
    }
    if any(len(studies) < 2 or test.assigned_fold.eq(fold).sum() == 0 for fold, studies in test_studies.items()):
        return None
    fold_train = {}
    fold_caps = []
    for fold in range(n_folds):
        test_part = test[test.assigned_fold.eq(fold)]
        blocked_sc = set(test_part.Scaffold_connectivity_group)
        blocked_conn = set(test_part.Connectivity_SMILES)
        blocked_graphs = set(test_part.Graph_ID)
        train = dev[
            ~dev.held_sc
            & dev.assigned_fold.ne(fold)
            & ~dev.Graph_ID.isin(blocked_graphs)
            & ~dev.Connectivity_SMILES.isin(blocked_conn)
            & ~dev.Scaffold_connectivity_group.isin(blocked_sc)
        ].copy()
        fold_train[fold] = train
        cap = pool_caps(train)
        cap["fold"] = fold
        cap["n_test_graphs"] = int(test_part.Graph_ID.nunique())
        cap["n_test_labels"] = int(len(test_part))
        cap["n_test_studies"] = int(test_part.Study_key.nunique())
        cap["n_train_labels"] = int(len(train))
        cap["n_train_graphs"] = int(train.Graph_ID.nunique())
        fold_caps.append(cap)
    return {
        "assignment": assignment,
        "study_fold": study_fold,
        "held_sc": held,
        "dev": dev,
        "test": test,
        "fold_train": fold_train,
        "fold_caps": fold_caps,
    }


def graph_table(frame: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "Graph_ID",
        "Connectivity_SMILES",
        "Scaffold_connectivity_group",
        "Scaffold_FP06_component",
        "MW",
        "cohort",
        "bin50",
    ]
    return frame[cols].drop_duplicates("Graph_ID")


def connectivity_of(pool: pd.DataFrame, ids: list[str]) -> set[str]:
    if not ids:
        return set()
    values = pool.set_index("Graph_ID").loc[list(ids), "Connectivity_SMILES"]
    if len(ids) == 1:
        return {str(values)}
    return set(values.tolist())


def sample_components(train: pd.DataFrame, b: int, k: int, fold: int, seed: int, block: str) -> dict:
    small = representatives(train[train.cohort.eq("nonmacro_le500")])
    macro = representatives(train[train.cohort.eq("macro_500_700")])
    donor = representatives(train[train.cohort.eq("nonmacro_500_700")])
    quotas = largest_remainder(
        k,
        [int(min((macro.bin50 == name).sum(), (donor.bin50 == name).sum())) for name in BINS],
    )
    if quotas is None:
        raise RuntimeError(f"quotas disappeared for fold {fold}")
    chosen_b = sample_ids(small.Graph_ID.tolist(), b, block, fold, seed, "B", b)
    used = connectivity_of(small, chosen_b)
    s_pool = small.loc[~small.Connectivity_SMILES.isin(used), "Graph_ID"].tolist()
    chosen_s = sample_ids(s_pool, k, block, fold, seed, "S", k)
    used.update(connectivity_of(small, chosen_s))

    def draw(pool: pd.DataFrame, prefix: str) -> list[str]:
        picked: list[str] = []
        for name, quota in zip(BINS, quotas):
            cand = pool.loc[(pool.bin50 == name) & ~pool.Connectivity_SMILES.isin(used), "Graph_ID"].tolist()
            draw_ids = sample_ids(cand, quota, block, fold, seed, prefix, name, quota)
            picked.extend(draw_ids)
            used.update(connectivity_of(pool, draw_ids))
        return picked

    chosen_m = draw(macro, "M")
    chosen_n = draw(donor, "N")
    return {"B": chosen_b, "S": chosen_s, "M": chosen_m, "N": chosen_n, "quotas": quotas}


def nested_subset(train: pd.DataFrame, selected: dict, k_small: int, fold: int, seed: int) -> dict[str, list[str]]:
    info = graph_table(train).set_index("Graph_ID")
    counts = [sum(info.loc[gid, "bin50"] == name for gid in selected["M"]) for name in BINS]
    quotas = largest_remainder(k_small, counts)
    if quotas is None:
        raise RuntimeError("nested quota is not contained in the frozen 50")
    out = {"B": list(selected["B"]), "quotas": quotas}
    out["S"] = sample_ids(selected["S"], k_small, "B02NEST", fold, seed, "S", k_small)
    for key in ("M", "N"):
        picked: list[str] = []
        for name, quota in zip(BINS, quotas):
            cand = [gid for gid in selected[key] if info.loc[gid, "bin50"] == name]
            picked.extend(sample_ids(cand, quota, "B02NEST", fold, seed, key, name, quota))
        out[key] = picked
    return out


def labels_for(train: pd.DataFrame, graph_id: str, banned_study: str | None = None) -> pd.DataFrame:
    part = train[train.Graph_ID.eq(graph_id)]
    if banned_study is not None:
        part = part[part.Study_key.ne(banned_study)]
    return part


def refill(selected: list[str], pool: pd.DataFrame, n: int, banned_conn: set[str], *seed_parts: object) -> list[str] | None:
    if pool.empty and n > 0:
        return None
    info = pool.drop_duplicates("Graph_ID").set_index("Graph_ID")
    kept = [gid for gid in selected if gid in info.index]
    used = connectivity_of(pool, kept)
    used.update(banned_conn)
    if len(kept) > n:
        return sample_ids(kept, n, *seed_parts, "trim")
    short = n - len(kept)
    if short == 0:
        return kept
    cand = [gid for gid in info.index if info.loc[gid, "Connectivity_SMILES"] not in used and gid not in kept]
    if len(cand) < short:
        return None
    kept.extend(sample_ids(cand, short, *seed_parts, "refill"))
    return kept


def member_rows(block: str, fold: int, seed: int, arm: str, component: str, ids: list[str], info: pd.DataFrame) -> list[dict]:
    rows = []
    for gid in ids:
        rec = info.loc[gid]
        rows.append(
            {
                "block": block,
                "fold": fold,
                "seed": seed,
                "arm": arm,
                "component": component,
                "Graph_ID": gid,
                "Connectivity_SMILES": rec.Connectivity_SMILES,
                "sc_id": rec.Scaffold_connectivity_group,
                "fp06": rec.Scaffold_FP06_component,
                "bin50": rec.bin50,
                "MW": rec.MW,
                "cohort": rec.cohort,
            }
        )
    return rows


def expand_arms(block: str, fold: int, seed: int, parts: dict, info: pd.DataFrame) -> list[dict]:
    rows = []
    mapping = {"B": parts["B"], "S": parts["S"], "M": parts["M"], "N": parts["N"]}
    for arm, components in ARMS.items():
        for component in components:
            rows.extend(member_rows(block, fold, seed, arm, component, mapping[component], info))
    return rows
