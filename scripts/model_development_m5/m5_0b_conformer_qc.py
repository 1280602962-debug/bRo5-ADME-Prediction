"""Compare strict MMFF convergence with the original acceptance rule on one shared embedding."""

from __future__ import annotations

import json
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

import sys

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "model_development_m4_cpu"))
from m4_chem import bond_geometry_ok, conformer_descriptors, graph_seed, heavy_rmsd, parse_mol  # noqa: E402

OUT = REPO / "runs" / "model_development_m5"
M4 = REPO / "runs" / "model_development_m4_cpu"


def relative_change(left: float, right: float, floor: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), floor)


def choose_molecules() -> pd.DataFrame:
    eligibility = pd.read_csv(M4 / "conformer_eligibility.csv")
    eligible = eligibility[eligibility.eligible_main.astype(str).str.lower().eq("true")].copy()
    stability = pd.read_csv(M4 / "conformer_pilot.csv")
    stability_ids = set(stability.Graph_ID.astype(str))
    base = eligible[eligible.Graph_ID.astype(str).isin(stability_ids)].copy()
    base["qc_role"] = "m4_stability_20"
    groups = pd.read_csv(REPO / "data" / "modeling_v1" / "molecule_groups.csv.gz")
    groups["Graph_ID"] = groups.Graph_ID.astype(str)
    pool = eligible[~eligible.Graph_ID.astype(str).isin(stability_ids)].merge(
        groups[["Graph_ID", "SMILES"]], on="Graph_ID", how="left"
    )
    chiral = []
    for row in pool.itertuples(index=False):
        mol = parse_mol(row.SMILES) if isinstance(row.SMILES, str) else None
        centers = 0 if mol is None else sum(atom.HasProp("_CIPCode") for atom in mol.GetAtoms())
        doubles = 0 if mol is None else sum(bond.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ) for bond in mol.GetBonds())
        chiral.append(centers + doubles > 0)
    pool["has_specified_stereo"] = chiral
    pool["rot_bin"] = pd.qcut(pool.n_rotatable.rank(method="first"), 3, labels=["rot_low", "rot_mid", "rot_high"])
    pool["stratum"] = pool.cohort.astype(str) + "|" + pool.rot_bin.astype(str) + "|" + pool.has_specified_stereo.astype(str)
    rng = np.random.default_rng(np.random.SeedSequence([20261126, 5]))
    strata = list(pool.stratum.unique())
    rng.shuffle(strata)
    picked = []
    while len(picked) < 20 and strata:
        for stratum in list(strata):
            part = pool[pool.stratum.eq(stratum) & ~pool.Graph_ID.isin(picked)]
            if part.empty:
                strata.remove(stratum)
                continue
            picked.append(part.iloc[int(rng.integers(0, len(part)))].Graph_ID)
            if len(picked) == 20:
                break
    extra = pool[pool.Graph_ID.isin(picked)].copy()
    extra["qc_role"] = "eligible_draw_20"
    smiles = groups.set_index("Graph_ID").SMILES
    selected = pd.concat([base, extra[base.columns.intersection(extra.columns)]], ignore_index=True)
    if "SMILES" not in selected.columns:
        selected = selected.merge(groups[["Graph_ID", "SMILES"]], on="Graph_ID", how="left")
    else:
        selected["SMILES"] = selected.SMILES.fillna(selected.Graph_ID.map(smiles))
    selected = selected.drop_duplicates("Graph_ID")
    return selected


def cip_codes(mol) -> dict[int, str]:
    codes = {}
    for atom in mol.GetAtoms():
        if atom.HasProp("_CIPCode"):
            codes[atom.GetIdx()] = atom.GetProp("_CIPCode")
    return codes


def bond_stereo(mol) -> dict[tuple[int, int], str]:
    found = {}
    for bond in mol.GetBonds():
        stereo = bond.GetStereo()
        if stereo in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ):
            key = tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
            found[key] = str(stereo)
    return found


def stereo_changes(reference, probed, conf_id: int) -> dict:
    probe = Chem.Mol(probed)
    Chem.AssignStereochemistryFrom3D(probe, confId=conf_id, replaceExistingTags=True)
    Chem.AssignStereochemistry(probe, force=True, cleanIt=True)
    ref_atom = cip_codes(reference)
    new_atom = cip_codes(probe)
    atom_checked = 0
    atom_changed = 0
    for index, code in ref_atom.items():
        atom_checked += 1
        if new_atom.get(index) != code:
            atom_changed += 1
    ref_bond = bond_stereo(reference)
    new_bond = bond_stereo(probe)
    bond_checked = 0
    bond_changed = 0
    for key, code in ref_bond.items():
        bond_checked += 1
        if new_bond.get(key) != code:
            bond_changed += 1
    return {
        "stereo_atoms_checked": atom_checked,
        "stereo_atoms_changed": atom_changed,
        "stereo_bonds_checked": bond_checked,
        "stereo_bonds_changed": bond_changed,
    }


def dedup_and_window(mol_h, records: list[dict]) -> list[dict]:
    ordered = sorted(records, key=lambda item: item["energy"])
    kept = []
    for record in ordered:
        if any(heavy_rmsd(mol_h, record["conf_id"], previous["conf_id"]) < 0.5 for previous in kept):
            continue
        kept.append(record)
    if not kept:
        return []
    minimum = kept[0]["energy"]
    window = []
    for record in kept:
        relative = record["energy"] - minimum
        if relative <= 5.0:
            item = dict(record)
            item["relative_energy_kcal"] = relative
            window.append(item)
    return window


def ensemble_status(n_window: int, requested: int) -> str:
    flexible = requested == 100
    if n_window <= 0:
        return "failed"
    if n_window == 1:
        return "single_conformer"
    reliable = n_window >= (5 if flexible else 3)
    if reliable:
        return "ensemble"
    if n_window == 2 and not flexible:
        return "few_conformers"
    return "undersampled"


def stats(mol_h, window: list[dict]) -> dict:
    if not window:
        return {"n_window": 0, "Rgyr_median": np.nan, "PolarSASA_median": np.nan, "PolarSASA_p10": np.nan, "IMHB_mean": np.nan, "IMHB_fraction": np.nan, "lowest_energy": np.nan, "lowest_conf": -1}
    rows = []
    for record in window:
        desc = conformer_descriptors(mol_h, record["conf_id"])
        rows.append(desc)
    frame = pd.DataFrame(rows)
    return {
        "n_window": int(len(window)),
        "Rgyr_median": float(frame.Rgyr.median()),
        "PolarSASA_median": float(frame.PolarSASA_NOH.median()),
        "PolarSASA_p10": float(frame.PolarSASA_NOH.quantile(0.1)),
        "IMHB_mean": float(frame.IMHB_count.mean()),
        "IMHB_fraction": float((frame.IMHB_count > 0).mean()),
        "lowest_energy": float(window[0]["energy"]),
        "lowest_conf": int(window[0]["conf_id"]),
    }


def run_one(task: tuple) -> dict:
    graph_id, smiles, requested, role = task
    mol = parse_mol(smiles)
    mol_h = Chem.AddHs(mol)
    properties = AllChem.MMFFGetMoleculeProperties(mol_h, mmffVariant="MMFF94s")
    params = AllChem.ETKDGv3()
    params.randomSeed = graph_seed(graph_id, requested, "a")
    params.enforceChirality = True
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    params.useMacrocycle14config = True
    params.numThreads = 1
    params.pruneRmsThresh = 0.5
    ids = AllChem.EmbedMultipleConfs(mol_h, numConfs=int(requested), params=params)
    records = []
    nonconverged = 0
    for conf_id in range(mol_h.GetNumConformers()):
        status1 = AllChem.MMFFOptimizeMolecule(mol_h, mmffVariant="MMFF94s", maxIters=500, confId=conf_id)
        status2 = np.nan
        final = status1
        if status1 == -1:
            continue
        if status1 == 1:
            status2 = AllChem.MMFFOptimizeMolecule(mol_h, mmffVariant="MMFF94s", maxIters=1000, confId=conf_id)
            final = status2
            if status2 == -1:
                continue
        field = AllChem.MMFFGetMoleculeForceField(mol_h, properties, confId=conf_id)
        if field is None:
            continue
        energy = float(field.CalcEnergy())
        if not np.isfinite(energy) or not bond_geometry_ok(mol_h, conf_id):
            continue
        if final != 0:
            nonconverged += 1
        records.append({"conf_id": conf_id, "energy": energy, "status1": int(status1), "status2": status2, "converged": bool(final == 0)})
    original = dedup_and_window(mol_h, records)
    strict = dedup_and_window(mol_h, [record for record in records if record["converged"]])
    original_stats = stats(mol_h, original)
    strict_stats = stats(mol_h, strict)
    stereo = stereo_changes(mol, mol_h, original[0]["conf_id"]) if original else {"stereo_atoms_checked": 0, "stereo_atoms_changed": 0, "stereo_bonds_checked": 0, "stereo_bonds_changed": 0}
    original_member = ensemble_status(original_stats["n_window"], requested)
    strict_member = ensemble_status(strict_stats["n_window"], requested)
    row = {
        "Graph_ID": graph_id,
        "qc_role": role,
        "n_requested": requested,
        "n_embedded": len(ids),
        "n_energy_ok": len(records),
        "n_nonconverged_energy_ok": nonconverged,
        "nonconverged_fraction": nonconverged / len(records) if records else np.nan,
        "original_status": original_member,
        "strict_status": strict_member,
        "membership_flipped": original_member != strict_member,
        "lowest_conf_changed": original_stats["lowest_conf"] != strict_stats["lowest_conf"],
    }
    for name in ("n_window", "Rgyr_median", "PolarSASA_median", "PolarSASA_p10", "IMHB_mean", "IMHB_fraction", "lowest_energy"):
        row[f"original_{name}"] = original_stats[name]
        row[f"strict_{name}"] = strict_stats[name]
    row["Rgyr_exceeds"] = relative_change(original_stats["Rgyr_median"], strict_stats["Rgyr_median"], 0.5) > 0.05
    row["PolarSASA_exceeds"] = relative_change(original_stats["PolarSASA_median"], strict_stats["PolarSASA_median"], 5.0) > 0.05
    row["PolarSASA_p10_exceeds"] = relative_change(original_stats["PolarSASA_p10"], strict_stats["PolarSASA_p10"], 5.0) > 0.10
    row["IMHB_exceeds"] = abs(original_stats["IMHB_mean"] - strict_stats["IMHB_mean"]) > 0.5
    row["IMHB_fraction_exceeds"] = abs(original_stats["IMHB_fraction"] - strict_stats["IMHB_fraction"]) > 0.15
    row.update(stereo)
    return row


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    selected = choose_molecules()
    selected.to_csv(OUT / "conformer_qc_members.csv", index=False)
    tasks = [
        (row.Graph_ID, row.SMILES, int(row.requested_conformers), row.qc_role)
        for row in selected.itertuples(index=False)
    ]
    rows = []
    with ProcessPoolExecutor(max_workers=2, mp_context=mp.get_context("spawn")) as pool:
        futures = [pool.submit(run_one, task) for task in tasks]
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            print(json.dumps({"Graph_ID": row["Graph_ID"], "flip": row["membership_flipped"], "nonconverged": row["nonconverged_fraction"]}), flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "conformer_qc_sensitivity.csv", index=False)
    descriptor_any = frame[["Rgyr_exceeds", "PolarSASA_exceeds", "PolarSASA_p10_exceeds", "IMHB_exceeds", "IMHB_fraction_exceeds"]].any(axis=1)
    frame["descriptor_exceeds"] = descriptor_any
    stereo_molecules = int(((frame.stereo_atoms_changed > 0) | (frame.stereo_bonds_changed > 0)).sum())
    rates = frame.groupby("n_requested").descriptor_exceeds.mean()
    categorical = bool(((frame.groupby("n_requested").size() >= 5) & (rates > 0.5)).any())
    summary = {
        "molecules": int(len(frame)),
        "mean_nonconverged_fraction": float(frame.nonconverged_fraction.mean()),
        "membership_flips": int(frame.membership_flipped.sum()),
        "descriptor_threshold_molecules": int(descriptor_any.sum()),
        "descriptor_threshold_fraction": float(descriptor_any.mean()),
        "stereo_molecules": stereo_molecules,
        "stereo_atom_changes": int(frame.stereo_atoms_changed.sum()),
        "stereo_bond_changes": int(frame.stereo_bonds_changed.sum()),
        "categorical_descriptor_bias": categorical,
        "pause_interpretation": bool(frame.membership_flipped.sum() >= 10 or stereo_molecules >= 5 or categorical),
    }
    (OUT / "m5_0b_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
