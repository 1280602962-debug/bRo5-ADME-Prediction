"""One-molecule xTB probe. A failed or slow calculation does not expand the panel."""

from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from m5_lib import M4, OUT

CHEM = Path(__file__).resolve().parents[1] / "model_development_m4_cpu"
import sys

sys.path.insert(0, str(CHEM))
from m4_chem import bond_geometry_ok, conformer_descriptors, graph_seed, optimize_and_collect, parse_mol  # noqa: E402

XTB = Path("/home/gwj/miniconda3/envs/m5_xtb/bin/xtb")
WORK = OUT / "solvent_pilot"
EH_TO_KCAL = 627.509474


def choose_molecule() -> pd.Series:
    frame = pd.read_csv(M4 / "conformer_eligibility.csv")
    eligible = frame.eligible_main.astype(str).str.lower().eq("true")
    pool = frame[eligible & frame.Topology.eq("NONMACRO_SCREEN") & frame.MW.between(500, 700) & frame.formal_charge.eq(0) & frame.unspecified_tetrahedral.eq(0) & frame.unspecified_double_bond.eq(0)].copy()
    pool = pool[pool.n_rotatable.ge(pool.n_rotatable.median())].sort_values("Graph_ID")
    rng = np.random.default_rng(np.random.SeedSequence([20261126, 53]))
    return pool.iloc[int(rng.integers(0, len(pool)))]


def embed(smiles: str, requested: int, seed: int):
    mol = parse_mol(smiles)
    mol_h = Chem.AddHs(mol)
    properties = AllChem.MMFFGetMoleculeProperties(mol_h, mmffVariant="MMFF94s")
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.enforceChirality = True
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    params.useMacrocycle14config = True
    params.numThreads = 1
    params.pruneRmsThresh = 0.5
    AllChem.EmbedMultipleConfs(mol_h, numConfs=int(requested), params=params)
    kept = []
    for conf_id in range(mol_h.GetNumConformers()):
        result = optimize_and_collect(mol_h, properties, conf_id)
        if result is not None:
            kept.append((conf_id, result[0]))
    kept.sort(key=lambda item: item[1])
    minimum = kept[0][1]
    window = [(conf_id, energy) for conf_id, energy in kept if energy - minimum <= 5.0]
    return mol_h, window


def write_xyz(path: Path, mol, conf_id: int) -> None:
    conf = mol.GetConformer(conf_id)
    lines = [str(mol.GetNumAtoms()), f"conf {conf_id}"]
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol()} {pos.x:.8f} {pos.y:.8f} {pos.z:.8f}")
    path.write_text("\n".join(lines) + "\n")


def total_energy(text: str) -> float:
    found = []
    for line in text.splitlines():
        if "| TOTAL ENERGY" not in line or "Eh" not in line:
            continue
        tokens = line.replace("|", " ").split()
        found.append(float(tokens[tokens.index("Eh") - 1]))
    if not found:
        raise RuntimeError("xTB total energy was not found")
    return found[-1]


def run_xtb(directory: Path, xyz: Path, solvent: str, charge: int, uhf: int, optimize: bool) -> dict:
    directory.mkdir(parents=True, exist_ok=True)
    command = [str(XTB), str(xyz), "--gfn", "2", "--alpb", solvent, "--chrg", str(charge), "--uhf", str(uhf), "--norestart"]
    if optimize:
        command.append("--opt")
    started = time.perf_counter()
    env = dict(**{key: value for key, value in __import__("os").environ.items()})
    env.update({"OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"})
    completed = subprocess.run(command, cwd=directory, capture_output=True, text=True, timeout=1200, env=env)
    elapsed = time.perf_counter() - started
    output = completed.stdout + "\n" + completed.stderr
    (directory / "xtb.out").write_text(output)
    energy = total_energy(output) if completed.returncode == 0 else float("nan")
    return {"returncode": int(completed.returncode), "energy_eh": energy, "seconds": elapsed, "core_hours": elapsed / 3600}


def load_optimized(mol, xyz: Path):
    lines = [line for line in xyz.read_text().splitlines()[2:] if line.strip()]
    if len(lines) != mol.GetNumAtoms():
        return None, "atom_count_changed"
    copy = Chem.Mol(mol)
    copy.RemoveAllConformers()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for index, line in enumerate(lines):
        symbol, x, y, z = line.split()[:4]
        if symbol != mol.GetAtomWithIdx(index).GetSymbol():
            return None, "atom_order_changed"
        conf.SetAtomPosition(index, (float(x), float(y), float(z)))
    conf_id = copy.AddConformer(conf, assignId=True)
    if not bond_geometry_ok(copy, conf_id):
        return copy, "connectivity_or_geometry_changed"
    return copy, "geometry_retained"


def main() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    selected = choose_molecule()
    groups = pd.read_csv(OUT.parents[1] / "data" / "modeling_v1" / "molecule_groups.csv.gz", usecols=["Graph_ID", "SMILES"])
    smiles = groups.set_index("Graph_ID").loc[selected.Graph_ID, "SMILES"]
    requested = int(selected.requested_conformers)
    mol, window = embed(smiles, requested, graph_seed(selected.Graph_ID, requested, "a"))
    stored = pd.read_csv(M4 / "conformer_parts" / f"{selected.Graph_ID}_{requested}_a.csv")
    stored_window = stored[stored.within_5kcal.astype(str).str.lower().eq("true")].sort_values("mmff94s_energy")
    lowest_id, lowest_energy = window[0]
    fresh = conformer_descriptors(mol, lowest_id)
    stored_row = stored_window.iloc[0]
    reproduced = (
        abs(float(fresh["Rgyr"]) - float(stored_row.Rgyr)) <= 0.02
        and abs(float(fresh["PolarSASA_NOH"]) - float(stored_row.PolarSASA_NOH)) <= 1.0
        and int(fresh["IMHB_count"]) == int(stored_row.IMHB_count)
        and abs(float(lowest_energy) - float(stored_row.mmff94s_energy)) <= 1e-3
    )
    geometry_status = "matches_stored_m4_conformer" if reproduced else "independent_pilot_geometry"
    charge = int(Chem.GetFormalCharge(mol))
    electrons = int(sum(atom.GetAtomicNum() for atom in mol.GetAtoms()) - charge)
    uhf = electrons % 2
    rows = []
    spent = 0.0
    for rank, (conf_id, energy) in enumerate(window[:3]):
        if spent > 45 * 60:
            break
        start = WORK / f"start_{rank}"
        start.mkdir(parents=True, exist_ok=True)
        xyz = start / "start.xyz"
        write_xyz(xyz, mol, conf_id)
        start_desc = conformer_descriptors(mol, conf_id)
        for solvent in ("water", "chcl3"):
            single = run_xtb(start / f"{solvent}_sp", xyz, solvent, charge, uhf, False)
            spent += single["seconds"]
            optimized = run_xtb(start / f"{solvent}_opt", xyz, solvent, charge, uhf, True)
            spent += optimized["seconds"]
            opt_xyz = start / f"{solvent}_opt" / "xtbopt.xyz"
            opt_mol, retention = load_optimized(mol, opt_xyz) if opt_xyz.exists() else (None, "missing_output")
            opt_desc = conformer_descriptors(opt_mol, 0) if opt_mol is not None and retention == "geometry_retained" else {}
            rows.append(
                {
                    "Graph_ID": selected.Graph_ID,
                    "geometry_status": geometry_status,
                    "start_rank": rank,
                    "solvent": solvent,
                    "charge": charge,
                    "uhf": uhf,
                    "single_point_eh": single["energy_eh"],
                    "optimized_eh": optimized["energy_eh"],
                    "seconds": single["seconds"] + optimized["seconds"],
                    "structure_retention": retention,
                    "delta_Rgyr": opt_desc.get("Rgyr", np.nan) - start_desc["Rgyr"] if opt_desc else np.nan,
                    "delta_PolarSASA": opt_desc.get("PolarSASA_NOH", np.nan) - start_desc["PolarSASA_NOH"] if opt_desc else np.nan,
                    "delta_IMHB": opt_desc.get("IMHB_count", np.nan) - start_desc["IMHB_count"] if opt_desc else np.nan,
                }
            )
            print(json.dumps({"rank": rank, "solvent": solvent, "seconds": rows[-1]["seconds"], "retention": retention}), flush=True)
    frame = pd.DataFrame(rows)
    if len(frame):
        wide = frame.pivot(index="start_rank", columns="solvent", values="single_point_eh")
        if {"water", "chcl3"}.issubset(wide.columns):
            proxy = (wide["chcl3"] - wide["water"]) * EH_TO_KCAL
            frame["deltaE_transfer_proxy_kcal"] = frame.start_rank.map(proxy)
            frame["deltaE_name"] = "E_chcl3_minus_E_water_on_the_same_geometry_not_a_membrane_free_energy"
    frame.to_csv(OUT / "solvent_pilot_summary.csv", index=False)
    summary = {
        "xtb_version": "6.8.0",
        "solvation": "ALPB water and chcl3; GBSA names H2O and CHCl3 were not used",
        "molecules": 1,
        "geometry_status": geometry_status,
        "starts": int(frame.start_rank.nunique()) if len(frame) else 0,
        "cpu_core_hours": float(frame.seconds.sum() / 3600) if len(frame) else 0.0,
        "crest": "not_installed",
        "l2c": "not_executed",
        "formal_model": False,
        "expanded_panel": False,
    }
    (OUT / "m5_3_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
