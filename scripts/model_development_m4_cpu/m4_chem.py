"""Shared structure checks for M4. Stereochemistry and conformer descriptors."""

from __future__ import annotations

import hashlib
import os

import numpy as np
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, ChemicalFeatures, rdMolAlign, rdFreeSASA, rdMolDescriptors

FEATURE_FACTORY = ChemicalFeatures.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef"))
PERIODIC = Chem.GetPeriodicTable()
ACID_SMARTS = {
    "carboxylic_acid": "[CX3](=O)[OX2H1]",
    "phenol": "[OX2H][c]",
    "tetrazole": "[nH]1nnnc1",
    "sulfonamide": "[SX4](=O)(=O)[NX3;H1,H2]",
}
BASE_SMARTS = {
    "aliphatic_amine": "[NX3;H2,H1;!$(NC=O);!$(NS=O)]",
    "aromatic_nitrogen": "[nX2;H0]",
    "guanidine": "[NX3][CX3](=[NX2])[NX3]",
}
ACID_PATTERNS = {name: Chem.MolFromSmarts(smarts) for name, smarts in ACID_SMARTS.items()}
BASE_PATTERNS = {name: Chem.MolFromSmarts(smarts) for name, smarts in BASE_SMARTS.items()}
PERMANENT_CATION = Chem.MolFromSmarts("[N+,n+;!$([N+](=O)[O-])]")
PERMANENT_ANION = Chem.MolFromSmarts("[-;!$([O-][N+](=O))]")


def graph_seed(*parts: object) -> int:
    digest = hashlib.sha256("|".join(map(str, ("bRo5-Perm-m4", *parts))).encode()).digest()
    return int.from_bytes(digest[:4], "little") % (2**31 - 1)


def parse_mol(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    return mol


def stereo_counts(mol) -> tuple[int, int]:
    unspecified_tet = 0
    unspecified_double = 0
    for item in Chem.FindPotentialStereo(mol):
        if item.specified != Chem.StereoSpecified.Unspecified:
            continue
        if item.type == Chem.StereoType.Atom_Tetrahedral:
            unspecified_tet += 1
        elif item.type == Chem.StereoType.Bond_Double:
            unspecified_double += 1
    return unspecified_tet, unspecified_double


def site_counts(mol) -> dict:
    acids = {name: len(mol.GetSubstructMatches(pattern)) for name, pattern in ACID_PATTERNS.items()}
    bases = {name: len(mol.GetSubstructMatches(pattern)) for name, pattern in BASE_PATTERNS.items()}
    n_acid = int(sum(acids.values()))
    n_base = int(sum(bases.values()))
    n_cation = len(mol.GetSubstructMatches(PERMANENT_CATION))
    n_anion = len(mol.GetSubstructMatches(PERMANENT_ANION))
    if n_cation or n_anion:
        kind = "permanent_charge"
    elif n_acid == 0 and n_base == 0:
        kind = "no_potential_site"
    elif n_acid == 1 and n_base == 0:
        kind = "single_potential_acid"
    elif n_base == 1 and n_acid == 0:
        kind = "single_potential_base"
    elif n_acid > 0 and n_base > 0:
        kind = "potential_zwitterion"
    else:
        kind = "multi_potential_site"
    return {
        "potential_acid_sites": n_acid,
        "potential_base_sites": n_base,
        "permanent_cation_atoms": int(n_cation),
        "permanent_anion_atoms": int(n_anion),
        "potential_ionization_class": kind,
        **{f"potential_{name}": count for name, count in {**acids, **bases}.items()},
    }


def request_conformers(mw: float, topology: str, n_rotatable: int) -> int:
    name = str(topology).upper()
    macro = name.startswith("MACRO")
    if macro or n_rotatable >= 10:
        return 100
    if mw > 500:
        return 50
    return 20


def mmff_properties(mol):
    return AllChem.MMFFGetMoleculeProperties(Chem.AddHs(mol), mmffVariant="MMFF94s")


def vdw_radii(mol) -> list[float]:
    return [float(PERIODIC.GetRvdw(atom.GetAtomicNum())) for atom in mol.GetAtoms()]


def atom_sasa(mol, conf_id: int) -> np.ndarray:
    options = rdFreeSASA.SASAOpts()
    options.probeRadius = 1.4
    rdFreeSASA.CalcSASA(mol, vdw_radii(mol), conf_id, None, options)
    return np.array([float(atom.GetProp("SASA")) for atom in mol.GetAtoms()], dtype=np.float64)


def polar_mask(mol) -> np.ndarray:
    mask = np.zeros(mol.GetNumAtoms(), dtype=bool)
    for atom in mol.GetAtoms():
        number = atom.GetAtomicNum()
        if number in {7, 8}:
            mask[atom.GetIdx()] = True
        elif number == 1 and atom.GetNeighbors() and atom.GetNeighbors()[0].GetAtomicNum() in {7, 8}:
            mask[atom.GetIdx()] = True
    return mask


def imhb_count(mol, conf_id: int) -> tuple[int, int, int]:
    conf = mol.GetConformer(conf_id)
    donors = []
    for feature in FEATURE_FACTORY.GetFeaturesForMol(mol, includeOnly="Donor"):
        for atom_id in feature.GetAtomIds():
            atom = mol.GetAtomWithIdx(int(atom_id))
            hydrogens = [neighbor.GetIdx() for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 1]
            if hydrogens:
                donors.append((int(atom_id), hydrogens))
    acceptors = []
    for feature in FEATURE_FACTORY.GetFeaturesForMol(mol, includeOnly="Acceptor"):
        acceptors.extend(int(atom_id) for atom_id in feature.GetAtomIds())
    pairs = set()
    donor_ids = set()
    acceptor_ids = set()
    positions = conf.GetPositions()
    for donor, hydrogens in donors:
        for acceptor in acceptors:
            if donor == acceptor:
                continue
            if len(Chem.GetShortestPath(mol, donor, acceptor)) - 1 < 3:
                continue
            for hydrogen in hydrogens:
                d_h = positions[hydrogen]
                d_a = positions[acceptor]
                d_d = positions[donor]
                ha = float(np.linalg.norm(d_h - d_a))
                da = float(np.linalg.norm(d_d - d_a))
                cosine = np.dot(d_d - d_h, d_a - d_h) / max(np.linalg.norm(d_d - d_h) * ha, 1e-8)
                angle = float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))
                if da <= 3.5 and ha <= 2.5 and angle >= 120.0:
                    pairs.add((donor, acceptor))
                    donor_ids.add(donor)
                    acceptor_ids.add(acceptor)
                    break
    return len(pairs), len(donor_ids), len(acceptor_ids)


def max_heavy_distance(mol, conf_id: int) -> float:
    conf = mol.GetConformer(conf_id)
    coords = np.array([conf.GetAtomPosition(atom.GetIdx()) for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1], dtype=np.float64)
    if len(coords) < 2:
        return 0.0
    delta = coords[:, None, :] - coords[None, :, :]
    return float(np.sqrt((delta * delta).sum(-1).max()))


def conformer_descriptors(mol, conf_id: int) -> dict:
    sasa = atom_sasa(mol, conf_id)
    polar = sasa[polar_mask(mol)].sum()
    total = float(sasa.sum())
    pmi1 = float(rdMolDescriptors.CalcPMI1(mol, confId=conf_id))
    pmi2 = float(rdMolDescriptors.CalcPMI2(mol, confId=conf_id))
    pmi3 = float(rdMolDescriptors.CalcPMI3(mol, confId=conf_id))
    imhb, donors, acceptors = imhb_count(mol, conf_id)
    return {
        "Rgyr": float(rdMolDescriptors.CalcRadiusOfGyration(mol, confId=conf_id)),
        "PMI1": pmi1,
        "PMI2": pmi2,
        "PMI3": pmi3,
        "NPR1": pmi1 / pmi3 if pmi3 else np.nan,
        "NPR2": pmi2 / pmi3 if pmi3 else np.nan,
        "Asphericity": float(rdMolDescriptors.CalcAsphericity(mol, confId=conf_id)),
        "Eccentricity": float(rdMolDescriptors.CalcEccentricity(mol, confId=conf_id)),
        "SpherocityIndex": float(rdMolDescriptors.CalcSpherocityIndex(mol, confId=conf_id)),
        "MolecularVolume": float(AllChem.ComputeMolVolume(mol, confId=conf_id, gridSpacing=0.5)),
        "MaximumDimension": max_heavy_distance(mol, conf_id),
        "TotalSASA": total,
        "PolarSASA_NOH": float(polar),
        "NonpolarSASA_proxy": float(total - polar),
        "PolarFraction": float(polar / total) if total else np.nan,
        "IMHB_count": int(imhb),
        "IMHB_donor_count": int(donors),
        "IMHB_acceptor_count": int(acceptors),
    }


def bond_geometry_ok(mol, conf_id: int) -> bool:
    conf = mol.GetConformer(conf_id)
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetAtomicNum() == 1 or bond.GetEndAtom().GetAtomicNum() == 1:
            continue
        length = conf.GetAtomPosition(bond.GetBeginAtomIdx()).Distance(conf.GetAtomPosition(bond.GetEndAtomIdx()))
        limit = PERIODIC.GetRcovalent(bond.GetBeginAtom().GetAtomicNum()) + PERIODIC.GetRcovalent(bond.GetEndAtom().GetAtomicNum())
        if length > 2.2 * limit or length < 0.4 * limit:
            return False
    return True


def optimize_and_collect(mol_h, properties, conf_id: int) -> tuple[float, int] | None:
    status = AllChem.MMFFOptimizeMolecule(mol_h, mmffVariant="MMFF94s", maxIters=500, confId=conf_id)
    steps = 500
    if status == -1:
        return None
    if status == 1:
        status = AllChem.MMFFOptimizeMolecule(mol_h, mmffVariant="MMFF94s", maxIters=1000, confId=conf_id)
        steps = 1000
        if status == -1:
            return None
    field = AllChem.MMFFGetMoleculeForceField(mol_h, properties, confId=conf_id)
    if field is None:
        return None
    energy = float(field.CalcEnergy())
    if not np.isfinite(energy) or not bond_geometry_ok(mol_h, conf_id):
        return None
    return energy, steps


def heavy_rmsd(mol_h, left: int, right: int) -> float:
    heavy = Chem.RemoveHs(mol_h)
    return float(rdMolAlign.GetBestRMS(heavy, heavy, left, right))


def generate_conformers(smiles: str, n_requested: int, seed: int, sdf_path: str | None = None) -> tuple[list[dict], dict]:
    mol = parse_mol(smiles)
    if mol is None:
        return [], {"status": "smiles_failed", "n_embedded": 0, "n_optimized": 0, "n_kept": 0}
    mol_h = Chem.AddHs(mol)
    properties = AllChem.MMFFGetMoleculeProperties(mol_h, mmffVariant="MMFF94s")
    if properties is None:
        return [], {"status": "mmff94s_failed", "n_embedded": 0, "n_optimized": 0, "n_kept": 0}
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.enforceChirality = True
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    params.useMacrocycle14config = True
    params.numThreads = 1
    params.pruneRmsThresh = 0.5
    ids = AllChem.EmbedMultipleConfs(mol_h, numConfs=int(n_requested), params=params)
    embedded = len(ids)
    optimized = []
    for conf_id in range(mol_h.GetNumConformers()):
        result = optimize_and_collect(mol_h, properties, conf_id)
        if result is not None:
            optimized.append((conf_id, result[0], result[1]))
    optimized.sort(key=lambda item: item[1])
    kept = []
    for conf_id, energy, steps in optimized:
        if any(heavy_rmsd(mol_h, conf_id, previous) < 0.5 for previous, _, _ in kept):
            continue
        kept.append((conf_id, energy, steps))
    if not kept:
        return [], {"status": "no_valid_conformer", "n_embedded": embedded, "n_optimized": len(optimized), "n_kept": 0}
    minimum = kept[0][1]
    rows = []
    for rank, (conf_id, energy, steps) in enumerate(kept):
        relative = energy - minimum
        if relative > 10.0:
            continue
        row = conformer_descriptors(mol_h, conf_id)
        row.update({"conformer_rank": rank, "mmff94s_energy": energy, "relative_energy_kcal": relative, "optimize_steps": steps, "within_5kcal": relative <= 5.0})
        rows.append(row)
    status = "ok" if any(row["within_5kcal"] for row in rows) else "no_window_conformer"
    if sdf_path and rows:
        writer = Chem.SDWriter(sdf_path)
        for conf_id, energy, _steps in kept:
            if energy - minimum <= 5.0:
                writer.write(mol_h, confId=conf_id)
        writer.close()
    return rows, {"status": status, "n_embedded": embedded, "n_optimized": len(optimized), "n_kept": int(sum(row["within_5kcal"] for row in rows)), "n_within_10kcal": len(rows)}
