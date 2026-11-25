"""Compute the full RDKit 2D descriptor matrix for the frozen development graphs."""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "runs" / "model_development_v1"


def main() -> None:
    rdBase.DisableLog("rdApp.*")
    dev = pd.read_csv(OUT / "development_pool.csv.gz")
    smiles = pd.read_csv(REPO / "data" / "splits_v1" / "test_labels.csv.gz", usecols=["Graph_ID", "smiles_std"])
    groups = pd.read_csv(REPO / "data" / "modeling_v1" / "molecule_groups.csv.gz", usecols=["Graph_ID", "SMILES"])
    lookup = smiles.drop_duplicates("Graph_ID").set_index("Graph_ID").smiles_std.to_dict()
    for row in groups.itertuples(index=False):
        lookup.setdefault(row.Graph_ID, row.SMILES)
    graphs = dev[["Graph_ID"]].drop_duplicates().sort_values("Graph_ID")
    names = [name for name, _ in Descriptors.descList]
    calculator = MolecularDescriptorCalculator(names)
    rows = []
    failed = []
    started = time.perf_counter()
    for graph_id in graphs.Graph_ID:
        smiles_value = lookup.get(graph_id)
        mol = Chem.MolFromSmiles(smiles_value) if isinstance(smiles_value, str) else None
        if mol is None:
            failed.append({"Graph_ID": graph_id, "reason": "mol_from_smiles"})
            rows.append([np.nan] * len(names))
            continue
        try:
            values = calculator.CalcDescriptors(mol)
        except Exception as exc:  # descriptor implementation failure is recorded, not hidden
            failed.append({"Graph_ID": graph_id, "reason": type(exc).__name__})
            rows.append([np.nan] * len(names))
            continue
        rows.append(values)
    elapsed = time.perf_counter() - started
    matrix = pd.DataFrame(rows, columns=names)
    matrix.insert(0, "Graph_ID", graphs.Graph_ID.to_numpy())
    finite = np.isfinite(matrix[names].to_numpy(np.float64))
    dictionary = pd.DataFrame(
        {
            "name": names,
            "source": "rdkit.Chem.Descriptors.descList",
            "function": "MolecularDescriptorCalculator.CalcDescriptors",
            "rdkit_version": rdBase.rdkitVersion,
            "n_nonfinite": [int((~finite[:, i]).sum()) for i in range(len(names))],
        }
    )
    matrix.to_csv(OUT / "descriptors_rdkit2d.csv.gz", index=False)
    dictionary.to_csv(OUT / "descriptor_dictionary.csv", index=False)
    manifest = {
        "n_graphs": int(len(graphs)),
        "n_descriptors": int(len(names)),
        "n_failed_molecules": int(len(failed)),
        "failed": failed,
        "seconds": round(elapsed, 3),
        "rdkit_version": rdBase.rdkitVersion,
        "identity_note": "Descriptor calculation does not change Graph_ID, stereo records, or frozen MW.",
    }
    (OUT / "feature_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps({k: manifest[k] for k in ["n_graphs", "n_descriptors", "n_failed_molecules", "seconds", "rdkit_version"]}))


if __name__ == "__main__":
    main()
