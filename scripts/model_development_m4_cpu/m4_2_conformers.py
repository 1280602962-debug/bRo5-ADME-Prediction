"""RDKit conformer eligibility, pilot, and development-pool generation."""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from m4_chem import generate_conformers, graph_seed, mmff_properties, parse_mol, request_conformers, site_counts, stereo_counts

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "runs" / "model_development_m4_cpu"
PARTS = OUT / "conformer_parts"
STATIC_COLUMNS = ["Rgyr", "NPR1", "NPR2", "MolecularVolume", "MaximumDimension", "TotalSASA", "PolarSASA_NOH", "PolarFraction", "IMHB_count"]


def molecule_table() -> pd.DataFrame:
    dev = pd.read_csv(REPO / "runs" / "model_development_v1" / "development_pool.csv.gz")
    holdout = set(pd.read_csv(REPO / "data" / "splits_v1" / "holdout_members.csv").Graph_ID.astype(str))
    groups = pd.read_csv(REPO / "data" / "modeling_v1" / "molecule_groups.csv.gz")
    frame = dev.drop_duplicates("Graph_ID")[["Graph_ID", "cohort", "MW"]].copy()
    frame["Graph_ID"] = frame.Graph_ID.astype(str)
    if holdout & set(frame.Graph_ID):
        raise RuntimeError("holdout graph is in the conformer universe")
    groups["Graph_ID"] = groups.Graph_ID.astype(str)
    frame = frame.merge(groups[["Graph_ID", "SMILES", "Topology"]], on="Graph_ID", how="left")
    return frame


def eligibility(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in frame.itertuples(index=False):
        mol = parse_mol(row.SMILES) if isinstance(row.SMILES, str) else None
        if mol is None:
            rows.append({"Graph_ID": row.Graph_ID, "eligible_main": False, "failure": "smiles_failed"})
            continue
        tet, double = stereo_counts(mol)
        sites = site_counts(mol)
        rotatable = int(mol.GetNumAtoms() and __import__("rdkit").Chem.rdMolDescriptors.CalcNumRotatableBonds(mol))
        tpsa = float(__import__("rdkit").Chem.rdMolDescriptors.CalcTPSA(mol))
        fragments = len(__import__("rdkit").Chem.GetMolFrags(mol))
        mmff_ok = mmff_properties(mol) is not None
        eligible = fragments == 1 and tet == 0 and double == 0 and mmff_ok
        reason = ""
        if fragments != 1:
            reason = "multi_fragment"
        elif tet or double:
            reason = "unspecified_stereo"
        elif not mmff_ok:
            reason = "mmff94s_incomplete"
        rows.append(
            {
                "Graph_ID": row.Graph_ID,
                "cohort": row.cohort,
                "MW": row.MW,
                "Topology": row.Topology,
                "n_fragments": fragments,
                "unspecified_tetrahedral": tet,
                "unspecified_double_bond": double,
                "mmff94s_ok": mmff_ok,
                "n_rotatable": rotatable,
                "TPSA": tpsa,
                "formal_charge": int(__import__("rdkit").Chem.GetFormalCharge(mol)),
                "eligible_main": eligible,
                "failure": reason,
                "requested_conformers": request_conformers(float(row.MW), str(row.Topology), rotatable),
                **sites,
            }
        )
    return pd.DataFrame(rows)


def choose_pilot(eligible: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(np.random.SeedSequence([20261126, 2]))
    work = eligible[eligible.eligible_main].copy()
    work["rot_bin"] = pd.qcut(work.n_rotatable.rank(method="first"), 3, labels=["rot_low", "rot_mid", "rot_high"])
    work["tpsa_bin"] = pd.qcut(work.TPSA.rank(method="first"), 3, labels=["tpsa_low", "tpsa_mid", "tpsa_high"])
    work["stratum"] = work.cohort.astype(str) + "|" + work.rot_bin.astype(str) + "|" + work.tpsa_bin.astype(str)
    picked = []
    strata = list(work.stratum.unique())
    rng.shuffle(strata)
    while len(picked) < 100 and strata:
        for stratum in list(strata):
            pool = work[work.stratum.eq(stratum) & ~work.Graph_ID.isin(picked)]
            if pool.empty:
                strata.remove(stratum)
                continue
            choice = pool.iloc[int(rng.integers(0, len(pool)))]
            picked.append(choice.Graph_ID)
            if len(picked) == 100:
                break
    pilot = work[work.Graph_ID.isin(picked)].copy()
    pilot["pilot_set"] = "conformer_100"
    stability = []
    for count, quota in [(20, 7), (50, 7), (100, 6)]:
        part = pilot[pilot.requested_conformers.eq(count) & ~pilot.Graph_ID.isin(stability)]
        take = part.sample(n=min(quota, len(part)), random_state=20261126) if len(part) else part
        stability.extend(take.Graph_ID.tolist())
    if len(stability) < 20:
        rest = pilot[~pilot.Graph_ID.isin(stability)]
        stability.extend(rest.sample(n=min(20 - len(stability), len(rest)), random_state=20261126).Graph_ID.tolist())
    pilot["stability_set"] = pilot.Graph_ID.isin(stability)
    return pilot


def part_path(graph_id: str, n_requested: int, seed_tag: str) -> Path:
    return PARTS / f"{graph_id}_{n_requested}_{seed_tag}.csv"


def run_one(task: tuple) -> dict:
    graph_id, smiles, n_requested, seed_tag = task
    path = part_path(graph_id, n_requested, seed_tag)
    if path.exists():
        return {"Graph_ID": graph_id, "resumed": True, "file": path.name}
    started = time.perf_counter()
    rows, info = generate_conformers(smiles, n_requested, graph_seed(graph_id, n_requested, seed_tag))
    frame = pd.DataFrame(rows)
    frame["Graph_ID"] = graph_id
    frame["n_requested"] = n_requested
    frame["seed_tag"] = seed_tag
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)
    info.update({"Graph_ID": graph_id, "n_requested": n_requested, "seed_tag": seed_tag, "seconds": round(time.perf_counter() - started, 3), "resumed": False})
    return info


def jobs_for(pilot: pd.DataFrame, smiles: pd.Series) -> list[tuple]:
    tasks = []
    for row in pilot.itertuples(index=False):
        tasks.append((row.Graph_ID, smiles.loc[row.Graph_ID], int(row.requested_conformers), "a"))
        if row.stability_set:
            for count in (20, 50, 100):
                for seed_tag in ("a", "b"):
                    task = (row.Graph_ID, smiles.loc[row.Graph_ID], count, seed_tag)
                    if task not in tasks:
                        tasks.append(task)
    return tasks


def relative_change(left: float, right: float, floor: float) -> float:
    denom = max(abs(left), abs(right), floor)
    return abs(left - right) / denom


def aggregate_stability(pilot: pd.DataFrame) -> pd.DataFrame:
    rows = []
    focus = pilot[pilot.stability_set]
    for row in focus.itertuples(index=False):
        summaries = {}
        for count in (20, 50, 100):
            for seed_tag in ("a", "b"):
                path = part_path(row.Graph_ID, count, seed_tag)
                frame = pd.read_csv(path)
                kept = frame[frame.within_5kcal.astype(str).str.lower().eq("true")] if len(frame) else frame
                if kept.empty:
                    summaries[(count, seed_tag)] = None
                    continue
                summaries[(count, seed_tag)] = {
                    "Rgyr": float(kept.Rgyr.median()),
                    "PolarSASA": float(kept.PolarSASA_NOH.median()),
                    "PolarSASA_p10": float(kept.PolarSASA_NOH.quantile(0.1)),
                    "IMHB_mean": float(kept.IMHB_count.mean()),
                    "IMHB_fraction": float((kept.IMHB_count > 0).mean()),
                    "n_kept": int(len(kept)),
                }
        base = summaries.get((int(row.requested_conformers), "a"))
        other = summaries.get((int(row.requested_conformers), "b"))
        high = summaries.get((100, "a"))
        mid = summaries.get((50, "a"))
        def flag(left, right):
            if left is None or right is None:
                return True
            checks = [
                relative_change(left["Rgyr"], right["Rgyr"], 0.5) > 0.05,
                relative_change(left["PolarSASA"], right["PolarSASA"], 5.0) > 0.05,
                relative_change(left["PolarSASA_p10"], right["PolarSASA_p10"], 5.0) > 0.10,
                abs(left["IMHB_mean"] - right["IMHB_mean"]) > 0.5,
                abs(left["IMHB_fraction"] - right["IMHB_fraction"]) > 0.15,
            ]
            return any(checks)
        rows.append(
            {
                "Graph_ID": row.Graph_ID,
                "requested_conformers": int(row.requested_conformers),
                "seed_unstable": flag(base, other),
                "count_50_vs_100_unstable": flag(mid, high),
                "n_kept_default": None if base is None else base["n_kept"],
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--pool", action="store_true")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    PARTS.mkdir(parents=True, exist_ok=True)
    molecules = molecule_table()
    table = eligibility(molecules)
    table.to_csv(OUT / "conformer_eligibility.csv", index=False)
    print("eligible", int(table.eligible_main.sum()), "of", len(table), flush=True)
    if not args.pilot and not args.pool:
        return
    if args.pool and not args.pilot:
        tasks = []
        smiles = molecules.set_index("Graph_ID").SMILES
        for row in table[table.eligible_main].itertuples(index=False):
            tasks.append((row.Graph_ID, smiles.loc[row.Graph_ID], int(row.requested_conformers), "a"))
        print("pool jobs", len(tasks), flush=True)
        done = 0
        with ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context("spawn")) as pool:
            futures = [pool.submit(run_one, task) for task in tasks]
            for future in as_completed(futures):
                info = future.result()
                done += 1
                if done % 50 == 0 or done == len(tasks):
                    print(json.dumps({"done": done, "of": len(tasks), **info}), flush=True)
        return
    pilot = choose_pilot(table)
    pilot.to_csv(OUT / "conformer_pilot_members.csv", index=False)
    smiles = molecules.set_index("Graph_ID").SMILES
    tasks = jobs_for(pilot, smiles)
    print("pilot jobs", len(tasks), flush=True)
    records = []
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context("spawn")) as pool:
        futures = [pool.submit(run_one, task) for task in tasks]
        for future in as_completed(futures):
            info = future.result()
            records.append(info)
            print(json.dumps(info), flush=True)
    pd.DataFrame(records).to_csv(OUT / "conformer_pilot_jobs.csv", index=False)
    stability = aggregate_stability(pilot)
    stability.to_csv(OUT / "conformer_pilot.csv", index=False)
    print("stability unstable seeds", int(stability.seed_unstable.sum()), "count", int(stability.count_50_vs_100_unstable.sum()), flush=True)


if __name__ == "__main__":
    main()
