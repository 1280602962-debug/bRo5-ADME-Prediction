"""Build evaluation A/B and the shared evaluation-C inner splits.

G1 graph membership is the Stage B ALL_DEV membership. G2 and G3 change only
which of those training labels are allowed, and they share validation labels.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "stage_b"))
from stageb_lib import assign_studies  # noqa: E402

SPLIT = REPO / "data" / "splits_v1"
OUT = REPO / "runs" / "model_development_v1"
RANDOM_SEED = 20261125


def development_table() -> pd.DataFrame:
    outer = pd.read_csv(SPLIT / "outer_members.csv.gz")
    dev = outer[outer.role.isin(["test", "train_eligible", "held_crossfold_sc"])].copy()
    dev["g1"] = dev.role.isin(["test", "train_eligible"])
    dev["g2"] = dev.g1 & dev.cohort.astype(str).str.startswith("nonmacro")
    dev["g3"] = dev.g1 & dev.cohort.eq("nonmacro_500_700")
    return dev


def random_group_folds(frame: pd.DataFrame, group_col: str, evaluation: str) -> pd.DataFrame:
    groups = np.array(sorted(frame[group_col].unique()))
    rng = np.random.default_rng(RANDOM_SEED)
    order = rng.permutation(len(groups))
    assignment = {groups[i]: int(rank % 5) for rank, i in enumerate(order)}
    out = frame[["label_id", "Graph_ID", "Study_key", "Scaffold_connectivity_group", "Connectivity_SMILES", "cohort", "MW", "y"]].copy()
    out["evaluation"] = evaluation
    out["outer_fold"] = out[group_col].map(assignment).astype(int)
    out["inner_fold"] = -1
    out["role"] = "pool"
    out["group_id"] = out[group_col]
    return out


def inner_state(pool: pd.DataFrame, n_folds: int) -> tuple[pd.DataFrame, dict] | None:
    rows = []
    for study, part in pool.groupby("Study_key"):
        target = set(part.loc[part.cohort.eq("nonmacro_500_700"), "Graph_ID"])
        rows.append({"study": study, "n_target": len(target), "target_graphs": target, "graphs": set(part.Graph_ID)})
    assignment = assign_studies(rows, n_folds)
    study_fold = {study: fold for fold, studies in enumerate(assignment) for study in studies}
    marked = pool.copy()
    marked["inner_fold"] = marked.Study_key.map(study_fold).astype(int)
    held = set()
    for sc, part in marked.groupby("Scaffold_connectivity_group"):
        if part.inner_fold.nunique() > 1:
            held.add(sc)
    marked["held_inner_sc"] = marked.Scaffold_connectivity_group.isin(held)
    validation = marked[marked.cohort.eq("nonmacro_500_700") & ~marked.held_inner_sc]
    for fold in range(n_folds):
        studies = set(validation.loc[validation.inner_fold.eq(fold), "Study_key"])
        if len(studies) < 2 or validation.inner_fold.eq(fold).sum() == 0:
            return None
    info = {
        "n_folds": n_folds,
        "n_held_sc": int(len(held)),
        "n_validation_labels": int(len(validation)),
        "validation_by_fold": {str(fold): int(validation.inner_fold.eq(fold).sum()) for fold in range(n_folds)},
    }
    return marked, info


def build_inner(dev: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    frames = []
    notes = []
    g1 = dev[dev.g1].copy()
    for outer in range(5):
        pool = g1[g1.assigned_fold.ne(outer)].copy()
        built = None
        for n_folds in (3, 2):
            built = inner_state(pool, n_folds)
            if built is not None:
                break
        if built is None:
            raise RuntimeError(f"outer fold {outer} has no isolated inner split")
        marked, info = built
        info["outer_fold"] = outer
        notes.append(info)
        validation = marked[marked.cohort.eq("nonmacro_500_700") & ~marked.held_inner_sc]
        for inner in range(info["n_folds"]):
            val_ids = set(validation.loc[validation.inner_fold.eq(inner), "label_id"])
            val_studies = set(validation.loc[validation.inner_fold.eq(inner), "Study_key"])
            val_sc = set(validation.loc[validation.inner_fold.eq(inner), "Scaffold_connectivity_group"])
            blocked = marked.Study_key.isin(val_studies) | marked.Scaffold_connectivity_group.isin(val_sc) | marked.held_inner_sc
            part = marked[["label_id", "Graph_ID", "Study_key", "Scaffold_connectivity_group", "Connectivity_SMILES", "cohort", "MW", "y", "g2", "g3"]].copy()
            part["evaluation"] = "C_inner"
            part["outer_fold"] = outer
            part["inner_fold"] = inner
            part["is_validation"] = part.label_id.isin(val_ids)
            part["train_blocked"] = blocked.to_numpy()
            part = part[part.is_validation | ~part.train_blocked]
            part["role"] = np.where(part.is_validation, "validation", "train_candidate")
            frames.append(part)
    inner = pd.concat(frames, ignore_index=True)
    return inner, notes


def compare_alldev(dev: pd.DataFrame) -> pd.DataFrame:
    members = pd.read_csv(SPLIT / "sampling_members.csv.gz")
    alldev = members[members.block.eq("B01") & members.arm.eq("ALL_DEV")]
    rows = []
    g1 = dev[dev.g1]
    for fold in range(5):
        g1_ids = set(g1.loc[g1.assigned_fold.ne(fold), "Graph_ID"])
        ad_ids = set(alldev.loc[alldev.fold.eq(fold), "Graph_ID"])
        rows.append(
            {
                "fold": fold,
                "g1_graphs": len(g1_ids),
                "alldev_graphs": len(ad_ids),
                "only_g1": len(g1_ids - ad_ids),
                "only_alldev": len(ad_ids - g1_ids),
                "identical": g1_ids == ad_ids,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    dev = development_table()
    dev.to_csv(OUT / "development_pool.csv.gz", index=False)
    eval_a = random_group_folds(dev, "Connectivity_SMILES", "A")
    eval_b = random_group_folds(dev, "Scaffold_connectivity_group", "B")
    inner, notes = build_inner(dev)
    keep = [
        "evaluation",
        "outer_fold",
        "inner_fold",
        "role",
        "label_id",
        "Graph_ID",
        "Study_key",
        "Scaffold_connectivity_group",
        "Connectivity_SMILES",
        "cohort",
        "MW",
        "y",
        "group_id",
    ]
    a_out = eval_a.copy()
    b_out = eval_b.copy()
    flags = dev.set_index("label_id")[["g2", "g3"]]
    for frame in (a_out, b_out):
        frame["g2"] = flags.g2.reindex(frame.label_id).to_numpy()
        frame["g3"] = flags.g3.reindex(frame.label_id).to_numpy()
        frame["is_validation"] = False
        frame["train_blocked"] = False
    common = [
        "evaluation",
        "outer_fold",
        "inner_fold",
        "role",
        "label_id",
        "Graph_ID",
        "Study_key",
        "Scaffold_connectivity_group",
        "Connectivity_SMILES",
        "cohort",
        "MW",
        "y",
        "g2",
        "g3",
        "is_validation",
        "train_blocked",
    ]
    splits = pd.concat([a_out[common], b_out[common], inner[common]], ignore_index=True)
    splits.to_csv(OUT / "evaluation_splits.csv.gz", index=False)
    identity = compare_alldev(dev)
    identity.to_csv(OUT / "g1_alldev_identity.csv", index=False)
    summary = {
        "development_labels": int(len(dev)),
        "development_graphs": int(dev.Graph_ID.nunique()),
        "g1_labels": int(dev.g1.sum()),
        "holdout_excluded": True,
        "random_seed": RANDOM_SEED,
        "eval_a_group": "Connectivity_SMILES",
        "eval_b_group": "Scaffold_connectivity_group",
        "inner": notes,
        "g1_equals_alldev_graphs": bool(identity.identical.all()),
    }
    (OUT / "split_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(identity.to_string(index=False))


if __name__ == "__main__":
    main()
