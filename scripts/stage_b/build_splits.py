"""Freeze the Stage B holdout, source folds, features, and sampling members."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stageb_lib import (  # noqa: E402
    BINS,
    BUDGETS,
    DESC_NAMES,
    SEEDS,
    assert_quota_examples,
    bin_of,
    budget_fits,
    build_fold_state,
    connectivity_of,
    expand_arms,
    featurize,
    graph_table,
    largest_remainder,
    load_main,
    member_rows,
    nested_subset,
    pool_caps,
    refill,
    representatives,
    sample_components,
    tie_key,
)


def reason_frame(frame: pd.DataFrame, hold: pd.DataFrame) -> pd.Series:
    study = set(hold.Study_key)
    conn = set(hold.Connectivity_SMILES)
    sc = set(hold.Scaffold_connectivity_group)
    fp06 = set(hold.Scaffold_FP06_component)
    reasons = []
    for row in frame.itertuples(index=False):
        found = []
        if row.Study_key in study:
            found.append("LINKED_STUDY")
        if row.Connectivity_SMILES in conn:
            found.append("LINKED_CONNECTIVITY")
        if row.Scaffold_connectivity_group in sc:
            found.append("LINKED_SC")
        if row.Scaffold_FP06_component in fp06:
            found.append("LINKED_FP06")
        reasons.append("|".join(found))
    return pd.Series(reasons, index=frame.index)


def adapt_members(parts: dict, train: pd.DataFrame, b: int, k: int, quotas: list[int], fold: int, seed: int, tag: str) -> dict | None:
    info = graph_table(train).set_index("Graph_ID")
    small = representatives(train[train.cohort.eq("nonmacro_le500")])
    macro = representatives(train[train.cohort.eq("macro_500_700")])
    donor = representatives(train[train.cohort.eq("nonmacro_500_700")])

    def still(ids: list[str]) -> list[str]:
        return [gid for gid in ids if gid in info.index]

    chosen_b = refill(still(parts["B"]), small, b, set(), tag, fold, seed, "B")
    if chosen_b is None:
        return None
    used = connectivity_of(small, chosen_b)
    chosen_s = refill(still(parts["S"]), small, k, used, tag, fold, seed, "S")
    if chosen_s is None:
        return None
    used |= connectivity_of(small, chosen_s)

    def binned(pool: pd.DataFrame, ids: list[str], prefix: str) -> list[str] | None:
        picked: list[str] = []
        local = set(used)
        for name, quota in zip(BINS, quotas):
            kept = [gid for gid in still(ids) if info.loc[gid, "bin50"] == name]
            sub = pool[pool.bin50.eq(name)]
            got = refill(kept, sub, quota, local, tag, fold, seed, prefix, name)
            if got is None:
                return None
            picked.extend(got)
            local |= connectivity_of(sub, got)
        return picked

    chosen_m = binned(macro, parts["M"], "M")
    if chosen_m is None:
        return None
    used |= connectivity_of(macro, chosen_m)
    chosen_n = binned(donor, parts["N"], "N")
    if chosen_n is None:
        return None
    return {"B": chosen_b, "S": chosen_s, "M": chosen_m, "N": chosen_n, "quotas": quotas}


def largest_macro_study(train: pd.DataFrame) -> str:
    macro = train[train.cohort.eq("macro_500_700")]
    reps = set(representatives(macro).Graph_ID)
    counts = macro[macro.Graph_ID.isin(reps)].groupby("Study_key").Graph_ID.nunique()
    if counts.empty:
        raise RuntimeError("M pool is empty while choosing a source to remove")
    ranked = sorted(counts.items(), key=lambda item: (-int(item[1]), tie_key(str(item[0]))))
    return str(ranked[0][0])


def make_pilot(state: dict, outer_fold: int) -> dict | None:
    train = state["fold_train"][outer_fold]
    inner = build_fold_state(train, 2)
    if inner is None:
        return None
    validation = inner["test"][inner["test"].assigned_fold.eq(0)]
    internal_train = inner["fold_train"][0]
    if validation.Study_key.nunique() < 2 or validation.Graph_ID.nunique() == 0:
        return None
    caps = [pool_from_train(internal_train)]
    chosen = next(((b, k) for b, k in BUDGETS if budget_fits(caps, b, k)), None)
    if chosen is None:
        return None
    b, k = chosen
    parts = sample_components(internal_train, b, k, outer_fold, 17, "PILOT")
    return {
        "outer_fold": outer_fold,
        "b": b,
        "k": k,
        "parts": parts,
        "train": internal_train,
        "test": validation,
        "n_components_later": True,
    }


def pool_from_train(train: pd.DataFrame) -> dict:
    return pool_caps(train)


def main() -> None:
    assert_quota_examples()
    repo = Path(__file__).resolve().parents[2]
    out = repo / "data" / "splits_v1"
    run = repo / "runs" / "stage_b"
    out.mkdir(parents=True, exist_ok=True)
    run.mkdir(parents=True, exist_ok=True)
    frame = load_main(repo)
    frame["bin50"] = frame.MW.map(bin_of)
    hold = frame[frame.cohort.eq("nonmacro_gt700")].copy()
    dev_mw = frame[frame.MW.le(700) & frame.cohort.ne("nonmacro_gt700")].copy()
    dev_mw["link_reason"] = reason_frame(dev_mw, hold).to_numpy()
    linked = dev_mw[dev_mw.link_reason.ne("")].copy()
    dev = dev_mw[dev_mw.link_reason.eq("")].copy()
    graphs = dev[["Graph_ID", "smiles_std"]].drop_duplicates("Graph_ID")
    desc, fp, failed = featurize(graphs)
    failed_set = set(failed)
    dev = dev[~dev.Graph_ID.isin(failed_set)].copy()
    from rdkit import rdBase

    meta = {
        "rdkit": rdBase.rdkitVersion,
        "descriptor_names": DESC_NAMES,
        "ecfp4": {"radius": 2, "n_bits": 2048, "useChirality": False, "useFeatures": False},
        "n_graphs_featurized": int(len(desc)),
        "n_graphs_failed": int(len(failed_set)),
        "failed_graph_id": failed,
        "classification_mw": "frozen Data V1.0 MW; recomputed descriptor MW is not used for cohorts or bins",
    }
    (run / "feature_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    np.savez_compressed(run / "features_dev.npz", graph_id=desc.Graph_ID.to_numpy(), desc=desc[DESC_NAMES].to_numpy(np.float32), fp=fp)
    desc.to_csv(run / "descriptors_dev.csv.gz", index=False)

    attempts = []
    locked = None
    for n_folds in (5, 4, 3):
        state = build_fold_state(dev, n_folds)
        if state is None:
            attempts.append({"n_folds": n_folds, "status": "test_fold_below_two_studies_or_empty"})
            continue
        for b, k in BUDGETS:
            fits = budget_fits(state["fold_caps"], b, k)
            attempts.append(
                {
                    "n_folds": n_folds,
                    "b": b,
                    "k": k,
                    "status": "feasible" if fits else "capacity_short",
                    "per_fold": state["fold_caps"],
                }
            )
            if fits and locked is None:
                locked = {"n_folds": n_folds, "b": b, "k": k, "state": state}
                break
        if locked is not None:
            break
    if locked is None:
        (out / "budget_decision.json").write_text(
            json.dumps({"status": "NOT_ESTIMABLE", "attempts": attempts}, indent=2, default=str),
            encoding="utf-8",
        )
        raise SystemExit("no fold/budget combination meets the prespecified capacity rules")

    state = locked["state"]
    b, k = locked["b"], locked["k"]
    n_folds = locked["n_folds"]
    primary_quotas = [largest_remainder(k, cap["common_bins"]) for cap in state["fold_caps"]]
    rows: list[dict] = []
    primary_parts: dict[tuple[int, int], dict] = {}
    for fold in range(n_folds):
        train = state["fold_train"][fold]
        info = graph_table(train).set_index("Graph_ID")
        all_ids = train.loc[train.MW.le(700), "Graph_ID"].drop_duplicates().tolist()
        rows.extend(member_rows("B01", fold, -1, "ALL_DEV", "ALL", all_ids, info))
        for seed in SEEDS:
            parts = sample_components(train, b, k, fold, seed, "B01")
            if parts["quotas"] != primary_quotas[fold]:
                raise RuntimeError("sampled quotas do not match the frozen allocation")
            primary_parts[(fold, seed)] = parts
            rows.extend(expand_arms("B01", fold, seed, parts, info))
            if k == 50:
                nested = nested_subset(train, parts, 25, fold, seed)
                for gid in nested["S"] + nested["M"] + nested["N"]:
                    if gid not in set(parts["S"] + parts["M"] + parts["N"] + parts["B"]):
                        raise RuntimeError("nested 25 is not inside the frozen 50")
                rows.extend(expand_arms("B02", fold, seed, nested, info))
                primary_parts[(fold, seed, 25)] = nested

    def add_block(block: str, fold: int, seed: int, parts: dict, train: pd.DataFrame) -> None:
        info = graph_table(train).set_index("Graph_ID")
        rows.extend(expand_arms(block, fold, seed, parts, info))

    def b03_parts_for(mode: str) -> list[tuple] | None:
        collected = []
        for fold in range(n_folds):
            train = state["fold_train"][fold]
            test_fp = set(state["test"].loc[state["test"].assigned_fold.eq(fold), "Scaffold_FP06_component"])
            restricted = train[~train.Scaffold_FP06_component.isin(test_fp)].copy()
            for seed in SEEDS:
                source = primary_parts[(fold, seed)] if mode == "primary" else primary_parts[(fold, seed, 25)]
                use_k = k if mode == "primary" else 25
                adapted = adapt_members(source, restricted, b, use_k, source["quotas"], fold, seed, "B03")
                if adapted is None:
                    return None
                collected.append((fold, seed, adapted, restricted))
        return collected

    b03_mode = "primary" if b03_parts_for("primary") is not None else ("nested25" if k == 50 and b03_parts_for("nested25") is not None else "")
    if b03_mode:
        for fold, seed, adapted, restricted in b03_parts_for(b03_mode):
            add_block("B03", fold, seed, adapted, restricted)
        b03_status = {"run": True, "budget_b": b, "budget_k": 25 if b03_mode == "nested25" else k, "mode": b03_mode}
    else:
        b03_status = {"run": False, "reason": "FP06-restricted pool cannot keep a matched primary or nested quota"}

    def b04_parts_for(mode: str) -> list[tuple] | None:
        collected = []
        for fold in range(n_folds):
            train = state["fold_train"][fold]
            dropped = largest_macro_study(train)
            reduced = train[~((train.cohort.eq("macro_500_700")) & (train.Study_key.eq(dropped)))].copy()
            for seed in SEEDS:
                source = primary_parts[(fold, seed)] if mode == "primary" else primary_parts[(fold, seed, 25)]
                use_k = k if mode == "primary" else 25
                adapted = adapt_members(source, reduced, b, use_k, source["quotas"], fold, seed, "B04")
                if adapted is None:
                    return None
                collected.append((fold, seed, adapted, train, reduced, dropped))
        return collected

    b04_mode = "primary" if b04_parts_for("primary") is not None else ("nested25" if k == 50 and b04_parts_for("nested25") is not None else "")
    if b04_mode:
        b04_sources = {}
        for fold, seed, adapted, train, reduced, dropped in b04_parts_for(b04_mode):
            b04_sources[str(fold)] = dropped
            info_full = graph_table(train).set_index("Graph_ID")
            info_m = graph_table(reduced).set_index("Graph_ID")
            rows.extend(member_rows("B04", fold, seed, "B", "B", adapted["B"], info_full))
            rows.extend(member_rows("B04", fold, seed, "BM", "B", adapted["B"], info_full))
            rows.extend(member_rows("B04", fold, seed, "BM", "M", adapted["M"], info_m))
            rows.extend(member_rows("B04", fold, seed, "BN", "B", adapted["B"], info_full))
            rows.extend(member_rows("B04", fold, seed, "BN", "N", adapted["N"], info_full))
        b04_status = {
            "run": True,
            "budget_b": b,
            "budget_k": 25 if b04_mode == "nested25" else k,
            "mode": b04_mode,
            "dropped_study_by_fold": b04_sources,
            "reference": "B02" if b04_mode == "nested25" else "B01",
        }
    else:
        b04_status = {"run": False, "reason": "removing the largest M source breaks the matched quota"}

    pilot = []
    for outer_fold in (0, 1):
        if outer_fold >= n_folds:
            continue
        made = make_pilot(state, outer_fold)
        if made is None:
            pilot.append({"outer_fold": outer_fold, "run": False, "reason": "internal source split or capacity unavailable"})
            continue
        info = graph_table(made["train"]).set_index("Graph_ID")
        rows.extend(expand_arms("PILOT", outer_fold, 17, made["parts"], info))
        pilot.append({"outer_fold": outer_fold, "run": True, "b": made["b"], "k": made["k"], "n_test_graphs": int(made["test"].Graph_ID.nunique()), "n_test_studies": int(made["test"].Study_key.nunique())})
        made["test"].assign(outer_fold=outer_fold).to_csv(out / f"pilot_test_fold{outer_fold}.csv.gz", index=False)

    members = pd.DataFrame(rows)
    members.to_csv(out / "sampling_members.csv.gz", index=False)
    state["test"].to_csv(out / "test_labels.csv.gz", index=False)

    dev_marked = state["dev"].copy()
    dev_marked["role"] = np.where(dev_marked.held_sc, "held_crossfold_sc", "train_eligible")
    test_ids = set(state["test"].label_id)
    dev_marked.loc[dev_marked.label_id.isin(test_ids), "role"] = "test"
    hold_out = hold.copy()
    hold_out["role"] = "holdout"
    hold_out["assigned_fold"] = np.nan
    hold_out["held_sc"] = False
    hold_out["link_reason"] = "HOLDOUT_GT700"
    linked_out = linked.copy()
    linked_out["role"] = "excluded_link"
    linked_out["assigned_fold"] = np.nan
    linked_out["held_sc"] = False
    failed_out = dev_mw[dev_mw.Graph_ID.isin(failed_set) & dev_mw.link_reason.eq("")].copy()
    failed_out["role"] = "excluded_feature"
    failed_out["assigned_fold"] = np.nan
    failed_out["held_sc"] = False
    failed_out["link_reason"] = "FEATURE_UNAVAILABLE"
    columns = [
        "label_id",
        "Graph_ID",
        "Study_key",
        "condition_id",
        "y",
        "MW",
        "cohort",
        "bin50",
        "Connectivity_SMILES",
        "Scaffold_connectivity_group",
        "Scaffold_FP06_component",
        "assigned_fold",
        "role",
        "link_reason",
    ]
    pieces = []
    for part in (hold_out, linked_out, failed_out, dev_marked):
        use = part.copy()
        if "link_reason" not in use.columns:
            use["link_reason"] = ""
        pieces.append(use.reindex(columns=columns))
    outer = pd.concat(pieces, ignore_index=True)
    missing = set(frame.label_id) - set(outer.label_id)
    if missing:
        rest = frame[frame.label_id.isin(missing)].copy()
        rest["role"] = "excluded_mw_gt700"
        rest["assigned_fold"] = np.nan
        rest["link_reason"] = "MW_GT_700_NOT_NONMACRO_HOLDOUT"
        outer = pd.concat([outer, rest.reindex(columns=columns)], ignore_index=True)
    if set(outer.label_id) != set(frame.label_id) or outer.label_id.duplicated().any():
        raise RuntimeError("outer member file does not cover each main label once")
    outer.to_csv(out / "outer_members.csv.gz", index=False)
    outer[outer.role.ne("train_eligible")].to_csv(out / "split_exclusions.csv.gz", index=False)
    hold_out[columns].to_csv(out / "holdout_members.csv", index=False)

    capacity_rows = []
    for attempt in attempts:
        if "per_fold" not in attempt:
            capacity_rows.append({"record": "attempt", **{key: attempt[key] for key in ("n_folds", "status")}})
            continue
        for cap in attempt["per_fold"]:
            capacity_rows.append(
                {
                    "record": "fold_pool",
                    "n_folds": attempt["n_folds"],
                    "b": attempt["b"],
                    "k": attempt["k"],
                    "status": attempt["status"],
                    "fold": cap["fold"],
                    "n_small_connectivity": cap["n_small"],
                    "common_bins": "|".join(map(str, cap["common_bins"])),
                    "macro_bins": "|".join(map(str, cap["macro_bins"])),
                    "donor_bins": "|".join(map(str, cap["donor_bins"])),
                    "common_total": cap["common_total"],
                    "n_test_graphs": cap["n_test_graphs"],
                    "n_test_labels": cap["n_test_labels"],
                    "n_test_studies": cap["n_test_studies"],
                    "n_train_graphs": cap["n_train_graphs"],
                    "n_train_labels": cap["n_train_labels"],
                }
            )
    pd.DataFrame(capacity_rows).to_csv(out / "fold_capacity.csv", index=False)
    decision = {
        "status": "LOCKED",
        "n_folds": n_folds,
        "b": b,
        "k": k,
        "seeds": SEEDS,
        "quotas_by_fold": primary_quotas,
        "b02": {"run": k == 50, "k": 25 if k == 50 else None, "nested_in_primary": True},
        "b03": b03_status,
        "b04": b04_status,
        "b05": {"run": False, "reason": "Data V1.0 has no comparable PAMPA auxiliary task"},
        "pilot": pilot,
        "holdout_graphs": int(hold.Graph_ID.nunique()),
        "holdout_labels": int(len(hold)),
        "linked_labels": int(len(linked)),
        "feature_failed_graphs": int(len(failed_set)),
        "held_crossfold_sc": int(len(state["held_sc"])),
        "held_crossfold_labels": int(dev_marked.role.eq("held_crossfold_sc").sum()) if False else int((dev_marked.role == "held_crossfold_sc").sum()),
        "attempts": [{key: value for key, value in attempt.items() if key != "per_fold"} | {"per_fold": attempt.get("per_fold")} for attempt in attempts],
        "group_semantics": {
            "Graph_ID": "standardized molecular graph",
            "Connectivity_SMILES": "full connectivity; one budget slot",
            "connectivity_id": "SC group inherited from Scaffold_connectivity_group",
            "fp06_component": "inherited Scaffold_FP06_component; not recomputed",
            "budget_unit": "unique Graph_ID with one representative per Connectivity_SMILES",
            "representative": "lexicographically smallest Graph_ID; labels are not used",
            "same_sc_policy": "different connectivities in one SC may co-occur in training",
        },
        "weight": "each qualified label of a training graph receives sample_weight 1/n",
    }
    (out / "budget_decision.json").write_text(json.dumps(decision, indent=2, default=str), encoding="utf-8")
    print(json.dumps({key: decision[key] for key in ("status", "n_folds", "b", "k", "b02", "b03", "b04", "pilot", "holdout_graphs", "linked_labels", "feature_failed_graphs", "held_crossfold_sc")}, indent=2, default=str))


if __name__ == "__main__":
    main()
