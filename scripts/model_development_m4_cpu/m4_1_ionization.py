"""MolGpKa site pKa in a separate environment, plus SMARTS potential-site flags.

pH 6.5 and 7.4 are reference calculations. They are not the recorded assay pH.
Henderson–Hasselbalch fractions are written only for a single predicted site
with no permanent charge and a neutral recorded structure.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, rdMolDescriptors

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "model_development_m4_cpu"))
from m4_chem import parse_mol, site_counts  # noqa: E402

OUT = REPO / "runs" / "model_development_m4_cpu"
MOLGPKA = Path("/home/gwj/m4_tools/MolGpKa")
REFERENCE_PH = (6.5, 7.4)
MODEL_COLUMNS = [
    "molgpka_n_acid",
    "molgpka_n_base",
    "molgpka_has_acid",
    "molgpka_has_base",
    "molgpka_permanent",
    "molgpka_min_acid_pka_x_has",
    "molgpka_max_base_pka_x_has",
]


def molecule_table() -> pd.DataFrame:
    dev = pd.read_csv(REPO / "runs" / "model_development_v1" / "development_pool.csv.gz")
    holdout = set(pd.read_csv(REPO / "data" / "splits_v1" / "holdout_members.csv").Graph_ID.astype(str))
    groups = pd.read_csv(REPO / "data" / "modeling_v1" / "molecule_groups.csv.gz")
    frame = dev.drop_duplicates("Graph_ID")[["Graph_ID", "cohort", "MW"]].copy()
    frame["Graph_ID"] = frame.Graph_ID.astype(str)
    if holdout & set(frame.Graph_ID):
        raise RuntimeError("holdout graph entered the ionization universe")
    groups["Graph_ID"] = groups.Graph_ID.astype(str)
    frame = frame.merge(groups[["Graph_ID", "SMILES", "Topology", "Max_ring"]], on="Graph_ID", how="left")
    return frame


def load_models():
    import utils.ionization_group as groups
    from predict_pka import load_model

    groups.smarts_file = str(MOLGPKA / "src" / "utils" / "smarts_pattern.tsv")
    acid = load_model(str(MOLGPKA / "models" / "weight_acid.pth"), device="cpu")
    base = load_model(str(MOLGPKA / "models" / "weight_base.pth"), device="cpu")
    return acid, base


def unique_ids(values) -> list[int]:
    seen = []
    for value in values:
        item = int(value)
        if item not in seen:
            seen.append(item)
    return seen


def predict_sites(mol, acid_model, base_model) -> tuple[list[float], list[float]]:
    from rdkit.Chem.MolStandardize import rdMolStandardize
    from predict_pka import model_pred
    from utils.ionization_group import get_ionization_aid

    uncharger = rdMolStandardize.Uncharger()
    neutral = Chem.MolFromSmiles(Chem.MolToSmiles(uncharger.uncharge(mol)))
    if neutral is None:
        raise RuntimeError("uncharge_failed")
    explicit = AllChem.AddHs(neutral)
    acids = unique_ids(get_ionization_aid(explicit, acid_or_base="acid"))
    bases = unique_ids(get_ionization_aid(explicit, acid_or_base="base"))
    acid_pka = [float(model_pred(explicit, aid, acid_model, device="cpu")) for aid in acids]
    base_pka = [float(model_pred(explicit, aid, base_model, device="cpu")) for aid in bases]
    return acid_pka, base_pka


def fraction_neutral(kind: str, pka: float, ph: float) -> float:
    if kind == "single_acid":
        return float(1.0 / (1.0 + 10 ** (ph - pka)))
    if kind == "single_base":
        return float(1.0 / (1.0 + 10 ** (pka - ph)))
    raise KeyError(kind)


def classify_prediction(sites: dict, formal_charge: int, acid_pka: list[float], base_pka: list[float], status: str) -> dict:
    permanent = int(sites["permanent_cation_atoms"] + sites["permanent_anion_atoms"])
    row = {
        "molgpka_status": status,
        "molgpka_n_acid": np.nan,
        "molgpka_n_base": np.nan,
        "molgpka_acid_pka": "",
        "molgpka_base_pka": "",
        "molgpka_has_acid": np.nan,
        "molgpka_has_base": np.nan,
        "molgpka_permanent": permanent,
        "molgpka_min_acid_pka": np.nan,
        "molgpka_max_base_pka": np.nan,
        "molgpka_min_acid_pka_x_has": np.nan,
        "molgpka_max_base_pka_x_has": np.nan,
        "recorded_formal_charge": formal_charge,
        "ionization_class": "prediction_failure",
        "hh_applicable": False,
        "hh_kind": "",
        "f_neutral_ph6_5": np.nan,
        "f_neutral_ph7_4": np.nan,
        "logD_neutral_proxy_ph6_5": np.nan,
        "logD_neutral_proxy_ph7_4": np.nan,
    }
    if status != "ok":
        return row
    row["molgpka_n_acid"] = len(acid_pka)
    row["molgpka_n_base"] = len(base_pka)
    row["molgpka_acid_pka"] = ";".join(f"{value:.4f}" for value in sorted(acid_pka))
    row["molgpka_base_pka"] = ";".join(f"{value:.4f}" for value in sorted(base_pka))
    row["molgpka_has_acid"] = int(len(acid_pka) > 0)
    row["molgpka_has_base"] = int(len(base_pka) > 0)
    if acid_pka:
        row["molgpka_min_acid_pka"] = float(min(acid_pka))
        row["molgpka_min_acid_pka_x_has"] = float(min(acid_pka))
    else:
        row["molgpka_min_acid_pka_x_has"] = 0.0
    if base_pka:
        row["molgpka_max_base_pka"] = float(max(base_pka))
        row["molgpka_max_base_pka_x_has"] = float(max(base_pka))
    else:
        row["molgpka_max_base_pka_x_has"] = 0.0
    if permanent:
        kind = "permanent_charge"
    elif not acid_pka and not base_pka:
        kind = "no_ionizable_site"
    elif len(acid_pka) == 1 and not base_pka:
        kind = "single_acid"
    elif len(base_pka) == 1 and not acid_pka:
        kind = "single_base"
    elif acid_pka and base_pka:
        kind = "multi_site_or_zwitterion"
    else:
        kind = "multi_site"
    row["ionization_class"] = kind
    single = kind in {"single_acid", "single_base"} and formal_charge == 0 and permanent == 0
    row["hh_applicable"] = bool(single)
    row["hh_kind"] = kind if single else ""
    if single:
        pka = acid_pka[0] if kind == "single_acid" else base_pka[0]
        for ph, column in ((6.5, "f_neutral_ph6_5"), (7.4, "f_neutral_ph7_4")):
            row[column] = fraction_neutral(kind, pka, ph)
    return row


def annotate_logd(row: dict, mol) -> dict:
    if not row["hh_applicable"] or mol is None:
        return row
    logp = float(Crippen.MolLogP(mol))
    for ph, source, target in (
        (6.5, "f_neutral_ph6_5", "logD_neutral_proxy_ph6_5"),
        (7.4, "f_neutral_ph7_4", "logD_neutral_proxy_ph7_4"),
    ):
        fraction = float(row[source])
        row[target] = logp + float(np.log10(max(fraction, 1e-12)))
    return row


def one_molecule(mol, sites: dict, acid_model, base_model) -> dict:
    formal = int(Chem.GetFormalCharge(mol)) if mol is not None else 0
    if mol is None:
        row = classify_prediction(sites, formal, [], [], "smiles_failed")
        return row
    try:
        acid_pka, base_pka = predict_sites(mol, acid_model, base_model)
    except Exception as exc:  # noqa: BLE001
        row = classify_prediction(sites, formal, [], [], "prediction_failure")
        row["failure_message"] = f"{type(exc).__name__}: {exc}"
        return row
    if not all(np.isfinite(acid_pka + base_pka)):
        return classify_prediction(sites, formal, [], [], "nonfinite_pka")
    row = classify_prediction(sites, formal, acid_pka, base_pka, "ok")
    return annotate_logd(row, mol)


def potential_row(graph_id: str, meta: pd.Series, mol) -> dict:
    sites = site_counts(mol) if mol is not None else {
        "potential_acid_sites": np.nan,
        "potential_base_sites": np.nan,
        "permanent_cation_atoms": np.nan,
        "permanent_anion_atoms": np.nan,
        "potential_ionization_class": "smiles_failed",
    }
    rotatable = int(rdMolDescriptors.CalcNumRotatableBonds(mol)) if mol is not None else np.nan
    return {
        "Graph_ID": graph_id,
        "cohort": meta.cohort,
        "MW": meta.MW,
        "Topology": meta.Topology,
        "n_rotatable": rotatable,
        **sites,
    }


def choose_pilot(frame: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(np.random.SeedSequence([20261126, 1]))
    quotas = [("nonmacro_le500", 20), ("nonmacro_500_700", 20), ("macro_500_700", 20)]
    classes = [
        "no_potential_site",
        "single_potential_acid",
        "single_potential_base",
        "potential_zwitterion",
        "multi_potential_site",
        "permanent_charge",
    ]
    picked = []
    shortfall = []
    for cohort, quota in quotas:
        pool = frame[frame.cohort.eq(cohort) & ~frame.Graph_ID.isin(picked)].copy()
        chosen = []
        while len(chosen) < quota and not pool.empty:
            progressed = False
            for kind in classes:
                part = pool[pool.potential_ionization_class.eq(kind) & ~pool.Graph_ID.isin(chosen)]
                if part.empty:
                    continue
                choice = part.iloc[int(rng.integers(0, len(part)))]
                chosen.append(choice.Graph_ID)
                progressed = True
                if len(chosen) == quota:
                    break
            if not progressed:
                break
        if len(chosen) < quota:
            shortfall.append({"cohort": cohort, "requested": quota, "selected": len(chosen)})
        picked.extend(chosen)
    pilot = frame[frame.Graph_ID.isin(picked)].copy()
    pilot["pilot_role"] = "technical_pilot_60"
    pilot.attrs["shortfall"] = shortfall
    return pilot


def macrocycle(mol) -> tuple[int, bool]:
    if mol is None:
        return 0, False
    rings = [len(ring) for ring in mol.GetRingInfo().AtomRings()]
    largest = max(rings) if rings else 0
    return largest, largest >= 12


def assignment_errors(experimental: list[float], predicted: list[float]) -> tuple[list[float], int, str]:
    import itertools

    if not experimental:
        return [], 0, "empty"
    if not predicted:
        return [], len(experimental), "no_prediction"
    k = min(len(experimental), len(predicted))
    if k > 5:
        remaining = predicted.copy()
        errors = []
        for value in experimental:
            if not remaining:
                break
            index = int(np.argmin([abs(value - item) for item in remaining]))
            errors.append(abs(value - remaining.pop(index)))
        return errors, len(experimental) - len(errors), "greedy"
    best = None
    if len(experimental) <= len(predicted):
        for chosen in itertools.permutations(predicted, k):
            error = [abs(left - right) for left, right in zip(experimental, chosen)]
            if best is None or sum(error) < sum(best):
                best = error
        return best, 0, "exact"
    best = None
    for subset in itertools.combinations(experimental, k):
        for chosen in itertools.permutations(predicted, k):
            error = [abs(left - right) for left, right in zip(subset, chosen)]
            if best is None or sum(error) < sum(best):
                best = error
    return best, len(experimental) - k, "exact"


def validate_sampl(acid_model, base_model) -> pd.DataFrame:
    source = OUT / "external" / "SAMPL_combined_final_processed.tsv"
    source.parent.mkdir(parents=True, exist_ok=True)
    raw = Path("/tmp/pkahub/SAMPL_combined_final_processed.tsv")
    if not source.exists():
        source.write_bytes(raw.read_bytes())
    table = pd.read_csv(source, sep="\t")
    rows = []
    for molid, part in table.groupby("molid", sort=False):
        smiles = str(part.reprSMILES_processed.iloc[0])
        mol = parse_mol(smiles)
        largest, macro = macrocycle(mol)
        experimental = [float(value) for value in part.Exp_Macro_pKa_processed]
        status = "ok"
        message = ""
        acid_pka: list[float] = []
        base_pka: list[float] = []
        if mol is None:
            status = "smiles_failed"
        else:
            try:
                acid_pka, base_pka = predict_sites(mol, acid_model, base_model)
            except Exception as exc:  # noqa: BLE001
                status = "prediction_failure"
                message = f"{type(exc).__name__}: {exc}"
        predicted = acid_pka + base_pka
        errors, unmatched, method = assignment_errors(experimental, predicted)
        rows.append(
            {
                "molid": molid,
                "dataset": part.dataset.iloc[0],
                "smiles": smiles,
                "n_experimental": len(experimental),
                "experimental_pka": ";".join(f"{value:.4f}" for value in experimental),
                "MW": float(Descriptors.MolWt(mol)) if mol is not None else np.nan,
                "heavy_atoms": int(mol.GetNumHeavyAtoms()) if mol is not None else np.nan,
                "largest_ring": largest,
                "macrocycle_ring_ge_12": macro,
                "n_predicted_sites": len(predicted),
                "predicted_pka": ";".join(f"{value:.4f}" for value in sorted(predicted)),
                "status": status,
                "failure_message": message,
                "match_method": method,
                "n_matched": len(errors),
                "n_unmatched_experimental": unmatched,
                "matched_mae": float(np.mean(errors)) if errors else np.nan,
                "matched_abs_errors": ";".join(f"{value:.4f}" for value in errors),
            }
        )
    return pd.DataFrame(rows)


def run_graphs(frame: pd.DataFrame, acid_model, base_model) -> pd.DataFrame:
    rows = []
    for index, meta in enumerate(frame.itertuples(index=False), start=1):
        mol = parse_mol(meta.SMILES) if isinstance(meta.SMILES, str) else None
        base = potential_row(meta.Graph_ID, meta, mol)
        sites = {key: base[key] for key in base if key.startswith("potential") or key.startswith("permanent")}
        started = time.perf_counter()
        try:
            predicted = one_molecule(mol, sites, acid_model, base_model)
        except Exception:  # noqa: BLE001
            predicted = {"molgpka_status": "prediction_failure", "failure_message": traceback.format_exc(limit=1)}
        predicted["seconds"] = round(time.perf_counter() - started, 3)
        rows.append({**base, **predicted})
        if index % 100 == 0 or index == len(frame):
            print(json.dumps({"ionization_done": index, "of": len(frame)}), flush=True)
    return pd.DataFrame(rows)


def write_report(features: pd.DataFrame, pilot: pd.DataFrame, sampl: pd.DataFrame, shortfall: list[dict]) -> None:
    ok = features[features.molgpka_status.eq("ok")]
    matched = sampl[sampl.n_matched.gt(0)]
    errors = []
    for text in matched.matched_abs_errors:
        errors.extend(float(item) for item in str(text).split(";") if item)
    summary = {
        "tool": "MolGpKa",
        "repository": "https://github.com/Xundrug/MolGpKa",
        "commit": "4dc835219baad18c2e705274ac54d49597703966",
        "weights": ["models/weight_acid.pth", "models/weight_base.pth"],
        "environment": "/home/gwj/miniconda3/envs/m4_pka",
        "torch_scatter": "pure torch scatter_add shim; torch-scatter had no wheel for torch 2.14.1+cpu",
        "pka_definition": "microscopic site pKa emitted by MolGpKa after RDKit Uncharger",
        "reference_ph": [6.5, 7.4],
        "reference_ph_meaning": "calculation environment, not the assay pH",
        "development_graphs": int(features.Graph_ID.nunique()),
        "prediction_ok": int(features.molgpka_status.eq("ok").sum()),
        "prediction_failure": int(features.molgpka_status.ne("ok").sum()),
        "class_counts": features.ionization_class.value_counts(dropna=False).to_dict(),
        "hh_applicable": int(features.hh_applicable.astype(str).str.lower().eq("true").sum()) if "hh_applicable" in features else 0,
        "pilot_graphs": int(pilot.Graph_ID.nunique()),
        "pilot_shortfall": shortfall,
        "sampl_molecules": int(len(sampl)),
        "sampl_matched_mae": float(np.mean(errors)) if errors else None,
        "sampl_matched_pairs": int(len(errors)),
        "sampl_unmatched_experimental": int(sampl.n_unmatched_experimental.sum()),
        "sampl_mw_min": float(sampl.MW.min()),
        "sampl_mw_max": float(sampl.MW.max()),
        "sampl_macrocycles": int(sampl.macrocycle_ring_ge_12.sum()),
        "sampl_multi_experimental": int(sampl.n_experimental.gt(1).sum()),
        "model_columns": MODEL_COLUMNS,
        "fractions_not_model_inputs": True,
        "overlap_excluded_benchmark_files": "not published as separate files; not reconstructed",
    }
    (OUT / "m4_1_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    hh = int(features.hh_applicable.fillna(False).astype(bool).sum())
    lines = [
        "# M4.1 电离",
        "",
        "主工具是 MolGpKa（Xundrug/MolGpKa，提交 4dc8352，MIT）。酸性与碱性权重分开，输出的是微观位点 pKa。预测前按该工具的默认流程用 RDKit Uncharger 生成中性副本，原始记录结构的形式电荷被单独保留。计算在独立环境 `m4_pka` 中完成，没有改动 M3 使用的解释器。",
        "",
        "torch 2.14.1 CPU 没有可用的 torch-scatter 轮子。MolGpKa 的 GCN 只在度归一化处调用 `scatter_add`，这里用等价的 PyTorch `scatter_add` 代替。乙酸、苯甲酸、苯酚、吡啶和苯胺的已知 pKa 与预测相差约 0.3 以内；乙胺的预测明显偏离实验值。这是工具在脂肪胺上的失败，不是把 SMARTS 计数写成电离状态。",
        "",
        "pKaSolver 的可再分发包不是原文使用的 Epik 微调模型，本轮没有把它称作 pKaSolver。QupKake 需要本机没有的 xTB 6.4.1。Uni-pKa 权重不在本地 CPU 包中。工具选择发生在项目渗透标签计算之前。",
        "",
        "pKaHub 论文的单质子、两性离子和多质子评价子集（11780、927、2025）没有作为独立文件发布，本轮没有重建。公开实验核查使用 pKaHub 仓库中的 SAMPL6/7/8 与 euroSAMPL1 合并表。这些集合是常用外部集，不是论文里已经排除训练重叠的那一组。SAMPL 的宏观 pKa 与 MolGpKa 的微观 pKa 不是同一物理量，匹配 MAE 不能拿来和论文表格或本项目的渗透 MAE 排名。",
        "",
        f"SAMPL 分子 {len(sampl)} 个，分子量 {sampl.MW.min():.1f}–{sampl.MW.max():.1f}，大环（最大环至少 12）{int(sampl.macrocycle_ring_ge_12.sum())} 个，实验 pKa 多于一条的分子 {int(sampl.n_experimental.gt(1).sum())} 个。匹配位点 {len(errors)} 对，匹配 MAE {float(np.mean(errors)) if errors else float('nan'):.3f}，未匹配实验值 {int(sampl.n_unmatched_experimental.sum())} 条。高分子量、多电离和大环在这个公开集合里覆盖不足。",
        "",
        f"60 个项目分子只做技术先导，选择只用结构分层（种子 20261126 的独立随机流），不用预测误差或渗透标签。公开实验 pKa 与这 60 个分子不能互相替代。开发池成功预测 {int(features.molgpka_status.eq('ok').sum())} / {len(features)}。类别计数：{features.ionization_class.value_counts().to_dict()}。",
        "",
        f"单位点 Henderson–Hasselbalch 只在恰好一个预测位点、没有永久电荷、记录结构形式电荷为 0 时计算。满足条件的图有 {hh} 个。pH 6.5 和 7.4 是参考计算环境，不是原始实验的实际 pH。中性分数和 logD 代理没有进入后续模型矩阵，避免对多位点分子做中位数填补。模型只用完整观测的位点计数、有无标记，以及“无位点时记 0”的最酸/最碱 pKa 乘积编码。",
        "",
        "RDKit SMARTS 列以 potential_ 开头，表示潜在酸碱位点，不是已经验证的电离状态。预测失败的分子保留 Graph_ID，不从数据中删除。",
        "",
    ]
    if shortfall:
        lines.append(f"技术先导分层名额没有凑满：{shortfall}。")
        lines.append("")
    (OUT / "M4_1_IONIZATION_REPORT.md").write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if str(MOLGPKA / "src") not in sys.path:
        sys.path.insert(0, str(MOLGPKA / "src"))
    import torch

    torch.set_num_threads(1)
    molecules = molecule_table()
    base_rows = []
    for meta in molecules.itertuples(index=False):
        mol = parse_mol(meta.SMILES) if isinstance(meta.SMILES, str) else None
        base_rows.append({**potential_row(meta.Graph_ID, meta, mol), "SMILES": meta.SMILES})
    annotated = pd.DataFrame(base_rows)
    pilot = choose_pilot(annotated)
    shortfall = list(pilot.attrs.get("shortfall", []))
    print("pilot", len(pilot), "shortfall", shortfall, flush=True)
    acid_model, base_model = load_models()
    pilot_pred = run_graphs(pilot, acid_model, base_model)
    pilot_pred.to_csv(OUT / "ionization_pilot.csv", index=False)
    if int(pilot_pred.molgpka_status.eq("ok").sum()) < 30:
        pilot_pred[pilot_pred.molgpka_status.ne("ok")].to_csv(OUT / "ionization_failures.csv", index=False)
        raise RuntimeError("MolGpKa technical pilot did not clear 30 successful molecules")
    work = annotated if args.limit <= 0 else annotated.head(args.limit)
    features = run_graphs(work, acid_model, base_model)
    features.to_csv(OUT / "ionization_features.csv.gz", index=False)
    failures = features[features.molgpka_status.ne("ok")]
    failures.to_csv(OUT / "ionization_failures.csv", index=False)
    sampl = validate_sampl(acid_model, base_model)
    sampl.to_csv(OUT / "ionization_tool_validation.csv", index=False)
    write_report(features, pilot_pred, sampl, shortfall)
    print("ionization complete", features.molgpka_status.value_counts().to_dict(), flush=True)


if __name__ == "__main__":
    main()
