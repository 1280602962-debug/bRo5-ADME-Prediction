"""Summarize kept conformers into static and ensemble feature tables."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from m4_chem import generate_conformers, graph_seed

OUT = Path(__file__).resolve().parents[2] / "runs" / "model_development_m4_cpu"
PARTS = OUT / "conformer_parts"


def kept_frame(graph_id: str, n_requested: int, seed_tag: str = "a") -> pd.DataFrame:
    path = PARTS / f"{graph_id}_{n_requested}_{seed_tag}.csv"
    if not path.exists():
        return pd.DataFrame()
    frame = pd.read_csv(path)
    if frame.empty or "within_5kcal" not in frame.columns:
        return pd.DataFrame()
    kept = frame[frame.within_5kcal.astype(str).str.lower().eq("true")].copy()
    return kept


def summarize(graph_id: str, n_requested: int, flexible: bool) -> dict:
    kept = kept_frame(graph_id, n_requested)
    row = {"Graph_ID": graph_id, "n_requested": n_requested, "n_window": int(len(kept))}
    if kept.empty:
        row.update({"static_available": False, "ensemble_reliable": False, "sampling_status": "failed"})
        return row
    lowest = kept.sort_values("relative_energy_kcal").iloc[0]
    for name in ["Rgyr", "NPR1", "NPR2", "MolecularVolume", "MaximumDimension", "TotalSASA", "PolarSASA_NOH", "PolarFraction", "IMHB_count"]:
        row[f"static_{name}"] = float(lowest[name])
    row["static_available"] = True
    # A 100-request molecule needs at least five kept conformers before low percentiles are treated as a distribution.
    reliable = len(kept) >= (5 if flexible else 3)
    if len(kept) == 1:
        status = "single_conformer"
    elif reliable:
        status = "ensemble"
    elif len(kept) == 2 and not flexible:
        status = "few_conformers"
    else:
        status = "undersampled"
    row["sampling_status"] = status
    row["ensemble_reliable"] = reliable
    if not reliable:
        return row
    row.update(
        {
            "ens_Rgyr_median": float(kept.Rgyr.median()),
            "ens_Rgyr_p10": float(kept.Rgyr.quantile(0.1)),
            "ens_Rgyr_iqr": float(kept.Rgyr.quantile(0.75) - kept.Rgyr.quantile(0.25)),
            "ens_PolarSASA_median": float(kept.PolarSASA_NOH.median()),
            "ens_PolarSASA_p10": float(kept.PolarSASA_NOH.quantile(0.1)),
            "ens_PolarSASA_iqr": float(kept.PolarSASA_NOH.quantile(0.75) - kept.PolarSASA_NOH.quantile(0.25)),
            "ens_PolarFraction_median": float(kept.PolarFraction.median()),
            "ens_PolarFraction_p10": float(kept.PolarFraction.quantile(0.1)),
            "ens_IMHB_mean": float(kept.IMHB_count.mean()),
            "ens_IMHB_max": float(kept.IMHB_count.max()),
            "ens_IMHB_fraction": float((kept.IMHB_count > 0).mean()),
            "ens_NPR1_median": float(kept.NPR1.median()),
            "ens_NPR2_median": float(kept.NPR2.median()),
            "ens_Volume_median": float(kept.MolecularVolume.median()),
            "ens_MaximumDimension_median": float(kept.MaximumDimension.median()),
        }
    )
    return row


def main() -> None:
    eligibility = pd.read_csv(OUT / "conformer_eligibility.csv")
    rows = []
    level_paths = []
    for row in eligibility[eligibility.eligible_main].itertuples(index=False):
        flexible = int(row.requested_conformers) == 100
        rows.append(summarize(row.Graph_ID, int(row.requested_conformers), flexible))
        path = PARTS / f"{row.Graph_ID}_{int(row.requested_conformers)}_a.csv"
        if path.exists():
            level_paths.append(path)
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "conformer_generation_summary.csv", index=False)
    summary[summary.static_available].to_csv(OUT / "l1_static_features.csv.gz", index=False)
    summary[summary.ensemble_reliable].to_csv(OUT / "l1_ensemble_features.csv.gz", index=False)
    failed = summary[summary.sampling_status.eq("failed")][["Graph_ID", "n_requested", "sampling_status"]]
    ineligible = eligibility[~eligibility.eligible_main.astype(str).str.lower().eq("true")][
        ["Graph_ID", "failure", "unspecified_tetrahedral", "unspecified_double_bond", "MW", "Topology"]
    ]
    failed.to_csv(OUT / "conformer_failures.csv", index=False)
    ineligible.to_csv(OUT / "conformer_exclusions.csv", index=False)
    if level_paths:
        level = pd.concat([pd.read_csv(path) for path in level_paths], ignore_index=True)
        level = level[level.within_5kcal.astype(str).str.lower().eq("true")]
        level.to_csv(OUT / "conformer_level_descriptors.csv.gz", index=False)
    sdf_dir = OUT / "representative_conformers"
    sdf_dir.mkdir(parents=True, exist_ok=True)
    groups = pd.read_csv(Path(__file__).resolve().parents[2] / "data" / "modeling_v1" / "molecule_groups.csv.gz")
    smiles = groups.assign(Graph_ID=groups.Graph_ID.astype(str)).drop_duplicates("Graph_ID").set_index("Graph_ID").SMILES
    examples = []
    for status in ("single_conformer", "ensemble", "undersampled", "few_conformers"):
        part = summary[summary.sampling_status.eq(status)]
        if not part.empty:
            examples.append(part.iloc[0])
    for example in examples:
        destination = sdf_dir / f"{example.Graph_ID}_{example.sampling_status}.sdf"
        generate_conformers(smiles.loc[example.Graph_ID], int(example.n_requested), graph_seed(example.Graph_ID, int(example.n_requested), "a"), str(destination))
    counts = summary.sampling_status.value_counts(dropna=False).to_dict()
    excluded = ineligible.failure.value_counts(dropna=False).to_dict()
    lines = [
        "# M4.2 构象",
        "",
        "构象请求在稳定性先导之后保持原配置：大环或可旋转键不少于 10 的分子请求 100 个，其余分子量高于 500 Da 的非大环请求 50 个，其余低分子量非大环请求 20 个。20 个先导分子里，按工程阈值有 5 个对种子不稳定，2 个在 50 与 100 个请求之间不稳定。这些阈值是工程稳定性参考，不是热力学收敛标准。没有为了缩短时间而降低高柔性分子的请求数。",
        "",
        "嵌入使用 ETKDGv3，`pruneRmsThresh` 0.5 Å。优化后按重原子对称 RMSD 0.5 Å 去重，主集合保留相对 MMFF94s 能量不超过 5 kcal/mol 的构象。MMFF94s 失败时不改用 UFF。",
        "",
        "立体化学用 RDKit `FindPotentialStereo`。未指定四面体手性和未指定双键立体化学都不进入主三维计算，也不指定单一异构体。",
        "",
        f"开发池 2511 个图中，主计算合格 {int(eligibility.eligible_main.sum())} 个。排除原因：{excluded}。合格分子的采样状态：{counts}。至少有一个有效构象的图 {int(summary.static_available.sum())} 个；集合描述符只对可靠集合写入，共 {int(summary.ensemble_reliable.sum())} 个。请求 100 个构象时，5 kcal/mol 窗口内少于 5 个构象记为采样不足，其低分位数和离散度不进入模型。单个构象仍保留静态描述符。",
        "",
        "三维描述符来自记录的标准化结构，不声称是 pH 7.4 的优势微观物种。PolarSASA_NOH 是 N、O 以及直接连在 N 或 O 上的氢的溶剂可及表面积，探针半径 1.4 Å，使用元素范德华半径。它不是 EPSA。分子内氢键阈值是操作性几何规则。",
        "",
    ]
    (OUT / "M4_2_CONFORMER_REPORT.md").write_text("\n".join(lines))
    print(counts)
    print("static", int(summary.static_available.sum()), "ensemble", int(summary.ensemble_reliable.sum()))


if __name__ == "__main__":
    main()
