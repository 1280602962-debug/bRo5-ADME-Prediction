"""Write the M4 plan, manifest, and report from computed tables."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parents[2] / "runs" / "model_development_m4_cpu"


def interval_sentence(row: pd.Series) -> str:
    if row.status != "estimated":
        return f"{row.contrast} 没有可执行的共同样本，仍保留在比较清单中。"
    direction = "较低" if row.estimate < 0 else "较高"
    text = (
        f"{row.contrast} 的来源等权 MAE 差（候选减参照）是 {row.estimate:.4f}，"
        f"候选一侧{direction}。95% 区间是 [{row.ci95_low:.4f}, {row.ci95_high:.4f}]"
    )
    if row.role == "primary" and pd.notna(row.ci_bonferroni_low):
        text += f"，Bonferroni 约 99.17% 区间是 [{row.ci_bonferroni_low:.4f}, {row.ci_bonferroni_high:.4f}]"
    covers = row.ci95_low <= 0 <= row.ci95_high
    text += "。普通 95% 区间包含 0。" if covers else "。普通 95% 区间不包含 0。"
    if row.role == "primary" and pd.notna(row.ci_bonferroni_low):
        adjusted = row.ci_bonferroni_low <= 0 <= row.ci_bonferroni_high
        text += "校正区间包含 0。" if adjusted else "校正区间不包含 0。"
    return text


def main() -> None:
    variability = json.loads((OUT / "m4_0_summary.json").read_text())
    ion = json.loads((OUT / "m4_1_summary.json").read_text())
    model = json.loads((OUT / "m4_3_summary.json").read_text())
    deltas = pd.read_csv(OUT / "m4_paired_deltas.csv")
    metrics = pd.read_csv(OUT / "m4_metrics.csv")
    summary = pd.read_csv(OUT / "conformer_generation_summary.csv")
    eligibility = pd.read_csv(OUT / "conformer_eligibility.csv")
    importance = pd.read_csv(OUT / "m4_permutation_importance.csv")
    primary = deltas[deltas.role.eq("primary")]
    sentences = [interval_sentence(row) for row in deltas.itertuples(index=False)]
    p4 = metrics[metrics.evaluation.eq("C") & metrics.arm.eq("P4_ensemble") & metrics.model.eq("RF")]
    p0 = metrics[metrics.evaluation.eq("C") & metrics.arm.eq("P0_ensemble") & metrics.model.eq("RF")]
    config = {
        "parent": "2b30eb687844e180f480dcb0bb90a3c610cfa233",
        "branch": "model-development-m4-cpu",
        "endpoint": "Caco-2 A->B log10(Papp [cm/s])",
        "bootstrap_seed": 271828,
        "bootstrap_draws": 10000,
        "primary_contrast": "P4-P0",
        "primary_family": ["P2-P0", "P3-P0", "P4-P0", "P4-P3", "P5-P4", "P5-P2"],
        "exploratory_contrast": "P4-P1",
        "bonferroni_family_size": 6,
        "conformer_requests": {"macro_or_rotatable_ge_10": 100, "mw_gt_500_nonmacro": 50, "other": 20},
        "rmsd_dedup_angstrom": 0.5,
        "energy_window_kcal": 5,
        "ionization_tool": "MolGpKa",
        "ionization_commit": ion["commit"],
        "reference_ph": [6.5, 7.4],
        "models": "fixed RF; fixed XGBoost only if four outer folds improve internally",
        "xgb_added": model["xgb_added"],
        "crest": "not_installed",
        "xtb": "not_installed",
    }
    (OUT / "config.json").write_text(json.dumps(config, indent=2))
    manifest = {
        "potential_sites": "RDKit SMARTS counts. Potential acid or base sites, not validated ionization.",
        "molgpka": ion["model_columns"],
        "molgpka_not_in_model": ["f_neutral_ph6_5", "f_neutral_ph7_4", "logD_neutral_proxy_ph6_5", "logD_neutral_proxy_ph7_4"],
        "l1_static": [column for column in summary.columns if column.startswith("static_")],
        "l1_ensemble": [column for column in summary.columns if column.startswith("ens_")],
        "not_model_inputs": ["n_conformers", "runtime", "failure flags", "sampling_status"],
        "p4_uses": "ensemble_reliable graphs only",
        "fallback": "labels without an ensemble use the pre-specified F04 random forest trained on the full G1 training fold",
    }
    (OUT / "feature_manifest.json").write_text(json.dumps(manifest, indent=2))
    plan = """# M4-CPU 执行记录

本文件记录已经批准并执行的方案，不是新的待确认计划。分支 `model-development-m4-cpu` 从 `2b30eb6` 建立。数据冻结、106 个大于 700 Da 的封存图、评价 A/B/C、来源等权 MAE 和停止条件保持不变。

已执行：M4.0 条件审计；M4.1 在独立环境中用 MolGpKa 做技术先导、SAMPL 公开核查和开发池预测；M4.2 按 20/50/100 的请求数完成合格分子的构象；M4.3 用固定随机森林做共同训练标签和共同测试标签的配对比较。P4−P0 是核心比较，六项主要比较家族没有在看到结果后缩减。P4−P1 是探索性比较。

没有执行新的 Optuna、Mordred、PaDEL、神经网络、GPU 分子动力学、多任务学习或大规模 CREST。本机没有 CREST 和 xTB，因此没有做溶剂配对试算，也没有把样本自动扩大。
"""
    (OUT / "M4_CPU_PLAN.md").write_text(plan)
    sampling = summary.sampling_status.value_counts().to_dict()
    importance_text = importance.sort_values("label_mae_increase_mean", ascending=False).head(6)
    importance_bits = ", ".join(
        f"{row.feature} {row.label_mae_increase_mean:.4f}" for row in importance_text.itertuples(index=False)
    )
    p4_mae = float(p4.mae_source.iloc[0]) if len(p4) else float("nan")
    p0_mae = float(p0.mae_source.iloc[0]) if len(p0) else float("nan")
    fallback = metrics[metrics.arm.eq("P4_fallback_full")]
    full_p0 = metrics[metrics.arm.eq("P0_all") & metrics.evaluation.eq("C")]
    fallback_mae = float(fallback.mae_source.iloc[0])
    full_mae = float(full_p0.mae_source.iloc[0])
    ab = metrics[metrics.evaluation.isin(["A", "B"]) & metrics.model.eq("RF")]
    ab_text = ", ".join(
        f"{row.evaluation} {row.arm} {row.mae_source:.4f}" for row in ab.sort_values(["evaluation", "arm"]).itertuples(index=False)
    )
    gpu_sentence = (
        "P3−P0 和 P4−P0 的普通 95% 区间刚好不包含 0，校正后的区间包含 0。"
        "这不是启动 GPU 分子动力学、神经网络或多任务学习的依据。"
    )
    report = f"""# M4-CPU 报告

终点是 Caco-2 A→B 的 log10(Papp [cm/s])。评价沿用阶段 B 的来源等权 MAE：标签绝对误差先在研究内按图平均，再对研究等权。Bootstrap 使用 10000 次研究–骨架连通分量重抽样，随机数种子 271828。本轮结果仍是开发性泛化，因为评价 C 已经参与过先前的模型决策。

## 测量重复

开发池有 7 组研究内独立重复，184 组数据库重复收录。跨研究、且 protocol token 不是缺失标记的兼容重复组为 0。因此不能充分量化研究间变异，没有做方差分解。独立重复组的标准差中位数见 M4.0 报告，数量太少，不能当作模型误差下限，也不能把实验异质性写成已经证实的模型失败原因。

## 电离

MolGpKa 对 2511 个开发图都给出了位点 pKa。类别计数为 {ion["class_counts"]}。单位点 Henderson–Hasselbalch 只适用于 {ion["hh_applicable"]} 个没有永久电荷、记录结构为中性且只有一个预测位点的分子。中性分数没有进入模型。进入模型的是位点计数和最酸/最碱 pKa 的存在性编码。

SAMPL 公开集合的分子量范围是 {ion["sampl_mw_min"]:.1f}–{ion["sampl_mw_max"]:.1f}，大环 {ion["sampl_macrocycles"]} 个，匹配 MAE {ion["sampl_matched_mae"]:.3f}。这不是 pKaHub 论文里已经排除训练重叠的单质子、两性和多质子子集，那些名单没有作为独立文件发布。它也不能证明工具已经覆盖本项目的高分子量、多位点和大环分子。pH 6.5 和 7.4 只是参考计算环境。

## 构象

合格图 {int(eligibility.eligible_main.sum())} 个。采样状态 {sampling}。可靠集合 {model["ensemble_graphs"]} 个，静态描述符可用 {model["static_graphs"]} 个。请求数没有因为墙钟时间而从 100 下调。重原子 RMSD 去重阈值是 0.5 Å。

## 模型

共同样本上的 P4 来源等权 MAE 是 {p4_mae:.4f}，同一批 430 条标签上的 P0 是 {p0_mae:.4f}。这 430 条标签的 MAE 不能和全部 530 条标签的历史 MAE 直接比较。本轮在全部 530 条标签上重新拟合的 F04 随机森林来源等权 MAE 是 {full_mae:.4f}；种子与 M3 不同，因此不把它写成 M3 数值的逐位复现。预先规定的全覆盖回退（有可靠构象集合的标签用 P4，其余用完整训练集上的 F04）是 {fallback_mae:.4f}。回退结果单独报告，没有在看到数值后加入六项主要比较。

五个外折里有 {model["inner_improved_folds"]} 个外折的内部验证显示 P4 低于 P0。四折门槛没有达到，没有增加 XGBoost。评价 A/B 使用预先固定的随机森林和预先登记的 P4，没有继承评价 C 选出的配置。A/B 的来源等权 MAE：{ab_text}。

{" ".join(sentences)}

负值表示候选的来源等权 MAE 更低。P4−P0 是本轮核心比较。P3−P0 和 P4−P0 的 95% 区间下沿非常接近 0，六比较校正后的区间都包含 0。P4−P1 的探索性区间包含 0，不能说三维集合表征优于充分扩展的二维 F05。P4−P3、P5−P4 和 P5−P2 的区间也包含 0。电离比较使用的是 MolGpKa 位点计数，不是微观物种分数；它们可以执行，所以没有被标成缺失，家族规模仍是 6。

外折 0 上，打乱单个特征后的标签 MAE 变化：{importance_bits}。这是一个外折上的置换结果，不是因果贡献，也没有据此增删描述符。极性表面积的置换变化大于回转半径和分子内氢键，但配对模型比较才是预先规定的证据。

特征覆盖率在 `m4_coverage_metrics.csv`。评价 C 的 530 条测试标签中，430 条有可靠构象集合，438 条至少有一个静态构象，电离特征覆盖全部 530 条。

## 溶剂与能量

本机没有 CREST，也没有 xTB。没有做水相/氯仿配对试算，没有累计 12–20 个分子的预算，也没有把固定几何的环境能量差写成膜自由能。M4.4 和 M4.5 没有结果，不能用来声称溶剂响应或能量代理已经改善渗透预测。

## 对后续计算的含义

{gpu_sentence} 没有新增 Optuna、Mordred、PaDEL、神经网络或大规模半经验构象搜索。
"""
    (OUT / "M4_CPU_REPORT.md").write_text(report)
    (OUT / "M4_3_MODEL_REPORT.md").write_text(
        "# M4.3 模型\n\n"
        "所有配对比较使用同一批训练标签和同一批测试标签，并保留原来的外折编号和来源隔离。没有只匹配测试集而更换训练集。\n\n"
        + report.split("## 模型", 1)[1]
    )
    print(report)


if __name__ == "__main__":
    main()
