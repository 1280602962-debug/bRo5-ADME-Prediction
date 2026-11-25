"""Assemble the model-development report from computed tables."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parents[2] / "runs" / "model_development_v1"


def csv_block(frame: pd.DataFrame) -> str:
    return "```\n" + frame.to_csv(index=False) + "```"


def fmt(value: float) -> str:
    return f"{value:.3f}"


def main() -> None:
    m0_delta = pd.read_csv(OUT / "m0_paired_deltas.csv")
    screen = pd.read_csv(OUT / "model_screening.csv")
    metrics = pd.read_csv(OUT / "metrics_pooled.csv")
    deltas = pd.read_csv(OUT / "paired_deltas.csv")
    identity = pd.read_csv(OUT / "g1_alldev_identity.csv")
    strategy = pd.read_csv(OUT / "strategy_inner_summary.csv")
    train = pd.read_csv(OUT / "training_summary.csv")
    target = metrics[metrics.subset.eq("nonmacro_500_700")].copy()
    full_c = target[target.evaluation.eq("C") & target.n_labels.eq(530)]
    show_cols = ["evaluation", "model", "feature", "parameter_source", "mae_source", "rmse", "r2", "spearman", "n_labels", "n_studies"]
    rf = screen[screen.model.eq("RF") & screen.feature.isin(["F01", "F02", "F04", "F05"])]
    rf_mean = rf.groupby("feature").mae_source.mean()
    focus = strategy[strategy.model.eq("RF") & strategy.feature.isin(["F04", "F05"])]
    c_train = train[train.evaluation.eq("C") & train.model.eq("RF") & train.feature.eq("F05") & train.parameter_source.eq("fixed_selected")]
    favorable = deltas[(deltas.estimate < 0) & deltas.interval_excludes_zero]
    median_gap = float(m0_delta.loc[m0_delta.left.eq("N-only Median") & m0_delta.right.eq("TL-N0"), "estimate"].iloc[0])
    ridge_gap = float(m0_delta.loc[m0_delta.left.eq("N-only Ridge") & m0_delta.right.eq("TL-N0"), "estimate"].iloc[0])
    text = f"""# bRo5-Perm 模型开发报告（M0–M3）

终点是 Caco-2 A→B 的 log10(Papp [cm/s])。评价 C 使用阶段 B 冻结的 515 个图、530 条标签和 153 个研究。来源等权 MAE 先对重复种子取标签绝对误差的平均，再在研究内按图平均，然后对研究等权。这些函数直接调用阶段 B 已复核的实现。本轮结果是开发性泛化：阶段 B/C 的结果已经参与过研究决策，不再称为从未接触的独立测试。

106 个 MW>700 Da 非大环图及其协议规定的关联隔离数据没有进入训练、内部选择、阈值或评价。预测表中也没有这些图。

G1 与阶段 B ALL_DEV 的训练图在五个外折上完全一致。`g1_alldev_identity.csv` 的 `identical` 全部为真，`only_g1` 和 `only_alldev` 都是 0。五个外折的 G1 图数分别是 {", ".join(str(int(v)) for v in identity.g1_graphs)}。

## 阶段 C 的完成状态

C02 的八个策略全部保留，预期并核对为 42400 行预测、每臂 5300 行。独立重算的来源等权 MAE 与已发布数值的最大绝对差为 0。C01 只有被 TL-B 复用的 D-MPNN B 臂预测记录；C01 的 XGBoost 以及 S、M、N、ALL_DEV 没有补跑。仓库没有阶段 C 训练源码，因此没有把这些网络写成已经独立复现的训练。TL-N0 的支持证据是现有训练记录：编码器可训练、训练步数为 60、训练前后哈希都发生变化，且没有 `reused_from`。

## M0

N-only 模型只使用阶段 B 已经冻结的每个种子 50 个 N 供体。Median 是图总权重为 1 的加权中位数。

{csv_block(pd.read_csv(OUT / "m0_metrics.csv"))}

{csv_block(m0_delta)}

N-only Median 的来源等权 MAE 比 TL-N0 低 {abs(median_gap):.3f}，95% 区间不包含 0。因此没有观察到 TL-N0 相对于这个常量预测的明确优势；常量预测更好。N-only Ridge 比 TL-N0 高 {ridge_gap:.3f}，区间也不包含 0。Ridge 与 TL-N0 接近，不能解释成二者都只学习了标签中心。TL-N0 的预测标准差约为实验值标准差的 0.490，Spearman 为 0.074；它有一定预测范围，但排序很弱。N-only Median 的预测标准差比约为 0.135，因为每一折内部几乎是常数，合并后的 Spearman 可以为负。N-only RF 的来源等权 MAE 为 0.613，Spearman 为 0.049，预测标准差比为 0.260，低于 TL-N0，但排序仍然弱。

## 特征、预处理和内部筛选

F01 是阶段 B 的 13 个描述符，训练侧完全无变异的 FormalCharge 被删除后保留 12 列。F02 是 RDKit 2025.09.6 的 217 个二维描述符。F03 是阶段 B 的 ECFP4。F04 是 F01 加 ECFP4。F05 是 F02 加 ECFP4。F06 没有计算。

预处理只在对应的内部训练分区拟合。非有限值先视为缺失；训练侧缺失超过 20% 的列删除；其余用训练侧中位数填补。训练侧标准差低于 1e-6 的列删除。Ridge 和 SVR 只标准化描述符，标准差低于 1e-6 的缩放因子设为 1，指纹保持 0/1。预处理后检查训练和验证矩阵都是有限值。LightGBM 在设置 subsample 时使用 subsample_freq=1。没有使用早停，因为内部验证集就是选择目标。

12 个进入筛选的组合是 RF、XGBoost、LightGBM、SVR 与 F02、F03、F05。F01、F04、Median 和 Ridge 是基线，不因为看到结果后再送入 Optuna。每个外层训练池使用同一个来源和骨架连接隔离的三折内部划分。G1、G2、G3 的验证标签相同，只改变训练侧分子。选择指标是三折预测合并后的一次来源等权 MAE，不是大小不同的验证折 MAE 的平均。

RF 在五个外折上的内部来源等权 MAE：

{csv_block(rf.pivot_table(index="outer_fold", columns="feature", values="mae_source").reset_index())}

F04 在五个外折都低于 F02 和 F05。五个外折的平均分别是 F01 {fmt(rf_mean["F01"])}、F02 {fmt(rf_mean["F02"])}、F04 {fmt(rf_mean["F04"])}、F05 {fmt(rf_mean["F05"])}。对 RF，完整 RDKit 二维描述符相对 F01 的平均差只有约 0.006，而且第 3 折 F02 高于 F01。F04 与 F05 的比较使用同一算法和同一内部划分。Top50/Top100 没有执行。预先规定的筛选对象是 RDKit 描述符列，不包含指纹位；本轮没有把整个组合特征空间拿去排序。LightGBM、XGBoost 和 SVR 的 F02 相对 F01 大约低 0.02–0.03，但五个外折的第一名都不是这些组合，因此没有追加描述符筛选。

12 个组合里，每个外折的前两名是：

{csv_block(pd.read_csv(OUT / "selected_pairs.csv"))}

五个外折的第一名都是 RF-F05。第二名只覆盖选出它的外折，不能与 153 个研究的完整基线配对。

G1 的固定筛选是 12 个组合加 7 个基线，共 19 个配置、5 个外折、3 个内折，285 次内部拟合。其中 12 个组合的 180 次加上随后 30 次试验、2 个组合、5 折、3 内折的 900 次，构成 1080 次主要内部拟合上限。7 个基线的 105 次另计。Optuna 使用 TPE，采样种子 20261125，每个研究 30 次试验，记录在 `hyperparameter_trials.csv`，共 300 行。

固定参数 RF-F05 的外层训练标签 MAE 平均为 {c_train.train_mae_label.mean():.3f}。同一模型的内部来源等权 MAE 约为 0.567，外层为 0.569。训练误差远低于验证和测试，说明模型记住了训练标签，但目标域误差没有跟着下降。

## G1、G2、G3

三种策略先在相同固定参数下比较。G2 只保留 G1 训练池中的非大环分子，G3 只保留其他外折中的 500–700 Da 非大环分子。验证标签不变。G2 和 G3 的 Optuna 试验数是 0，G1 是 30。这个限制是明确的：G1 的 30 次试验没有降低外层来源等权 MAE，资源优先留给已经完成的 G1，不再扩展算法。登记过的 G2/G3 上限是每种策略 15 次试验、2 个组合、5 折、3 内折，各 450 次内部拟合；这些拟合没有执行。G2/G3 额外完成的是固定参数内部拟合，包括四个已选组合和 RF-F04。

下表是五个外折各自合并内部预测后的来源等权 MAE，再对五个外折取平均。这不是把五个外层训练池的验证标签合成一个样本。

{csv_block(focus.sort_values(["feature", "strategy"]))}

G3 的内部和外部误差都高于 G1。G2 与 G1 的差异随外折改变符号，平均差小于 0.003。固定参数比较没有支持用更窄的训练集替代 G1。外层固定参数结果见下表，区间都包含 0。

## 外层结果

500–700 Da 非大环子集：

{csv_block(target[show_cols].sort_values(["evaluation", "mae_source"]))}

评价 A 按 Connectivity_SMILES、评价 B 按 Scaffold_connectivity_group 做随机五折，随机种子都是 20261125。二者只使用预先固定的 RF、XGBoost、LightGBM 和完整 F05。没有把评价 C 选出的配置拿到 A/B 上，因此没有探索性 A/B 结果。也没有复制完整超参数搜索。

评价 C 的完整覆盖模型相对阶段 B 的配对差等于新模型减参考模型。负值表示新模型的来源等权 MAE 更低。置信区间来自阶段 B 的研究–骨架连接成分重抽样，种子 271828，10000 次。重复种子先平均，不当成独立样本。

{csv_block(deltas)}

完整覆盖模型中，估计值为负且区间不包含 0 的比较有 {len(favorable)} 项。固定 RF-F05 相对阶段 B RF B 的差是 -0.015，相对 B+M 是 -0.007，相对 ALL_DEV 是 -0.014，相对本轮 G1 RF-F04 是 -0.013；四个区间都包含 0。G1 RF-F04 相对历史 ALL_DEV 的差是 -0.001，区间也包含 0。成员相同、表征相同、参数相同，剩余差别来自随机种子和实现，不是训练成员不同。Optuna RF-F05 的外层来源等权 MAE 是 0.580，高于固定 RF-F05 的 0.569。内部验证上的小幅下降没有转化成外层改进。

评价 A 的 500–700 Da 非大环子集上，XGBoost、LightGBM 和 RF 的来源等权 MAE 约为 0.494、0.509 和 0.513，R² 约为 0.379、0.340 和 0.347，Spearman 约为 0.623、0.601 和 0.592。评价 B 同一子集的 R² 约为 0.259、0.222 和 0.252。评价 C 固定 RF-F05 的 R² 是 0.009，Spearman 是 0.219。随机划分和骨架划分都比严格来源隔离容易，两者都需要报告。MAE 不高于 0.45 和 R² 不低于 0.4 是原先的期望，不是本轮的通过门槛。

## 直接回答

1. G1 与历史 ALL_DEV 的训练图一致。扩大到全部合格开发数据后，RF-F04 的外层误差与阶段 B 的 500+50 RF 以及 ALL_DEV 的配对区间都包含 0。全数据本身没有稳定改善严格来源隔离下的误差。
2. 对 RF，13 个描述符加 ECFP4 在五个外折的内部误差都低于完整 RDKit 二维描述符，也低于完整 RDKit 加 ECFP4。完整二维描述符不是当前数据上 RF 的更好默认表征。
3. 在预先登记的 12 个组合中，RF-F05 是五个外折的第一名。F04 基线更好，但它不在这 12 个组合里，看到结果后没有补进 Optuna。
4. 30 次试验没有降低 RF-F05 的外层来源等权 MAE。
5. 评价 A/B 的目标子集 Spearman 约为 0.47–0.62、R² 约为 0.22–0.38；评价 C 的 Spearman 约为 0.22、R² 约为 0.01。A/B 使用固定配置，不是评价 C 的优选模型。
6. 没有任何完整覆盖的新模型在相对阶段 B RF B、B+M 或 ALL_DEV 时给出不含 0 的有利区间。不能写成稳定优于阶段 B 的随机森林。
7. G1、G2、G3 使用相同的内部验证标签。固定参数下 G3 更差，G1 与 G2 没有稳定差别。G2/G3 的优化预算是 0 次试验，少于 G1 的 30 次。
8. 固定树模型和 SVR 的训练标签 MAE 约为 0.09–0.18，内部和外层来源等权 MAE 约为 0.55–0.60。继续扩大同一批二维特征的超参数搜索没有优先依据。构象、极性和电离仍是后续问题，本轮没有计算。

## 文献记录

下列数字摘自各文自己的数据、划分和指标。它们不是与本项目来源等权 MAE 的排名，也不能换算成谁的模型更好。

CaliciBoost（Journal of Cheminformatics, 2025, DOI 10.1186/s13321-025-01137-7）使用 TDC Caco-2 的 906 个化合物和官方骨架划分，训练 728、测试 182；另有 OCHEM 分层划分，训练 4377、测试 1095。表征包括 Morgan、Avalon、ErG、RDKit、MACCS、PaDEL、Mordred 和 CDDD。模型是 AutoGluon 0.7.0 的自动集成，种子 42，再用排列重要性和 SHAP 选特征，并用贝叶斯优化重训。TDC 上未筛选的 RDKit、Morgan、Avalon MAE 约为 0.2898、0.2995、0.3009。PaDEL 顶部特征模型在该文 TDC 划分上的 MAE 为 0.2525，RMSE 为 0.3216，R² 为 0.7805，Pearson 为 0.8835。按 TDC 官方协议、五个种子报告的 CaliciBoost MAE 为 0.2560 ± 0.006。本轮没有使用 PaDEL、Mordred、CDDD 或 AutoGluon。

Steinbauer 等（JCIM, 2024, DOI 10.1021/acs.jcim.4c01083）使用拜耳内部 Caco-2 数据，全文 t-SNE 图注明 n=33398。随机 10% 内部测试的训练数是 30058、测试数是 3340。ECFP4 上 LightGBM 的 MAE/AAFE 为 0.357/2.27，RMSE 为 0.479，R² 为 0.620。描述符比较中 RDKit 最好，MAE 0.318，R² 0.683；PaDEL 为 MAE 0.327，R² 0.668。项目相关测试集 n=189，只含 eRo5/bRo5。局部相似模型相对全局模型只有边际变化。算法比较使用 LightGBM 4.1.0 和 scikit-learn 1.3.2。

Wang 等（Journal of Cheminformatics, 2025, DOI 10.1186/s13321-025-00947-z）合并公开数据 7861 条，整理后 5654 个化合物。标签是 log10(Papp，单位为 10^-6 cm/s)，范围 -8.778 到 -3.510，均值 -5.340。表征包括 Morgan 半径 2、1024 位，RDKit2D，以及二者组合；D-MPNN 和 CombinedNet 使用分子图。经典模型用五折交叉验证和网格搜索，神经网络用贝叶斯优化。文中给出的测试 R² 大约在 0.4 到 0.65，RMSE 大约在 0.5 到 0.6。适用域过滤后，XGB_combined 的 R² 为 0.663、RMSE 为 0.479，SVM_combined 的 R² 为 0.667、RMSE 为 0.475。上海齐鲁内部数据上，最高 Tanimoto 相似度为 0.45 以下；XGB_fp 的 Pearson 为 0.702，Spearman 为 0.704，MAE 为 0.629，RMSE 为 0.771。

Doi 与 Esaki（JCIM, 2026, DOI 10.1021/acs.jcim.6c00597）从 ChEMBL 34 整理五个表观渗透率任务。排除 MW>700 Da 后，Caco-2、MDCK、RRCK、LLC-PK1 和 PAMPA 分别为 5480、1276、895、451 和 4378 个化合物。模型包括指纹随机森林、单任务 GCN、迁移 GCN 和多任务 GCN。随机森林的超参数在内部三折里选择。摘要报告的 10 折交叉验证 R² 为 Caco-2 0.503、MDCK 0.500、LLC-PK1 0.558、RRCK 0.449；PAMPA 上 RF 的 R² 为 0.532。评价同时使用随机划分和骨架划分。本轮没有做多任务学习或误差模型。

Kim 等（JCIM, 2025, DOI 10.1021/acs.jcim.5c01600）使用 32 个异双功能降解剂。终点是 Caco-2 A→B 与 B→A 的几何平均被动渗透率。模型为 RF、PLS 和线性 SVM，特征为二维、WT-MetaD 氯仿系综的三维描述符，以及二者合并。评价是 100 次 50/50 划分。摘要报告 PLS 的交叉验证 r² 从只使用二维特征的 0.29 提高到加入三维特征后的 0.48。本轮没有运行 CREST、MD 或 WT-MetaD。

## 复现和停止点

推理由保存的特征表、`selected_parameters.json`、固定参数和 `model_lib.py` 的预处理重新拟合得到。最终模型权重和缩放器没有单独序列化。图在 `runs/model_development_v1/figures/`。

本轮到结果、报告和复核后停止。没有新增 Mordred、PaDEL、GPU 神经网络、CREST、分子动力学、多任务学习或新一轮参数搜索。
"""
    (OUT / "MODEL_DEVELOPMENT_REPORT.md").write_text(text)
    print("report", len(text), "favorable", len(favorable))


if __name__ == "__main__":
    main()
