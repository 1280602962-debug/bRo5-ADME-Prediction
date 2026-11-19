# bRo5 permeability prediction

高分子量非大环化合物的膜渗透性预测及训练数据组成研究。数据快照日期为 2026-11-18。 最近一次项目核对增量日期为 2026-11-19。

## 当前进度

数据索引、审核报告、复现脚本、可直接读取的工作集及全量原始材料均已归档。完整性检查结果见 `verification/remote_verification.json` 和 Releases 中的 `UPLOAD_VERIFIED.json`。模型尚待训练，正式终点、标签聚合和验证划分仍需冻结。

## 当前正式执行协议（V2.3 立项审查后有限修订，2026-11-19）

[正式V2.3执行协议](protocols/v2.3/BRO5_EXECUTION_PROTOCOL_V2.3.md) 为当前执行入口，规范性继承[数据迁移V2.2](protocols/v2.2/BRO5_EXECUTION_PROTOCOL_V2.2.md)、[审稿后V2.1](protocols/v2.1/BRO5_EXECUTION_PROTOCOL_V2.1.md)、[V2.0](protocols/BRO5_EXECUTION_PROTOCOL_V2_20261119.md)和[V1.0](protocols/BRO5_EXECUTION_PROTOCOL_V1_20261119.md)未变更条款。原始数据及历史协议保留。

V2.3 保留A–F六阶段、Caco-2 A→B主终点、PAMPA独立复现、500<MW≤700 Da非大环主评价、>700 Da严格全局封存。H1/P01唯一主要比较、H2/H3、匹配H4、TL01及DS01不改变。新增纯ECFP4强二维供体选择基线D07、锁定外推F06、无标签供体选择重叠及覆盖诊断、先导估计精度分析、L1/L2和机制关联诊断；保留条件性多任务，不新增复杂神经网络或全排列。创新主张须对标Steinbauer等2024及相关研究，尚未确认首创。

| V2.3文件 | 用途 |
| --- | --- |
| [closest_work_register.md](protocols/v2.3/closest_work_register.md) | 近邻论文、实际读取的第三方skills及创新边界 |
| [donor_diagnostics_addendum.md](protocols/v2.3/donor_diagnostics_addendum.md) | ECFP-only对照、覆盖、选择重叠及构象诊断 |
| [precision_plan.md](protocols/v2.3/precision_plan.md) | 训练侧先导精度、来源/组件级统计 |
| [pilot_gate_plan.md](protocols/v2.3/pilot_gate_plan.md) | 有限证据恢复、pilot及正式运行准入 |
| [comparison_registry.csv](protocols/v2.3/comparison_registry.csv) | 完整继承统计比较及新增ECFP辅助比较 |
| [run_matrix.csv](protocols/v2.3/run_matrix.csv) | 原21类运行块及新增D07/F06 |
| [compute_budget.json](protocols/v2.3/compute_budget.json) | 原计算上限、新基线及pilot上限 |

专项方案继续使用V2.2的[迁移学习](protocols/v2.2/transfer_learning_plan.md)、[供体选择](protocols/v2.2/donor_selection_plan.md)、[多任务](protocols/v2.2/multitask_plan.md)、[追加验收](protocols/v2.2/acceptance_addendum.md)及V2.1的[三维资格](protocols/v2.1/3d_eligibility_rules.json)、[构象稳定性](protocols/v2.1/descriptor_stability_plan.md)、[MD选择](protocols/v2.1/md_selection_protocol.md)、[阶段验收](protocols/v2.1/stage_acceptance_checklist.md)。

当前仅完成研究设计与模板，真实逐行数据、划分容量、3D覆盖、预测、模拟及硬件成本仍未验收。先A/B数据核查与受限pilot；通过对应门槛且实际资源签署后才能启动批量。扩展失败不更改H1或封存边界；全部拟评价模型与配置锁定后才开展F。V2.3不等于已经执行实验或确认JCIM录用。

## 2026-11-19 最新数据核对（候选，未冻结）

项目内最新核对已记录：PD-L1 来源 17 条记录由 Caco-2 改归 MDCKII-WT；CB30、CB31 补入 Caco-2 A→B 独立测量，分别为 25.8、21.7 × 10⁻⁶ cm/s；Rzepiela CSV 与官方记录匹配 3,348 条，确认 log10(Peff)（cm/s）。逐行修订尚未落实到下列 2026-11-18 压缩工作集。

| 终点 | 500 < MW ≤ 700 Da 非大环候选 | MW > 700 Da 非大环候选 | 500 < MW ≤ 700 Da 大环供体 |
| --- | ---: | ---: | ---: |
| PAMPA | 764 | 144 | 1,481 |
| Caco-2 A→B | 600 | 88 | 207 |

PAMPA 的 144 个高分子量非大环候选含 20 个无环肽来源结构，排除该类后为 124（敏感性口径）。Rzepiela 的 120 组 BAS 身份歧义、54 条 −8 值的检测限语义、Kim 协议及其他结构冲突仍待核实。上述为项目内报告的最新候选数量，并非本次重新计算的训练集或正式冻结人数。

详见 [11 月 19 日增量核对记录](data/current_design_audit_20261119/AUDIT_UPDATE_20261119.md) 和 [候选分组统计 CSV](data/current_design_audit_20261119/candidate_group_counts.csv)。原始快照及校验记录继续保留，未伪造新的逐行数据或哈希。

## 数据入口

| 位置 | 内容 |
| --- | --- |
| `data/curated/working_graph_study_labels.csv.gz` | 当前保守工作集的结构、来源、任务和暂定渗透性标签 |
| `data/curated/split_eligible_working_labels.csv.gz` | 来源折处理后的划分工作集 |
| `data/curated/molecule_groups.csv.gz` | 分子图、SMILES、分子量、环结构筛选和化学分组 |
| `data/current_design_audit_20261119/` | 最新候选统计及增量核对记录（未物化为逐行数据） |
| `data/current_design_audit_20261118/` | 当前方案的数量、划分容量、冲突清单及独立核查结果 |
| `data/bro5_stage2/` | 来源、终点和化学分组审核记录 |
| `data/bro5_stage3/` | 原文证据、结构修订、检测限及身份冲突记录 |
| `data/literature_inventory/` | 分子类别、实验终点和文献核实清单 |
| `snapshot/snapshot_manifest.csv` | 全量文件路径、字节数、SHA-256及去重映射 |
| `scripts/restore_snapshot.py` | 从 Releases 还原文件并逐个核对 SHA-256 |

Pandas 可直接读取 `.csv.gz`。同一分子可能有多个来源标签，记录条数、分子图数量及独立实验次数应按各阶段报告的定义区分。

完整归档下载见 [数据快照 Release](https://github.com/1280602962-debug/bRo5-ADME-Prediction/releases/tag/data-snapshot-20261118)。独立检查已确认 176 个数据分卷的服务器端 SHA-256，以及三个直接读取工作集文件的 Git 内容哈希。旧传输包的 1,511 个文件与已归档原文件逐一 SHA-256 一致。

## 全量快照

全量清单覆盖 1,686 个文件，逻辑大小 15,015,151,924 字节。按相同文件内容去重后为 1,579 个文件，其中 1,577 个随项目压缩包保存，另两个为 ChEMBL 34、37 官方 SQLite 压缩包。官方数据库按原始下载地址获取，核对与本地一致的 SHA-256 后保存到 Releases。

快照路径 `collected/` 保存原始数据库、文献及补充附件、全量提取表、标准化数据、审核与早期脚本，`early_audit/` 保存初始数据资格审核，`current/` 保存第二、三阶段和当前方案的工作集及复算记录。

归档采用内容去重。还原脚本依据清单恢复同内容文件的原有路径。历史传输缓存和 Python 字节码缓存在快照清单范围之外，科研原文件的范围见文件清单。

完整下载和还原需要 GitHub CLI 与 Python 3.12 或以上版本。

```bash
python scripts/restore_snapshot.py --destination restored_bro5
```

## 2026-11-18 基线快照的数据范围（未按 11 月 19 日修订重算）

| 终点 | 保守工作集中的不同结构 | 分子–来源标签 | MW>500 Da 非大环结构 |
| --- | ---: | ---: | ---: |
| PAMPA | 6,872 | 7,423 | 924 |
| Caco-2 A→B | 2,902 | 3,040 | 718 |

上述工作集采用暂缓身份冲突、异常上限、研究内数值大幅分歧及底层研究来源待明确记录的规则。原始记录和暂缓原因完整保留在快照中。具体规则及其限制见当前方案核查报告。

## 研究问题

在共同基础训练集上，分别加入等量的常规分子量非大环、高分子量大环和高分子量非大环数据，比较其预测增广贡献。**Caco-2 A→B 为主终点，PAMPA 独立重复，500<MW≤700 Da 非大环为主评价空间，MW>700 Da 非大环为预先封存的严格外推留出集。** 先纠错与冻结数据，再完成来源隔离及预算复算，随后在固定样本上开展等预算增广、二维/三维表征、深度学习、等样本微调及训练侧构象感知供体选择，以代表分子的溶剂依赖构象动力学辅助解释，最终进行锁定外推评价。

上游数据库和文献文件的权利及使用条件继续适用。各文件原有来源和许可说明随材料保留。
