# bRo5 permeability prediction

高分子量非大环化合物的膜渗透性预测及训练数据组成研究。数据快照日期为 2026-11-18。 最近一次项目核对增量日期为 2026-11-19。

## 当前进度

数据索引、审核报告、复现脚本、可直接读取的工作集及全量原始材料均已归档。完整性检查结果见 `verification/remote_verification.json` 和 Releases 中的 `UPLOAD_VERIFIED.json`。模型尚待训练，正式终点、标签聚合和验证划分仍需冻结。

## 当前正式执行协议（V2.2 数据高效迁移与构象感知选择，2026-11-19）

[正式V2.2执行协议](protocols/v2.2/BRO5_EXECUTION_PROTOCOL_V2.2.md) 为当前入口，规范性继承[审稿后V2.1](protocols/v2.1/BRO5_EXECUTION_PROTOCOL_V2.1.md)、[V2.0](protocols/BRO5_EXECUTION_PROTOCOL_V2_20261119.md)、[V1.0](protocols/BRO5_EXECUTION_PROTOCOL_V1_20261119.md) 未变更条款。历史协议及原始快照保留。

V2.2 保留A–F六阶段、Caco-2 A→B主终点、PAMPA独立复现、500<MW≤700 Da非大环主评价及>700 Da严格全局封存；H1/P01为唯一主要比较，H2/H3为D3D同样本3D消融及交互，H4为匹配机制病例。新增C02等样本D-MPNN联合训练/微调和D05构象感知供体选择（固定RF-2D，训练侧独立结构参考集R及候选池，等预算、同MW配额）；DS01和TL01为预登记扩展比较FTRANSFER。C03 Caco-2/PAMPA多任务为条件性，构象加权、全网络微调与联合表征不默认全排列。

| V2.2执行文件 | 用途 |
| --- | --- |
| [transfer_learning_plan.md](protocols/v2.2/transfer_learning_plan.md) | 同样本联合/微调、总优化步数及负迁移 |
| [donor_selection_plan.md](protocols/v2.2/donor_selection_plan.md) | 参考结构成本、固定池、构象距离与供体选择 |
| [multitask_plan.md](protocols/v2.2/multitask_plan.md) | 条件性多任务、跨任务泄漏和辅助标签预算 |
| [comparison_registry.csv](protocols/v2.2/comparison_registry.csv) | 完整继承比较及新增比较家族 |
| [run_matrix.csv](protocols/v2.2/run_matrix.csv) | 完整计划运行块与启动门槛 |
| [compute_budget.json](protocols/v2.2/compute_budget.json) | 原资源边界及新增任务上限 |
| [acceptance_addendum.md](protocols/v2.2/acceptance_addendum.md) | 新模块验收与停止规则 |

三维资格、构象稳定性与MD病例规则继续使用V2.1的 [3D资格](protocols/v2.1/3d_eligibility_rules.json)、[构象稳定性](protocols/v2.1/descriptor_stability_plan.md)、[MD选择](protocols/v2.1/md_selection_protocol.md)及[阶段验收](protocols/v2.1/stage_acceptance_checklist.md)。

当前只归档协议和模板，未完成数据逐行修订、样本与划分容量复算、正式训练、构象或MD。GPU/CPU和存储需实测后锁定。立即执行A/B，后续开展C/D；E按质量及资源验收，F前锁定全部拟评价方法。扩展模块失败不更改H1与封存边界。

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
