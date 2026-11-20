# Data V1.0 数据整理报告

日期：2026-11-20。本文件是阶段A的正式数据结果，不是训练、划分或三维计算结果。

## 纳入标准

Caco-2 主标签只保留能明确识别的 A→B Papp，数值为精确点值，单位已是 log10(cm/s)。进入主集还要求 stage2 已把来源解析为单一 DOI 或专利号。MolMeDB 中由多条参考文献拼出的 DOI 不进入主集。数据库存储 DOI（10.6019/）不作为原始研究。已知外加抑制剂或诱导的记录单独成集。温度或浓度缺失本身不删除终点、方向、单位、结构和来源已经明确的记录，但缺失条件的记录不会自动并入同一实验组。

聚合只发生在同一分子图、同一原始研究、同一终点、同一方向和同一已确认实验协议内。主标签是组内 log10 值的中位数。跨研究的数值不合并。同一数据库记录被两个库收录时只计一次测量。极差超过 0.5 log 单位时先检查单位错位；已确认属于同一协议的测量保留并报告极差，无法解释的冲突暂缓。0.5 是本项目的操作性核查规则，不是文献里的标准差阈值。

## 已知问题

- PD-L1：DOI 10.1016/j.ejmech.2024.116730、ChEMBL assay 2489294 的 17 条活性记录（activity 26123123–26123139，化合物 CB13–CB31）从 Caco-2:AB 改归 MDCKII-WT A→B。17 是出现记录数。CB30 与 CB31 的 ChEMBL 原值分别是 6.9 和 3.9 × 10^-6 cm/s，属于这 17 条之中，已随整组改归，没有重复补录。
- CB30/CB31 的 25.8 与 21.7 × 10^-6 cm/s：底表中没有这两条数值。本次没有重新读到原文表格，结构、条件和独立 Caco-2 测量不能同时复核，因此不新增记录。
- Rzepiela（DOI 10.1021/acs.jmedchem.1c02090）：匹配 CSV 3,348 条，官方库与归档库各收录一次，共 6696 条出现。数值字段已是 log10(Peff)，单位 cm/s，没有再次取对数。其中 48 个编号的记录值为 −8；检测限语义仍未由原文闭合，不进入精确回归。归档说明中的 54 条与这次匹配到的 48 个编号不一致，没有补造差额记录。
- BAS：120 个未写明立体化学的结构各对应 2 个 BAS 编号。其中 118 个结构的两个编号数值不同，相关出现记录暂缓，不把编号当成两个分子，也不把两个数值合成一个标签。2 个结构的两个编号数值相同，合并为一次测量，手性按未知处理，可以进入二维标签。
- Kim（DOI 10.1021/acs.jcim.5c01600）：来源列单位为 nm/s，换算是 log10(数值)−7，与列单位一致；原文实验条件仍未复核，协议标记为未知。方向为 A→B 的精确记录保留在主集，条件缺失的记录不互相合并。B→A 与主集分开。
- 同研究同名样品冲突沿用 stage3 的 10 组。连接关系冲突和无法确定立体归属的记录暂缓。仅仅是手性中心未指定、且记录能对应到唯一结构的分子，保留在二维数据中，并标为三维资格受限。
- log10(cm/s) > −3 的记录暂缓。log10(cm/s) < −9 的记录保留原值并标记，不截断。限定符不是等号的记录不进入精确回归。

## 终点

PAMPA 按系数分开。辅助任务选择看的是系数是否有明确实验定义、500<MW≤700 非大环是否至少有 30 个图和 5 个原始研究，以及最大来源占比是否不超过 80%。未指定系数的混合 PAMPA 不因为行数多而被选中。本次辅助任务：未指定，各系数保留为探索性任务。Caco-2 主集不依赖这个选择。

## 正式规模

下表是精确标签。数据库出现次数、测量记录数、聚合标签数、分子图数和原始研究数分开计数。测量记录数不是独立实验重复次数。

| 集合 | 分组 | 标签 | 分子图 | 连接组 | 原始研究 | 数据库出现 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| caco2_ab_main | macro_500_700 | 238 | 222 | 133 | 35 | 437 |
| caco2_ab_main | nonmacro_500_700 | 644 | 608 | 408 | 185 | 652 |
| caco2_ab_main | nonmacro_gt700 | 167 | 106 | 93 | 54 | 167 |
| caco2_ab_main | nonmacro_le500 | 1879 | 1740 | 1037 | 465 | 1891 |
| caco2_ab_main | other | 296 | 240 | 126 | 42 | 488 |
| caco2_known_intervention | macro_500_700 | 48 | 48 | 15 | 2 | 94 |
| caco2_known_intervention | nonmacro_500_700 | 1 | 1 | 1 | 1 | 1 |
| caco2_known_intervention | nonmacro_le500 | 2 | 2 | 2 | 1 | 2 |
| caco2_known_intervention | other | 56 | 55 | 26 | 5 | 102 |
| pampa_PAMPA_APPARENT_exploratory | nonmacro_500_700 | 4 | 3 | 3 | 2 | 4 |
| pampa_PAMPA_APPARENT_exploratory | nonmacro_le500 | 115 | 72 | 52 | 3 | 115 |
| pampa_PAMPA_COEFFICIENT_UNSPECIFIED_exploratory | macro_500_700 | 211 | 127 | 81 | 18 | 211 |
| pampa_PAMPA_COEFFICIENT_UNSPECIFIED_exploratory | nonmacro_500_700 | 864 | 769 | 481 | 169 | 879 |
| pampa_PAMPA_COEFFICIENT_UNSPECIFIED_exploratory | nonmacro_gt700 | 155 | 125 | 96 | 37 | 169 |
| pampa_PAMPA_COEFFICIENT_UNSPECIFIED_exploratory | nonmacro_le500 | 3047 | 2183 | 1210 | 380 | 3064 |
| pampa_PAMPA_COEFFICIENT_UNSPECIFIED_exploratory | other | 251 | 226 | 129 | 21 | 257 |
| pampa_PAMPA_COEFFICIENT_UNSPECIFIED_exploratory_nonstrict | macro_500_700 | 1399 | 710 | 90 | 28 | 1399 |
| pampa_PAMPA_COEFFICIENT_UNSPECIFIED_exploratory_nonstrict | nonmacro_500_700 | 1876 | 1219 | 717 | 280 | 1890 |
| pampa_PAMPA_COEFFICIENT_UNSPECIFIED_exploratory_nonstrict | nonmacro_gt700 | 275 | 172 | 113 | 46 | 275 |
| pampa_PAMPA_COEFFICIENT_UNSPECIFIED_exploratory_nonstrict | nonmacro_le500 | 8820 | 4454 | 2162 | 737 | 8973 |
| pampa_PAMPA_COEFFICIENT_UNSPECIFIED_exploratory_nonstrict | other | 12442 | 6347 | 1056 | 61 | 12446 |
| pampa_PAMPA_INTRINSIC_exploratory | macro_500_700 | 1 | 1 | 1 | 1 | 1 |
| pampa_PAMPA_INTRINSIC_exploratory | nonmacro_500_700 | 11 | 10 | 9 | 4 | 11 |
| pampa_PAMPA_INTRINSIC_exploratory | nonmacro_gt700 | 5 | 5 | 4 | 1 | 5 |
| pampa_PAMPA_INTRINSIC_exploratory | nonmacro_le500 | 143 | 114 | 65 | 10 | 143 |
| pampa_PAMPA_INTRINSIC_exploratory_nonstrict | nonmacro_le500 | 9 | 9 | 8 | 1 | 9 |
| pampa_PAMPA_PAPP_exploratory | macro_500_700 | 98 | 88 | 59 | 15 | 102 |
| pampa_PAMPA_PAPP_exploratory | nonmacro_500_700 | 30 | 27 | 17 | 13 | 30 |
| pampa_PAMPA_PAPP_exploratory | nonmacro_gt700 | 45 | 35 | 22 | 6 | 45 |
| pampa_PAMPA_PAPP_exploratory | nonmacro_le500 | 332 | 260 | 151 | 52 | 333 |
| pampa_PAMPA_PAPP_exploratory | other | 155 | 148 | 101 | 17 | 171 |
| pampa_PAMPA_PEFF_exploratory | macro_500_700 | 1256 | 1256 | 808 | 3 | 2506 |
| pampa_PAMPA_PEFF_exploratory | other | 1859 | 1859 | 912 | 3 | 3675 |

500<MW≤700 的 50 Da 箱按分子图计数：

| 集合 | 分组 | 箱 | 分子图 | 标签 |
| --- | --- | --- | ---: | ---: |
| caco2_ab_main | macro_500_700 | 500_550 | 82 | 87 |
| caco2_ab_main | macro_500_700 | 550_600 | 64 | 71 |
| caco2_ab_main | macro_500_700 | 600_650 | 45 | 49 |
| caco2_ab_main | macro_500_700 | 650_700 | 31 | 31 |
| caco2_ab_main | nonmacro_500_700 | 500_550 | 329 | 343 |
| caco2_ab_main | nonmacro_500_700 | 550_600 | 137 | 138 |
| caco2_ab_main | nonmacro_500_700 | 600_650 | 76 | 85 |
| caco2_ab_main | nonmacro_500_700 | 650_700 | 66 | 78 |

## 来源与旧候选数

2026-11-19 的候选数（Caco-2 非大环 500–700 为 600，大环 207，>700 为 88；PAMPA 为 764、1,481、144）是修订前的工作集汇总。当时把组内极差超过 0.5 的记录一律暂缓，并把不同 PAMPA 系数放在一起计数。本版只在无法确认同一实验协议时暂缓冲突，并把 Peff、intrinsic、apparent 和未指明系数分开。因此正式分子图数不要求与那些数字相同。

Caco-2 主集、500<MW≤700 非大环最大原始研究是 DOI:10.1021/acs.jmedchem.3c00795，33 个分子图，占该组按研究计数的分子图合计 5.4%。同一分子图可以出现在多个研究中，因此这个合计可以大于唯一分子图数。

该组按分子图数排列的前五个原始研究：

- DOI:10.1021/acs.jmedchem.3c00795：33 个分子图
- DOI:10.1021/acsmedchemlett.2c00502：31 个分子图
- DOI:10.1016/j.bmcl.2009.03.124：18 个分子图
- DOI:10.1021/acs.jmedchem.4c02062：15 个分子图
- DOI:10.1016/j.bmc.2021.116163：14 个分子图

该组精确标签 644 个，不同测量记录 647 条，数据库出现 652 次，分子图 608 个，连接组 408 个，Murcko 骨架 408 个，FP06 分量 185 个，原始研究 185 个。连接组、Murcko 骨架和 FP06 分量分别保存。

MW>700 非大环精确标签只标记为外推候选，本阶段不封存、不划分、不抽样。

## 结构与三维资格

未改结构继承 2026-11-18 的 Graph_ID、Connectivity_ID（Scaffold_connectivity_group）、Murcko 和 FP06。FP06 是 Morgan 半径 2、2048 位、不含手性、Tanimoto≥0.6 的单链接分量，并包含骨架/连接约束；当时 RDKit 为 2026.03.3。本次环境的 RDKit 是 2025.09.6，只用于确认缓存与 34,150 个 SMILES 一致，没有重算分组。stage3 将含 Pb 的错误图改指到已有 MZ3 图，没有新增立体异构体。

初步三维资格只表示：结构对应关系明确，且继承的未指定手性中心数为 0。它不是构象生成或力场检查。手性未知的记录仍可留在二维主集。
Caco-2 主集 500<MW≤700 非大环中，手性标记完整的分子图 517 个，手性未指定的 91 个。初步三维资格按手性标记完整计数，为 517 个图。

## 描述符

下表使用冻结分子表中已有的 MW、cLogP、HBD、HBA，以及标签 y。TPSA 和可旋转键不在冻结表中，留到阶段B。

| 集合 | 分组 | n | MW中位 | cLogP中位 | TPSA | HBD中位 | HBA中位 | RotB | y中位 |
| --- | --- | ---: | ---: | ---: | --- | ---: | ---: | --- | ---: |
| caco2_ab_main | macro_500_700 | 222 | 568.2 | 3.73 | 未在冻结表 | 3.0 | 7.0 | 未在冻结表 | -5.824 |
| caco2_ab_main | nonmacro_500_700 | 608 | 544.9 | 4.14 | 未在冻结表 | 2.0 | 7.0 | 未在冻结表 | -5.506 |
| caco2_ab_main | nonmacro_gt700 | 106 | 820.4 | 5.40 | 未在冻结表 | 3.0 | 11.0 | 未在冻结表 | -6.432 |
| caco2_ab_main | nonmacro_le500 | 1740 | 417.5 | 3.32 | 未在冻结表 | 2.0 | 5.0 | 未在冻结表 | -5.036 |

## 暂缓与排除

| 决定 | 出现记录 |
| --- | ---: |
| IN_EXACT_LABEL | 65710 |
| HELD | 1059 |
| HELD_AFTER_AGGREGATION | 2 |

输入出现记录 66771 条，决定表 66771 条，两者相同。没有为 CB30/CB31 新增行。

暂缓原因按出现记录计数。一条记录可以同时有多个原因，所以各项之和可以大于暂缓记录数。

| 原因 | 出现记录 |
| --- | ---: |
| BAS_SHARED_STRUCTURE_DISCORDANT_VALUES | 472 |
| LOG10_CM_S_GT_MINUS3 | 446 |
| RZEPIELA_MINUS8_CENSORING_UNRESOLVED | 96 |
| NAMED_SAMPLE_CONNECTIVITY_CONFLICT | 23 |
| NAMED_SAMPLE_STEREO_ASSIGNMENT_AMBIGUOUS | 14 |
| STAGE2_HOLD:ORIGINAL_CENSORED_POINT_OR_POTENTIAL_CROSS_SOURCE_COPY_REVIEW | 5 |
| STAGE2_HOLD:PRIMARY_LOW_RECOVERY_AND_CENSORED_RESULT_REVIEW | 3 |
| AGGREGATION_CONFLICT | 2 |
| POINT_TO_PRIMARY_BOUND | 2 |
| QUALIFIER_NOT_EXACT | 2 |

## 进入阶段B需要的信息

阶段B可以从 `label_table.csv.gz` 的 `analysis_set=caco2_ab_main` 取开发标签，用 `cohort=nonmacro_gt700` 标记外推候选，用 `connectivity_id` 和 `fp06_component` 做隔离。PAMPA 若辅助任务为空，则各 `pampa_*_exploratory` 集合都不要并入 Caco-2。本阶段没有生成折、预算或封存成员。

## 参考文献与本项目规则的差别

Esaki 等，Journal of Cheminformatics，2024，DOI 10.1186/s13321-024-00826-z：本项目沿用其按终点、单位、方向和细胞体系保留已审核记录的做法。Feng 等，Scientific Data，2025，DOI 10.1038/s41597-024-04302-z：本项目参考其对大环、终点类别、单位和限定符的区分；大环阈值使用本项目已有的最大计算环≥12。Wang 等，Journal of Cheminformatics，2025，DOI 10.1186/s13321-025-00947-z：参考重复测量和 RDKit 标准化的必要性；组内中位数是本项目预先规定的标签，不是该文的原句。Pharmaceutics，2022，DOI 10.3390/pharmaceutics14101998：该文讨论测量变异和标准差，本项目的 0.5 是 log 极差核查规则，二者不是同一个阈值。Fourches 等，JCIM，2010，DOI 10.1021/ci100176x：结构和重复记录需要可追溯整理。CycPeptMPDB，JCIM，2023，DOI 10.1021/acs.jcim.2c01573：立体化学、实验条件和检测限要分开保存。Steinbauer 等，JCIM，2024，DOI 10.1021/acs.jcim.4c01083：相同实验条件才适合放在一起比较；本项目据此不把缺失条件自动视为兼容。

## 复核

独立复核脚本是 `scripts/data_v1/verify_data_v1.py`。它重新检查出现记录守恒、17 条不再进入 Caco-2、单位没有二次取对数、不一致的 BAS 副本没有进入精确标签、测量没有按数据库副本重复计数，以及分组标识与分子表一致。结果写在 `verification_result.json`。本次结果为 PASS：Caco-2 主集精确标签 3,224 条，分子图 2,916 个。
