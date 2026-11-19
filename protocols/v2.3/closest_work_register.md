# 立项审查依据与近邻研究登记（2026-11-19）

对应V2.2固定提交 585dd140640589f8b21001421c22b171f903314f，本版根据前轮四套科研skills审查意见修订。以下为定向核查而非系统综述、竞争论文全量复现或独立同行评审。实际使用的是指导文件中的方法；第三方skills不是期刊官方标准，未执行第三方自动脚本。

## Skills证据

| 来源及读取内容 | Git blob SHA（此前审查实际读取） | 应用 |
| --- | --- | --- |
| [Academic-Research-Agent-Skill](https://github.com/ngtiendong/Academic-Research-Agent-Skill/blob/main/SKILL.md) | 707ea86e5c06c8c4d06e06ff94ae018c1776c0d4 | 区分创新和执行准入 |
| [Reality Gate](https://github.com/ngtiendong/Academic-Research-Agent-Skill/blob/main/references/reality_gate.md) | 668f7c2de08d3ea999256ea86983523668399246 | 数据、有效产出、最小试算、终止 |
| [Novelty Gate](https://github.com/ngtiendong/Academic-Research-Agent-Skill/blob/main/references/novelty_gate.md) | 6707cdcb2c501376c68ee37b665c397f75fdeca0 | 识别最近竞争研究与剩余贡献 |
| [Novelty Assessment](https://github.com/lingzhi227/agent-research-skills/blob/main/skills/novelty-assessment/SKILL.md) | 1a6577bbeafd3fb4698502bda124c1c4339b2586 | 多轮检索、重叠审查 |
| [Assessment prompts](https://github.com/lingzhi227/agent-research-skills/blob/main/skills/novelty-assessment/references/assessment-prompts.md) | b61063dcff24fb5ccc04d9b94a8b3b5a547c3da8 | 方法及应用创新区分 |
| [Scientific Critical Thinking](https://github.com/K-Dense-AI/scientific-agent-skills/blob/main/skills/scientific-critical-thinking/SKILL.md) v1.5 | 2a2fc151b5b13e4e11ece7410f21a14d5b25e30d | 独立单位、偏差、统计和解释 |
| [Experimental Design](https://github.com/K-Dense-AI/scientific-agent-skills/blob/main/skills/scientific-critical-thinking/references/experimental_design.md) | 4be76ae49eb6339e6a72c8d3ee69338c4ea3e347 | 对照、可证伪、精度 |
| [Research Idea Rubric](https://github.com/UCL-ERL/skills/blob/main/skills/ideation/research-idea-rubric/SKILL.md) | 1d2bb4012cad925910c1bf45757c05fef85c7c89 | 改用于计算药化的压力测试 |
| [Idea Stress Tests](https://github.com/UCL-ERL/skills/blob/main/skills/ideation/research-idea-rubric/references/idea-stress-tests.md) | b1500ab6ef9c63a1d956ee3086c2efe628d8b351 | 最便宜的证伪及停止条件 |

## 核心近邻研究（本次重新核对可访问出版页面/摘要和相关正文）

| 来源 | 已有贡献/重叠 | 本项目处理 |
| --- | --- | --- |
| Steinbauer, Lehr, Reichel. JCIM 2024, 64, 9163–9172. [DOI](https://doi.org/10.1021/acs.jcim.4c01083) | 针对Caco-2测试分子以Tanimoto选择训练集，用LightGBM及RDKit性质预测；局部改善有限，数据效率值得研究。 | 最高优先DS对标，新增同池ECFP-only。本文R为训练侧参考而非逐测试分子局部模型，终点口径亦不同；不可直接比较绝对误差/宣称复现。 |
| Hudec & Mercado. [arXiv:2610.07998](https://arxiv.org/abs/2610.07998), 2026-10-06，预印本 | 扩展PROTAC PAMPA数据并显示跨来源转移不对称，指出数据组成的重要性。 | H1不主张首次；具体差异为跨拓扑等预算、MW配平、来源隔离和构象交互。 |
| Kim et al. JCIM 2025, 65, 12563–12578. [DOI](https://doi.org/10.1021/acs.jcim.5c01600) | 异双功能降解剂构象3D描述符及MetaD渗透性研究。 | H2/H4不能单独主张新颖；检验选择用途与预测用途独立贡献。SI包含结构/渗透性相关字段，但需既有来源查重和实验语义审核。 |
| Peteani et al. Nat Commun 2024, 15, 5764. [DOI](https://doi.org/10.1038/s41467-024-49979-3) | 面向靶向降解剂的性质预测及模型适配研究。 | TL仅同数据控制辅助，不声称首次对降解剂迁移。 |
| Doi & Esaki. JCIM 2026, 66, 10326–10338. [DOI](https://doi.org/10.1021/acs.jcim.6c00597) | 五类渗透性实验多任务GCN及RF/ST/TL对照。 | C03条件性，非主要方法创新；SI数据需审计原始来源、标签及重叠。 |
| Yao et al. Commun Chem 2024, 7, 85. [DOI](https://doi.org/10.1038/s42004-024-01169-4) | 利用梯度量度源任务可迁移性、指导源数据集选择。 | 本研究为同终点跨拓扑分子供体排序而非源任务梯度图；不声称首次设计可迁移性度量。 |
| Brenner et al. Mol Pharm 2026, 23, 3608–3620. [DOI](https://doi.org/10.1021/acs.molpharmaceut.5c01908) | 3271个PROTAC涉及多种渗透/外排实验；条件影响可比性。 | 仅数据候选线索，3271不等于公开Caco-2 AB精确标签数。逐分子结构、AB/BA、单位/条件及许可可获得性尚未确认。 |

补充证据：Steinbauer代码链接 https://github.com/ffstghc/caco2ml （来自前轮核查，未在本次运行）；MoTSE 2022 https://pmc.ncbi.nlm.nih.gov/articles/PMC9579493/ 为源任务相似性背景，前轮已登记，本次未重新全文核查。

## 论文差异与失效条件

P01：受控跨拓扑数据组成差异，若来源依赖、误差不精确，仅报告适用范围；不能解释独立拓扑因果。H3：同样本3D交互，缺3D共同支撑则未估计。DS01：同池相同RF-2D条件的构象引导选择对照，若HYB不优于强2D对照或没有差异，不能宣称改善预测；若只有某供体获益，限定范围。TL：训练策略对照，若无改善报告负/不明确结果。所有模块均须有明确证据，不能用大量模拟代替核心预测对照。
投稿前补充最近竞争全文、SI、实验终点/划分/基线差异表及相关开放数据许可与来源去重。数据扩充只在正式冻结前依现有统一规则纳入；冻结后新来源独立另版，不事后改变留出。审查结论非首创性证明。
