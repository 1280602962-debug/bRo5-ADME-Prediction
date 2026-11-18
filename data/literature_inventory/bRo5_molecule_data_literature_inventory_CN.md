# bRo5分子类别、数据覆盖与待核实原文清单

日期：2026-11-18。按最新第三阶段审核明细重新计数。

## 统计口径

当前有34,140个规范结构图表示、66,761条点值来源出现记录。它们是技术候选，不是全量原文审核通过的最终化学实体；同论文同名样品尚有结构冲突。下面的分子数量均沿用这个口径。仅有界限的数据不在本次点值分母中；不能将本表当作全部结构资源或全部删失数据普查。

## 1. 按分子量与环型精确分成8组

| 类别 | 结构数 | Caco-2全部方向 | 其中AB候选 | PAMPA一般型 | MDCK |
| --- | --- | --- | --- | --- | --- |
| MW_500_700__MACRO_SCREEN | 2641 | 373 | 271 | 2208 | 84 |
| MW_500_700__NONMACRO_SCREEN | 4883 | 2168 | 1372 | 1893 | 875 |
| MW_700_1000__MACRO_SCREEN | 4267 | 317 | 153 | 4004 | 59 |
| MW_700_1000__NONMACRO_SCREEN | 662 | 314 | 204 | 273 | 97 |
| MW_GT1000__MACRO_SCREEN | 3105 | 1033 | 10 | 2402 | 43 |
| MW_GT1000__NONMACRO_SCREEN | 129 | 59 | 32 | 66 | 5 |
| MW_LE500__MACRO_SCREEN | 2434 | 317 | 169 | 1995 | 244 |
| MW_LE500__NONMACRO_SCREEN | 16019 | 5906 | 3934 | 6465 | 3999 |

分类互斥，8行结构数可相加得到34,140。每行终点人数可能重叠，不能横向相加。MW边界：≤500、500<MW≤700、700<MW≤1000、>1000 Da。MACRO_SCREEN表示已有RDKit最大环≥12的算法筛查，不是人工大环归类；NONMACRO_SCREEN也不等于非肽或PROTAC。

MW>500 Da共15687个：大环筛查10013个、非大环筛查5674个。MW>700 Da非大环为791个，其中662个在700–1000 Da、129个>1000 Da。MW>500只是尺寸分层，不是完整bRo5定义。

## 2. 药化类型目前可以给出的来源注释

| 来源注释集合 | 结构数 | Caco-2 | PAMPA一般型 | MDCK | RRCK |
| --- | --- | --- | --- | --- | --- |
| PROTAC_SOURCE | 126 | 61 | 68 | 1 | 0 |
| CYCPEPT_SOURCE | 7613 | 1242 | 6859 | 68 | 165 |
| NPMMP_SOURCE | 4329 | 551 | 3543 | 291 | 2 |
| OTHER_SOURCES_ONLY | 22122 | 8670 | 8886 | 5046 | 0 |

PROTAC_SOURCE来自PROTACDB4、Hudec新PROTAC整理或Kim数据的来源注释；CYCPEPT_SOURCE来自CycPept当前及历史快照；NPMMP_SOURCE来自NPMMP当前及存档。OTHER_SOURCES_ONLY表示没有前三种来源标签，不是已证实全部属于常规小分子。三个来源集合有重叠；目前也可能有仅在ChEMBL或MolMeDB中出现但未补药化类型注释的PROTAC/肽，故这些数字不是全项目最终互斥药化分类数。每组的终点覆盖使用合并表中该组结构的全部保留标签。

目前准确的类别是第1表的8组。药化上包含PROTAC、环肽、非肽/半肽大环及其他结构，但还不能提供每种药化类别的全量最终合格数量。

## 3. 每个结构实际有哪些数据

- 结构层：规范异构SMILES、MW、环型筛查、Murcko骨架/连接分组及未指定手性中心标记。9,459个结构含未指定手性中心，不能直接当成唯一立体构型的3D输入。
- 标签层：实验终点家族、方向/任务、原始报告值、换算后的log10(cm/s)、关系符号、来源数据库及来源行、论文/专利标识。每个结构至少有一个保留点值；不保证具有多个实验终点。
- 条件层：可获取的实验描述和pH、转运体/抑制剂等提示；缺失仍保留UNKNOWN。不同MDCK体系、方向、PAMPA膜/协议以及几何平均不能自动合并。
- 本轮核实的统一任务是膜渗透性。溶解度、代谢稳定性、血浆蛋白结合率及口服生物利用度尚未形成全体统一审核的标签集，不能将这些分子称作完整ADME数据集。部分原文含额外性质，不等于项目已经统一提取通过。

| 终点 | 全部MW结构数 | MW>500结构数 |
| --- | --- | --- |
| Caco-2 | 10487 | 4264 |
| PAMPA | 19306 | 10846 |
| MDCK | 5406 | 1163 |
| RRCK | 167 | 165 |
| LLC-PK1 | 405 | 81 |

PAMPA一般型为19,306个；若连同单列的intrinsic/apparent子类，结构并集为19,375个。方向可识别的Caco-2 AB候选为6,145个，含根据实验描述自动识别的记录，尚非全部原文核实。所有终点均可能交叉。

## 4. 原文到底还需下载、阅读哪些

下表35篇是当前优先清单：31篇关联高分子量非大环测试/校准候选；另4篇关联已定位的同名结构冲突。五篇已有局部原文证据，仍需补SI/结构或表格，不等于全文和全部记录均审核完成。其余多数需完整正文加全部相关SI；下载到文件和阅读审核是两个状态。

| 编号 | 原文 | 当前数据类型 | 还要核实 |
| --- | --- | --- | --- |
| 1 | Hepatitis C virus NS5A replication complex inhibitors: the discovery of daclatasvir.；[10.1021/jm401836p](https://doi.org/10.1021/jm401836p) | Caco-2 / PAMPA | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 2 | Understanding and Improving the Membrane Permeability of VH032-Based PROTACs.；[10.1021/acsmedchemlett.0c00265](https://doi.org/10.1021/acsmedchemlett.0c00265) | PAMPA | 补SI方法、原始结构图/作者SMILES和PAMPA完整数值；核对CM9、CM11连接长度及CMP98立体信息；正文已有定点取证。 |
| 3 | The Antioxidant Additive Approach for Alzheimer's Disease Therapy: New Ferulic (Lipoic) Acid Plus Melatonin Modified Tacrines as Cholinesterases Inhibitors, Direct Antioxidants, and Nuclear Factor (Erythroid-Derived 2)-Like 2 Activators.；[10.1021/acs.jmedchem.6b01178](https://doi.org/10.1021/acs.jmedchem.6b01178) | PAMPA | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 4 | Discovery and Optimization of Thiazolidinyl and Pyrrolidinyl Derivatives as Inhaled PDE4 Inhibitors for Respiratory Diseases.；[10.1021/acs.jmedchem.7b01044](https://doi.org/10.1021/acs.jmedchem.7b01044) | Caco-2 | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 5 | In Vitro ADME Profiling of PROTACs: Successes, Challenges, and Lessons Learned from Analysis of Clinical PROTACs from a Diverse Physicochemical Space.；[10.1021/acs.jmedchem.5c00358](https://doi.org/10.1021/acs.jmedchem.5c00358) | Caco-2 | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 6 | Development of an Orally Bioavailable LCK PROTAC Degrader as a Potential Therapeutic Approach to T-Cell Acute Lymphoblastic Leukemia.；[10.1021/acs.jmedchem.4c00481](https://doi.org/10.1021/acs.jmedchem.4c00481) | Caco-2 | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 7 | The cyclimids: Degron-inspired cereblon binders for targeted protein degradation；[10.1016/j.chembiol.2024.01.003](https://doi.org/10.1016/j.chembiol.2024.01.003) | Caco-2 | 补Supplementary Figure10g及其source data：逐项数值、单位、回收率、作者化合物名；核实PROTACDB记录级引文与具体测定来源。 |
| 8 | Linker-Dependent Folding Rationalizes PROTAC Cell Permeability.；[10.1021/acs.jmedchem.2c00877](https://doi.org/10.1021/acs.jmedchem.2c00877) | Caco-2 / PAMPA | 补SI的Caco-2双向原始表和作者结构CSV；核对PAMPA化合物3上界、回收率及同一分子的立体表示；Table2和几何平均定义已取证。 |
| 9 | Structure-Guided Discovery of Selective USP7 Inhibitors with In Vivo Activity.；[10.1021/acs.jmedchem.4c01472](https://doi.org/10.1021/acs.jmedchem.4c01472) | Caco-2 | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 10 | Diastereomeric Branched-Ester dBET1 Analogs Exhibit Conformation-Dependent Differences in Passive Membrane Permeability；[10.1021/acs.jmedchem.5c02791](https://doi.org/10.1021/acs.jmedchem.5c02791) | PAMPA | 补作者结构CSV/SI，逐原子核对dBET1、1、2a、2b及CRBN位点样品状态；Table1四个PAMPA均值和协议已确认。 |
| 11 | Neurotensin(8-13) analogs as dual NTS1 and NTS2 receptor ligands with enhanced effects on a mouse model of Parkinson's disease.；[10.1016/j.ejmech.2023.115386](https://doi.org/10.1016/j.ejmech.2023.115386) | PAMPA | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 12 | Design, Synthesis, and Biological Evaluations of DOT1L Peptide Mimetics Targeting the Protein-Protein Interactions between DOT1L and MLL-AF9/MLL-ENL.；[10.1021/acs.jmedchem.2c00083](https://doi.org/10.1021/acs.jmedchem.2c00083) | Caco-2 | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 13 | Discovery and Optimization of the First Highly Effective and Orally Available Galectin-3 Inhibitors for Treatment of Fibrotic Disease.；[10.1021/acs.jmedchem.2c00660](https://doi.org/10.1021/acs.jmedchem.2c00660) | Caco-2 | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 14 | Discovery of an orally active VHL-recruiting PROTAC that achieves robust HMGCR degradation and potent hypolipidemic activity in vivo；[10.1016/j.apsb.2020.11.001](https://doi.org/10.1016/j.apsb.2020.11.001) | Caco-2 | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 15 | Dual antagonists of α5β1/αvβ1 integrin for airway hyperresponsiveness.；[10.1016/j.bmcl.2020.127578](https://doi.org/10.1016/j.bmcl.2020.127578) | Caco-2 | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 16 | Discovery of the thieno[2,3-d]pyrimidine-2,4-dione derivative 21a: A potent and orally bioavailable gonadotropin-releasing hormone receptor antagonist.；[10.1016/j.ejmech.2022.114679](https://doi.org/10.1016/j.ejmech.2022.114679) | Caco-2 | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 17 | Discovery of a highly potent, selective, and stable d-amino acid-containing peptide inhibitor of CDK9/cyclin T1 interaction for the treatment of prostate cancer.；[10.1016/j.ejmech.2025.117248](https://doi.org/10.1016/j.ejmech.2025.117248) | Caco-2 | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 18 | Design and characterization of a heterobifunctional degrader of KEAP1；[10.1016/j.redox.2022.102552](https://doi.org/10.1016/j.redox.2022.102552) | PAMPA | 补SI TableS2：化合物8、14的数值、关系符号、单位、回收率；正文PAMPA方法已有定点取证。 |
| 19 | Rational Design of Sulfonyl-γ-AApeptides as Highly Potent HIV-1 Fusion Inhibitors with Broad-Spectrum Activity.；[10.1021/acs.jmedchem.3c01412](https://doi.org/10.1021/acs.jmedchem.3c01412) | PAMPA | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 20 | Rational In Silico Design of Selective TMPRSS6 Peptidomimetic Inhibitors via Exploitation of the S2 Subpocket.；[10.1021/acs.jmedchem.4c00922](https://doi.org/10.1021/acs.jmedchem.4c00922) | PAMPA | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 21 | Identification of a Potent and Selective CDK9 Degrader as a Targeted Therapeutic Option for the Treatment of Small-Cell Lung Cancer.；[10.1021/acs.jmedchem.4c01621](https://doi.org/10.1021/acs.jmedchem.4c01621) | Caco-2 | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 22 | Potent, Efficacious, and Stable Cyclic Opioid Peptides with Long Lasting Antinociceptive Effect after Peripheral Administration.；[10.1021/acs.jmedchem.9b01963](https://doi.org/10.1021/acs.jmedchem.9b01963) | PAMPA | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 23 | Synthesis and anti-tumor activity of menthylamine derivatives: Focusing on the potent inhibitory effect of compound W8 on glioma growth.；[10.1016/j.bmc.2024.118024](https://doi.org/10.1016/j.bmc.2024.118024) | Caco-2 | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 24 | Calcitonin gene-related peptide (CGRP) receptor antagonists: Heterocyclic modification of a novel azepinone lead.；[10.1016/j.bmcl.2021.128077](https://doi.org/10.1016/j.bmcl.2021.128077) | Caco-2 / PAMPA | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 25 | Identification of substituted 5-membered heterocyclic compounds as potential anti-leukemic agents.；[10.1016/j.ejmech.2018.12.059](https://doi.org/10.1016/j.ejmech.2018.12.059) | Caco-2 | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 26 | Rational design and biological evaluation of a new class of thiazolopyridyl tetrahydroacridines as cholinesterase and GSK-3 dual inhibitors for Alzheimer's disease.；[10.1016/j.ejmech.2020.112751](https://doi.org/10.1016/j.ejmech.2020.112751) | PAMPA | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 27 | Identification of novel and potent PROTACs targeting FAK for non-small cell lung cancer: Design, synthesis, and biological study.；[10.1016/j.ejmech.2022.114373](https://doi.org/10.1016/j.ejmech.2022.114373) | Caco-2 | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 28 | Discovery and characterization of novel potent BCR-ABL degraders by conjugating allosteric inhibitor.；[10.1016/j.ejmech.2022.114810](https://doi.org/10.1016/j.ejmech.2022.114810) | Caco-2 | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 29 | Novel gut-restricted bivalent agonists targeting mucosal 5-HT4R: Design, synthesis, and biological evaluation.；[10.1016/j.ejmech.2025.117425](https://doi.org/10.1016/j.ejmech.2025.117425) | PAMPA | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 30 | Discovery of the First Potent, Selective, and In Vivo Efficacious Polo-like Kinase 4 Proteolysis Targeting Chimera Degrader for the Treatment of TRIM37-Amplified Breast Cancer.；[10.1021/acs.jmedchem.3c00505](https://doi.org/10.1021/acs.jmedchem.3c00505) | Caco-2 | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 31 | Discovery of APH02174 as a Highly Selective and Orally Bioavailable IRAK4 Degrader for the Treatment of Inflammatory Diseases.；[10.1021/acs.jmedchem.5c00711](https://doi.org/10.1021/acs.jmedchem.5c00711) | Caco-2 | 完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。 |
| 32 | Discovery of macrocyclic covalent inhibitors for severe acute respiratory syndrome coronavirus 2 3CL protease；[10.1016/j.bmc.2024.117846](https://doi.org/10.1016/j.bmc.2024.117846) | Caco-2 | 完整正文+SI/作者SMILES，裁决同论文同编号样品的结构差异：25a；27g；27j |
| 33 | Structure-Permeability Relationship of Semipeptidic Macrocycles—Understanding and Optimizing Passive Permeability and Efflux Ratio；[10.1021/acs.jmedchem.0c00013](https://doi.org/10.1021/acs.jmedchem.0c00013) | Caco-2 / PAMPA | 完整正文+SI/作者SMILES，裁决同论文同编号样品的结构差异：39 |
| 34 | Cell Permeability of Isomeric Macrocycles: Predictions and NMR Studies；[10.1021/acsmedchemlett.1c00126](https://doi.org/10.1021/acsmedchemlett.1c00126) | Caco-2 | 完整正文+SI/作者SMILES，裁决同论文同编号样品的结构差异：6 |
| 35 | Property-Based Optimization of Hydroxamate-Based gamma-Lactam HDAC Inhibitors to Improve Their Metabolic Stability and Pharmacokinetic Profiles；[10.1021/jm3009376](https://doi.org/10.1021/jm3009376) | Caco-2 / PAMPA | 完整正文+SI/作者SMILES，裁决同论文同编号样品的结构差异：8f；8d |


## 5. 优先清单不是全部文献

整个当前点值候选层登记3607个来源键：DOI来源键3513、专利号12、CycPept作者年份别名52、文献引文别名29、未解析URI1。这些是来源标识，不是已经裁决的独立研究数，也不能直接相加为独立实验数。

all_source_identifiers_and_paper_links.csv完整列出来源标识、DOI/专利链接、涉及结构数、终点和优先级。35篇优先原文以外，其他未核实来源仍需按正式研究范围安排阅读；不能误读为只补完这35篇就足以确认全量34,140个结构全部合格。对于没有DOI的引文/作者年份别名，先恢复原始参考文献，再定位正文/SI。

## 文件说明

- molecular_category_endpoint_counts.csv：8组互斥结构与终点计数。
- source_tag_endpoint_counts.csv：来源注释集合及其终点覆盖。
- molecule_categories_and_available_data.csv：34,140个结构的SMILES、类别、来源和可用终点逐项清单。
- priority35_papers_and_exact_evidence_gaps.csv：35篇优先论文完整英文标题、DOI、原文链接和具体核查内容。
- all_source_identifiers_and_paper_links.csv：全项目3,607个来源键及优先级。
- inventory_summary.json：分母、集合重叠及手性缺失统计。

只做现有数据清点与待核实项整理，未扩大合格标准、未训练模型、未将局部原文核实当成全部通过。现有结构冲突仍保留待裁决。
