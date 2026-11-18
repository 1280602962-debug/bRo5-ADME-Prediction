import csv,json,collections,re,html,zipfile,hashlib
from pathlib import Path
P=Path(__file__).resolve().parent;D=P/'data';O=P/'deliverables';O.mkdir(exist_ok=True)
def rows(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(n,rr):
 rr=list(rr)
 with (O/n).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
def md(h,rr):return '| '+' | '.join(h)+' |\n| '+' | '.join(['---']*len(h))+' |\n'+''.join('| '+' | '.join(str(v).replace('|',' / ') for v in r)+' |\n' for r in rr)
L=rows(D/'stage3/point_candidates_after_primary_overlay.csv');M={r['Graph_ID']:r for r in rows(D/'stage2/molecule_groups.csv')};gids={r['Graph_ID'] for r in L}
assert len(L)==66761 and len(gids)==34140
bygraph=collections.defaultdict(list)
for r in L:bygraph[r['Graph_ID']].append(r)
tags={'PROTAC_SOURCE':{'PROTACDB4','Hudec_new_parsed','Kim2025'},'CYCPEPT_SOURCE':{'CycPept_current','CycPeptMPDB_historical'},'NPMMP_SOURCE':{'NPMMP_archive','NPMMPD_official_current'}}
tagsets={t:{r['Graph_ID'] for r in L if r['Source'] in ss} for t,ss in tags.items()};covered=set().union(*tagsets.values());tagsets['OTHER_SOURCES_ONLY']=gids-covered
def category(g):
 m=M[g];mw=float(m['MW']);b='MW_LE500' if mw<=500 else 'MW_500_700' if mw<=700 else 'MW_700_1000' if mw<=1000 else 'MW_GT1000'
 return b+'__'+m['Topology']
def stats(gg):
 return dict(Structures=len(gg),Point_occurrences=sum(len(bygraph[g]) for g in gg),Caco2=len({g for g in gg if any(r['Endpoint']=='Caco-2' for r in bygraph[g])}),Caco2_AB_candidate=len({g for g in gg if any(r['Task']=='Caco-2:AB' for r in bygraph[g])}),PAMPA_general=len({g for g in gg if any(r['Endpoint']=='PAMPA' for r in bygraph[g])}),PAMPA_any=len({g for g in gg if any(r['Endpoint'].startswith('PAMPA') for r in bygraph[g])}),MDCK=len({g for g in gg if any(r['Endpoint']=='MDCK' for r in bygraph[g])}),RRCK=len({g for g in gg if any(r['Endpoint']=='RRCK' for r in bygraph[g])}),LLC_PK1=len({g for g in gg if any(r['Endpoint']=='LLC-PK1' for r in bygraph[g])}))
cells=collections.defaultdict(set)
for g in gids:cells[category(g)].add(g)
matrix=[dict(Category=k,**stats(v)) for k,v in sorted(cells.items())]
write('molecular_category_endpoint_counts.csv',matrix)
write('source_tag_endpoint_counts.csv',[dict(Source_tag=k,**stats(v),Classification_status='DATABASE_SOURCE_TAG_NOT_FINAL_MEDCHEM_MODALITY') for k,v in tagsets.items()])
write('molecule_categories_and_available_data.csv',[dict(Graph_ID=g,Canonical_isomeric_SMILES=M[g]['SMILES'],MW=M[g]['MW'],Category=category(g),Source_tags='|'.join(t for t,gs in tagsets.items() if g in gs),Sources='|'.join(sorted({r['Source'] for r in bygraph[g]})),Endpoint_families='|'.join(sorted({r['Endpoint'] for r in bygraph[g]})),Tasks='|'.join(sorted({r['Task'] for r in bygraph[g]})),Point_occurrences=len(bygraph[g]),Unassigned_chiral_centers=M[g]['Unassigned_chiral_centers'],Identity_status='GRAPH_REPRESENTATION_NOT_FINAL_CHEMICAL_ENTITY',Primary_numeric_confirmed=any(r['Primary_point_check']=='PRIMARY_TABLE_VALUE_AND_SCALE_CONFIRMED' for r in bygraph[g])) for g in sorted(gids)])
pri=rows(D/'stage3/priority1_evidence_gaps.csv');ev=json.loads((D/'stage3/primary_evidence.json').read_text());evmap={r['DOI']:r for r in ev}
fill={
'10.1016/j.chembiol.2024.01.003':'The cyclimids: Degron-inspired cereblon binders for targeted protein degradation',
'10.1021/acs.jmedchem.5c02791':'Diastereomeric Branched-Ester dBET1 Analogs Exhibit Conformation-Dependent Differences in Passive Membrane Permeability',
'10.1016/j.apsb.2020.11.001':'Discovery of an orally active VHL-recruiting PROTAC that achieves robust HMGCR degradation and potent hypolipidemic activity in vivo',
'10.1016/j.redox.2022.102552':'Design and characterization of a heterobifunctional degrader of KEAP1',
'10.1016/j.bmc.2024.117846':'Discovery of macrocyclic covalent inhibitors for severe acute respiratory syndrome coronavirus 2 3CL protease',
'10.1021/acs.jmedchem.0c00013':'Structure-Permeability Relationship of Semipeptidic Macrocycles—Understanding and Optimizing Passive Permeability and Efflux Ratio',
'10.1021/acsmedchemlett.1c00126':'Cell Permeability of Isomeric Macrocycles: Predictions and NMR Studies',
'10.1021/jm3009376':'Property-Based Optimization of Hydroxamate-Based gamma-Lactam HDAC Inhibitors to Improve Their Metabolic Stability and Pharmacokinetic Profiles'}
titles={r['DOI']:re.sub('<[^>]+>','',html.unescape(html.unescape(r['Title']))) for r in pri};titles.update(fill)
conf=rows(D/'stage3/source_sample_identity_conflicts.csv');confbydoi=collections.defaultdict(list)
for r in conf:confbydoi[r['Study_key'][4:]].append(r)
prioritydois=[r['DOI'] for r in pri];extra=sorted(set(confbydoi)-set(prioritydois))
assert len(prioritydois)==31 and len(extra)==4
bystudy=collections.defaultdict(list)
for r in L:bystudy[r['Study_key']].append(r)
def original_need(doi):
 special={
 '10.1021/acsmedchemlett.0c00265':'补SI方法、原始结构图/作者SMILES和PAMPA完整数值；核对CM9、CM11连接长度及CMP98立体信息；正文已有定点取证。',
 '10.1021/acs.jmedchem.2c00877':'补SI的Caco-2双向原始表和作者结构CSV；核对PAMPA化合物3上界、回收率及同一分子的立体表示；Table2和几何平均定义已取证。',
 '10.1021/acs.jmedchem.5c02791':'补作者结构CSV/SI，逐原子核对dBET1、1、2a、2b及CRBN位点样品状态；Table1四个PAMPA均值和协议已确认。',
 '10.1016/j.redox.2022.102552':'补SI TableS2：化合物8、14的数值、关系符号、单位、回收率；正文PAMPA方法已有定点取证。',
 '10.1016/j.chembiol.2024.01.003':'补Supplementary Figure10g及其source data：逐项数值、单位、回收率、作者化合物名；核实PROTACDB记录级引文与具体测定来源。'}
 if doi in special:return special[doi]
 if doi in confbydoi:return '完整正文+SI/作者SMILES，裁决同论文同编号样品的结构差异：'+'；'.join(r['Original_names'] for r in confbydoi[doi])
 return '完整正文+SI（PDF和作者CSV/XLSX/SMILES如有）：作者化合物编号、原始数值及不等号、单位/尺度、实验方向/膜或细胞体系、pH/抑制剂、回收率、系列和重复来源。'
papers=[]
for i,doi in enumerate(prioritydois+extra,1):
 rs=bystudy.get('DOI:'+doi,[]);high={r['Graph_ID'] for r in rs if float(M[r['Graph_ID']]['MW'])>700 and M[r['Graph_ID']]['Topology']=='NONMACRO_SCREEN'}
 papers.append(dict(No=i,Priority_scope='TEST_CALIBRATION_PRIORITY31' if doi in prioritydois else 'EXPLICIT_IDENTITY_CONFLICT_ADDITIONAL',DOI=doi,Title=titles[doi],Landing_page='https://doi.org/'+doi,Known_readable_primary=evmap.get(doi,{}).get('URL',''),Current_primary_review='PARTIAL_FIELDS_ONLY' if doi in evmap else 'FULL_PRIMARY_REVIEW_PENDING',All_MW_graphs=len({r['Graph_ID'] for r in rs}),HighMW_nonmacro_graphs=len(high),Current_endpoint_families='|'.join(sorted({r['Endpoint'] for r in rs})),Current_tasks='|'.join(sorted({r['Task'] for r in rs})),Need_to_obtain_and_check=original_need(doi)))
write('priority35_papers_and_exact_evidence_gaps.csv',papers)
catalog=[]
for key,rs in sorted(bystudy.items()):
 kind=key.split(':',1)[0];ident=key.split(':',1)[1] if ':' in key else key;gs={r['Graph_ID'] for r in rs};high={g for g in gs if float(M[g]['MW'])>700 and M[g]['Topology']=='NONMACRO_SCREEN'}
 catalog.append(dict(Source_key=key,Source_key_kind=kind,DOI_or_ID=ident,Title_if_known=titles.get(ident,''),Landing_page='https://doi.org/'+ident if kind=='DOI' else 'https://patents.google.com/patent/'+ident+'/en' if kind=='PATENT' else '',Priority='P1_TEST_CALIBRATION' if ident in prioritydois else 'P1_IDENTITY_CONFLICT' if ident in extra else 'P2_HIGHMW_NONMACRO' if high else 'OTHER_DOMAIN',Point_graphs=len(gs),HighMW_nonmacro_graphs=len(high),Point_occurrences=len(rs),Endpoints='|'.join(sorted({r['Endpoint'] for r in rs})),Primary_evidence='PARTIAL_FIELDS_ONLY' if ident in evmap else 'NOT_ESTABLISHED',Interpretation='NORMALIZED_SOURCE_IDENTIFIER_NOT_CONFIRMED_INDEPENDENT_STUDY'))
write('all_source_identifiers_and_paper_links.csv',catalog)
stage3=json.loads((D/'stage3/stage3_summary.json').read_text());coverage=collections.Counter(r['Source_key_kind'] for r in catalog)
outstats=dict(point_graphs=len(gids),point_occurrences=len(L),MW_GT500=sum(float(M[g]['MW'])>500 for g in gids),MW_GT500_NONMACRO=sum(float(M[g]['MW'])>500 and M[g]['Topology']=='NONMACRO_SCREEN' for g in gids),MW_GT500_MACRO=sum(float(M[g]['MW'])>500 and M[g]['Topology']=='MACRO_SCREEN' for g in gids),all_source_keys=len(catalog),source_key_kinds=dict(coverage),priority31=31,additional_conflict_papers=4,source_tag_sets_overlap={a+'__'+b:len(tagsets[a]&tagsets[b]) for a in tags for b in tags if a<b},stereo_unspecified_graphs=sum(int(M[g]['Unassigned_chiral_centers'])>0 for g in gids),limits='Not final chemical-entity or model-ready molecule counts; only the existing point-candidate layer; bound-only structures not included.')
(O/'inventory_summary.json').write_text(json.dumps(outstats,ensure_ascii=False,indent=2),encoding='utf-8')
assert sum(r['Structures'] for r in matrix)==34140
paperlines=md(['编号','原文','当前数据类型','还要核实'],[(r['No'],r['Title']+'；['+r['DOI']+']('+r['Landing_page']+')',r['Current_endpoint_families'],r['Need_to_obtain_and_check']) for r in papers])
report=f'''# bRo5分子类别、数据覆盖与待核实原文清单

日期：2026-11-18。按最新第三阶段审核明细重新计数。

## 统计口径

当前有34,140个规范结构图表示、66,761条点值来源出现记录。它们是技术候选，不是全量原文审核通过的最终化学实体；同论文同名样品尚有结构冲突。下面的分子数量均沿用这个口径。仅有界限的数据不在本次点值分母中；不能将本表当作全部结构资源或全部删失数据普查。

## 1. 按分子量与环型精确分成8组

{md(['类别','结构数','Caco-2全部方向','其中AB候选','PAMPA一般型','MDCK'],[(r['Category'],r['Structures'],r['Caco2'],r['Caco2_AB_candidate'],r['PAMPA_general'],r['MDCK']) for r in matrix])}
分类互斥，8行结构数可相加得到34,140。每行终点人数可能重叠，不能横向相加。MW边界：≤500、500<MW≤700、700<MW≤1000、>1000 Da。MACRO_SCREEN表示已有RDKit最大环≥12的算法筛查，不是人工大环归类；NONMACRO_SCREEN也不等于非肽或PROTAC。

MW>500 Da共{outstats['MW_GT500']}个：大环筛查{outstats['MW_GT500_MACRO']}个、非大环筛查{outstats['MW_GT500_NONMACRO']}个。MW>700 Da非大环为791个，其中662个在700–1000 Da、129个>1000 Da。MW>500只是尺寸分层，不是完整bRo5定义。

## 2. 药化类型目前可以给出的来源注释

{md(['来源注释集合','结构数','Caco-2','PAMPA一般型','MDCK','RRCK'],[(k,stats(v)['Structures'],stats(v)['Caco2'],stats(v)['PAMPA_general'],stats(v)['MDCK'],stats(v)['RRCK']) for k,v in tagsets.items()])}
PROTAC_SOURCE来自PROTACDB4、Hudec新PROTAC整理或Kim数据的来源注释；CYCPEPT_SOURCE来自CycPept当前及历史快照；NPMMP_SOURCE来自NPMMP当前及存档。OTHER_SOURCES_ONLY表示没有前三种来源标签，不是已证实全部属于常规小分子。三个来源集合有重叠；目前也可能有仅在ChEMBL或MolMeDB中出现但未补药化类型注释的PROTAC/肽，故这些数字不是全项目最终互斥药化分类数。每组的终点覆盖使用合并表中该组结构的全部保留标签。

目前准确的类别是第1表的8组。药化上包含PROTAC、环肽、非肽/半肽大环及其他结构，但还不能提供每种药化类别的全量最终合格数量。

## 3. 每个结构实际有哪些数据

- 结构层：规范异构SMILES、MW、环型筛查、Murcko骨架/连接分组及未指定手性中心标记。9,459个结构含未指定手性中心，不能直接当成唯一立体构型的3D输入。
- 标签层：实验终点家族、方向/任务、原始报告值、换算后的log10(cm/s)、关系符号、来源数据库及来源行、论文/专利标识。每个结构至少有一个保留点值；不保证具有多个实验终点。
- 条件层：可获取的实验描述和pH、转运体/抑制剂等提示；缺失仍保留UNKNOWN。不同MDCK体系、方向、PAMPA膜/协议以及几何平均不能自动合并。
- 本轮核实的统一任务是膜渗透性。溶解度、代谢稳定性、血浆蛋白结合率及口服生物利用度尚未形成全体统一审核的标签集，不能将这些分子称作完整ADME数据集。部分原文含额外性质，不等于项目已经统一提取通过。

{md(['终点','全部MW结构数','MW>500结构数'],[(ep,len({r['Graph_ID'] for r in L if r['Endpoint']==ep}),len({r['Graph_ID'] for r in L if r['Endpoint']==ep and float(M[r['Graph_ID']]['MW'])>500})) for ep in ['Caco-2','PAMPA','MDCK','RRCK','LLC-PK1']])}
PAMPA一般型为19,306个；若连同单列的intrinsic/apparent子类，结构并集为19,375个。方向可识别的Caco-2 AB候选为6,145个，含根据实验描述自动识别的记录，尚非全部原文核实。所有终点均可能交叉。

## 4. 原文到底还需下载、阅读哪些

下表35篇是当前优先清单：31篇关联高分子量非大环测试/校准候选；另4篇关联已定位的同名结构冲突。五篇已有局部原文证据，仍需补SI/结构或表格，不等于全文和全部记录均审核完成。其余多数需完整正文加全部相关SI；下载到文件和阅读审核是两个状态。

{paperlines}

## 5. 优先清单不是全部文献

整个当前点值候选层登记{len(catalog)}个来源键：DOI来源键{coverage['DOI']}、专利号{coverage['PATENT']}、CycPept作者年份别名{coverage['CYC_SOURCE']}、文献引文别名{coverage['CITATION']}、未解析URI{coverage['UNRESOLVED_URI']}。这些是来源标识，不是已经裁决的独立研究数，也不能直接相加为独立实验数。

all_source_identifiers_and_paper_links.csv完整列出来源标识、DOI/专利链接、涉及结构数、终点和优先级。35篇优先原文以外，其他未核实来源仍需按正式研究范围安排阅读；不能误读为只补完这35篇就足以确认全量34,140个结构全部合格。对于没有DOI的引文/作者年份别名，先恢复原始参考文献，再定位正文/SI。

## 文件说明

- molecular_category_endpoint_counts.csv：8组互斥结构与终点计数。
- source_tag_endpoint_counts.csv：来源注释集合及其终点覆盖。
- molecule_categories_and_available_data.csv：34,140个结构的SMILES、类别、来源和可用终点逐项清单。
- priority35_papers_and_exact_evidence_gaps.csv：35篇优先论文完整英文标题、DOI、原文链接和具体核查内容。
- all_source_identifiers_and_paper_links.csv：全项目3,607个来源键及优先级。
- inventory_summary.json：分母、集合重叠及手性缺失统计。

只做现有数据清点与待核实项整理，未扩大合格标准、未训练模型、未将局部原文核实当成全部通过。现有结构冲突仍保留待裁决。
'''
(O/'bRo5_molecule_data_literature_inventory_CN.md').write_text(report,encoding='utf-8')
write('input_sha256.csv',[dict(Path=str(p),Bytes=p.stat().st_size,SHA256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in [D/'stage3/point_candidates_after_primary_overlay.csv',D/'stage2/molecule_groups.csv',D/'stage3/priority1_evidence_gaps.csv',D/'stage3/primary_evidence.json',D/'stage3/source_sample_identity_conflicts.csv']])
zpath=O/'bRo5_molecule_data_literature_inventory_20261118.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(O.iterdir()):
  if p.suffix in ['.csv','.md','.json']:z.write(p,p.name)
 z.write(Path(__file__),'summarize_inventory.py')
with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
print(json.dumps(dict(summary=outstats,category_counts=matrix,zip=str(zpath),zip_bytes=zpath.stat().st_size),ensure_ascii=False,indent=2))
