import csv,json,re,collections,hashlib,zipfile,math,sys
from pathlib import Path
from rdkit import Chem
sys.stdout.reconfigure(encoding='utf-8')
P=Path(__file__).resolve().parent;S=P.parent/'bro5_stage2';B=Path(r'C:\Users\guwenjing\bRo5 ADME Prediction\bRo5_Permeability_Benchmark')
def rows(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(n,rr):
 with (P/n).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
def md(h,rr):return '| '+' | '.join(h)+' |\n| '+' | '.join(['---']*len(h))+' |\n'+''.join('| '+' | '.join(map(str,r))+' |\n' for r in rr)
L=rows(P/'labels_primary_overlay.csv');before=rows(S/'labels_task_source_audit.csv');M={r['Graph_ID']:r for r in rows(S/'molecule_groups.csv')};beforeid={r['Occurrence_ID']:r for r in before};counts=rows(P/'stage3_counts.csv')
assert len(L)==len(before)==66771 and len({r['Occurrence_ID'] for r in L})==66771
assert {r['Occurrence_ID'] for r in L}==set(beforeid)
points=[r for r in L if r['Stage3_status']=='POINT_CANDIDATE'];other=[r for r in L if r['Stage3_status']!='POINT_CANDIDATE']
assert len(points)+len(other)==66771 and len(points)==66761
assert len({r['Graph_ID'] for r in points})==34140
bound=rows(P/'new_primary_bound_occurrences.csv');assert len(bound)==2 and len({r['Reported_result_ID'] for r in bound})==1
assert all(r['Relation']=='<' and float(r['Log10_bound_cm_s'])==-7.37 for r in bound)
assert len(rows(P/'primary_numeric_protocol_confirmed_means.csv'))==4
changed=[]
for r in L:
 old=beforeid[r['Occurrence_ID']]
 if (r['Graph_ID'],r['Task'],r['Log10_cm_s'],r['Relation'])!=(old['Graph_ID'],old['Task'],old['Log10_cm_s'],old['Relation']):changed.append(r['Occurrence_ID'])
checks=[]
for mode,field in [('SC','Scaffold_connectivity_group'),('FP06','Scaffold_FP06_component')]:
 for identity in ['Graph_ID','Study_key',field]:
  ps=collections.defaultdict(set)
  for r in points:
   part=r['Partition_'+mode]
   if part not in ['TRAIN','CALIBRATION','TEST']:continue
   val=M[r['Graph_ID']][field] if identity==field else r[identity]
   ps[val].add(part)
  n=sum(len(v)>1 for v in ps.values())
  checks.append(dict(Split_mode=mode,Identity=identity,Cross_partition_conflicts=n,Status='PASS' if n==0 else 'NEEDS_REPARTITION_REVIEW'))
write('post_correction_leakage_checks.csv',checks)
assert all(r['Cross_partition_conflicts']==0 for r in checks)
for row in counts:
 def within(g):
  m=M[g];mw=float(m['MW']);non=m['Topology']=='NONMACRO_SCREEN';scope=row['Scope']
  return scope=='ALL_MW' or scope=='MW_GT500' and mw>500 or scope=='MW_GT700_NONMACRO' and mw>700 and non or scope=='MW_700_1000_NONMACRO' and 700<mw<=1000 and non
 rr=[r for r in points if within(r['Graph_ID']) and (row['Task']=='ALL_TASKS' or r['Task']==row['Task'])]
 assert len(rr)==int(row['Point_occurrences']) and len({r['Graph_ID'] for r in rr})==int(row['Graph_representations'])
ver=dict(status='PASS',point_occurrences=len(points),point_graphs=len({r['Graph_ID'] for r in points}),nonpoint_or_held_occurrences=len(other),modified_occurrences=len(changed),checks='Occurrence conservation, exact recount, bound inequality reversal, point removal, atom-fragment correction, post-correction graph/source/group partition separation',limits='Primary numeric verification is distinct from full atom identity, replicate-event reconstruction and full scheme completion.')
(P/'stage3_independent_verification.json').write_text(json.dumps(ver,ensure_ascii=False,indent=2),encoding='utf-8')
studies=rows(P/'study_review_inventory.csv');docs={r['doi'].lower():r for r in rows(B/'03_extraction/chembl/chembl_37/DOCUMENTS_FULL.csv') if r['doi']};evidence=json.loads((P/'primary_evidence.json').read_text(encoding='utf-8'));ev={r['DOI']:r for r in evidence};access=[]
for r in studies:
 if r['Priority']!='1':continue
 doi=r['Study_key'][4:];doc=docs.get(doi,{});e=ev.get(doi,{})
 access.append(dict(DOI=doi,Title=doc.get('title','METADATA_FROM_WEB_SEARCH_OR_DATABASE_ALIAS'),Structures=r['Structures'],Tasks=r['Tasks'],Primary_URL=e.get('URL','https://doi.org/'+doi),Targeted_content_review='YES_LIMITED_FIELDS' if e else 'NOT_COMPLETED',Evidence_ID=e.get('Evidence_ID',''),Required_next='Complete source sample and atom identity, primary numerical table, qualifiers, protocol, recovery and measurement-level source lineage',Blocking_status='SEE_EVIDENCE_UNCONFIRMED_FIELDS' if e else 'FULLTEXT_SI_NOT_REVIEWED_DO_NOT_UPGRADE',Access_note='Publisher/SI retrieval returned403 or unavailable for several tested ACS URLs; PMC pages sometimes returned verification pages. This is not evidence that all31 studies are inaccessible. All31 had metadata searches; no full-review claim.'))
write('priority1_evidence_gaps.csv',access)
split=[]
for task in ['Caco-2:AB','PAMPA:DIRECTION_NOT_APPLICABLE']:
 for part in ['TRAIN','CALIBRATION','TEST']:
  rs=[r for r in points if r['Task']==task and r['Partition_SC']==part and float(M[r['Graph_ID']]['MW'])>700 and M[r['Graph_ID']]['Topology']=='NONMACRO_SCREEN']
  split.append(dict(Task=task,Partition=part,Graph_representations=len({r['Graph_ID'] for r in rs}),Source_studies=len({r['Study_key'] for r in rs}),Status='CANDIDATE_NOT_FROZEN'))
write('updated_highMW_split_counts.csv',split)
summary=json.loads((P/'stage3_summary.json').read_text(encoding='utf-8'));flags=json.loads((P/'unusual_structure_flags.json').read_text(encoding='utf-8'));conf=rows(P/'source_sample_identity_conflicts.csv');proto=rows(P/'protocol_subtype_review.csv')
snapshot=[(r['Scope'],r['Graph_representations'],r['Point_occurrences']) for r in counts if r['Task']=='ALL_TASKS']
report=f'''# bRo5项目原文审核：第三阶段执行记录

日期：2026-11-18。**本阶段继续执行原文审核，原方案尚未全部完成。** 当前结果是已取得证据支持的更正和可追溯的待审队列，不是最终合格分子总数。

## 本次实际完成

1. 从上一阶段66,771条来源出现记录继续处理；锁定31篇优先研究、116个原候选测试/校准结构。31篇均执行了文献元数据检索，不能据此声称31篇全文审核完成。
2. 对{len(evidence)}篇研究取得定点原文证据，逐字段记录已确认及未确认项。包括PAMPA原始数值/关系、Caco-2几何平均定义、PAMPA协议、Caco-2双向协议和结构身份错误。
3. 全体34,150个原始技术结构做了元素及自由基筛查。共22个带元素提示：21个含Se、1个含Pb。**Se不是自动排除理由**。含Pb的MZ3记录经原文Phe说明纠正，其他提示不被当作错误。
4. 全体出现记录补回作者原始化合物名称，用同论文同名样品检测表示冲突。发现{len(conf)}个待核实组，其中{summary['connectivity_conflict_groups']}个涉及连接关系差异；没有把普通“混合非对映异构体”注释当作化合物别名来合并不同编号。
5. 保留全量修订覆盖表、关系更正、结构别名更正、原文协议证据、待审论文表以及独立复算。未训练模型或执行3D/MD，也未按性能重新挑划分。

## 已落地的更正

| 问题 | 原记录 | 本阶段处理 | 证据 |
| --- | --- | --- | --- |
| MZ3结构被误录 | 侧链含[CH2][Pb]，且带自由基 | 按原文Phe侧链替换为苄基；所得规范异构图与已有MZ3图完全相同，保留旧ID到新ID映射 | E3_001，原文PDF第5页Phe讨论 |
| PAMPA检测界限丢失 | 化合物3的−logPe为7.37点值 | 原文为>7.37，转为log10Pe<−7.37；两条来源出现对应同一个原文汇总结果，不再作为点值 | E3_002，Table2 |
| Caco-2终点未判明 | 30、11、6nm/s被列为方向待定 | 明确为sqrt(PappAB×PappBA)，独立命名Caco-2:GEOMETRIC_MEAN_AB_BA；不放入AB任务 | E3_002，Table2脚注c |

几何平均的3条来源数值换算正确，本次改的是终点含义。上述更正涉及{summary['changed_occurrences']}条出现记录。额外把一个PAMPA点值规范到原文精度−6.56，但未把数值精度修订当成新增实验。

## 已确认数值及协议的四条PAMPA汇总结果

{md(['原始化合物','Papp（10^-6cm/s）','log10Papp（cm/s）'],[(r['Original_compound'],r['Papp_1e6_cm_s'],r['Log10_cm_s']) for r in rows(P/'primary_numeric_protocol_confirmed_means.csv')])}
来自10.1021/acs.jmedchem.5c02791的Table1和PAMPA方法。协议：Corning Gentest预包被板，供体pH6.4/受体pH7.4，10µM，5%DMSO，37°C、4h、无振荡。作者报告四次技术重复；表注报告两次独立实验的均值。**未获取逐次实验数据，不能将这四条均值扩成八个或更多独立事件。** 回收率仍为UNKNOWN。

2a/2b的新引入连接链手性分别与R/S一致，分子式、合成转化和记录相符；CRBN位点未指定，文中NMR也描述epimer。它们不是已经锁定唯一全立体构型的3D输入。四条记录已确认的是数值和协议，完整逐原子身份审核仍列待审。

## 当前数量：仍是候选图表示

{md(['范围','修订后点值图表示数','点值出现记录数'],snapshot)}
全量点值候选由34,142改为34,140：含Pb错误图映射回已有正确图，另一个仅有点值的图转为界限；其他原有点值仍保留。前阶段8条暂缓继续保留，本阶段新增2条界限出现记录；66,761条点值+10条非点值/暂缓=66,771条输入守恒。

**34,140仍不能称作最终独立化学实体数。** 已证实同论文同名样品存在连接、手性或数据库表示差异；这些问题不是全局去掉手性或取相同SMILES就能解决。化学实体和图表示两个层次必须继续分开。

修订后的MW>700非大环候选划分：

{md(['任务','训练候选','校准候选','测试候选'],[(t,*[next(r['Graph_representations'] for r in split if r['Task']==t and r['Partition']==p) for p in ['TRAIN','CALIBRATION','TEST']]) for t in ['Caco-2:AB','PAMPA:DIRECTION_NOT_APPLICABLE']])}
本次只投影更正后的身份与点值资格，未重新分配研究、改变种子或恢复先前purge的记录。更正后图、来源标识和化学分组跨分区检查均通过，**不意味着原文谱系、论文/专利同族及外部验证资格已完成**。

## 协议与独立系列的进展

KEAP1研究的PAMPA协议已确认：pH7.4，10µM，1%DMSO，1%lecithin/dodecane，室温24h；TableS2数值和回收率仍未核实。Cyclimids研究的Caco-2协议已确认：5µM，pH7.4，0.5%DMSO，37°C、2h，双向取样并由LC-MS/MS测定；具体数值位于SupplementaryFigure10g，尚未逐条读出，不能仅凭PROTACDB记录级引文认定测定事件来自该论文。

原文可支持部分系列标签（MZ、AT、MZP、CM/CMP、SL-X，以及特定BRD4/CRBN连接链系列），但这不是全项目独立药化系列数量。相同VHL或CRBN配体使不同命名系列仍可能有强相似性，应同时受全局化学分组约束。

另有{len(proto)}条协议提示，含PAMPA-BBB/皮肤/胃肠体系及Caco-2有机溶剂取样条件。提示来自规则筛查，不是原文确认，也没有据此把全部PAMPA直接合成一个同协议任务。

## 剩余工作与具体障碍

| 剩余项 | 当前状态与影响 |
| --- | --- |
| 31篇优先研究的逐条数值、原图结构和SI审核 | 未完成。见priority1_evidence_gaps.csv逐篇清单；元数据检索不算全文审核 |
| CM9/CM11/CMP98等同名来源结构冲突 | 已定位，待原图或作者结构文件裁决；不按来源多数投票选择 |
| 独立原始实验事件数 | 未恢复原始重复数据；已将出现、原文汇总结果和独立事件分开 |
| 独立药化系列与论文/专利同族 | 仅部分原文系列标签明确，全项目数量未定 |
| 正式训练、校准和外部测试集冻结 | 未完成；候选检查通过不能替代原文审核 |

访问障碍已实际发生：多条ACS原文/SI链接返回403或不可访问；若干PMC页面返回验证页；本地网络请求无法解析开放文献服务地址。合法作者存档已补充了两篇可读PDF。**这些障碍不证明所有缺失文献都不可获得，也不能把未审核条目标为通过。** 当前交付列明已核实证据和缺失内容，未将方案标记为全部完成。

下一项应先完成尚未核实的测试/校准来源表格、原始化合物编号和结构，解决明确的连接/手性冲突，再更新合格实体与实验事件计数。没有这些原始证据时，继续重复运行同一统计脚本不会得到最终合格数。

## 复现与交付

labels_primary_overlay.csv保留全部66,771条出现记录及原字段，point_candidates_after_primary_overlay.csv与nonpoint_or_held_occurrences.csv守恒；没有覆盖原始文件。primary_evidence.json给出URL、原文位置和未确认项；报告不附带未经授权再分发的原文全文。新包与先前stage2包配套，脚本顶部路径可按本机位置修改。stage3_independent_verification.json结果PASS。输入指纹见input_manifest.csv。
'''
(P/'bRo5_stage3_primary_review_CN.md').write_text(report,encoding='utf-8')
input_paths=[S/'labels_task_source_audit.csv',S/'molecule_groups.csv',B/'03_extraction/NONCHEMBL_ASSAY_OCCURRENCES.csv',B/'03_extraction/chembl/chembl_37/DOCUMENTS_FULL.csv',P/'primary_evidence.json']
manifest=[]
for p in input_paths:manifest.append(dict(Path=str(p),Bytes=p.stat().st_size,SHA256=hashlib.sha256(p.read_bytes()).hexdigest()))
write('input_manifest.csv',manifest)
zip_path=P/'bRo5_stage3_primary_review_20261118.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(P.iterdir()):
  if p.is_file() and p.suffix in ['.csv','.json','.md','.py'] and p.name not in ['local_primary_inventory.csv','highMW_review_occurrences.csv']:z.write(p,p.name)
with zipfile.ZipFile(zip_path) as z:assert z.testzip() is None
print(json.dumps(dict(verification=ver,zip=str(zip_path),bytes=zip_path.stat().st_size,split=split,conflict_groups=len(conf)),ensure_ascii=False,indent=2))
