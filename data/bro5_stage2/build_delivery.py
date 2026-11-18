import csv,json,zipfile,collections,hashlib,sys
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
P=Path(__file__).resolve().parent
def read(n):
 with (P/n).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
S=json.loads((P/'stage2_summary.json').read_text(encoding='utf-8'));V=json.loads((P/'independent_stage2_verification.json').read_text(encoding='utf-8'))
F=read('endpoint_feasibility.csv');D=read('split_feasibility.csv');L=read('labels_task_source_audit.csv');M={r['Graph_ID']:r for r in read('molecule_groups.csv')}
def md(headers,rows):return '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+''.join('| '+' | '.join(str(v) for v in r)+' |\n' for r in rows)
targets=['Caco-2:AB','PAMPA:DIRECTION_NOT_APPLICABLE','MDCK_MDR1_MENTION:AB','MDCK_WT:AB','MDCK_LE:AB']
def table(scope):
 rows=[]
 for r in F:
  if r['Scope']!=scope or r['Task'] not in targets:continue
  task=r['Task'];split={x['Partition']:x for x in D if x['Scope']==scope and x['Task']==task and x['Split_mode']=='SC'}
  rows.append([task,r['Stage2_point_structures'],r['Task_defined_traceable_structures'],r['Scaffold_groups'],r['FP06_components'],split['TRAIN']['Structures'],split['CALIBRATION']['Structures'],split['TEST']['Structures']])
 return md(['终点任务','修正后点值结构','有DOI/专利且任务可定义','骨架/连接分组¹','相似性连通分组¹','训练候选²','校准候选²','测试候选²'],rows)
pointids={r['Graph_ID'] for r in L if r['Stage2_point_eligible']=='True'}
scopes=[]
for scope in ['ALL_MW','MW_GT500','MW_GT700_NONMACRO','MW_700_1000_NONMACRO']:
 def ok(g):
  m=M[g];mw=float(m['MW']);non=m['Topology']=='NONMACRO_SCREEN'
  return scope=='ALL_MW' or scope=='MW_GT500' and mw>500 or scope=='MW_GT700_NONMACRO' and mw>700 and non or scope=='MW_700_1000_NONMACRO' and 700<mw<=1000 and non
 pp={g for g in pointids if ok(g)};eligible={r['Graph_ID'] for r in L if r['Split_eligible']=='True' and ok(r['Graph_ID'])}
 scopes.append([scope,len(pp),len(eligible)])
copies=read('confirmed_database_record_copies.csv');held=read('stage2_held_point_labels.csv')
counts=collections.Counter(r['Stage2_hold_reason'] for r in held)
report=f'''# bRo5渗透性项目：终点、来源与验证供给审核

日期：2026-11-18。上一轮34,150个技术合格结构是本轮输入。本轮完成已有点值集的来源补全、方向/细胞体系细分、已知原文冲突排查、全量化学分组、候选划分和守恒复算。未训练模型，未执行3D或MD。

## 1. 当前可以确定的数量

修正后仍有**{S['stage2_point_structures']:,}个结构**至少具备一条未触发本轮暂缓规则的点值，共**{S['stage2_point_fields']:,}条来源出现记录**。其中**{S['task_defined_traceable_structures']:,}个结构**至少具有一条DOI/专利标识可解析、方向/终点可定义且未触发已知冲突的标签，可进入二维基线的有限元数据候选层。

这两个数量都不代表全部逐篇原文审核通过。后一个数量也不能跨终点合并成一个回归任务。

{md(['范围','修正后点值结构并集','至少一条来源/任务明确的结构'],scopes)}
MW>500只是尺寸分层；大环依据既有最大环≥12筛查，不代表人工modality裁决。各行有包含关系。

上一轮共66,771条点值出现记录，本轮暂缓{S['held_point_fields']}条；有{S['held_only_structures']}个结构因没有其他保留点值而退出点值并集。未因某条测定失败而删除同一分子的其他测定。

{md(['暂缓原因','出现记录数'],sorted(counts.items()))}
已直接重读2020年原始SI工作簿（Summary行15/19/23，SHA256 db4dbb36823fc5c6cbb638a3c4854606257522aae45d7e5bdf5b00eac02b024b）：RC-1、RNC-1、IRC-1回收率24.01%、29.48%、18.06%，原始结果为>8.16、>8.68、>7.61的负对数界限。它们不能作为普通点值，也不能仅因有界限就认定为可靠删失标签。原始original_pampa含不等号而整理值丢失符号者转入暂缓；同分子/同论文/PAMPA的跨库潜在副本一并待核实，其他实验条件的结果可经原文复核恢复。

## 2. 高分子量非大环的实际任务供给

### MW>700 Da非大环

{table('MW_GT700_NONMACRO')}
### 700<MW≤1000 Da非大环

{table('MW_700_1000_NONMACRO')}
### 全部分子量：用于低分子量基线与迁移研究

{table('ALL_MW')}
¹ 分组数针对该任务技术候选结构，是计算代理，不是确认的独立药化系列。分组在全体结构上建立后再对子集计数，避免子集之外的桥接分子被忽略。

² 训练/校准/测试是同一套全局候选划分，尚非正式外部验证集。只纳入单一DOI/专利可解析且任务可定义记录；pH、抑制剂及回收率未报告仍明确保留UNKNOWN，未据此认定有效或无效。Caco-2未知方向、MDCK未明细胞体系、多来源归属未解者留在补充层。PAMPA一般型、intrinsic和apparent分开，未据此断言膜组成相同。MDCK_MDR1_MENTION只代表描述中提及转运体，抑制剂条件仍须核验。

## 3. 独立实验与来源去重

- 统一DOI大小写、URL前缀与专利号；解析MolMeDB快照文献元数据。ref1030、ref1031在该本地快照没有可用书目字段，含这两者的多参考记录没有被任意指定成单一原始研究。
- 数据库记录ID、结构、终点、方向和数值一致且跨快照重复的组为{len(copies):,}组，保留完整出现关系。此证据支持“同一数据库记录重复出现”，不支持声称原始实验事件已全部确认。
- 共有{S['candidate_event_signatures']:,}个协议/结构/来源/任务/数值签名候选；**不能把这个数写成独立实验数**。不同协议不能因数值相等直接合并，相同协议中的重复测定也可能无法从现有字段辨识。
- 同结构/研究/任务的值域>0.5 log10单位有{S['same_study_task_value_review_groups']:,}组，逐组列为不同条件或冲突候选，未自动取均值或剔除。

**独立实验事件数和原文确认的药化系列数，目前仍不可确认为单一最终值。** 未知不是0，也不应以骨架数或数据库行数替代。本次已产出带来源锚点的逐条复核队列。

## 4. 化学分组与划分规则

RDKit {S['rdkit']}；先把相同无立体连接表示及相同非空Murcko骨架连接，得到{S['scaffold_connectivity_groups']:,}组。空骨架仅按连接表示归组，不把所有无环分子并成一组。此处不跨电荷/互变异构体自动统一。

补充敏感性分组采用Morgan半径2、2048位、不含手性、Tanimoto≥0.60的单链接连通分量，并合入前述骨架/连接约束。共精确计算{S['fingerprint_pairs']:,}对，阈值边{S['fingerprint_edges_ge06']:,}条，得到{S['fp06_components']:,}组；无近邻截断。单链接可能经桥接造成巨大分量，不能将一个分量当作一个真实系列。

候选研究标识按SHA256('20261118|Study_key')前16个十六进制数转整数、模100划分：<15测试、15–29校准、其余训练。未使用标签值、预测性能或反复挑种子。任何化学组跨研究分区时保留优先级测试>校准>训练，低优先级出现记录移入PURGED_CHEMICAL_OVERLAP。相同分子跨终点采用全局约束。

SC为骨架/连接分组隔离，FP06为更严格敏感性方案。两者均通过研究标识、精确结构、无立体连接表示和所属化学组的跨分区零重叠检查。**研究标识零重叠不等于论文/专利家族及所有原始实验谱系已完全核清；这仍是候选验证设计。** SC不能声称全部近邻相似度<0.60；FP06才提供该阈值下的全图分离。

全部分区和流失见split_feasibility.csv，不能把被purge的标签视作实验无效。测试后才发现的同源关系应按版本规则重新审核，不能为了保留规模而隐藏。

## 5. 推进结论

项目可以推进到**限定终点、明确补充层的二维基线准备与原文审核**。当前不宜把34,150或修正后总数写成最终训练样本，更不能据此宣称已经具备同等规模的高分子量非大环测试集。

优先路线是PAMPA与Caco-2 A→B分别开展；高分子量非大环作为明确的目标域，低分子量数据作为迁移/多任务的候选支持域。是否采用多任务，应在主终点和验证队列锁定后决定。MDCK按细胞体系与方向分开，不能拼成单一MDCK任务。

正式建模前必须核实拟保留测试与校准记录的原始研究、方向、细胞/膜体系、单位、检测界限、回收率与原始化合物身份；确认原文药化系列，并检查论文/专利同族。未完成这些项时，当前表只能支持技术可行性判断。3D/构象增益研究还须单独解决立体身份与构象预算。

## 6. 文件与复现

- labels_task_source_audit.csv：66,771条来源出现及任务、来源、暂缓、候选分区。
- stage2_point_labels.csv / stage2_held_point_labels.csv：本轮保留/暂缓点值，守恒对应上一轮点值集；不覆盖既有审核文件。
- molecule_groups.csv：全局结构与两套计算分组。
- event_candidates.csv / confirmed_database_record_copies.csv / same_study_value_review.csv：事件候选、已识别数据库副本与冲突队列。
- endpoint_feasibility.csv / split_feasibility.csv / study_concentration.csv：分母、验证预算与来源集中度。
- highMW_nonmacro_primary_review_queue.csv：MW>700非大环的逐条原文复核任务。
- reference_resolution.csv / rdf_reference_metadata.json：本地文献编号解析证据。
- leakage_checks.csv / independent_stage2_verification.json：跨分区检查与独立代码复算；结果{V['status']}。
- run_stage2.py / verify_stage2.py / input_manifest.csv：可重跑脚本及输入SHA256。脚本输入沿用本机原目录，迁移时改顶部路径。全量指纹已计算；复跑允许在核对图结构后复用缓存分组，指纹边没有再独立计算一遍。

完整数据包仅供本项目内部研究审核；未上传GitHub，不能据此替各原始来源作开放再分发授权。
'''
(P/'bRo5_stage2_report_CN.md').write_text(report,encoding='utf-8')
with (P/'scope_summary.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.writer(f);w.writerow(['Scope','Stage2_point_structures','Source_task_defined_structures']);w.writerows(scopes)
zip_path=P/'bRo5_stage2_audit_20261118.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(P.iterdir()):
  if p.is_file() and p.suffix in ['.csv','.json','.md','.py'] and p.name!='inspect_inputs.py':z.write(p,p.name)
with zipfile.ZipFile(zip_path) as z:assert z.testzip() is None
print(json.dumps({'summary':S,'scopes':scopes,'zip_bytes':zip_path.stat().st_size,'zip':str(zip_path),'highMW_task_table':table('MW_GT700_NONMACRO')},ensure_ascii=False,indent=2))
