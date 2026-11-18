import csv,json,collections,hashlib
from pathlib import Path
P=Path(__file__).resolve().parent
def read(n):
 with (P/n).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
L=read('labels_task_source_audit.csv');M={r['Graph_ID']:r for r in read('molecule_groups.csv')}
points=read('stage2_point_labels.csv');held=read('stage2_held_point_labels.csv')
assert len(points)+len(held)==len(L)==66771
assert len({r['Occurrence_ID'] for r in L})==len(L)
assert {r['Occurrence_ID'] for r in points}.isdisjoint(r['Occurrence_ID'] for r in held)
assert all(r['Stage2_point_eligible']=='True' and not r['Stage2_hold_reason'] for r in points)
assert all(r['Stage2_point_eligible']=='False' and r['Stage2_hold_reason'] for r in held)
for mode,field in [('SC','Scaffold_connectivity_group'),('FP06','Scaffold_FP06_component')]:
 graphparts=collections.defaultdict(set);studyparts=collections.defaultdict(set);groupparts=collections.defaultdict(set)
 for r in L:
  p=r['Partition_'+mode]
  if p not in ['TRAIN','CALIBRATION','TEST']:continue
  assert r['Split_eligible']=='True' and r['Stage2_point_eligible']=='True'
  assert r['Study_key'].startswith(('DOI:','PATENT:'))
  graphparts[r['Graph_ID']].add(p);studyparts[r['Study_key']].add(p);groupparts[M[r['Graph_ID']][field]].add(p)
 assert all(len(v)==1 for v in graphparts.values())
 assert all(len(v)==1 for v in studyparts.values())
 assert all(len(v)==1 for v in groupparts.values())
for r in read('endpoint_feasibility.csv'):
 def scope(g):
  m=M[g];mw=float(m['MW']);non=m['Topology']=='NONMACRO_SCREEN';s=r['Scope']
  return s=='ALL_MW' or s=='MW_GT500' and mw>500 or s=='MW_GT700_NONMACRO' and mw>700 and non or s=='MW_700_1000_NONMACRO' and 700<mw<=1000 and non
 rr=[x for x in L if x['Task']==r['Task'] and scope(x['Graph_ID'])]
 assert len({x['Graph_ID'] for x in rr})==int(r['Technical_structures'])
 assert len({x['Graph_ID'] for x in rr if x['Stage2_point_eligible']=='True'})==int(r['Stage2_point_structures'])
 assert len({x['Graph_ID'] for x in rr if x['Split_eligible']=='True'})==int(r['Task_defined_traceable_structures'])
out=dict(status='PASS',source_occurrences=len(L),stage2_points=len(points),held_fields=len(held),stage2_structures=len({r['Graph_ID'] for r in points}),verified='Conservation; occurrence uniqueness; held-point exclusion; DOI/patent, graph and chemical-group partition separation; endpoint/scope independent recount',limits='No primary-study full verification; no confirmed medicinal chemistry series; FP0.6 edges not independently recomputed')
(P/'independent_stage2_verification.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(out,ensure_ascii=False,indent=2))
