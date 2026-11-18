import csv,json,collections,re,hashlib
from pathlib import Path
from rdkit import Chem,RDLogger
RDLogger.DisableLog('rdApp.*')
P=Path(__file__).resolve().parent;S=P.parent/'bro5_stage2';B=Path(r'C:\Users\guwenjing\bRo5 ADME Prediction\bRo5_Permeability_Benchmark')
def rows(p):
 with p.open(encoding='utf-8-sig',newline='') as f:yield from csv.DictReader(f)
M={r['Graph_ID']:r for r in rows(S/'molecule_groups.csv')}
L=[r for r in rows(S/'labels_task_source_audit.csv') if float(M[r['Graph_ID']]['MW'])>700 and M[r['Graph_ID']]['Topology']=='NONMACRO_SCREEN']
for r in L:
 r['Review_priority']=1 if r['Partition_SC'] in ['TEST','CALIBRATION'] and r['Task'] in ['Caco-2:AB','PAMPA:DIRECTION_NOT_APPLICABLE'] else 2 if r['Split_eligible']=='True' else 3
 Lmw=M[r['Graph_ID']]['MW'];r['MW']=Lmw
def write(n,rr):
 with (P/n).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
write('highMW_review_occurrences.csv',L)
G=collections.defaultdict(list)
for r in L:G[r['Study_key']].append(r)
studies=[]
for k,v in G.items():
 studies.append(dict(Study_key=k,Priority=min(r['Review_priority'] for r in v),Structures=len({r['Graph_ID'] for r in v}),Occurrences=len(v),Tasks='|'.join(sorted({r['Task'] for r in v})),Sources='|'.join(sorted({r['Source'] for r in v})),Original_reference=v[0]['Reference_enriched']))
studies.sort(key=lambda r:(r['Priority'],-r['Structures'],r['Study_key']))
write('study_review_inventory.csv',studies)
print(json.dumps([r for r in studies if r['Priority']==1],ensure_ascii=False,indent=2))
files=[]
for p in (B/'01_literature').rglob('*'):
 if p.is_file() and p.suffix.lower() in ['.pdf','.xml','.txt','.xlsx','.csv','.docx']:
  files.append(dict(Path=str(p),Bytes=p.stat().st_size))
write('local_primary_inventory.csv',files)
print('STUDIES',len(studies),'PRIORITY1',sum(r['Priority']==1 for r in studies),'P1_GRAPHS',len({r['Graph_ID'] for r in L if r['Review_priority']==1}))
flags=[]
for r in M.values():
 m=Chem.MolFromSmiles(r['SMILES']);bad=sorted({a.GetSymbol() for a in m.GetAtoms() if a.GetAtomicNum() not in [1,5,6,7,8,9,14,15,16,17,35,53]});rad=sum(a.GetNumRadicalElectrons() for a in m.GetAtoms())
 if bad or rad:flags.append(dict(Graph_ID=r['Graph_ID'],Elements=','.join(bad),Radical_electrons=rad,MW=r['MW'],SMILES=r['SMILES']))
(P/'unusual_structure_flags.json').write_text(json.dumps(flags,ensure_ascii=False,indent=2),encoding='utf-8')
write('unusual_structure_flags.csv',flags)
