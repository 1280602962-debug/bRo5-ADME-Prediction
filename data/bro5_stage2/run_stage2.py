"""Source/event candidate audit and pre-model leakage-controlled feasibility census.
No label averaging, ML fitting, or assertions of primary-source validation.
"""
import csv,json,re,hashlib,gzip,collections,math,statistics,time,zipfile,io
from pathlib import Path
from rdkit import Chem,DataStructs,RDLogger,rdBase
from rdkit.Chem import rdFingerprintGenerator
RDLogger.DisableLog('rdApp.*')
B=Path(r'C:\Users\guwenjing\bRo5 ADME Prediction\bRo5_Permeability_Benchmark')
A=Path(r'C:\Users\guwenjing\.codex\.chatgpt-projects\g-p-6a9522ad97f0819199b8f5924494f4c6\bro5_audit_20261117')
OUT=Path(__file__).resolve().parent
INPUTS=set()
def rows(p):
 INPUTS.add(p)
 with p.open(encoding='utf-8-sig',newline='') as f:yield from csv.DictReader(f)
def write(name,rs,fields=None):
 rs=list(rs)
 if not rs and not fields:return
 with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rs[0]));w.writeheader();w.writerows(rs)
def hid(tag,x):return tag+'_'+hashlib.sha256(str(x).encode()).hexdigest()[:24]
print('LOAD',flush=True)
M={r['Graph_ID']:r for r in rows(A/'technical_qualified_molecules_enriched.csv')}
L=list(rows(A/'technical_point_labels.csv'))
assert len(M)==34150 and len(L)==66771
legacy={r['AssayOccurrence_ID']:r for r in rows(B/'03_extraction/NONCHEMBL_ASSAY_OCCURRENCES.csv')}
rawlegacy={k:json.loads(v['Source_raw_JSON']) for k,v in legacy.items()}
cited={r['AssayOccurrence_ID']:r['Reference'] for r in rows(B/'06_deduplication/batch02/CITED_REFERENCE_ASSAY_INDEX.csv')}
chemneed={r['Source_row_ID'] for r in L if r['Source']=='ChEMBL37'}
chem={r['a__activity_id']:r for r in rows(B/'03_extraction/chembl/chembl_37/CHEMBL_CURRENT_RAW_PERMEABILITY.csv') if r['a__activity_id'] in chemneed}
assert len(chem)==len(chemneed)
params=collections.defaultdict(list)
for r in rows(B/'03_extraction/chembl/chembl_37/ASSAY_PARAMETERS.csv'):params[r['assay_id']].append(r)
bad_primary=list(rows(B/'05_assay_audit/batch02/HUDEC_PRIMARY_PAMPA_BOUND_ADJUDICATION.csv'))
INPUTS.add(B/'01_literature/original_sources_batch02/PMC7450057_supplementary/41467_2020_17997_MOESM6_ESM.xlsx')
bad_so={r['SourceOccurrence_ID'] for r in bad_primary}
cyc_current={str(i):r for i,r in enumerate(rows(B.parent/'CycPeptMPDB_Peptide_All.csv'),2)}
npmzip=B.parent/'SI_NPMMP_DB.zip';INPUTS.add(npmzip)
with zipfile.ZipFile(npmzip) as z:npm_current={str(i):r for i,r in enumerate(csv.DictReader(io.StringIO(z.read('SI_NPMMP_DB/Data/csv/Overall.csv').decode('utf-8-sig'))),2)}
bound_graph_studies=set();failure_graph_studies=set()
for r in L:
 raw=rawlegacy.get(r['Source_row_ID'],{})
 if r['Source']=='Hudec_new_parsed' and re.search(r'[<>≤≥]',str(raw.get('original_pampa',''))):
  bound_graph_studies.add((r['Graph_ID'],r['Reference'].lower()))
  if legacy[r['Source_row_ID']]['SourceOccurrence_ID'] in bad_so:failure_graph_studies.add((r['Graph_ID'],r['Reference'].lower()))
print('RESOLVE_RDF_REFERENCES',flush=True)
rdfpath=B/'02_raw_data/browser_supplements/molmedb/molmedb.nt.gz'
cache=OUT/'rdf_reference_metadata.json'
if cache.exists():rdf=json.loads(cache.read_text(encoding='utf-8'))
else:
 rdf=collections.defaultdict(dict)
 INPUTS.add(rdfpath)
 with gzip.open(rdfpath,'rt',encoding='utf-8') as f:
  for line in f:
   if not line.startswith('<https://rdf.molmedb.upol.cz/reference/'):continue
   if 'bibliographicCitation>' not in line and not re.search(r'/(?:identifier|title|sameAs|doi)>|#sameAs>',line):continue
   mt=re.match(r'<([^>]+)> <([^>]+)> (.+?)\s*\.\s*$',line)
   if mt:rdf[mt[1]][mt[2]]=mt[3]
 cache.write_text(json.dumps(rdf,ensure_ascii=False,indent=2),encoding='utf-8')
INPUTS.add(rdfpath)
def refkey(text):
 text=text.strip();d=re.findall(r'10\.\d{4,9}/[^\s<>"|]+',text,re.I)
 if len(set(x.lower().rstrip('.,;') for x in d))==1:return 'DOI:'+d[0].lower().rstrip('.,;')
 pat=re.fullmatch(r'(?:https?://[^ ]*/)?((?:WO|EP|US|CN)\d{6,}[A-Z]\d?)',text,re.I)
 if pat:return 'PATENT:'+pat[1].upper()
 if not text:return ''
 if 'rdf.molmedb' in text:return 'UNRESOLVED_URI:'+text
 if text.startswith('CHEMBL_DOC:'):return 'DATABASE_DOC:'+text
 if re.fullmatch(r'\d{4}_[A-Za-z]+(?:_\d+)?',text):return 'CYC_SOURCE:'+text.lower()
 return 'CITATION:'+re.sub(r'\W+','',text.lower())
def direction(text,existing):
 t=text.lower().replace('→',' to ')
 ab=bool(re.search(r'\ba\s*(?:-|to|2|=>)\s*b\b|apical\s*(?:-|to)\s*(?:basolateral|basal)',t))
 ba=bool(re.search(r'\bb\s*(?:-|to|2|=>)\s*a\b|(?:basolateral|basal)\s*(?:-|to)\s*apical',t))
 if ab and ba:return existing if existing in ['AB','BA'] else 'BIDIRECTIONAL_UNRESOLVED'
 if ab:return 'AB'
 if ba:return 'BA'
 return 'AB' if existing.startswith('AB') else 'BA' if existing.startswith('BA') else 'UNKNOWN'
def cellmodel(ep,txt):
 t=txt.lower()
 if ep!='MDCK':return ep
 if re.search(r'mdr\s*[- ]?1|mdck\s*[- ]?ii\s*[- ]?mdr|p[- ]?gp',t):return 'MDCK_MDR1_MENTION'
 if re.search(r'low[- ]?efflux|mdck[- ]?le',t):return 'MDCK_LE'
 if re.search(r'bcrp',t):return 'MDCK_BCRP_MENTION'
 if re.search(r'wild[- ]?type|mdck[- ]?wt',t):return 'MDCK_WT'
 return 'MDCK_UNSPECIFIED'
enh=[];refmaps=[]
for i,r0 in enumerate(L):
 r=dict(r0);src=r['Source'];rid=r['Source_row_ID'];raw=rawlegacy.get(rid,{})
 text=r['Protocol'];ref=r['Reference'];compound='';locator='';snapshot='';desc=''
 if src=='ChEMBL37':
  q=chem[rid];desc=q['s__description'];text=desc+' '+q['a__standard_type'];compound=q['cr__compound_key'];locator='ChEMBL37 activity '+rid+' / assay '+q['a__assay_id'];snapshot='CHEMBL37:'+rid
 elif rid in legacy:
  q=legacy[rid];text+=' '+str(raw.get('Description',raw.get('note','')));compound=raw.get('Original_ID',raw.get('Compound',raw.get('Original_Name_in_Source_Literature','')));locator=q['Source_file']+':row:'+q['Source_row'];snapshot=src+':'+q['Source_record_ID']+':'+q['Endpoint_original']
  ref=ref or cited.get(rid,'') or raw.get('Source','')
  if src=='Kim2025':ref='10.1021/acs.jcim.5c01600'
 elif src in ['CycPept_current','NPMMP_archive']:locator=src+':CSV_row:'+rid
 if src in ['CycPept_current','CycPeptMPDB_historical']:
  cr=cyc_current.get(rid.split(':')[0],{}) if src=='CycPept_current' else raw
  record=cr.get('ID','');snapshot=hid('CYC_RECORD',(record,r['Graph_ID'],r['Endpoint'],round(float(r['Log10_cm_s']),6))) if record else ''
 if src in ['NPMMP_archive','NPMMPD_official_current']:
  nr=npm_current.get(rid,{}) if src=='NPMMP_archive' else raw
  record=nr.get('ID','');snapshot=hid('NPM_RECORD',(record,r['Graph_ID'],r['Endpoint'],r['Direction'],round(float(r['Log10_cm_s']),6))) if record else ''
 if src=='MolMeDB_v2':
  uris=[x.strip() for x in ref.split('|')];pieces=[]
  for u in uris:
   attrs=rdf.get(u,{})
   pieces.append(' '.join(attrs.values()) or u)
   refmaps.append(dict(URI=u,Resolved_metadata=' '.join(attrs.values()),Status='RESOLVED_CITATION' if attrs else 'UNRESOLVED'))
  ref=' | '.join(pieces)
 d=direction(text,r['Direction']);ep=r['Endpoint'];model=cellmodel(ep,text)
 task=ep+':'+d if ep in ['Caco-2','MDCK','RRCK','LLC-PK1'] else ep+':DIRECTION_NOT_APPLICABLE'
 if ep=='MDCK':task=model+':'+d
 rk=refkey(ref)
 ph='|'.join(sorted(set(re.findall(r'\bpH\s*[=:]?\s*([0-9]+(?:\.[0-9]+)?)',text,re.I))))
 inhibitor='MENTION_REQUIRES_REVIEW' if re.search(r'elacridar|tariquidar|zosuquidar|cyclosporin|verapamil|gf120918|inhibitor',text,re.I) else 'NOT_REPORTED_NOT_ABSENT'
 transporter='|'.join(sorted(set(re.findall(r'MDR[- ]?1|BCRP|P[- ]?gp',text,re.I))))
 sourceok=bool(rk) and not rk.startswith(('UNRESOLVED_URI:','DATABASE_DOC:')) and ' | ' not in ref
 semantic='TASK_DEFINED_LIMITED_METADATA' if (ep.startswith('PAMPA') or d in ['AB','BA']) else 'DIRECTION_REVIEW'
 if ep=='MDCK' and model=='MDCK_UNSPECIFIED':semantic='CELL_MODEL_REVIEW'
 review='';study_doi=rk[4:] if rk.startswith('DOI:') else ''
 if ep=='PAMPA' and (r['Graph_ID'],study_doi) in bound_graph_studies:review='ORIGINAL_CENSORED_POINT_OR_POTENTIAL_CROSS_SOURCE_COPY_REVIEW'
 if ep=='PAMPA' and (r['Graph_ID'],study_doi) in failure_graph_studies:review='PRIMARY_LOW_RECOVERY_AND_CENSORED_RESULT_REVIEW'
 if src=='Hudec_new_parsed' and re.search(r'[<>≤≥]',str(raw.get('original_pampa',''))):review=review or 'ORIGINAL_CENSORED_NOT_EXACT_POINT'
 primaryid=sourceok and rk.startswith(('DOI:','PATENT:'))
 r.update(Occurrence_ID=hid('O',src+'|'+rid),Reference_original=r['Reference'],Reference_enriched=ref,Study_key=rk,Source_link_status='TRACEABLE_REFERENCE' if sourceok else 'REFERENCE_REVIEW',Primary_identifier_status='DOI_OR_PATENT' if primaryid else 'CITATION_OR_DATABASE_ALIAS_REVIEW',Direction_refined=d,Task=task,Cell_model=model,pH_reported=ph or 'UNKNOWN',Inhibitor_status=inhibitor,Transporter_mention=transporter or 'UNKNOWN',Compound_original_ID=compound,Evidence_locator=locator or src+':'+rid,Assay_description=desc or text,Snapshot_record_key=snapshot,Semantic_status=semantic,Primary_validity='NOT_PRIMARY_VALIDATED',Stage2_point_eligible=not review,Stage2_hold_reason=review,Split_eligible=primaryid and semantic=='TASK_DEFINED_LIMITED_METADATA' and not review)
 enh.append(r)
print('TASKS_AND_CANDIDATE_EVENTS',flush=True)
write('reference_resolution.csv',{r['URI']:r for r in refmaps}.values())
lineages=collections.defaultdict(list)
for r in enh:
 if r['Snapshot_record_key'].startswith(('CYC_RECORD_','NPM_RECORD_')):lineages[r['Snapshot_record_key']].append(r)
write('confirmed_database_record_copies.csv',[dict(Database_record_lineage=k,Graph_ID=v[0]['Graph_ID'],Endpoint=v[0]['Endpoint'],Occurrences=len(v),Sources='|'.join(sorted({r['Source'] for r in v})),Occurrence_IDs='|'.join(r['Occurrence_ID'] for r in v),Evidence='SAME_DATABASE_RECORD_ID_GRAPH_ENDPOINT_DIRECTION_AND_VALUE_ACROSS_SNAPSHOTS_NOT_INDEPENDENT_EVENT_CONFIRMATION') for k,v in sorted(lineages.items()) if len({r['Source'] for r in v})>1])
# Never equate identical values with independently confirmed events.
# Exact protocol-specific signatures are only event candidates.
groups=collections.defaultdict(list);broad=collections.defaultdict(list)
for r in enh:
 key=(r['Graph_ID'],r['Study_key'] or 'UNKNOWN:'+r['Source'],r['Task'],r['Protocol'],round(float(r['Log10_cm_s']),6))
 r['Candidate_event_ID']=hid('CE',key);groups[r['Candidate_event_ID']].append(r)
 broad[(r['Graph_ID'],r['Study_key'] or 'UNKNOWN:'+r['Source'],r['Task'])].append(r)
write('event_candidates.csv',[dict(Candidate_event_ID=k,Graph_ID=v[0]['Graph_ID'],Study_key=v[0]['Study_key'],Task=v[0]['Task'],Protocol=v[0]['Protocol'],Log10_cm_s=v[0]['Log10_cm_s'],Occurrences=len(v),Sources='|'.join(sorted({r['Source'] for r in v})),Occurrence_IDs='|'.join(r['Occurrence_ID'] for r in v),Resolution='SAME_SIGNATURE_NOT_CONFIRMED_EVENT') for k,v in sorted(groups.items())])
conflicts=[]
for k,v in broad.items():
 vals=[float(r['Log10_cm_s']) for r in v];span=max(vals)-min(vals)
 if span>0.5:conflicts.append(dict(Graph_ID=k[0],Study_key=k[1],Task=k[2],Fields=len(v),Min_log10=min(vals),Max_log10=max(vals),Range_log10=span,Status='DIFFERENT_CONDITIONS_OR_CONFLICT_REVIEW_NOT_AUTOMATIC_REJECTION'))
write('same_study_value_review.csv',conflicts)
class DSU:
 def __init__(self,keys):self.p={k:k for k in keys}
 def find(self,k):
  while self.p[k]!=k:self.p[k]=self.p[self.p[k]];k=self.p[k]
  return k
 def union(self,a,b):
  a=self.find(a);b=self.find(b)
  if a!=b:self.p[max(a,b)]=min(a,b)
def group_values(dsu,field):
 seen={}
 for k,m in sorted(M.items()):
  v=m[field]
  if v:
   if v in seen:dsu.union(k,seen[v])
   else:seen[v]=k
sc=DSU(M);group_values(sc,'Connectivity_SMILES');group_values(sc,'MurckoScaffold')
for k,m in M.items():m['Scaffold_connectivity_group']=hid('SC',sc.find(k))
print('FULL_EXHAUSTIVE_FINGERPRINT_PAIRS',len(M),flush=True)
# Full exact all-pairs bulk comparison, no nearest-neighbor truncation.
keys=sorted(M)
if (OUT/'molecule_groups.csv').exists() and (OUT/'stage2_summary.json').exists():
 cached={r['Graph_ID']:r for r in rows(OUT/'molecule_groups.csv')}
 assert set(cached)==set(M) and all(cached[k]['SMILES']==M[k]['SMILES'] for k in M)
 old=json.loads((OUT/'stage2_summary.json').read_text(encoding='utf-8'))
 near_edges=old['fingerprint_edges_ge06'];pairs=old['fingerprint_pairs']
 for k,m in M.items():m['Scaffold_FP06_component']=cached[k]['Scaffold_FP06_component']
 print('VERIFIED_FP_CACHE_REUSED',flush=True)
else:
 gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=False)
 fps=[gen.GetFingerprint(Chem.MolFromSmiles(M[k]['SMILES'])) for k in keys]
 fpdsu=DSU(M)
 for k in keys:fpdsu.union(k,sc.find(k))
 near_edges=0;pairs=0
 for i,k in enumerate(keys):
  sims=DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i]);pairs+=i
  for j,s in enumerate(sims):
   if s>=0.6:fpdsu.union(k,keys[j]);near_edges+=1
  if i%5000==0:print('FP_PROGRESS',i,flush=True)
 for k,m in M.items():m['Scaffold_FP06_component']=hid('FP',fpdsu.find(k))
assert pairs==len(M)*(len(M)-1)//2
write('molecule_groups.csv',[dict(Graph_ID=k,SMILES=m['SMILES'],Connectivity_SMILES=m['Connectivity_SMILES'],MurckoScaffold=m['MurckoScaffold'],Scaffold_connectivity_group=m['Scaffold_connectivity_group'],Scaffold_FP06_component=m['Scaffold_FP06_component'],MW=m['MW'],Topology=m['Topology'],Unassigned_chiral_centers=m['Unassigned_chiral_centers'],Series_status='COMPUTATIONAL_PROXY_NOT_CONFIRMED_MEDCHEM_SERIES') for k,m in sorted(M.items())])
# Source split is selected without labels. Source hashes frozen in advance.
# Prioritize TEST over CALIBRATION over TRAIN; purge any lower-priority chemical group.
def fold(study):
 q=int(hashlib.sha256(('20261118|'+study).encode()).hexdigest()[:16],16)%100
 return 'TEST' if q<15 else 'CALIBRATION' if q<30 else 'TRAIN'
priority={'TRAIN':0,'CALIBRATION':1,'TEST':2}
studyfold={r['Study_key']:fold(r['Study_key']) for r in enh if r['Split_eligible']}
write('source_split_assignments.csv',[dict(Study_key=k,Assigned_partition=v,Rule='SHA256(20261118|Study_key) mod 100; test<15; calibration<30; otherwise train') for k,v in sorted(studyfold.items())])
for mode,field in [('SC','Scaffold_connectivity_group'),('FP06','Scaffold_FP06_component')]:
 gp={}
 for r in enh:
  if not r['Split_eligible']:continue
  g=M[r['Graph_ID']][field];p=studyfold[r['Study_key']]
  if g not in gp or priority[p]>priority[gp[g]]:gp[g]=p
 for r in enh:
  if not r['Split_eligible']:part='SUPPLEMENT_ONLY'
  else:part=studyfold[r['Study_key']] if gp[M[r['Graph_ID']][field]]==studyfold[r['Study_key']] else 'PURGED_CHEMICAL_OVERLAP'
  r['Partition_'+mode]=part
write('labels_task_source_audit.csv',enh)
write('stage2_point_labels.csv',[r for r in enh if r['Stage2_point_eligible']])
write('stage2_held_point_labels.csv',[r for r in enh if not r['Stage2_point_eligible']])
def scopes(g):
 m=M[g];mw=float(m['MW']);non=m['Topology']=='NONMACRO_SCREEN'
 return ['ALL_MW']+(['MW_GT500'] if mw>500 else [])+(['MW_GT700_NONMACRO'] if mw>700 and non else [])+(['MW_700_1000_NONMACRO'] if 700<mw<=1000 and non else [])
idx=collections.defaultdict(list)
for r in enh:
 for s in scopes(r['Graph_ID']):idx[(s,r['Task'])].append(r)
stats=[];splitstats=[];concentration=[]
for (scope,task),rs in sorted(idx.items()):
 gids={r['Graph_ID'] for r in rs};known={r['Study_key'] for r in rs if r['Source_link_status']=='TRACEABLE_REFERENCE'}
 eligible=[r for r in rs if r['Split_eligible']];eg={r['Graph_ID'] for r in eligible}
 bystudy=collections.defaultdict(set)
 for r in rs:
  if r['Source_link_status']=='TRACEABLE_REFERENCE':bystudy[r['Study_key']].add(r['Graph_ID'])
 largest=max((len(v) for v in bystudy.values()),default=0)
 stats.append(dict(Scope=scope,Task=task,Technical_structures=len(gids),Stage2_point_structures=len({r['Graph_ID'] for r in rs if r['Stage2_point_eligible']}),Traceable_studies=len(known),Task_defined_traceable_structures=len(eg),Scaffold_groups=len({M[g]['Scaffold_connectivity_group'] for g in gids}),FP06_components=len({M[g]['Scaffold_FP06_component'] for g in gids}),Largest_study_structures=largest,Largest_study_fraction=largest/len(gids),Primary_validated_structures='NOT_ESTABLISHED',Confirmed_independent_series='NOT_ESTABLISHED'))
 for study,gs in bystudy.items():concentration.append(dict(Scope=scope,Task=task,Study_key=study,Structures=len(gs)))
 for mode,field in [('SC','Scaffold_connectivity_group'),('FP06','Scaffold_FP06_component')]:
  for part in ['TRAIN','CALIBRATION','TEST','PURGED_CHEMICAL_OVERLAP','SUPPLEMENT_ONLY']:
   rr=[r for r in rs if r['Partition_'+mode]==part];gg={r['Graph_ID'] for r in rr}
   splitstats.append(dict(Scope=scope,Task=task,Split_mode=mode,Partition=part,Label_fields=len(rr),Structures=len(gg),Studies=len({r['Study_key'] for r in rr if r['Study_key']}),Chemical_groups=len({M[g][field] for g in gg})))
write('endpoint_feasibility.csv',stats);write('split_feasibility.csv',splitstats);write('study_concentration.csv',sorted(concentration,key=lambda r:(r['Scope'],r['Task'],-r['Structures'])))
queue=[]
for r in enh:
 m=M[r['Graph_ID']]
 if float(m['MW'])<=700 or m['Topology']!='NONMACRO_SCREEN':continue
 reasons=['PRIMARY_PROTOCOL_AND_VALIDITY_CONFIRMATION']
 if not r['Split_eligible']:reasons.append(r['Semantic_status'] if r['Source_link_status']=='TRACEABLE_REFERENCE' else 'REFERENCE_REVIEW')
 if int(m['Unassigned_chiral_centers']):reasons.append('STEREO_IDENTITY_REVIEW')
 queue.append(dict(Graph_ID=r['Graph_ID'],Occurrence_ID=r['Occurrence_ID'],Task=r['Task'],MW=m['MW'],Study_key=r['Study_key'],Source=r['Source'],Evidence_locator=r['Evidence_locator'],Original_compound_ID=r['Compound_original_ID'],Protocol=r['Assay_description'],Review_actions='|'.join(reasons),Priority=1 if r['Task'] in ['Caco-2:AB','PAMPA:DIRECTION_NOT_APPLICABLE'] else 2))
write('highMW_nonmacro_primary_review_queue.csv',sorted(queue,key=lambda r:(r['Priority'],r['Study_key'],r['Graph_ID'])))
# Verify global source/structure/scaffold/FP component separation and conservation.
checks=[]
for mode,field in [('SC','Scaffold_connectivity_group'),('FP06','Scaffold_FP06_component')]:
 sets={p:dict(graph=set(),connection=set(),study=set(),group=set()) for p in priority}
 for r in enh:
  p=r['Partition_'+mode]
  if p not in sets:continue
  g=r['Graph_ID'];sets[p]['graph'].add(g);sets[p]['connection'].add(M[g]['Connectivity_SMILES']);sets[p]['group'].add(M[g][field]);sets[p]['study'].add(r['Study_key'])
 for x,y in [('TRAIN','CALIBRATION'),('TRAIN','TEST'),('CALIBRATION','TEST')]:
  for metric in sets[x]:
   overlap=len(sets[x][metric]&sets[y][metric]);assert overlap==0
   checks.append(dict(Split_mode=mode,A=x,B=y,Identity=metric,Overlap=overlap))
 assert sum(collections.Counter(r['Partition_'+mode] for r in enh).values())==len(L)
write('leakage_checks.csv',checks)
summary=dict(date='2026-11-18',rdkit=rdBase.rdkitVersion,technical_structures=len(M),point_fields=len(L),candidate_event_signatures=len(groups),confirmed_independent_events='NOT_ESTABLISHED',confirmed_independent_series='NOT_ESTABLISHED',reference_traceable_fields=sum(r['Source_link_status']=='TRACEABLE_REFERENCE' for r in enh),task_defined_traceable_fields=sum(r['Split_eligible'] for r in enh),task_defined_traceable_structures=len({r['Graph_ID'] for r in enh if r['Split_eligible']}),same_study_task_value_review_groups=len(conflicts),fingerprint_pairs=pairs,fingerprint_edges_ge06=near_edges,scaffold_connectivity_groups=len({m['Scaffold_connectivity_group'] for m in M.values()}),fp06_components=len({m['Scaffold_FP06_component'] for m in M.values()}),leakage_assertions='PASS',partitions={mode:dict(collections.Counter(r['Partition_'+mode] for r in enh)) for mode in ['SC','FP06']},interpretation='Pre-model technical/task/source census and proposed partitions, not final primary-validated model-ready dataset.')
(OUT/'stage2_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
summary['stage2_point_fields']=sum(r['Stage2_point_eligible'] for r in enh)
summary['stage2_point_structures']=len({r['Graph_ID'] for r in enh if r['Stage2_point_eligible']})
summary['held_point_fields']=sum(not r['Stage2_point_eligible'] for r in enh)
summary['held_only_structures']=len(M)-summary['stage2_point_structures']
(OUT/'stage2_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
manifest=[]
for p in sorted(INPUTS):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(2**20),b''):h.update(block)
 manifest.append(dict(Path=str(p),Bytes=p.stat().st_size,SHA256=h.hexdigest()))
write('input_manifest.csv',manifest)
print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)
