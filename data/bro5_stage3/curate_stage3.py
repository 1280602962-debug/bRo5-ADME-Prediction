import csv,json,re,hashlib,math,collections,zipfile,sys
from pathlib import Path
from rdkit import Chem,RDLogger
from rdkit.Chem import Descriptors,rdMolDescriptors
sys.stdout.reconfigure(encoding='utf-8');RDLogger.DisableLog('rdApp.*')
P=Path(__file__).resolve().parent;S=P.parent/'bro5_stage2';B=Path(r'C:\Users\guwenjing\bRo5 ADME Prediction\bRo5_Permeability_Benchmark')
def rows(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(n,rs):
 rs=list(rs)
 if not rs:return
 with (P/n).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]));w.writeheader();w.writerows(rs)
def hid(prefix,value):return prefix+'_'+hashlib.sha256(str(value).encode()).hexdigest()[:24]
M={r['Graph_ID']:r for r in rows(S/'molecule_groups.csv')};L=rows(S/'labels_task_source_audit.csv');legacy={r['AssayOccurrence_ID']:r for r in rows(B/'03_extraction/NONCHEMBL_ASSAY_OCCURRENCES.csv')}
E=json.loads((P/'primary_evidence.json').read_text(encoding='utf-8'))
changes=[];bounds=[];primary=[];lineages=[]
oldlead='G_cacbf9504e1f94e95bdacec8';newphe='G_e9503599b2a73e145750fbf8'
fixed=Chem.MolToSmiles(Chem.MolFromSmiles(M[oldlead]['SMILES'].replace('[CH2][Pb]','Cc9ccccc9')),isomericSmiles=True)
assert fixed==M[newphe]['SMILES']
value5={'dBET1':0.149,'1':0.437,'2a':0.348,'2b':1.2}
allnamed=collections.defaultdict(list)
for r in L:
 r['Original_Graph_ID']=r['Graph_ID'];r['Original_Task']=r['Task'];r['Original_Log10_cm_s']=r['Log10_cm_s']
 r['Stage3_status']='POINT_CANDIDATE' if r['Stage2_point_eligible']=='True' else 'STAGE2_HOLD'
 r['Primary_evidence_ID']='';r['Primary_point_check']='NOT_REVIEWED';r['Primary_protocol_check']='NOT_REVIEWED';r['Primary_identity_check']='NOT_REVIEWED';r['Reported_result_ID']='';r['Reported_series']='NOT_CONFIRMED';r['Recovery_primary']='UNKNOWN';r['Stage3_note']=''
 raw=json.loads(legacy[r['Source_row_ID']]['Source_raw_JSON']) if r['Source_row_ID'] in legacy else {}
 compound=r['Compound_original_ID'] or raw.get('compound_name','') or raw.get('Compound','') or raw.get('Original_ID','')
 r['Compound_original_ID_enriched']=compound
 study=r['Study_key'];graph=r['Graph_ID']
 if graph==oldlead:
  r['Graph_ID']=newphe;r['Primary_evidence_ID']='E3_001';r['Primary_identity_check']='PHE_FRAGMENT_CORRECTION_CONFIRMED_OTHER_STEREO_NOT_REVIEWED';r['Stage3_note']='Replaced impossible Pb-containing side chain with original Phe benzyl group; exact resulting canonical isomeric graph equals independently stored MZ3 graph.'
  changes.append(dict(Occurrence_ID=r['Occurrence_ID'],Action='CORRECT_GRAPH_ALIAS',Old=oldlead,New=newphe,Evidence_ID='E3_001'))
 if study=='DOI:10.1021/acs.jmedchem.2c00877':
  r['Primary_evidence_ID']='E3_002'
  if r['Endpoint']=='Caco-2' and r['Source']=='ChEMBL37' and compound in ['1','2','3']:
   expected={'1':30,'2':11,'3':6}[compound]
   assert abs(float(r['Log10_cm_s'])-(math.log10(expected)-7))<1e-6
   r['Task']='Caco-2:GEOMETRIC_MEAN_AB_BA';r['Direction_refined']='GEOMETRIC_MEAN_AB_BA';r['Primary_point_check']='PRIMARY_TABLE_VALUE_AND_SCALE_CONFIRMED';r['Primary_protocol_check']='PH74_AND_GEOMETRIC_MEAN_DEFINITION_CONFIRMED_SI_PENDING';r['Reported_result_ID']=hid('RME',(study,'Table2',compound,r['Task']));r['Reported_series']='PRIMARY_BR D4_CRBN_LINKER_SERIES'.replace('BR D4','BRD4')
   changes.append(dict(Occurrence_ID=r['Occurrence_ID'],Action='RESOLVE_GEOMETRIC_MEAN_TASK',Old=r['Original_Task'],New=r['Task'],Evidence_ID='E3_002'))
  if r['Endpoint']=='PAMPA':
   v=float(r['Log10_cm_s'])
   if abs(v+7.37)<1e-5:
    r['Log10_cm_s']='-7.37';r['Relation']='<';r['Stage3_status']='PRIMARY_BOUND';r['Primary_point_check']='PRIMARY_CENSORED_NOT_POINT';r['Reported_result_ID']=hid('RME',(study,'Table2','3','PAMPA'));r['Stage3_note']='Original minus-log Pe>7.37 becomes log10 Pe<-7.37. No raw repeat count inferred.'
    bounds.append(dict(Occurrence_ID=r['Occurrence_ID'],Graph_ID=r['Graph_ID'],Study_key=study,Compound='3',Endpoint='PAMPA',Relation='<',Log10_bound_cm_s=-7.37,Reported_result_ID=r['Reported_result_ID'],Recovery_status='UNKNOWN',Reliability='PRIMARY_BOUND_SEMANTICS_CONFIRMED_FULL_VALIDITY_PENDING',Evidence_ID='E3_002'))
    changes.append(dict(Occurrence_ID=r['Occurrence_ID'],Action='POINT_TO_PRIMARY_BOUND',Old=r['Original_Log10_cm_s'],New='<-7.37',Evidence_ID='E3_002'))
   elif abs(v+6.56)<1e-5:
    if r['Log10_cm_s']!='-6.56':changes.append(dict(Occurrence_ID=r['Occurrence_ID'],Action='NORMALIZE_PRIMARY_REPORTED_PRECISION',Old=r['Log10_cm_s'],New='-6.56',Evidence_ID='E3_002'))
    r['Log10_cm_s']='-6.56';r['Primary_point_check']='PRIMARY_TABLE_VALUE_AND_SCALE_CONFIRMED';r['Reported_result_ID']=hid('RME',(study,'Table2','1','PAMPA'));r['Stage3_note']='One reported mean; source indicates three repeats, not three recoverable raw events.'
   r['Reported_series']='PRIMARY_BRD4_CRBN_LINKER_SERIES'
 if study=='DOI:10.1021/acs.jmedchem.5c02791' and r['Source']=='Hudec_new_parsed':
  assert compound in value5
  assert abs(float(r['Log10_cm_s'])-(math.log10(value5[compound])-6))<1e-6
  r['Primary_evidence_ID']='E3_003';r['Primary_point_check']='PRIMARY_TABLE_VALUE_AND_SCALE_CONFIRMED';r['Primary_protocol_check']='PRIMARY_METHOD_CONFIRMED';r['Primary_identity_check']='SOURCE_NAME_FORMULA_TRANSFORMATION_AND_LINKER_CIP_CONSISTENT_FULL_ATOM_REVIEW_PENDING';r['Reported_result_ID']=hid('RME',(study,'Table1',compound,'PAMPA'));r['Reported_series']='PRIMARY_DBET1_ESTER_METHYL_LINKER_SERIES';r['Stage3_note']='Four distinct source samples. Preserve new R/S linker center; unspecified CRBN center is not resolved or enumerated. Two independent experiments and four technical replicates do not expose individual raw-event rows.'
  primary.append(dict(Reported_result_ID=r['Reported_result_ID'],Source_occurrence=r['Occurrence_ID'],DOI=study[4:],Original_compound=compound,Graph_ID=r['Graph_ID'],Papp_1e6_cm_s=value5[compound],Log10_cm_s=r['Log10_cm_s'],Donor_pH=6.4,Acceptor_pH=7.4,Temperature_C=37,Incubation_h=4,Concentration_uM=10,DMSO_percent=5,Membrane='Corning_Gentest_precoated',Agitation='NO',Technical_replicates=4,Reported_experiments_for_mean=2,Recovery='UNKNOWN',Raw_independent_event_count='NOT_RECOVERED',Source_series='dBET1_ester_methyl_linker',Full_atom_identity='PENDING',ThreeD_unique_stereo='NO',Evidence_ID='E3_003'))
 if study=='DOI:10.1016/j.redox.2022.102552':
  r['Primary_evidence_ID']='E3_004';r['Primary_protocol_check']='PRIMARY_METHOD_CONFIRMED';r['Stage3_note']='Methods verified; Table S2 values, qualifiers and per-compound recovery unavailable. Keep numerical eligibility unvalidated.';r['Reported_series']='PRIMARY_KEAP1_CRBN_LINKER_SERIES'
 if study=='DOI:10.1016/j.chembiol.2024.01.003':
  r['Primary_evidence_ID']='E3_005';r['Primary_protocol_check']='PRIMARY_METHOD_CONFIRMED';r['Stage3_note']='Caco-2 protocol verified. Supplementary Figure10g numerical rows and recovery not read; record-level database citation alone cannot confirm original measurement provenance.'
 if study=='DOI:10.1021/acsmedchemlett.0c00265':
  r['Primary_evidence_ID']='E3_001'
  aliases=[x.strip().upper() for x in re.split(r'[;|]',compound) if x.strip()]
  family=next((f for f in ['MZP','MZ','AT','CMP','CM','SL-X'] if any(a.startswith(f) for a in aliases)),None)
  if family:r['Reported_series']='PRIMARY_VHL_'+('CM_CMP' if family in ['CM','CMP'] else family)
  if compound in ['MZ3','9; MZ3']:r['Primary_identity_check']='PHE_SIDECHAIN_CONFIRMED_REMAINDER_PENDING'
 # Distinct naming in the same study is preserved. Aliases only connect exact tokens.
 aliases=[re.sub(r'\s+',' ',x.strip().lower()) for x in re.split(r'[;|]',compound) if x.strip()]
 aliases=[re.sub(r'^protac\s+','',x) for x in aliases if not re.search(r'mixture|racem|diastereo|enantiomer',x)]
 if aliases and study:allnamed[study].append((r,aliases))
# Named sample identity reconciliation: a source sample is not an independent assay event.
identity=[]
for study,rr in sorted(allnamed.items()):
 parent={a:a for _,aa in rr for a in aa}
 def find(a):
  while parent[a]!=a:parent[a]=parent[parent[a]];a=parent[a]
  return a
 for _,aa in rr:
  for a in aa[1:]:parent[find(a)]=find(aa[0])
 groups=collections.defaultdict(list)
 for r,aa in rr:groups[find(aa[0])].append(r)
 for _,rs in groups.items():
  gs={r['Graph_ID'] for r in rs};conns={M[g]['Connectivity_SMILES'] for g in gs};sources={r['Source'] for r in rs};names=sorted({r['Compound_original_ID_enriched'] for r in rs})
  candidate=hid('NS',(study,sorted(names)))
  for r in rs:r['Named_sample_candidate_ID']=candidate
  if len(gs)>1:
   identity.append(dict(Named_sample_candidate_ID=candidate,Study_key=study,Original_names='|'.join(names),Graph_count=len(gs),Connectivity_count=len(conns),Graph_IDs='|'.join(sorted(gs)),Sources='|'.join(sorted(sources)),Resolution='CONNECTIVITY_CONFLICT_REVIEW' if len(conns)>1 else 'STEREO_REPRESENTATION_REVIEW',Do_not_interpret_as='CONFIRMED_DISTINCT_MOLECULES_OR_CONFIRMED_SAME_ENTITY'))
for r in L:r.setdefault('Named_sample_candidate_ID','')
write('source_sample_identity_conflicts.csv',identity)
write('labels_primary_overlay.csv',L)
write('primary_numeric_protocol_confirmed_means.csv',primary)
write('new_primary_bound_occurrences.csv',bounds)
write('primary_change_log.csv',changes)
write('point_candidates_after_primary_overlay.csv',[r for r in L if r['Stage3_status']=='POINT_CANDIDATE'])
write('nonpoint_or_held_occurrences.csv',[r for r in L if r['Stage3_status']!='POINT_CANDIDATE'])
stats=[]
def scope(g,s):
 m=M[g];mw=float(m['MW']);non=m['Topology']=='NONMACRO_SCREEN'
 return s=='ALL_MW' or s=='MW_GT500' and mw>500 or s=='MW_GT700_NONMACRO' and mw>700 and non or s=='MW_700_1000_NONMACRO' and 700<mw<=1000 and non
for s in ['ALL_MW','MW_GT500','MW_GT700_NONMACRO','MW_700_1000_NONMACRO']:
 pool=[r for r in L if r['Stage3_status']=='POINT_CANDIDATE' and scope(r['Graph_ID'],s)]
 for task in ['ALL_TASKS','Caco-2:AB','PAMPA:DIRECTION_NOT_APPLICABLE','Caco-2:GEOMETRIC_MEAN_AB_BA']:
  rs=[r for r in pool if task=='ALL_TASKS' or r['Task']==task]
  stats.append(dict(Scope=s,Task=task,Point_occurrences=len(rs),Graph_representations=len({r['Graph_ID'] for r in rs}),Primary_numeric_confirmed_graphs=len({r['Graph_ID'] for r in rs if r['Primary_point_check']=='PRIMARY_TABLE_VALUE_AND_SCALE_CONFIRMED'}),Final_model_ready='NOT_FROZEN'))
write('stage3_counts.csv',stats)
# Refine protocols without reassigning partitions or selecting a new seed.
protocol_queues=[]
for r in L:
 if r['Stage3_status']!='POINT_CANDIDATE':continue
 text=r['Assay_description'];flags=[]
 if r['Task']=='PAMPA:DIRECTION_NOT_APPLICABLE':
  if re.search(r'\bBBB\b|blood[- ]brain|porcine brain|brain polar lipid',text,re.I):flags.append('PAMPA_BBB_MENTION')
  if re.search(r'\bskin\b|dermal',text,re.I):flags.append('PAMPA_SKIN_MENTION')
  if re.search(r'\bGIT\b|gastrointestinal',text,re.I):flags.append('PAMPA_GIT_MENTION')
 if r['Endpoint']=='Caco-2' and re.search(r'followed by ACN|acetonitrile',text,re.I):flags.append('ACN_EXTRACTION_TIMING_PROTOCOL_REVIEW')
 if flags:protocol_queues.append(dict(Occurrence_ID=r['Occurrence_ID'],Graph_ID=r['Graph_ID'],Study_key=r['Study_key'],Task=r['Task'],Protocol_flags='|'.join(flags),Assay_description=text,Status='RULE_FLAG_NOT_PRIMARY_PROTOCOL_APPROVAL'))
write('protocol_subtype_review.csv',protocol_queues)
assert len(L)==66771 and len({r['Occurrence_ID'] for r in L})==66771
assert all(not (r['Stage3_status']=='POINT_CANDIDATE' and r['Relation']!='=') for r in L)
assert all(r['Graph_ID']!=oldlead for r in L)
assert len(primary)==4 and len({r['Reported_result_ID'] for r in primary})==4
assert len(bounds)==2 and len({r['Reported_result_ID'] for r in bounds})==1
assert len({r['Task'] for r in L if r['Source']=='ChEMBL37' and r['Source_row_ID'] in ['24698084','24698085','24698086']})==1
summary=dict(status='PASS',input_occurrences=len(L),primary_studies_with_targeted_evidence=len(E),priority1_studies=31,priority1_graphs=116,changed_occurrences=len({r['Occurrence_ID'] for r in changes}),change_actions=dict(collections.Counter(r['Action'] for r in changes)),new_bound_occurrences=len(bounds),new_bound_reported_means=1,primary_numeric_protocol_confirmed_means=4,primary_numeric_confirmed_occurrences=sum(r['Primary_point_check']=='PRIMARY_TABLE_VALUE_AND_SCALE_CONFIRMED' for r in L),named_sample_identity_conflict_groups=len(identity),connectivity_conflict_groups=sum(r['Connectivity_count']>1 for r in identity),counts=stats,full_scheme_complete=False,final_qualified_molecules='NOT_ESTABLISHED',reason='Primary full-text/SI, individual source-sample identities and raw independent-event lineage remain unresolved; no model fit and no split re-optimization.')
(P/'stage3_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
