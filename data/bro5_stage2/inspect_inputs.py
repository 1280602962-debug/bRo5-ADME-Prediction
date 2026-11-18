import csv,json,gzip,re
from pathlib import Path
B=Path(r'C:\Users\guwenjing\bRo5 ADME Prediction\bRo5_Permeability_Benchmark')
def rows(p):
 with p.open(encoding='utf-8-sig',newline='') as f:yield from csv.DictReader(f)
for n in ['03_extraction/NONCHEMBL_ASSAY_OCCURRENCES.csv','03_extraction/batch02/ALL_SOURCE_ASSAY_GRAPH_INDEX.csv']:
 d={}
 for r in rows(B/n):
  if r['Source'] not in d:d[r['Source']]=r
 print(json.dumps(d,ensure_ascii=False))
