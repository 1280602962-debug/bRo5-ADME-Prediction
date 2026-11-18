from pathlib import Path
import csv, gzip, hashlib, io, json, tarfile

root=Path(__file__).resolve().parent
delivery=root/'github_delivery_20261118'
out=delivery/'portable'
out.mkdir(exist_ok=True)
entries=list(csv.DictReader((delivery/'snapshot_manifest.csv').open(encoding='utf-8-sig')))
roots={'collected':Path(r'C:\Users\guwenjing\bRo5 ADME Prediction'),
       'early_audit':Path(r'C:\Users\guwenjing\.codex\.chatgpt-projects\g-p-6a9522ad97f0819199b8f5924494f4c6'), 'current':root}
databases=[]
selected=[]
for row in entries:
    if row['Path'] != row['Canonical_path']:
        continue
    if row['Path'].endswith(('chembl_34_sqlite.tar.gz','chembl_37_sqlite.tar.gz')):
        version=34 if 'chembl_34_' in row['Path'] else 37
        databases.append({**row,'URL':f'https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_{version}/chembl_{version}_sqlite.tar.gz'})
    else:
        prefix,relative=row['Path'].split('/',1)
        selected.append((roots[prefix]/relative,row['Path']))
(delivery/'official_databases.json').write_text(json.dumps(databases,indent=2),encoding='utf-8')

class Writer(io.RawIOBase):
    def __init__(self):
        self.stream=None;self.index=0;self.size=0;self.entries=[];self.total=0
    def writable(self):return True
    def write(self,block):
        count=len(block)
        while block:
            if self.stream is None or self.size==16*1024*1024:
                self.finish();self.index+=1;self.path=out/f'project_data.tar.gz.part{self.index:03d}'
                self.stream=self.path.open('wb');self.sha=hashlib.sha256();self.size=0
            part=block[:16*1024*1024-self.size];self.stream.write(part);self.sha.update(part)
            self.size+=len(part);self.total+=len(part);block=block[len(part):]
        return count
    def finish(self):
        if self.stream:
            self.stream.close();self.entries.append(dict(Name=self.path.name,Bytes=self.size,SHA256=self.sha.hexdigest()));self.stream=None
            print(self.entries[-1],flush=True)
w=Writer()
with gzip.GzipFile(fileobj=w,mode='wb',compresslevel=1,mtime=0) as compressed:
    with tarfile.open(mode='w|',fileobj=compressed) as archive:
        for i,(path,name) in enumerate(selected):
            archive.add(path,arcname=name,recursive=False)
            if i%150==0:print('FILES',i,flush=True)
w.finish()
(delivery/'portable_checksums.json').write_text(json.dumps(w.entries,indent=2),encoding='utf-8')
print('COMPLETE',len(w.entries),w.total,flush=True)
