"""Prepare a content-deduplicated snapshot for the independent GitHub repository."""
from pathlib import Path
import csv
import gzip
import hashlib
import io
import json
import shutil
import tarfile
import time

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "github_delivery_20261118"
OUT.mkdir(exist_ok=True)
DATA = Path(r"C:\Users\guwenjing\bRo5 ADME Prediction")
EARLY = Path(r"C:\Users\guwenjing\.codex\.chatgpt-projects\g-p-6a9522ad97f0819199b8f5924494f4c6")
BLOCK = 4 * 1024 * 1024


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(BLOCK), b""):
            h.update(block)
    return h.hexdigest()


files = []
for child in DATA.iterdir():
    if child.name == "_transfer":
        continue
    if child.is_file():
        files.append((child, "collected/" + child.name))
    else:
        files.extend((p, "collected/" + str(p.relative_to(DATA)).replace("\\", "/")) for p in child.rglob("*") if p.is_file())
for name in ["bro5_audit_20261117", "sources"]:
    base = EARLY / name
    files.extend((p, "early_audit/" + str(p.relative_to(EARLY)).replace("\\", "/")) for p in base.rglob("*") if p.is_file())
for name in ["bro5_stage2", "bro5_stage3", "current_design_audit_20261118"]:
    base = ROOT / name
    files.extend((p, "current/" + str(p.relative_to(ROOT)).replace("\\", "/")) for p in base.rglob("*") if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc")
for name in ["audit_current_design.py", "verify_current_design_audit.py", "package_current_design_audit.py"]:
    files.append((ROOT / name, "current/" + name))
assert len({arc for _, arc in files}) == len(files)
manifest = []
seen = {}
unique = []
start = time.monotonic()
for i, (path, arc) in enumerate(sorted(files, key=lambda item: item[1])):
    sha = digest(path)
    canonical = seen.setdefault(sha, arc)
    manifest.append(dict(Path=arc, Bytes=path.stat().st_size, SHA256=sha, Canonical_path=canonical))
    if arc == canonical:
        unique.append((path, arc))
    if i % 100 == 0:
        print("HASH", i, "/", len(files), flush=True)
with (OUT / "snapshot_manifest.csv").open("w", encoding="utf-8-sig", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(manifest[0]))
    writer.writeheader()
    writer.writerows(manifest)
summary = dict(Files=len(files), Logical_bytes=sum(x["Bytes"] for x in manifest),
               Unique_files=len(unique), Unique_bytes=sum(p.stat().st_size for p, _ in unique),
               Excluded=["collected/_transfer (redundant transport copies)", "early workspace work/tmp (includes unrelated projects)", "__pycache__ and .pyc"],
               Snapshot_complete=False)
(OUT / "snapshot_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
print("MANIFEST", json.dumps(summary), flush=True)


class Volumes(io.RawIOBase):
    def __init__(self):
        self.stream = None
        self.index = 0
        self.size = 0
        self.position = 0
        self.metadata = []
        self.limit = 1024 ** 3

    def writable(self):
        return True

    def tell(self):
        return self.position

    def write(self, block):
        original = len(block)
        while block:
            if self.stream is None or self.size == self.limit:
                self.finish_part()
                self.index += 1
                self.path = OUT / f"bro5_snapshot_20261118.tar.gz.part{self.index:03d}"
                self.stream = self.path.open("wb")
                self.sha = hashlib.sha256()
                self.size = 0
            chunk = block[:self.limit-self.size]
            self.stream.write(chunk)
            self.sha.update(chunk)
            self.size += len(chunk)
            self.position += len(chunk)
            block = block[len(chunk):]
        return original

    def finish_part(self):
        if self.stream:
            self.stream.close()
            self.metadata.append(dict(Name=self.path.name, Bytes=self.size, SHA256=self.sha.hexdigest()))
            self.stream = None
            print("VOLUME", self.metadata[-1], flush=True)


volume = Volumes()
with gzip.GzipFile(fileobj=volume, mode="wb", compresslevel=1, mtime=0) as compressed:
    with tarfile.open(fileobj=compressed, mode="w|") as archive:
        for i, (path, arc) in enumerate(unique):
            archive.add(path, arcname=arc, recursive=False)
            if i % 75 == 0 or path.stat().st_size > 100_000_000:
                print("ARCHIVED", i, "/", len(unique), arc, "GB", round(volume.position/1e9, 3), flush=True)
volume.finish_part()
(OUT / "volume_checksums.json").write_text(json.dumps(volume.metadata, indent=2), encoding="utf-8")
summary.update(Snapshot_complete=True, Volumes=len(volume.metadata), Compressed_bytes=volume.position,
               Elapsed_seconds=round(time.monotonic()-start))
(OUT / "snapshot_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
print("DONE", json.dumps(summary), flush=True)
