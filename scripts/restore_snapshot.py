"""Restore a GitHub project snapshot, then verify every original file hash.

Usage: python restore_snapshot.py --destination restored_bro5
Requires GitHub CLI (gh) to download release assets. Local existing assets can be
used with --assets PATH.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import io
import json
import shutil
import subprocess
import tarfile

parser = argparse.ArgumentParser()
parser.add_argument("--destination", type=Path, required=True)
parser.add_argument("--assets", type=Path)
args = parser.parse_args()
target = args.destination.resolve()
target.mkdir(parents=True, exist_ok=True)
assets = args.assets or target / "_downloaded_assets"
assets.mkdir(parents=True, exist_ok=True)
if not args.assets:
    subprocess.run(["gh", "release", "download", "data-snapshot-20261118", "--repo",
                    "1280602962-debug/bRo5-ADME-Prediction", "--dir", str(assets), "--skip-existing"], check=True)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for data in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(data)
    return h.hexdigest()


class Joined(io.RawIOBase):
    def __init__(self, paths):
        self.paths = iter(paths)
        self.stream = None
    def readable(self):
        return True
    def readinto(self, buffer):
        while True:
            if self.stream is None:
                path = next(self.paths, None)
                if path is None:
                    return 0
                self.stream = path.open("rb")
            size = self.stream.readinto(buffer)
            if size:
                return size
            self.stream.close()
            self.stream = None


parts = json.loads((assets / "portable_checksums.json").read_text())
for item in parts:
    path = assets / item["Name"]
    assert path.stat().st_size == item["Bytes"] and digest(path) == item["SHA256"], path
with io.BufferedReader(Joined([assets / item["Name"] for item in parts])) as stream:
    with tarfile.open(fileobj=stream, mode="r|gz") as archive:
        archive.extractall(target, filter="data")
official = json.loads((assets / "official_asset_checksums.json").read_text())
for original in sorted({item["Original_path"] for item in official}):
    output = target / original
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as destination:
        for item in sorted((x for x in official if x["Original_path"] == original), key=lambda x: x["Name"]):
            path = assets / item["Name"]
            assert digest(path) == item["SHA256"], path
            with path.open("rb") as source:
                shutil.copyfileobj(source, destination)
rows = list(csv.DictReader((assets / "snapshot_manifest.csv").open(encoding="utf-8-sig")))
for row in rows:
    if row["Path"] != row["Canonical_path"]:
        destination = target / row["Path"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(target / row["Canonical_path"], destination)
for row in rows:
    path = target / row["Path"]
    assert path.stat().st_size == int(row["Bytes"]), path
    assert digest(path) == row["SHA256"], path
print("Verified original files", len(rows))
