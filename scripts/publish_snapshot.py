"""Verify uploaded project chunks and publish the full raw-data snapshot."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import urllib.request

TAG = "data-snapshot-20261118"
ROOT = Path.cwd()
SNAPSHOT = ROOT / "snapshot"


def gh(*args, check=True):
    return subprocess.run(["gh", *map(str, args)], check=check, capture_output=True, text=True)


def sha256(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for data in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            result.update(data)
    return result.hexdigest()


complete = json.loads((SNAPSHOT / "UPLOAD_COMPLETE.json").read_text())
parts = json.loads((SNAPSHOT / "portable_checksums.json").read_text())
assert complete["Portable_parts"] == len(parts)
for item in parts:
    path = SNAPSHOT / "parts" / item["Name"]
    assert path.stat().st_size == item["Bytes"], path.name
    assert sha256(path) == item["SHA256"], path.name
print("Project archive chunk checks passed", flush=True)
if gh("release", "view", TAG, check=False).returncode:
    gh("release", "create", TAG, "--target", os.environ["GITHUB_SHA"],
       "--title", "Complete bRo5 permeability data snapshot — 2026-11-18",
       "--notes", "Project data, literature attachments, extracted records, curation audits and scripts. Restore with scripts/restore_snapshot.py and verify snapshot_manifest.csv. Original ChEMBL archives are reproduced from their recorded official URLs and checked against local SHA-256 fingerprints.")
for item in parts:
    gh("release", "upload", TAG, SNAPSHOT / "parts" / item["Name"], "--clobber")
for name in ["snapshot_manifest.csv", "snapshot_summary.json", "portable_checksums.json", "official_databases.json"]:
    gh("release", "upload", TAG, SNAPSHOT / name, "--clobber")
gh("release", "upload", TAG, ROOT / "scripts" / "restore_snapshot.py", "--clobber")

official_assets = []
for record in json.loads((SNAPSHOT / "official_databases.json").read_text()):
    name = Path(record["Path"]).name
    download = ROOT / name
    print("Downloading official archive", name, flush=True)
    subprocess.run(["curl", "--fail", "--location", "--retry", "5", "--retry-all-errors",
                    "--output", str(download), record["URL"]], check=True)
    assert download.stat().st_size == int(record["Bytes"]), name
    assert sha256(download) == record["SHA256"], name
    print("Official archive matches local fingerprint", name, flush=True)
    with download.open("rb") as stream:
        for index in range(1, 100):
            data = stream.read(1900 * 1024 * 1024)
            if not data:
                break
            output = ROOT / f"{name}.part{index:03d}"
            output.write_bytes(data)
            item = dict(Name=output.name, Bytes=len(data), SHA256=hashlib.sha256(data).hexdigest(), Original_path=record["Path"])
            del data
            gh("release", "upload", TAG, output, "--clobber")
            official_assets.append(item)
            output.unlink()
    download.unlink()
(ROOT / "official_asset_checksums.json").write_text(json.dumps(official_assets, indent=2))
gh("release", "upload", TAG, ROOT / "official_asset_checksums.json", "--clobber")
(ROOT / "UPLOAD_VERIFIED.json").write_text(json.dumps(dict(Status="COMPLETE", Logical_files=complete["Logical_files"],
    Portable_parts=len(parts), Official_assets=len(official_assets), Local_fingerprints_matched=True), indent=2))
gh("release", "upload", TAG, ROOT / "UPLOAD_VERIFIED.json", "--clobber")
print("COMPLETE — all project data assets published and verified", flush=True)
