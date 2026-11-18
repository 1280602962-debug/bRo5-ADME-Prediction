"""Verify release asset counts, sizes, server SHA-256 and snapshot inventories."""
from pathlib import Path
import hashlib
import json
import requests

REPO = "1280602962-debug/bRo5-ADME-Prediction"
TAG = "data-snapshot-20261118"
API = "https://api.github.com/repos/" + REPO
session = requests.Session()
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent if SCRIPT_DIR.name == "scripts" else SCRIPT_DIR
CURATED = ROOT / "data" / "curated"
if not CURATED.exists():
    CURATED = ROOT / "curated"


def get_json(url):
    response = session.get(url, timeout=60)
    response.raise_for_status()
    return response.json()


release = get_json(API + "/releases/tags/" + TAG)
assets = []
for page in range(1, 100):
    batch = get_json(API + f"/releases/{release['id']}/assets?per_page=100&page={page}")
    assets.extend(batch)
    if len(batch) < 100:
        break
by_name = {asset["name"]: asset for asset in assets}


def release_json(name):
    response = session.get(by_name[name]["browser_download_url"], timeout=60)
    response.raise_for_status()
    return response.json()


complete = release_json("UPLOAD_VERIFIED.json")
assert complete["Status"] == "COMPLETE"
expected = release_json("portable_checksums.json") + release_json("official_asset_checksums.json")
for item in expected:
    asset = by_name[item["Name"]]
    assert asset["state"] == "uploaded", item["Name"]
    assert asset["size"] == item["Bytes"], item["Name"]
    assert asset.get("digest") == "sha256:" + item["SHA256"], item["Name"]
supplement = release_json("supplement_checksums.json")
for item in supplement:
    if item["Path"].startswith("archives/"):
        asset = by_name[Path(item["Path"]).name]
        assert asset["size"] == item["Bytes"]
        assert asset.get("digest") == "sha256:" + item["SHA256"]

head = get_json(API + "/git/ref/heads/main")["object"]["sha"]
tree = get_json(API + "/git/trees/" + head + "?recursive=1")
assert not tree.get("truncated")
paths = {item["path"]: item for item in tree["tree"] if item["type"] == "blob"}
assert len([p for p in paths if p.startswith("snapshot/parts/")]) == complete["Portable_parts"]
for path in paths:
    assert not any(s in path for s in ["Dual_Target_Docking", "NLRP3", "JNK1", "URAT1"]), path
for item in supplement:
    if item["Path"].startswith("data/curated/"):
        local = CURATED / Path(item["Path"]).name
        data = local.read_bytes()
        git_sha = hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()
        assert paths[item["Path"]]["sha"] == git_sha
result = dict(Status="PASS", Repository=REPO, Commit=head, Logical_source_files=complete["Logical_files"],
              Release_assets=len(assets), Verified_data_asset_digests=len(expected),
              Curated_git_files_verified=3, Portable_parts=complete["Portable_parts"],
              Official_parts=complete["Official_assets"], Unrelated_project_paths=0)
(ROOT / "remote_verification.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result))
