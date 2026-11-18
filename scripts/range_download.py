"""Download an official archive using validated independent byte ranges."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import json
import subprocess
import time
import urllib.request


def download(url, destination, size, tag):
    destination = Path(destination)
    block_size = 64 * 1024 * 1024
    intervals = [(start, min(start + block_size, size) - 1) for start in range(0, size, block_size)]
    with destination.open("wb") as stream:
        stream.truncate(size)

    def fetch(interval):
        start, end = interval
        for attempt in range(4):
            try:
                request = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}", "User-Agent": "bRo5-research-data-archive"})
                with urllib.request.urlopen(request, timeout=120) as response:
                    expected = f"bytes {start}-{end}/{size}"
                    if response.status != 206 or response.headers.get("Content-Range") != expected:
                        raise RuntimeError("Official server did not honor the requested byte range")
                    count = 0
                    with destination.open("r+b") as output:
                        output.seek(start)
                        while True:
                            data = response.read(min(8 * 1024 * 1024, end - start + 1 - count))
                            if not data:
                                break
                            output.write(data)
                            count += len(data)
                            if count == end - start + 1:
                                break
                    if count != end - start + 1:
                        raise RuntimeError("Incomplete range response")
                return count
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(2 ** attempt)

    finished = 0
    with ThreadPoolExecutor(max_workers=8) as pool:
        futures = [pool.submit(fetch, interval) for interval in intervals]
        for index, future in enumerate(as_completed(futures), start=1):
            finished += future.result()
            print(destination.name, finished, "/", size, flush=True)
            if index % 16 == 0 or index == len(intervals):
                progress = Path("official_download_progress.json")
                progress.write_text(json.dumps(dict(Archive=destination.name, Bytes_downloaded=finished, Total_bytes=size, Full_hash_verified=False)))
                subprocess.run(["gh", "release", "upload", tag, str(progress), "--clobber"], check=True, capture_output=True)
