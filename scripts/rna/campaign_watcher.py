#!/usr/bin/env python3
"""Detached RNA-seq campaign watcher (read-only).

Appends one JSON status snapshot per poll interval to
logs/campaign_monitor/<date>_snapshots.jsonl and writes
logs/campaign_monitor/PRODUCER_EXITED (then exits) when the producer
process is gone. Never signals, writes to, or otherwise touches the
live data root; the status script is invoked read-only.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN_PROJECT = PROJECT_ROOT / "projects" / "hymenoptera_amalgkit"
DATA_ROOT = Path("/Volumes/external_drive/Data/amalgkit")
PRODUCER_PATTERN = "[r]un_all_species.py"
POLL_SECONDS = 1800


def producer_alive() -> bool:
    probe = subprocess.run(
        ["pgrep", "-f", PRODUCER_PATTERN],
        capture_output=True,
        text=True,
        check=False,
    )
    if not probe.stdout.strip():
        return False
    pids = ",".join(probe.stdout.split())
    verify = subprocess.run(
        ["ps", "-p", pids, "-o", "command="],
        capture_output=True,
        text=True,
        check=False,
    )
    return str(DATA_ROOT) in verify.stdout


def snapshot(out_dir: Path, stamp: str) -> None:
    result = subprocess.run(
        [
            "env",
            "-u",
            "VIRTUAL_ENV",
            "uv",
            "run",
            "python",
            "scripts/report_campaign_status.py",
            "--data-root",
            str(DATA_ROOT),
            "--json",
        ],
        cwd=CAMPAIGN_PROJECT,
        capture_output=True,
        text=True,
        check=False,
        timeout=3600,
    )
    if result.returncode != 0:
        err = result.stderr.strip().splitlines()[-1:] or ["unknown error"]
        line = {"captured_at": stamp, "error": f"status rc={result.returncode}: {err[0]}"}
    else:
        payload = json.loads(result.stdout)
        counts = payload.get("state_counts") or {}
        line = {
            "captured_at": stamp,
            "state_counts": counts,
            "downloaded_fastq_gb": payload.get("downloaded_fastq_gb"),
            "derived_success_rate_per_hour": payload.get("derived_success_rate_per_hour"),
        }
    with (out_dir / f"{stamp[:10]}_snapshots.jsonl").open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(line, sort_keys=True) + "\n")


def main() -> int:
    mon = CAMPAIGN_PROJECT / "logs" / "campaign_monitor"
    mon.mkdir(parents=True, exist_ok=True)
    print(f"watcher started poll={POLL_SECONDS}s root={DATA_ROOT}", flush=True)
    while True:
        stamp = time.strftime("%Y-%m-%dT%H:%M:%S%z")
        if not producer_alive():
            with (mon / f"{stamp[:10]}_snapshots.jsonl").open("a", encoding="utf-8") as fh:
                fh.write(json.dumps({"captured_at": stamp, "producer_exited": True}) + "\n")
            (mon / "PRODUCER_EXITED").write_text(f"{stamp}\n", encoding="utf-8")
            print(f"producer exited at {stamp}", flush=True)
            return 0
        try:
            snapshot(mon, stamp)
        except Exception as exc:  # keep the loop alive on transient I/O errors
            with (mon / f"{stamp[:10]}_snapshots.jsonl").open("a", encoding="utf-8") as fh:
                fh.write(json.dumps({"captured_at": stamp, "error": str(exc)[:200]}) + "\n")
        time.sleep(POLL_SECONDS)


if __name__ == "__main__":
    sys.exit(main())
