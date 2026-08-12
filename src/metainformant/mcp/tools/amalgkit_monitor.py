#!/usr/bin/env python3
"""
Amalgkit Status Monitor (MCP Tool).

This script functions as a Model Context Protocol (MCP) tool.
It parses the Amalgkit pipeline logs and process table to return
structured JSON status information for AI agents.

Usage:
    python3 -m metainformant.mcp.tools.amalgkit_monitor

Output (JSON):
    {
        "status": "running" | "stalled" | "completed" | "error",
        "pid": 12345,
        "workers": 20,
        "progress": {
            "processed": 35,
            "total": 5615,
            "percent": 0.62
        },
        "last_log": "[35/5615] ✓ SRR26150181 (133.1s)...",
        "system": {
            "load_avg": [15.5, 12.0, 10.0],
            "free_disk_gb": 450.2
        }
    }
"""

import json
import os
import re
import shutil
import sqlite3
from argparse import ArgumentParser
from pathlib import Path

import psutil

# Configuration - could be dynamic in future
WORK_DIR = Path("output/amalgkit")
LOG_FILE = Path("output/amalgkit/run_all_species_incremental.log")


def is_pipeline_cmdline(cmdline):
    """Return True when a process command line looks like an Amalgkit workflow."""
    if any("amalgkit_monitor" in arg for arg in cmdline):
        return False
    if any("run_all_species" in arg for arg in cmdline):
        return True
    if any(
        marker in arg
        for arg in cmdline
        for marker in ("run_all_species.py", "process_species.py", "run_full_campaign.sh")
    ):
        return True
    return any(Path(arg).name == "amalgkit" for arg in cmdline)


def get_process_status():
    """Find the active Amalgkit process."""
    current_pid = os.getpid()
    for proc in psutil.process_iter(["pid", "name", "cmdline"]):
        try:
            if proc.info["pid"] == current_pid:
                continue
            cmdline = proc.info["cmdline"] or []
            if is_pipeline_cmdline(cmdline):
                return {"running": True, "pid": proc.info["pid"], "cmd": " ".join(cmdline)}
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue
    return {"running": False, "pid": None}


def parse_log_progress(log_file: Path = LOG_FILE):
    """Parse the last relevant line from the log file."""
    if not log_file.exists():
        return {"processed": 0, "total": 0, "last_line": "No log file found"}

    try:
        # Read last few lines to find progress
        # Using tail-like logic
        with open(log_file, "rb") as f:
            try:
                f.seek(-2000, os.SEEK_END)
            except OSError:
                f.seek(0)

            lines = f.readlines()

            # Decode and look for progress pattern: [35/5615]
            # Pattern: matches [123/456] key
            progress_re = re.compile(r"\[(\d+)/(\d+)\]")

            processed = 0
            total = 0
            last_line = ""

            for line in reversed(lines):
                line_str = line.decode("utf-8", errors="ignore").strip()
                if not line_str:
                    continue

                match = progress_re.search(line_str)
                if match:
                    processed = int(match.group(1))
                    total = int(match.group(2))
                    last_line = line_str
                    break

            return {
                "processed": processed,
                "total": total,
                "percent": round((processed / total) * 100, 2) if total > 0 else 0,
                "last_line": last_line,
            }

    except Exception as e:
        return {"error": str(e)}


def get_system_stats(work_dir: Path = WORK_DIR):
    """Get system load and disk space."""
    load = os.getloadavg()

    # Check disk space on work dir volume
    try:
        stat = os.statvfs(work_dir)
        free_gb = (stat.f_frsize * stat.f_bavail) / (1024**3)
    except OSError:
        free_gb = 0

    return {"load_avg": load, "free_disk_gb": round(free_gb, 1)}


def _readiness_snapshot(data_root: Path) -> dict[str, object]:
    """Return database-backed readiness without promoting scientific claims."""

    db_path = data_root / "pipeline_progress.db"
    executable_ready = shutil.which("amalgkit") is not None and shutil.which("kallisto") is not None
    result: dict[str, object] = {
        "executable": "ready" if executable_ready else "incomplete",
        "cohort": "unknown",
        "descriptive_analysis": "withheld",
        "biological_inference": "withheld",
        "database": str(db_path),
    }
    if not db_path.is_file():
        return result

    try:
        connection = sqlite3.connect(db_path)
        rows = connection.execute("SELECT species, state, COUNT(*) FROM samples GROUP BY species, state").fetchall()
        connection.close()
        counts: dict[str, dict[str, int]] = {}
        for species, state, count in rows:
            counts.setdefault(species, {})[state] = count
        total = sum(sum(values.values()) for values in counts.values())
        quantified = sum(values.get("quantified", 0) for values in counts.values())
        unresolved = sum(
            values.get(state, 0)
            for values in counts.values()
            for state in ("pending", "downloading", "downloaded", "quantifying")
        )
        failed = sum(values.get(state, 0) for values in counts.values() for state in ("failed", "quarantined"))
        cohort_ready = total > 0 and quantified == total and unresolved == 0 and failed == 0
        result["cohort"] = "ready" if cohort_ready else "partial_or_unresolved"

        required_steps = {"merge", "wsfilter", "finalize", "sanity"}
        receipt_count = 0
        for species_name in sorted(counts):
            receipt_path = data_root / species_name / "work" / ".metainformant_downstream_provenance.json"
            try:
                payload = json.loads(receipt_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            if (
                isinstance(payload, dict)
                and payload.get("schema") == "metainformant.rna.downstream.v2"
                and required_steps <= {str(step).strip().lower() for step in payload.get("steps", [])}
                and all(
                    isinstance(payload.get(key), str) and payload[key]
                    for key in ("amalgkit_version", "amalgkit_release_tag", "amalgkit_source_revision")
                )
            ):
                receipt_count += 1
        if counts and receipt_count == len(counts):
            result["descriptive_analysis"] = "receipt_present"
        elif receipt_count:
            result["descriptive_analysis"] = "partial_or_stale"
    except sqlite3.Error:
        result["cohort"] = "unreadable_database"
    return result


def build_status(
    *,
    data_root: Path,
    log_file: Path | None = None,
    inspect_processes: bool = True,
) -> dict[str, object]:
    """Build a machine-readable operational snapshot for future MCP adapters.

    Process and log observations are retained as diagnostics only. Readiness
    is derived from the campaign database and provenance-qualified downstream
    receipts, and biological inference is always withheld by this helper.
    """

    resolved_root = data_root.expanduser().resolve()
    resolved_log = (log_file or (resolved_root / "results" / "full_campaign.log")).expanduser().resolve()
    process = get_process_status() if inspect_processes else {"scanned": False, "running": None, "pid": None}
    readiness = _readiness_snapshot(resolved_root)
    writer_lock = any(
        (resolved_root / "results" / name).exists() for name in (".full_campaign.lock", ".finalization.lock")
    )
    log_info = parse_log_progress(resolved_log)
    return {
        "status": "running" if writer_lock or process.get("running") else "stopped",
        "process": process,
        "progress": {
            "processed": log_info.get("processed", 0),
            "total": log_info.get("total", 0),
            "percent": log_info.get("percent", 0),
        },
        "last_log": log_info.get("last_line", ""),
        "readiness": readiness,
        "evidence": {"data_root": str(resolved_root), "log_file": str(resolved_log), "writer_lock": writer_lock},
        "system": get_system_stats(resolved_root),
    }


def main():
    """Print an operational snapshot suitable for a terminal or adapter."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path(os.environ.get("AMALGKIT_DATA_ROOT", "output/amalgkit")))
    parser.add_argument("--log-file", type=Path)
    parser.add_argument(
        "--no-process-scan", action="store_true", help="Avoid process inspection; use locks and receipts only"
    )
    args = parser.parse_args()
    print(
        json.dumps(
            build_status(data_root=args.data_root, log_file=args.log_file, inspect_processes=not args.no_process_scan),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
