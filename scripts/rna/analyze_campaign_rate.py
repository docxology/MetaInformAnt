#!/usr/bin/env python3
"""Campaign rate analyzer for streaming_orchestrator logs.

Parses one or more logs/streaming_orchestrator_*.log files and reports:
- total tasks done (Done / Failed from "[N/M] <species> SRR...: <STATUS>" lines)
- per-hour task-completion counts
- cumulative finished-task rate over the log window
- download throughput from "Downloaded <SRR>.fastq.gz: <X> GB" lines (GB/h and MB/s)
- species mix over a trailing window (default 2h)
- optional ETA via --total N (with explicit caveats)

Read-only: never writes to the data root.

Usage:
    python scripts/rna/analyze_campaign_rate.py logs/streaming_orchestrator_20260831_181002.log
    python scripts/rna/analyze_campaign_rate.py --json <log> [<log> ...]
    python scripts/rna/analyze_campaign_rate.py --total 15620 --trailing-hours 2 <log>
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from collections.abc import Iterable
from datetime import datetime
from pathlib import Path

TASK_RE = re.compile(r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})[,\.]\d+ .*\[(\d+)/(\d+)\] " r"(\S+) (SRR\S+): (\S+)")
DOWNLOAD_RE = re.compile(r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})[,\.]\d+ .*Downloaded (\S+): ([0-9.]+) GB")


def parse_timestamp(raw: str) -> datetime:
    return datetime.strptime(raw, "%Y-%m-%d %H:%M:%S")


def analyze_log(path: Path, trailing_hours: float = 2.0) -> dict:
    """Parse one orchestrator log file and return summary statistics."""
    task_events: list[tuple[datetime, str, str, str]] = []
    download_events: list[tuple[datetime, float]] = []
    declared_total: int | None = None

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = TASK_RE.match(line)
            if m:
                ts = parse_timestamp(m.group(1))
                species, srr, status = m.group(4), m.group(5), m.group(6)
                task_events.append((ts, species, srr, status))
                declared_total = int(m.group(3))
                continue
            m = DOWNLOAD_RE.match(line)
            if m:
                ts = parse_timestamp(m.group(1))
                download_events.append((ts, float(m.group(3))))

    if not task_events and not download_events:
        return {"log": str(path), "error": "no parseable task or download events"}

    first_ts = min([ts for ts, *_ in task_events] + [ts for ts, _ in download_events])
    last_ts = max([ts for ts, *_ in task_events] + [ts for ts, _ in download_events])
    window_hours = max((last_ts - first_ts).total_seconds() / 3600.0, 1e-9)

    statuses = Counter(status for _, _, _, status in task_events)
    done = statuses.get("Done", 0)
    failed = statuses.get("Failed", 0)
    finished = done + failed

    # Per-hour task-completion counts (any terminal status).
    per_hour: Counter[str] = Counter()
    for ts, _, _, status in task_events:
        if status in ("Done", "Failed"):
            per_hour[ts.strftime("%Y-%m-%d %H:00")] += 1

    # Download throughput.
    total_gb = sum(gb for _, gb in download_events)
    gb_per_hour = total_gb / window_hours
    mb_per_s = total_gb * 1024 / (window_hours * 3600)

    # Species mix over trailing window.
    cutoff = last_ts.timestamp() - trailing_hours * 3600
    species_recent: Counter[str] = Counter()
    for ts, species, _, status in task_events:
        if status in ("Done", "Failed") and ts.timestamp() >= cutoff:
            species_recent[species] += 1

    return {
        "log": str(path),
        "window_start": first_ts.isoformat(),
        "window_end": last_ts.isoformat(),
        "window_hours": round(window_hours, 3),
        "tasks": {
            "done": done,
            "failed": failed,
            "finished_total": finished,
            "status_counts": dict(statuses),
            "declared_total_from_log": declared_total,
        },
        "per_hour_counts": dict(sorted(per_hour.items())),
        "cumulative_finished_rate_per_hour": round(finished / window_hours, 2),
        "downloads": {
            "files": len(download_events),
            "total_gb": round(total_gb, 2),
            "gb_per_hour": round(gb_per_hour, 2),
            "mb_per_s_equivalent": round(mb_per_s, 2),
        },
        "species_mix_trailing_hours": {
            "hours": trailing_hours,
            "counts": dict(species_recent.most_common()),
        },
    }


def eta_estimate(summary: dict, total: int | None) -> dict | None:
    """ETA from cumulative rate. Explicitly caveated: assumes constant rate,
    ignores overnight troughs, excludes still-queued re-quantification."""
    if total is None:
        declared = summary.get("tasks", {}).get("declared_total_from_log")
        if declared:
            total = declared
    if not total:
        return None
    finished = summary.get("tasks", {}).get("finished_total", 0)
    rate = summary.get("cumulative_finished_rate_per_hour", 0.0)
    if rate <= 0:
        return {"eta_hours": None, "caveat": "zero observed rate; ETA undefined"}
    remaining = max(total - finished, 0)
    return {
        "total": total,
        "finished": finished,
        "remaining": remaining,
        "rate_per_hour": rate,
        "eta_hours": round(remaining / rate, 2),
        "caveat": (
            "Linear extrapolation of cumulative average rate. Actual rate varies "
            "(overnight trough 35-65/h observed vs day peaks >150/h); ETA is a "
            "smoothed estimate, not a schedule."
        ),
    }


def aggregate(summaries: Iterable[dict]) -> dict:
    summaries = list(summaries)
    ok = [s for s in summaries if "error" not in s]
    if not ok:
        return {"logs": summaries, "error": "no parseable logs"}
    first = min(s["window_start"] for s in ok)
    last = max(s["window_end"] for s in ok)
    done = sum(s["tasks"]["done"] for s in ok)
    failed = sum(s["tasks"]["failed"] for s in ok)
    total_gb = sum(s["downloads"]["total_gb"] for s in ok)
    hours = sum(s["window_hours"] for s in ok)
    return {
        "logs": summaries,
        "aggregate": {
            "window_start": first,
            "window_end": last,
            "done": done,
            "failed": failed,
            "finished_total": done + failed,
            "downloads_total_gb": round(total_gb, 2),
            "window_hours_sum": round(hours, 3),
        },
    }


def _multi_log_eta(out: dict, total: int | None) -> dict | None:
    if total is None or "aggregate" not in out:
        return None
    agg = out["aggregate"]
    hours = agg["window_hours_sum"]
    if hours <= 0:
        return None
    rate = agg["finished_total"] / hours
    if rate <= 0:
        return None
    remaining = max(total - agg["finished_total"], 0)
    return {
        "total": total,
        "finished": agg["finished_total"],
        "remaining": remaining,
        "rate_per_hour": round(rate, 2),
        "eta_hours": round(remaining / rate, 2),
        "caveat": "Aggregate multi-log linear extrapolation; see single-log caveat.",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Campaign rate analyzer")
    parser.add_argument("logs", nargs="+", type=Path, help="orchestrator log file(s)")
    parser.add_argument("--json", action="store_true", help="JSON output mode")
    parser.add_argument("--trailing-hours", type=float, default=2.0, help="species-mix window (hours)")
    parser.add_argument("--total", type=int, default=None, help="explicit grand total for ETA")
    args = parser.parse_args(argv)

    for p in args.logs:
        if not p.is_file():
            print(f"error: log not found: {p}", file=sys.stderr)
            return 2

    summaries = [analyze_log(p, trailing_hours=args.trailing_hours) for p in args.logs]
    if len(args.logs) == 1:
        out = summaries[0]
        eta = eta_estimate(out, args.total) if "error" not in out else None
    else:
        out = aggregate(summaries)
        eta = _multi_log_eta(out, args.total)
    out["eta"] = eta

    if args.json:
        print(json.dumps(out, indent=2, sort_keys=True))
    else:
        _print_human(out)
    return 0


def _print_human(out: dict) -> None:
    def fmt(x: dict, indent: int = 0) -> None:
        pad = "  " * indent
        for key, val in x.items():
            if isinstance(val, dict):
                print(f"{pad}{key}:")
                fmt(val, indent + 1)
            else:
                print(f"{pad}{key}: {val}")

    fmt(out)


if __name__ == "__main__":
    sys.exit(main())
