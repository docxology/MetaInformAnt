#!/usr/bin/env python3
"""Monitor Apis mellifera processing progress.

Displays real-time statistics on:
- Number of quantified samples
- Processing rate (samples/hour)
- Estimated time to completion
- Disk space usage
- Recent activity

Usage:
    python scripts/rna/monitor_apis_mellifera.py [--watch]

With --watch flag, refreshes every 30 seconds.
"""

import argparse
import os
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

# Configuration
WORK_DIR = Path("output/amalgkit/apis_mellifera_all/work")
FASTQ_DIR = Path("output/amalgkit/apis_mellifera_all/fastq/getfastq")
QUANT_DIR = WORK_DIR / "quant"
METADATA_FILE = WORK_DIR / "metadata/metadata_selected.tsv"
TOTAL_SAMPLES = 7270  # Known total from metadata


def get_free_disk_gb() -> float:
    """Get free disk space in GB."""
    stat = os.statvfs(WORK_DIR)
    return (stat.f_frsize * stat.f_bavail) / (1024**3)


def get_disk_usage_gb(path: Path) -> float:
    """Get directory size in GB."""
    total = 0
    try:
        for entry in path.rglob("*"):
            if entry.is_file():
                total += entry.stat().st_size
    except Exception:
        pass
    return total / (1024**3)


def get_quantified_samples() -> list[dict]:
    """Get list of quantified samples with timestamps."""
    samples = []
    if not QUANT_DIR.exists():
        return samples

    for sample_dir in QUANT_DIR.iterdir():
        if sample_dir.is_dir():
            abundance_file = sample_dir / "abundance.tsv"
            if abundance_file.exists():
                mtime = abundance_file.stat().st_mtime
                samples.append({"id": sample_dir.name, "time": datetime.fromtimestamp(mtime)})

    return sorted(samples, key=lambda x: x["time"])


def get_processing_in_progress() -> list[str]:
    """Get samples currently being processed (have FASTQ but no quant)."""
    in_progress = []
    if not FASTQ_DIR.exists():
        return in_progress

    for sample_dir in FASTQ_DIR.iterdir():
        if sample_dir.is_dir() and sample_dir.name != "sra":
            quant_dir = QUANT_DIR / sample_dir.name
            abundance_file = quant_dir / "abundance.tsv"
            if not abundance_file.exists():
                in_progress.append(sample_dir.name)

    return in_progress


def calculate_rate(samples: list[dict], hours: float = 1.0) -> float:
    """Calculate processing rate over last N hours."""
    if not samples:
        return 0.0

    cutoff = datetime.now() - timedelta(hours=hours)
    recent = [s for s in samples if s["time"] > cutoff]

    if not recent:
        return 0.0

    return len(recent) / hours


def format_duration(seconds: float) -> str:
    """Format seconds as human-readable duration."""
    if seconds < 0:
        return "N/A"

    days = int(seconds // 86400)
    hours = int((seconds % 86400) // 3600)
    minutes = int((seconds % 3600) // 60)

    parts = []
    if days > 0:
        parts.append(f"{days}d")
    if hours > 0 or days > 0:
        parts.append(f"{hours}h")
    parts.append(f"{minutes}m")

    return " ".join(parts)


def print_progress():
    """Print current progress statistics."""
    print("=" * 60)
    print("🐝 Apis mellifera RNA-seq Processing Monitor")
    print("=" * 60)
    print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    # Get quantified samples
    samples = get_quantified_samples()
    quantified = len(samples)
    remaining = TOTAL_SAMPLES - quantified
    pct = (quantified / TOTAL_SAMPLES) * 100

    print(f"📊 Progress:")
    print(f"   Quantified: {quantified:,} / {TOTAL_SAMPLES:,} ({pct:.1f}%)")
    print(f"   Remaining:  {remaining:,}")

    # Progress bar
    bar_width = 40
    filled = int(bar_width * quantified / TOTAL_SAMPLES)
    bar = "█" * filled + "░" * (bar_width - filled)
    print(f"   [{bar}]")
    print()

    # Processing rate
    rate_1h = calculate_rate(samples, hours=1.0)
    rate_6h = calculate_rate(samples, hours=6.0)
    rate_24h = calculate_rate(samples, hours=24.0)

    print(f"⚡ Processing Rate:")
    print(f"   Last 1h:  {rate_1h:.1f} samples/hour")
    print(f"   Last 6h:  {rate_6h:.1f} samples/hour")
    print(f"   Last 24h: {rate_24h:.1f} samples/hour")
    print()

    # ETA
    if rate_24h > 0:
        eta_seconds = (remaining / rate_24h) * 3600
        eta = datetime.now() + timedelta(seconds=eta_seconds)
        print(f"⏱️ Estimated Time:")
        print(f"   Time remaining: {format_duration(eta_seconds)}")
        print(f"   ETA: {eta.strftime('%Y-%m-%d %H:%M')}")
    else:
        print(f"⏱️ Estimated Time: Calculating...")
    print()

    # Disk usage
    free_gb = get_free_disk_gb()
    fastq_gb = get_disk_usage_gb(FASTQ_DIR)
    quant_gb = get_disk_usage_gb(QUANT_DIR)

    print(f"💾 Disk Usage:")
    print(f"   Free:  {free_gb:.1f} GB")
    print(f"   FASTQ: {fastq_gb:.1f} GB (temporary)")
    print(f"   Quant: {quant_gb:.1f} GB (output)")
    print()

    # In-progress samples
    in_progress = get_processing_in_progress()
    print(f"🔄 In Progress: {len(in_progress)} samples")
    for sample in in_progress[:5]:
        print(f"   - {sample}")
    if len(in_progress) > 5:
        print(f"   ... and {len(in_progress) - 5} more")
    print()

    # Recent completions
    print(f"✅ Recent Completions:")
    for sample in samples[-5:]:
        elapsed = datetime.now() - sample["time"]
        elapsed_str = format_duration(elapsed.total_seconds())
        print(f"   - {sample['id']} ({elapsed_str} ago)")
    print()

    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(description="Monitor Apis mellifera processing")
    parser.add_argument("--watch", "-w", action="store_true", help="Refresh every 30 seconds")
    parser.add_argument("--interval", "-i", type=int, default=30, help="Refresh interval in seconds (default: 30)")
    args = parser.parse_args()

    if args.watch:
        try:
            while True:
                os.system("clear" if os.name == "posix" else "cls")
                print_progress()
                time.sleep(args.interval)
        except KeyboardInterrupt:
            print("\nMonitor stopped.")
    else:
        print_progress()


if __name__ == "__main__":
    main()
