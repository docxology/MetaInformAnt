#!/usr/bin/env python3
"""Sync amalgkit quant/curate/sanity/metadata results to a git-trackable location.

The pipeline stores results under blue/ which is a symlink to /Volumes/blue/data.
Git refuses to track files beyond symlinks, so we rsync the small, valuable
result files (quant abundances, curate tables, metadata) to output/amalgkit_results/
which is a real directory git can track.

Usage:
    uv run python scripts/rna/sync_quant_results.py [--dry-run]
"""

import argparse
import subprocess
import sys
from pathlib import Path

# Source: real data location (resolves the blue/ symlink)
SOURCE_BASE = Path("/Volumes/blue/data/amalgkit")

# Destination: git-trackable location inside the repo
DEST_BASE = Path("output/amalgkit_results")

# What to sync (small, valuable outputs only)
SYNC_DIRS = ["work/quant", "work/curate", "work/sanity", "work/metadata"]


def sync_species(species_dir: Path, dest_dir: Path, dry_run: bool = False) -> int:
    """Sync result directories for one species.

    Args:
        species_dir: Path to species data (e.g., /Volumes/blue/data/amalgkit/acromyrmex_echinatior)
        dest_dir: Destination path (e.g., output/amalgkit_results/acromyrmex_echinatior)
        dry_run: If True, only show what would be copied

    Returns:
        Number of files synced
    """
    total = 0
    for subdir in SYNC_DIRS:
        src = species_dir / subdir
        if not src.exists():
            continue

        dst = dest_dir / subdir
        dst.mkdir(parents=True, exist_ok=True)

        cmd = [
            "rsync", "-av", "--checksum",
            "--exclude", "*.fastq*",
            "--exclude", "*.sra",
            "--exclude", "*.fq*",
            f"{src}/",
            f"{dst}/",
        ]

        if dry_run:
            cmd.insert(2, "--dry-run")

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  ⚠ rsync failed for {src}: {result.stderr}", file=sys.stderr)
            continue

        # Count transferred files (lines that don't start with special chars)
        transferred = [
            line for line in result.stdout.strip().split("\n")
            if line and not line.startswith(("sending", "sent", "total", "created", "./"))
        ]
        total += len(transferred)

    return total


def main() -> int:
    parser = argparse.ArgumentParser(description="Sync amalgkit results to git-trackable location")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be synced")
    parser.add_argument("--species", type=str, help="Sync only this species")
    args = parser.parse_args()

    if not SOURCE_BASE.exists():
        print(f"Error: Source {SOURCE_BASE} not found", file=sys.stderr)
        return 1

    DEST_BASE.mkdir(parents=True, exist_ok=True)

    species_dirs = sorted(
        p for p in SOURCE_BASE.iterdir()
        if p.is_dir() and not p.name.startswith(".")
    )

    if args.species:
        species_dirs = [d for d in species_dirs if d.name == args.species]
        if not species_dirs:
            print(f"Error: Species '{args.species}' not found in {SOURCE_BASE}", file=sys.stderr)
            return 1

    total_files = 0
    mode = "DRY RUN" if args.dry_run else "SYNC"
    print(f"[{mode}] Syncing {len(species_dirs)} species from {SOURCE_BASE} → {DEST_BASE}")

    for species_dir in species_dirs:
        dest = DEST_BASE / species_dir.name
        count = sync_species(species_dir, dest, dry_run=args.dry_run)
        total_files += count
        if count > 0:
            print(f"  ✓ {species_dir.name}: {count} files")
        else:
            print(f"  · {species_dir.name}: up to date")

    print(f"\n{'Would sync' if args.dry_run else 'Synced'}: {total_files} files total")
    return 0


if __name__ == "__main__":
    sys.exit(main())
