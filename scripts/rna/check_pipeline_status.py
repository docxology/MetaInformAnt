#!/usr/bin/env python3
"""Pipeline status checker for amalgkit RNA-seq orchestration.

Shows per-species quantification progress, downstream step completion,
and overall pipeline health at a glance.

Uses the SQLite progress database for instant queries instead of
scanning potentially huge directories.

Usage:
    # Quick status overview
    python3 scripts/rna/check_pipeline_status.py

    # Verbose mode with sample-level detail
    python3 scripts/rna/check_pipeline_status.py -v

    # Check a single species
    python3 scripts/rna/check_pipeline_status.py --species solenopsis_invicta

    # Show failed samples
    python3 scripts/rna/check_pipeline_status.py --failed

    # Generate visual dashboard (PDF + PNG)
    python3 scripts/rna/check_pipeline_status.py --dashboard
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# Bootstrap src path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from metainformant.rna.engine.progress_db import ProgressDB
from metainformant.rna.engine.species import configured_data_root, discover_species_names

# ---------- constants ----------
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_ROOT = configured_data_root()
DB_PATH = DATA_ROOT / "pipeline_progress.db"

ALL_STATES = ["pending", "downloading", "downloaded", "quantifying", "quantified", "failed"]


# ---------- helpers ----------


def check_downstream(species: str, data_root: Path = DATA_ROOT) -> str:
    """Check current merge, filtering, finalization, and sanity evidence."""
    species_dir = data_root / species
    work_dir = species_dir / "work"

    def has_table(*candidates: Path) -> bool:
        return any(path.exists() and any(path.rglob("*.tsv")) for path in candidates)

    has_merge = has_table(work_dir / "merge", species_dir / "merged", species_dir / "merged" / "merge")
    has_wsfilter = has_table(work_dir / "wsfilter", species_dir / "wsfilter")
    has_finalize = has_table(work_dir / "finalize", species_dir / "finalize")
    has_sanity = any(
        marker.exists()
        for marker in (
            work_dir / "sanity" / "sanity_check.txt",
            work_dir / "sanity_check.txt",
            species_dir / "sanity" / "sanity_check.txt",
        )
    )

    if has_merge and has_finalize and has_sanity:
        return "✅ Complete"
    if has_finalize:
        return "⚠️  Finalized; sanity pending"
    if has_wsfilter:
        return "⚠️  Filtered; finalization pending"
    if has_merge:
        return "⚠️  Merge only"
    return "❌ Not run"


def check_process_running() -> bool:
    """Quick check if any pipeline processes are active."""
    try:
        import subprocess

        result = subprocess.run(
            ["pgrep", "-f", "run_all_species|streaming_orchestrator"], capture_output=True, text=True
        )
        return result.returncode == 0
    except Exception:
        return False


# ---------- main ----------


def main():
    parser = argparse.ArgumentParser(description="Amalgkit pipeline status checker")
    parser.add_argument("-v", "--verbose", action="store_true", help="Show per-state breakdown")
    parser.add_argument("--species", type=str, help="Check a single species only")
    parser.add_argument("--failed", action="store_true", help="Show failed samples")
    parser.add_argument("--dashboard", action="store_true", help="Generate visual dashboard PDF+PNG to output/")
    parser.add_argument("--data-root", type=Path, help="Data root containing pipeline_progress.db")
    args = parser.parse_args()

    # Change to project root so relative paths resolve
    os.chdir(PROJECT_ROOT)

    data_root = args.data_root.expanduser().resolve() if args.data_root else configured_data_root()
    db_path = data_root / "pipeline_progress.db"

    if not db_path.exists():
        print("  ⚠  No progress database found. Run the pipeline to initialize it.")
        print(f"  Expected at: {db_path}")
        return

    db = ProgressDB(db_path)
    counts = db.get_counts()

    config_species = discover_species_names(PROJECT_ROOT / "config" / "amalgkit")
    species_list = [args.species] if args.species else (config_species or sorted(counts))

    running = check_process_running()
    print("=" * 80)
    print(f"  Amalgkit Pipeline Status  │  Process: {'🟢 RUNNING' if running else '🔴 STOPPED'}")
    print("=" * 80)

    if args.verbose:
        header = f"{'Species':<28}"
        for s in ALL_STATES:
            header += f" {s[:7]:>8}"
        header += f" {'Total':>7}  {'Downstream':<16}"
        print(header)
    else:
        print(f"{'Species':<28} {'Quant':>8} {'Total':>8} {'%':>6}  {'Downstream':<16}")
    print("-" * 80)

    grand_quant = 0
    grand_total = 0

    for sp in species_list:
        sp_counts = counts.get(sp, {})
        total = sum(sp_counts.values())
        quant = sp_counts.get("quantified", 0)
        ds = check_downstream(sp, data_root)

        grand_quant += quant
        grand_total += total

        if args.verbose:
            row = f"  {sp:<26}"
            for s in ALL_STATES:
                row += f" {sp_counts.get(s, 0):>8}"
            row += f" {total:>7}  {ds}"
            print(row)
        else:
            pct = f"{100*quant/total:.0f}%" if total > 0 else "-"
            t_str = str(total) if total > 0 else "-"
            print(f"  {sp:<26} {quant:>8} {t_str:>8} {pct:>6}  {ds}")

    print("-" * 80)
    overall_pct = f"{100*grand_quant/grand_total:.1f}%" if grand_total > 0 else "?"
    if args.verbose:
        row = f"  {'TOTAL':<26}"
        totals_by_state = {}
        for sp in species_list:
            for s in ALL_STATES:
                totals_by_state[s] = totals_by_state.get(s, 0) + counts.get(sp, {}).get(s, 0)
        for s in ALL_STATES:
            row += f" {totals_by_state.get(s, 0):>8}"
        row += f" {grand_total:>7}"
        print(row)
    else:
        print(f"  {'TOTAL':<26} {grand_quant:>8} {grand_total:>8} {overall_pct:>6}")
    print("=" * 80)

    # Show failed samples if requested
    if args.failed:
        failed = db.get_failed()
        if failed:
            print(f"\n  Failed Samples ({len(failed)}):")
            print(f"  {'Species':<26} {'SRR ID':<16} {'Error':<30} {'When'}")
            print("  " + "-" * 76)
            for f in failed[:50]:
                print(f"  {f['species']:<26} {f['srr_id']:<16} {(f['error'] or '?'):<30} {f['updated_at']}")
            if len(failed) > 50:
                print(f"  ... and {len(failed) - 50} more")
        else:
            print("\n  No failed samples. 🎉")

    if not running:
        print("\n  ⚠  No active pipeline processes detected.")
        print("  To restart:  uv run python scripts/rna/run_all_species.py")
    print()

    db.close()

    # Generate dashboard if requested
    if args.dashboard:
        from metainformant.rna.engine.progress_dashboard import generate_dashboard

        generate_dashboard(db_path=db_path, output_dir=data_root)


if __name__ == "__main__":
    main()
