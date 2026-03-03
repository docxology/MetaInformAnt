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

# ---------- constants ----------
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_ROOT = Path("output/amalgkit")
DB_PATH = DATA_ROOT / "pipeline_progress.db"

SPECIES_ORDER = [
    "anoplolepis_gracilipes",
    "acromyrmex_echinatior",
    "dinoponera_quadriceps",
    "vollenhovia_emeryi",
    "odontomachus_brunneus",
    "formica_exsecta",
    "temnothorax_americanus",
    "wasmannia_auropunctata",
    "nylanderia_fulva",
    "temnothorax_curvispinosus",
    "pbarbatus",
    "cardiocondyla_obscurior",
    "temnothorax_nylanderi",
    "linepithema_humile",
    "atta_cephalotes",
    "ooceraea_biroi",
    "camponotus_floridanus",
    "solenopsis_invicta",
    "monomorium_pharaonis",
    "temnothorax_longispinosus",
    "harpegnathos_saltator",
    "amellifera",
]

ALL_STATES = ["pending", "downloading", "downloaded", "quantifying", "quantified", "failed"]


# ---------- helpers ----------

def check_downstream(species: str) -> str:
    """Check if merge/curate/sanity outputs exist."""
    work_dir = DATA_ROOT / species / "work"
    curate_dir = work_dir / "curate"
    merge_dir = DATA_ROOT / species / "merged"
    if not merge_dir.exists():
        merge_dir = work_dir / "merge"
    
    has_merge = merge_dir.exists() and any(merge_dir.rglob("*.tsv"))
    has_curate = curate_dir.exists() and any(curate_dir.rglob("*.tsv"))
    
    if has_merge and has_curate:
        return "✅ Complete"
    elif has_merge:
        return "⚠️  Merge only"
    else:
        return "❌ Not run"


def check_process_running() -> bool:
    """Quick check if any pipeline processes are active."""
    try:
        import subprocess
        result = subprocess.run(
            ["pgrep", "-f", "run_all_species|streaming_orchestrator|run_workflow"],
            capture_output=True, text=True
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
    parser.add_argument("--dashboard", action="store_true",
                        help="Generate visual dashboard PDF+PNG to output/")
    args = parser.parse_args()

    # Change to project root so relative paths resolve
    os.chdir(PROJECT_ROOT)

    if not DB_PATH.exists():
        print("  ⚠  No progress database found. Run the pipeline to initialize it.")
        print(f"  Expected at: {DB_PATH}")
        return

    db = ProgressDB(DB_PATH)
    counts = db.get_counts()

    species_list = [args.species] if args.species else SPECIES_ORDER

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
        ds = check_downstream(sp)

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
        print("  To restart:  cd /home/trim/Documents/Git/MetaInformAnt && "
              ".venv/bin/python scripts/rna/run_all_species.py")
    print()

    db.close()

    # Generate dashboard if requested
    if args.dashboard:
        from metainformant.rna.engine.progress_dashboard import generate_dashboard
        generate_dashboard(db_path=DB_PATH, output_dir=Path("output/amalgkit"))


if __name__ == "__main__":
    main()
