#!/usr/bin/env python3
"""Pipeline status checker for amalgkit RNA-seq orchestration.

Shows per-species quantification progress, downstream step completion,
and overall pipeline health at a glance.

Usage:
    # Quick status overview
    python3 scripts/rna/check_pipeline_status.py

    # Verbose mode with sample-level detail
    python3 scripts/rna/check_pipeline_status.py -v

    # Check a single species
    python3 scripts/rna/check_pipeline_status.py --species solenopsis_invicta
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# Bootstrap src path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

# ---------- constants ----------
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_ROOT = Path("blue/amalgkit")
LOG_DIR = PROJECT_ROOT / "output" / "amalgkit"

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

# ---------- helpers ----------

def count_quantified(species: str) -> int:
    """Count directories with an abundance file (quant completed)."""
    quant_dir = DATA_ROOT / species / "work" / "quant"
    if not quant_dir.exists():
        return 0
    return sum(
        1 for d in quant_dir.iterdir()
        if d.is_dir() and any(d.glob("*_abundance.tsv"))
    )


def count_total_samples(species: str) -> int | None:
    """Read metadata to count total samples (if metadata exists)."""
    meta_path = DATA_ROOT / species / "work" / "metadata" / "metadata.tsv"
    if not meta_path.exists():
        return None
    try:
        import pandas as pd
        df = pd.read_csv(meta_path, sep="\t", low_memory=False)
        return len(df)
    except Exception:
        # Fallback: line count minus header
        try:
            return sum(1 for _ in open(meta_path)) - 1
        except Exception:
            return None


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


def get_latest_log_line(species: str) -> str | None:
    """Get the last meaningful line from the workflow log."""
    log_path = LOG_DIR / f"{species}_workflow.log"
    if not log_path.exists():
        return None
    try:
        lines = log_path.read_text().strip().split("\n")
        for line in reversed(lines):
            if line.strip():
                return line.strip()[-100:]  # Last 100 chars
    except Exception:
        pass
    return None


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
    parser.add_argument("-v", "--verbose", action="store_true", help="Show sample-level detail")
    parser.add_argument("--species", type=str, help="Check a single species only")
    args = parser.parse_args()

    # Change to project root so relative paths resolve
    os.chdir(PROJECT_ROOT)

    species_list = [args.species] if args.species else SPECIES_ORDER

    running = check_process_running()
    print("=" * 72)
    print(f"  Amalgkit Pipeline Status  │  Process: {'🟢 RUNNING' if running else '🔴 STOPPED'}")
    print("=" * 72)
    print(f"{'Species':<30} {'Quant':>8} {'Total':>8} {'%':>6}  {'Downstream':<16}")
    print("-" * 72)

    total_quant = 0
    total_samples = 0
    for sp in species_list:
        q = count_quantified(sp)
        t = count_total_samples(sp)
        ds = check_downstream(sp)

        total_quant += q
        total_samples += t or 0

        pct = f"{100*q/t:.0f}%" if t and t > 0 else "?"
        t_str = str(t) if t is not None else "?"

        print(f"  {sp:<28} {q:>8} {t_str:>8} {pct:>6}  {ds}")

        if args.verbose:
            last_log = get_latest_log_line(sp)
            if last_log:
                print(f"    └─ {last_log}")

    print("-" * 72)
    overall_pct = f"{100*total_quant/total_samples:.1f}%" if total_samples > 0 else "?"
    print(f"  {'TOTAL':<28} {total_quant:>8} {total_samples:>8} {overall_pct:>6}")
    print("=" * 72)

    # Latest orchestrator log info
    logs = sorted(LOG_DIR.glob("streaming_orchestrator_*.log"))
    if logs:
        latest = logs[-1]
        mtime = os.path.getmtime(latest)
        from datetime import datetime
        age_mins = (datetime.now().timestamp() - mtime) / 60
        print(f"\n  Latest log: {latest.name}")
        print(f"  Last modified: {datetime.fromtimestamp(mtime):%Y-%m-%d %H:%M:%S} ({age_mins:.0f} min ago)")

    if not running:
        print("\n  ⚠  No active pipeline processes detected.")
        print("  To restart:  cd /home/trim/Documents/Git/MetaInformAnt && .venv/bin/python scripts/rna/run_all_species_parallel.py")
    print()


if __name__ == "__main__":
    main()
