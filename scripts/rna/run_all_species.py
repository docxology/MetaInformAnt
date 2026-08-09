#!/usr/bin/env python3
"""ENA-first sample-by-sample RNA-seq pipeline orchestrator.

This script is a thin wrapper around the StreamingPipelineOrchestrator.
It manages current metadata, acquisition, and quantification for multiple
species with robust error handling, retries, and resume capability. Downstream
merge, wsfilter, finalize, and sanity are owned by the external-volume
checkpoint runner after this producer stops.

Usage:
    python3 scripts/rna/run_all_species.py [--max-gb 50.0] [--workers 4] [--threads 8]
        [--discovery-workers N]
        [--quant-slots N] [--fastq-threads N] [--fastq-slots N]
        [--compression-threads N] [--validation-slots N] [--max-in-flight N]
        [--requantification-policy preserve|version-drift|all]
    python3 scripts/rna/run_all_species.py --species solenopsis_invicta --dry-run
"""

import argparse
import os
import sys
from pathlib import Path

# Add src to python path to allow importing metainformant modules
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from metainformant.rna.engine.species import configured_data_root, discover_species_config_names
from metainformant.rna.engine.streaming_orchestrator import StreamingPipelineOrchestrator

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = REPO_ROOT / "projects" / "hymenoptera_amalgkit" / "config" / "amalgkit"
if not CONFIG_DIR.is_dir():
    CONFIG_DIR = REPO_ROOT / "config" / "amalgkit"

DEFAULTS = {
    "max_gb": float(os.environ.get("PIPELINE_MAX_GB", 50.0)),
    "workers": int(os.environ.get("PIPELINE_WORKERS", 4)),
    "threads": int(os.environ.get("PIPELINE_THREADS", 8)),
}


def main() -> int:
    parser = argparse.ArgumentParser(description="ENA-first sample-by-sample RNA-seq pipeline")
    parser.add_argument(
        "--config-dir",
        type=Path,
        default=CONFIG_DIR,
        help=f"Directory containing amalgkit_*.yaml species configs (default: {CONFIG_DIR})",
    )
    parser.add_argument(
        "--species",
        action="append",
        help="Run only this species identifier; repeat for multiple species",
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        help="Data root for work files and the progress database (also sets AMALGKIT_DATA_ROOT)",
    )
    parser.add_argument("--dry-run", action="store_true", help="List the resolved species configs without running")
    parser.add_argument(
        "--max-gb",
        type=float,
        default=DEFAULTS["max_gb"],
        help=f"Max sample size in GB (default: {DEFAULTS['max_gb']})",
    )
    parser.add_argument(
        "--workers", type=int, default=DEFAULTS["workers"], help=f"Parallel workers (default: {DEFAULTS['workers']})"
    )
    parser.add_argument(
        "--threads", type=int, default=DEFAULTS["threads"], help=f"Total threads (default: {DEFAULTS['threads']})"
    )
    parser.add_argument(
        "--discovery-workers",
        type=int,
        default=None,
        help="Concurrent species metadata/reference discovery workers (default: 4; set 1 for serial diagnostics)",
    )
    parser.add_argument(
        "--quant-slots",
        type=int,
        default=None,
        help="Maximum concurrent Kallisto quantifications (default: bounded by workers and total threads)",
    )
    parser.add_argument(
        "--fastq-threads",
        type=int,
        default=None,
        help="Threads for fasterq-dump fallback (default: bounded from the quant budget)",
    )
    parser.add_argument(
        "--fastq-slots",
        type=int,
        default=None,
        help="Maximum concurrent fasterq-dump fallbacks (default: up to four; limits external-disk contention)",
    )
    parser.add_argument(
        "--compression-threads",
        type=int,
        default=None,
        help="Threads per pigz compression process (default: bounded from the quant budget)",
    )
    parser.add_argument(
        "--validation-slots",
        type=int,
        default=None,
        help="Concurrent full local FASTQ validations (default: up to four; limits volume contention)",
    )
    parser.add_argument(
        "--max-in-flight",
        type=int,
        default=None,
        help="Maximum submitted sample tasks (default: two times active workers)",
    )
    parser.add_argument(
        "--requantification-policy",
        choices=("preserve", "version-drift", "all"),
        default=os.environ.get("AMALGKIT_REQUANTIFICATION_POLICY", "preserve"),
        help="Handle contract-compatible version drift (default: preserve)",
    )
    args = parser.parse_args()

    config_dir = args.config_dir.expanduser().resolve()
    if args.data_root:
        os.environ["AMALGKIT_DATA_ROOT"] = str(args.data_root.expanduser().resolve())
    data_root = configured_data_root()

    config_names = discover_species_config_names(config_dir)
    if args.species:
        wanted = set(args.species)
        config_names = [name for name in config_names if name.removeprefix("amalgkit_").removesuffix(".yaml") in wanted]
    if not config_names:
        parser.error(f"no runnable species configurations found under {config_dir}")

    print("╔══════════════════════════════════════════════════════════╗")
    print("║  ENA-First Sample-by-Sample Pipeline (Current)          ║")
    print(
        f"║  Species: {len(config_names)} | Max: {args.max_gb}GB | Workers: {args.workers} | Threads: {args.threads}  ║"
    )
    print("╚══════════════════════════════════════════════════════════╝")

    print(f"Config directory: {config_dir}")
    print(f"Data root: {data_root}")
    print("Species configs:")
    for config_name in config_names:
        print(f"  - {config_name}")

    if args.dry_run:
        return 0

    orchestrator = StreamingPipelineOrchestrator(
        config_dir=config_dir,
        log_dir=data_root / "logs",
        db_path=data_root / "pipeline_progress.db",
    )
    orchestrator.run_all(
        config_names,
        args.max_gb,
        args.workers,
        args.threads,
        quant_slots=args.quant_slots,
        fasterq_threads=args.fastq_threads,
        fasterq_slots=args.fastq_slots,
        compression_threads=args.compression_threads,
        validation_slots=args.validation_slots,
        max_in_flight=args.max_in_flight,
        requantification_policy=args.requantification_policy,
        discovery_workers=args.discovery_workers,
    )
    return 0


if __name__ == "__main__":
    main()
