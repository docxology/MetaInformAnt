#!/usr/bin/env python3
"""ENA-first sample-by-sample RNA-seq pipeline orchestrator.

This script is a thin wrapper around the StreamingPipelineOrchestrator.
It manages the sequential execution of amalgkit pipeline steps for multiple species,
with robust error handling, retries, and resume capability.

Usage:
    python3 scripts/rna/run_all_species.py [--max-gb 5.0] [--workers 4] [--threads 12]
"""

import argparse
import sys
import os
from pathlib import Path

# Add src to python path to allow importing metainformant modules
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from metainformant.rna.engine.streaming_orchestrator import StreamingPipelineOrchestrator

# Species processing order: precisely the final 6 empty datasets
SPECIES_ORDER = [
    "amalgkit_athalia_rosae.yaml",
    "amalgkit_bombus_terrestris.yaml",
    "amalgkit_megachile_rotundata.yaml",
    "amalgkit_nasonia_vitripennis.yaml",
    "amalgkit_polistes_canadensis.yaml",
    "amalgkit_polistes_fuscatus.yaml",
]

DEFAULTS = {
    "max_gb": float(os.environ.get("PIPELINE_MAX_GB", 50.0)),
    "workers": int(os.environ.get("PIPELINE_WORKERS", 12)),
    "threads": int(os.environ.get("PIPELINE_THREADS", 6)),
}

def main():
    parser = argparse.ArgumentParser(
        description="ENA-first sample-by-sample RNA-seq pipeline"
    )
    parser.add_argument("--max-gb", type=float, default=DEFAULTS["max_gb"],
                        help=f"Max sample size in GB (default: {DEFAULTS['max_gb']})")
    parser.add_argument("--workers", type=int, default=DEFAULTS["workers"],
                        help=f"Parallel workers (default: {DEFAULTS['workers']})")
    parser.add_argument("--threads", type=int, default=DEFAULTS["threads"],
                        help=f"Total threads (default: {DEFAULTS['threads']})")
    args = parser.parse_args()

    print(f"╔══════════════════════════════════════════════════════════╗")
    print(f"║  ENA-First Sample-by-Sample Pipeline (Wrapped)          ║")
    print(f"║  Species: {len(SPECIES_ORDER)} | Max: {args.max_gb}GB | Workers: {args.workers} | Threads: {args.threads}  ║")
    print(f"╚══════════════════════════════════════════════════════════╝")

    orchestrator = StreamingPipelineOrchestrator()
    orchestrator.run_all(SPECIES_ORDER, args.max_gb, args.workers, args.threads)

if __name__ == "__main__":
    main()
