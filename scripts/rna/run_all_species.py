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
from pathlib import Path

# Add src to python path to allow importing metainformant modules
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from metainformant.rna.engine.streaming_orchestrator import StreamingPipelineOrchestrator

# Species processing order: all ants (smallest first), then bees
SPECIES_ORDER = [
    # ── Ants (21 species, sorted by sample count ascending) ──
    "amalgkit_anoplolepis_gracilipes.yaml",      #    2 samples
    "amalgkit_acromyrmex_echinatior.yaml",        #    8 samples
    "amalgkit_dinoponera_quadriceps.yaml",        #   13 samples
    "amalgkit_vollenhovia_emeryi.yaml",           #   15 samples
    "amalgkit_odontomachus_brunneus.yaml",        #   19 samples
    "amalgkit_formica_exsecta.yaml",              #   23 samples
    "amalgkit_temnothorax_americanus.yaml",       #   32 samples
    "amalgkit_wasmannia_auropunctata.yaml",       #   33 samples
    "amalgkit_nylanderia_fulva.yaml",             #   40 samples
    "amalgkit_temnothorax_curvispinosus.yaml",    #   43 samples
    "amalgkit_pbarbatus.yaml",                    #   95 samples
    "amalgkit_cardiocondyla_obscurior.yaml",      #  162 samples
    "amalgkit_temnothorax_nylanderi.yaml",        #  166 samples
    "amalgkit_linepithema_humile.yaml",           #  173 samples
    "amalgkit_atta_cephalotes.yaml",              #  220 samples
    "amalgkit_ooceraea_biroi.yaml",               #  237 samples
    "amalgkit_camponotus_floridanus.yaml",        #  304 samples
    "amalgkit_solenopsis_invicta.yaml",           #  349 samples
    "amalgkit_monomorium_pharaonis.yaml",         #  370 samples
    "amalgkit_temnothorax_longispinosus.yaml",    #  508 samples
    "amalgkit_harpegnathos_saltator.yaml",        #  689 samples
    # ── Bees ──
    "amalgkit_amellifera.yaml",                   # 3154 samples
]

DEFAULTS = {
    "max_gb": 5.0,
    "workers": 16,
    "threads": 12,
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
