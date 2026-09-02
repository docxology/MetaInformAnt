#!/usr/bin/env python3
"""Transcriptome SNP Calling Pipeline.

Thin orchestrator: re-downloads FASTQs for completed Amalgkit samples,
aligns to a reference genome with HISAT2, calls SNPs with bcftools, and
outputs per-biosample VCF files plus population genetics summaries.
All business logic lives in ``metainformant.eqtl``.

Usage:
    uv run python scripts/eqtl/rna_snp_pipeline.py --species amellifera --n-samples 3
    uv run python scripts/eqtl/rna_snp_pipeline.py --species amellifera --samples SRR21601882,SRR21601883
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

# ruff: noqa: E402 - this script prepends the local src tree before importing metainformant.
from metainformant.eqtl.pipeline import resolve_run_parameters, run_pipeline

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


def load_config(config_path: str) -> dict:
    """Load pipeline configuration from YAML file."""
    import yaml

    with open(config_path) as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(
        description="Transcriptome SNP Calling Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --config config/eqtl/eqtl_amellifera.yaml
  %(prog)s --species amellifera --n-samples 3
  %(prog)s --species amellifera --samples SRR21601882,SRR21601883
        """,
    )
    parser.add_argument("--config", help="Path to YAML config file")
    parser.add_argument("--species", help="Species name (e.g. amellifera)")
    parser.add_argument("--n-samples", type=int, default=3, help="Number of samples to process")
    parser.add_argument("--samples", help="Comma-separated list of SRR IDs to process")
    parser.add_argument("--threads", type=int, default=4, help="Threads for alignment")
    parser.add_argument("--no-cleanup", action="store_true", help="Keep FASTQ files after alignment")
    args = parser.parse_args()

    config = load_config(args.config) if args.config else None
    sample_ids_cli = args.samples.split(",") if args.samples else None

    try:
        params = resolve_run_parameters(
            species=args.species,
            sample_ids=sample_ids_cli,
            n_samples=args.n_samples,
            threads=args.threads,
            cleanup_fastq=not args.no_cleanup,
            config=config,
        )
    except ValueError as exc:
        parser.error(str(exc))
        return

    run_pipeline(
        species=params["species"],
        sample_ids=params["sample_ids"],
        n_samples=params["n_samples"],
        threads=params["threads"],
        cleanup_fastq=params["cleanup_fastq"],
    )


if __name__ == "__main__":
    main()
