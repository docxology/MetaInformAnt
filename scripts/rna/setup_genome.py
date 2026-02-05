#!/usr/bin/env python3
"""Genome setup orchestrator for RNA-seq workflows.

This is a thin wrapper that calls methods from metainformant.rna.genome_prep
to orchestrate complete genome setup pipelines.

Usage:
    # Full genome setup for a species
    python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

    # Verify status only
    python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --verify-only

    # Skip specific steps
    python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --skip-download --skip-prepare
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Import setup utilities (must be before other imports)
sys.path.insert(0, str(Path(__file__).parent))
from _setup_utils import ensure_venv_activated, check_environment_or_exit

# Auto-setup and activate venv
ensure_venv_activated(auto_setup=True)
check_environment_or_exit(auto_setup=True)

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna.genome_prep import orchestrate_genome_setup, verify_genome_status
from metainformant.core.utils.logging import get_logger

logger = get_logger("setup_genome")


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Orchestrate genome setup for a species",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Path to species workflow config file",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        help="Repository root directory (default: auto-detect)",
    )
    parser.add_argument(
        "--kmer-size",
        type=int,
        default=31,
        help="k-mer size for kallisto index (default: 31)",
    )
    parser.add_argument(
        "--skip-verify-initial",
        action="store_true",
        help="Skip initial verification step",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Skip genome download step",
    )
    parser.add_argument(
        "--skip-prepare",
        action="store_true",
        help="Skip transcriptome preparation step",
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Skip kallisto index building step",
    )
    parser.add_argument(
        "--skip-verify-final",
        action="store_true",
        help="Skip final verification step",
    )
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Only verify status, don't perform setup",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Dry run mode (don't actually perform operations)",
    )

    args = parser.parse_args()

    config_path = args.config.resolve()
    if not config_path.exists():
        logger.error(f"Config file not found: {config_path}")
        return 1

    # Determine repo root
    if args.repo_root:
        repo_root = args.repo_root.resolve()
    else:
        repo_root = Path(__file__).parent.parent.parent.resolve()

    # Verify only mode
    if args.verify_only:
        logger.info("Verifying genome status...")
        verification = verify_genome_status(config_path, repo_root=repo_root)
        logger.info(f"Genome downloaded: {verification.get('genome_downloaded', False)}")
        logger.info(f"RNA FASTA found: {verification.get('rna_fasta_found', False)}")
        logger.info(f"Kallisto index found: {verification.get('kallisto_index_found', False)}")
        if verification.get("error"):
            logger.warning(f"Error: {verification['error']}")
        return 0

    # Run orchestration
    logger.info(f"Orchestrating genome setup for {config_path.name}")
    results = orchestrate_genome_setup(
        config_path,
        repo_root=repo_root,
        skip_verify_initial=args.skip_verify_initial,
        skip_download=args.skip_download,
        skip_prepare=args.skip_prepare,
        skip_build=args.skip_build,
        skip_verify_final=args.skip_verify_final,
        kmer_size=args.kmer_size,
        dry_run=args.dry_run,
    )

    if results["success"]:
        logger.info(f"✅ Genome setup completed: {len(results.get('steps_completed', []))} steps")
        return 0
    else:
        logger.error(f"❌ Genome setup failed: {len(results.get('steps_failed', []))} steps failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())

