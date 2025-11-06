#!/usr/bin/env python3
"""Main workflow orchestrator for amalgkit RNA-seq analysis.

This is a thin wrapper that calls methods from metainformant.rna.orchestration
to run complete workflows for single species.

Usage:
    # Run full workflow for a species
    python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pbarbatus.yaml

    # Run specific steps
    python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pbarbatus.yaml --steps getfastq quant merge

    # Check status
    python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pbarbatus.yaml --status

    # Cleanup unquantified samples
    python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pbarbatus.yaml --cleanup-unquantified
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

from metainformant.rna.cleanup import cleanup_partial_downloads, fix_abundance_naming_for_species
from metainformant.rna.monitoring import analyze_species_status, check_workflow_progress
from metainformant.rna.orchestration import (
    cleanup_unquantified_samples,
    run_workflow_for_species,
)
from metainformant.core.logging import get_logger

logger = get_logger("run_workflow")


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Run amalgkit workflow for a single species",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Path to species workflow config file",
    )
    parser.add_argument(
        "--steps",
        nargs="+",
        help="Specific steps to run (default: all steps)",
    )
    parser.add_argument(
        "--status",
        action="store_true",
        help="Check workflow status instead of running",
    )
    parser.add_argument(
        "--detailed",
        action="store_true",
        help="Show detailed status (use with --status)",
    )
    parser.add_argument(
        "--cleanup-unquantified",
        action="store_true",
        help="Quantify downloaded samples and cleanup FASTQs",
    )
    parser.add_argument(
        "--cleanup-partial",
        action="store_true",
        help="Clean up partial downloads",
    )
    parser.add_argument(
        "--fix-abundance-naming",
        action="store_true",
        help="Fix abundance file naming for merge compatibility",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Stop on first failure",
    )

    args = parser.parse_args()

    config_path = args.config.resolve()
    if not config_path.exists():
        logger.error(f"Config file not found: {config_path}")
        return 1

    # Status check
    if args.status:
        if args.detailed:
            status = analyze_species_status(config_path)
        else:
            status = check_workflow_progress(config_path)
        logger.info(f"Status: {status}")
        return 0

    # Cleanup operations
    if args.cleanup_unquantified:
        logger.info("Cleaning up unquantified samples...")
        quantified, failed = cleanup_unquantified_samples(config_path)
        logger.info(f"Quantified: {quantified}, Failed: {failed}")
        return 0 if failed == 0 else 1

    if args.cleanup_partial:
        logger.info("Cleaning up partial downloads...")
        result = cleanup_partial_downloads(config_path, dry_run=False)
        logger.info(f"Deleted: {result['deleted']}, Freed: {result['freed_mb']}MB")
        return 0 if result["errors"] == 0 else 1

    if args.fix_abundance_naming:
        logger.info("Fixing abundance file naming...")
        created, exists = fix_abundance_naming_for_species(config_path)
        logger.info(f"Created: {created}, Already exists: {exists}")
        return 0

    # Run workflow
    logger.info(f"Running workflow for {config_path.name}")
    results = run_workflow_for_species(
        config_path,
        steps=args.steps,
        check=args.check,
    )

    if results["success"]:
        logger.info(f"✅ Workflow completed: {len(results.get('completed', []))} steps")
        return 0
    else:
        logger.error(f"❌ Workflow failed: {len(results.get('failed', []))} steps failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())

