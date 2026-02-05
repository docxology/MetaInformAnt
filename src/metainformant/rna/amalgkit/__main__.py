"""CLI entry point for amalgkit RNA-seq workflow.

Usage:
    python -m metainformant.rna.amalgkit --config path/to/config.yaml
    python -m metainformant.rna.amalgkit --config path/to/config.yaml --steps metadata select getfastq
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional

from metainformant.core import logging

logger = logging.get_logger(__name__)


def main(argv: Optional[List[str]] = None) -> int:
    """Main CLI entry point for amalgkit workflow execution.

    Args:
        argv: Command line arguments (defaults to sys.argv[1:])

    Returns:
        Exit code (0 for success, non-zero for failure)
    """
    parser = argparse.ArgumentParser(prog="metainformant.rna.amalgkit", description="Run amalgkit RNA-seq workflow")

    parser.add_argument("--config", "-c", type=str, required=True, help="Path to amalgkit YAML configuration file")

    parser.add_argument("--steps", "-s", nargs="*", help="Specific steps to run (default: all steps)")

    parser.add_argument("--species", type=str, help="Species name (overrides config if specified)")

    parser.add_argument("--dry-run", action="store_true", help="Print what would be executed without running")

    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose output")

    args = parser.parse_args(argv)

    # Validate config file exists
    config_path = Path(args.config)
    if not config_path.exists():
        logger.error(f"Configuration file not found: {config_path}")
        return 1

    if args.dry_run:
        logger.info("DRY RUN - would run amalgkit workflow with:")
        logger.info(f"  Config: {config_path}")
        logger.info(f"  Steps: {args.steps or 'all'}")
        logger.info(f"  Species: {args.species or 'from config'}")
        return 0

    # Import and run workflow
    try:
        from metainformant.rna.engine.orchestration import run_workflow_for_species

        logger.info(f"Starting amalgkit workflow with config: {config_path}")

        results = run_workflow_for_species(config_path=config_path, species=args.species, steps=args.steps)

        # Report results
        if results.get("success"):
            logger.info(f"Workflow completed successfully: {results.get('successful_steps', 0)} steps successful")
            return 0
        else:
            failed = results.get("failed", [])
            logger.error(f"Workflow failed: {len(failed)} steps failed: {', '.join(failed)}")
            return 1

    except Exception as e:
        logger.exception(f"Workflow execution failed: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
