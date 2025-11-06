#!/usr/bin/env python3
"""Pre-flight check script for RNA-seq workflows.

This is a thin wrapper that calls metainformant.rna.environment.validate_environment()
to verify all required tools and dependencies.

Usage:
    python3 scripts/rna/check_environment.py
"""

import sys
from pathlib import Path

# Import setup utilities (must be before other imports)
sys.path.insert(0, str(Path(__file__).parent))
from _setup_utils import ensure_venv_activated, check_environment_or_exit

# Auto-setup and activate venv
ensure_venv_activated(auto_setup=True)
check_environment_or_exit(auto_setup=True)

# Add repo root to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna import validate_environment
from metainformant.core.logging import get_logger

logger = get_logger("check_environment")


def main() -> int:
    """Run all environment checks and report results."""
    logger.info("=" * 80)
    logger.info("RNA-SEQ WORKFLOW ENVIRONMENT CHECK")
    logger.info("=" * 80)
    logger.info("")

    # Use metainformant function for comprehensive validation
    validation = validate_environment()
    deps = validation["dependencies"]

    checks = [
        ("Virtual Environment", deps["virtual_env"]),
        ("metainformant", deps["metainformant"]),
        ("amalgkit", deps["amalgkit"]),
        ("SRA Toolkit (fasterq-dump)", deps["sra_toolkit"]),
        ("kallisto", deps["kallisto"]),
    ]

    # Print results
    max_name_len = max(len(name) for name, _ in checks)
    for name, (passed, message) in checks:
        status = "✅" if passed else "❌"
        logger.info(f"{status} {name:<{max_name_len}} : {message}")

    logger.info("")
    logger.info("=" * 80)

    if validation["all_passed"]:
        logger.info("✅ ALL CHECKS PASSED - Ready to run RNA-seq workflows!")
        logger.info("")
        logger.info("To run workflow for a species:")
        logger.info("  python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_<species>.yaml")
        return 0
    else:
        logger.error("❌ SOME CHECKS FAILED - Please fix issues before running workflows")
        logger.info("")

        # Print recommendations
        if validation["recommendations"]:
            for rec in validation["recommendations"]:
                logger.info(f"  - {rec}")
            logger.info("")

        logger.info("For complete setup instructions, see:")
        logger.info("  docs/rna/GETTING_STARTED.md")
        return 1


if __name__ == "__main__":
    sys.exit(main())





