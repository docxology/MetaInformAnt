#!/usr/bin/env python3
"""
Pre-flight check script for RNA-seq workflows.

Verifies that all required tools and dependencies are properly installed
before running amalgkit workflows.

Usage:
    source .venv/bin/activate
    python3 scripts/rna/check_environment.py
"""

import sys
from pathlib import Path

# Import setup utilities (must be before other imports)
sys.path.insert(0, str(Path(__file__).parent))
from _setup_utils import check_environment

# Add repo root to path for imports
repo_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(repo_root / "src"))

from metainformant.rna import (
    check_amalgkit,
    check_dependencies,
    check_kallisto,
    check_metainformant,
    check_sra_toolkit,
    check_virtual_env,
    validate_environment,
)
from metainformant.core.logging import get_logger

logger = get_logger("check_environment")


def main():
    """Run all environment checks and report results."""
    logger.info("=" * 80)
    logger.info("RNA-SEQ WORKFLOW ENVIRONMENT CHECK")
    logger.info("=" * 80)
    logger.info("")
    
    # Use setup_utils for basic checks
    success, missing, warnings = check_environment()
    if not success:
        logger.warning("⚠️  Basic environment check failed:")
        for item in missing:
            logger.warning(f"  ❌ Missing: {item}")
    else:
        logger.info("✅ Basic environment check passed")
    
    if warnings:
        logger.info("⚠️  Optional tools not found:")
        for warning in warnings:
            logger.info(f"  ⚠️  {warning}")
    
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
        logger.info("To run multi-species workflow:")
        logger.info("  python3 scripts/rna/run_multi_species.py")
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





