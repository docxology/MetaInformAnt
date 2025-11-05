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

# Add repo root to path for imports
repo_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(repo_root))

from metainformant.rna import (
    check_amalgkit,
    check_dependencies,
    check_kallisto,
    check_metainformant,
    check_sra_toolkit,
    check_virtual_env,
    validate_environment,
)


def main():
    """Run all environment checks and report results."""
    print("=" * 80)
    print("RNA-SEQ WORKFLOW ENVIRONMENT CHECK")
    print("=" * 80)
    print()
    
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
        print(f"{status} {name:<{max_name_len}} : {message}")
    
    print()
    print("=" * 80)
    
    if validation["all_passed"]:
        print("✅ ALL CHECKS PASSED - Ready to run RNA-seq workflows!")
        print()
        print("To run multi-species workflow:")
        print("  python3 scripts/rna/run_multi_species.py")
        return 0
    else:
        print("❌ SOME CHECKS FAILED - Please fix issues before running workflows")
        print()
        
        # Print recommendations
        if validation["recommendations"]:
            for rec in validation["recommendations"]:
                print(f"  - {rec}")
            print()
        
        print("For complete setup instructions, see:")
        print("  docs/rna/SETUP.md")
        return 1


if __name__ == "__main__":
    sys.exit(main())





