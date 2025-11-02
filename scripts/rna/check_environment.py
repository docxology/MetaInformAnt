#!/usr/bin/env python3
"""
Pre-flight check script for RNA-seq workflows.

Verifies that all required tools and dependencies are properly installed
before running amalgkit workflows.

Usage:
    source .venv/bin/activate
    python3 scripts/rna/check_environment.py
"""

import shutil
import sys
from pathlib import Path


def check_amalgkit():
    """Check if amalgkit is available and get version."""
    amalgkit_path = shutil.which("amalgkit")
    if not amalgkit_path:
        return False, "Not found on PATH"
    
    try:
        import subprocess
        result = subprocess.run(
            ["amalgkit", "--version"],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0:
            # Extract version from output
            for line in result.stdout.split("\n"):
                if "version" in line.lower():
                    return True, line.strip()
            return True, "Found but version unknown"
        else:
            return False, f"Error: {result.stderr}"
    except Exception as e:
        return False, f"Error checking version: {e}"


def check_sra_toolkit():
    """Check if SRA Toolkit is installed."""
    fasterq_path = shutil.which("fasterq-dump")
    if not fasterq_path:
        return False, "fasterq-dump not found"
    
    try:
        import subprocess
        result = subprocess.run(
            ["fasterq-dump", "--version"],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0:
            for line in result.stdout.split("\n"):
                if "fasterq-dump" in line.lower():
                    return True, line.strip()
            return True, "Found but version unknown"
        else:
            return False, f"Error: {result.stderr}"
    except Exception as e:
        return False, f"Error checking version: {e}"


def check_kallisto():
    """Check if kallisto is installed."""
    kallisto_path = shutil.which("kallisto")
    if not kallisto_path:
        return False, "kallisto not found"
    
    try:
        import subprocess
        result = subprocess.run(
            ["kallisto", "version"],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0:
            return True, result.stdout.strip()
        else:
            return False, f"Error: {result.stderr}"
    except Exception as e:
        return False, f"Error checking version: {e}"


def check_metainformant():
    """Check if metainformant package is installed."""
    try:
        import metainformant
        version = getattr(metainformant, "__version__", "unknown")
        return True, f"v{version}"
    except ImportError as e:
        return False, f"Import error: {e}"


def check_virtual_env():
    """Check if running inside a virtual environment."""
    # Check for VIRTUAL_ENV environment variable
    import os
    venv_path = os.environ.get("VIRTUAL_ENV")
    
    if venv_path:
        return True, f"Active: {venv_path}"
    
    # Check if sys.prefix differs from sys.base_prefix
    if sys.prefix != sys.base_prefix:
        return True, f"Active: {sys.prefix}"
    
    return False, "Not in virtual environment"


def main():
    """Run all environment checks and report results."""
    print("=" * 80)
    print("RNA-SEQ WORKFLOW ENVIRONMENT CHECK")
    print("=" * 80)
    print()
    
    checks = [
        ("Virtual Environment", check_virtual_env),
        ("metainformant", check_metainformant),
        ("amalgkit", check_amalgkit),
        ("SRA Toolkit (fasterq-dump)", check_sra_toolkit),
        ("kallisto", check_kallisto),
    ]
    
    all_passed = True
    results = []
    
    for name, check_func in checks:
        try:
            passed, message = check_func()
            status = "✅" if passed else "❌"
            results.append((name, status, message))
            if not passed:
                all_passed = False
        except Exception as e:
            results.append((name, "❌", f"Check failed: {e}"))
            all_passed = False
    
    # Print results
    max_name_len = max(len(name) for name, _, _ in results)
    for name, status, message in results:
        print(f"{status} {name:<{max_name_len}} : {message}")
    
    print()
    print("=" * 80)
    
    if all_passed:
        print("✅ ALL CHECKS PASSED - Ready to run RNA-seq workflows!")
        print()
        print("To run multi-species workflow:")
        print("  python3 scripts/rna/run_multi_species.py")
        return 0
    else:
        print("❌ SOME CHECKS FAILED - Please fix issues before running workflows")
        print()
        
        # Provide specific guidance
        if not results[0][1] == "✅":  # Virtual environment
            print("Virtual Environment Issue:")
            print("  Activate the virtual environment first:")
            print("    source .venv/bin/activate")
            print()
        
        if not results[1][1] == "✅":  # metainformant
            print("metainformant Not Installed:")
            print("  Install the package:")
            print("    uv pip install -e .")
            print()
        
        if not results[2][1] == "✅":  # amalgkit
            print("amalgkit Not Installed:")
            print("  Install amalgkit:")
            print("    uv pip install git+https://github.com/kfuku52/amalgkit")
            print()
        
        if not results[3][1] == "✅":  # SRA Toolkit
            print("SRA Toolkit Not Installed:")
            print("  Install SRA Toolkit:")
            print("    sudo apt-get update")
            print("    sudo apt-get install -y sra-toolkit")
            print()
        
        if not results[4][1] == "✅":  # kallisto
            print("kallisto Not Installed:")
            print("  Install kallisto:")
            print("    sudo apt-get update")
            print("    sudo apt-get install -y kallisto")
            print()
        
        print("For complete setup instructions, see:")
        print("  docs/rna/SETUP.md")
        return 1


if __name__ == "__main__":
    sys.exit(main())

