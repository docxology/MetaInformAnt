#!/usr/bin/env python3
"""
Simple test runner for METAINFORMANT test suite.
"""

import sys
import os
import subprocess
from pathlib import Path

def main():
    # Add src to path
    repo_root = Path(__file__).resolve().parent
    src_dir = repo_root / "src"
    if str(src_dir) not in sys.path:
        sys.path.insert(0, str(src_dir))

    print("🚀 Running METAINFORMANT test suite...")

    # Change to repo root
    os.chdir(repo_root)

    # Run pytest with output
    cmd = [
        sys.executable, "-m", "pytest",
        "tests/",
        "-v",
        "--tb=short",
        "--durations=10",
        "--maxfail=5"
    ]

    print(f"Running: {' '.join(cmd)}")
    print("=" * 60)

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        print(result.stdout)
        if result.stderr:
            print("STDERR:")
            print(result.stderr)
        print("=" * 60)
        print(f"Return code: {result.returncode}")
        return result.returncode
    except subprocess.TimeoutExpired:
        print("❌ Test suite timed out after 5 minutes")
        return 1
    except Exception as e:
        print(f"❌ Error running tests: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
