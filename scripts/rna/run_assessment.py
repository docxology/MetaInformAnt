#!/usr/bin/env python3
"""Run comprehensive assessment using full_assessment.py."""
import subprocess
import sys
from pathlib import Path

script_path = Path(__file__).parent / 'full_assessment.py'
repo_root = Path(__file__).parent.parent.parent

print("Running comprehensive assessment...")
print("=" * 100)
print()

try:
    result = subprocess.run(
        [sys.executable, str(script_path)],
        cwd=str(repo_root),
        text=True,
        check=False
    )
    sys.exit(result.returncode)
except Exception as e:
    print(f"Error running assessment: {e}")
    sys.exit(1)



