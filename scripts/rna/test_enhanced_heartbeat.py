#!/usr/bin/env python3
"""Test heartbeat logging on a single sample."""

from __future__ import annotations

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

from metainformant.rna.amalgkit import getfastq

# Test with a single sample
test_sample = "SRR14740500"  # Change to any sample ID you want to test

# Create temporary metadata for single sample
metadata_content = f"run\n{test_sample}\n"
metadata_path = Path(f"output/amalgkit/pogonomyrmex_barbatus/work/metadata/metadata.single.{test_sample}.tsv")
metadata_path.parent.mkdir(parents=True, exist_ok=True)
metadata_path.write_text(metadata_content)

print(f"Testing heartbeat logging with sample: {test_sample}")
print("=" * 80)
print("Watch for heartbeat messages with size, rate, and current sample info")
print("=" * 80)
print()

# Run getfastq
result = getfastq(
    {
        "metadata": str(metadata_path.absolute()),
        "out_dir": "output/amalgkit/pogonomyrmex_barbatus/fastq",
        "threads": 24,
        "aws": True,
        "gcp": True,
        "ncbi": True,
        "fastp": True,
    },
    log_dir="output/amalgkit/pogonomyrmex_barbatus/logs",
    step_name="getfastq_test",
)

print()
print(f"Test completed with return code: {result.returncode}")

