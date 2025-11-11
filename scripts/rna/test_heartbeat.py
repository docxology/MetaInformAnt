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

from metainformant.rna.steps.getfastq import run as run_getfastq

# Use existing metadata file and filter to a single sample
# This avoids amalgkit bugs with minimal metadata
full_metadata = Path("output/amalgkit/pogonomyrmex_barbatus/work/metadata/metadata.tsv")
if not full_metadata.exists():
    print(f"Error: Full metadata file not found: {full_metadata}")
    sys.exit(1)

# Read full metadata and extract first sample
import csv
with open(full_metadata, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    rows = list(reader)
    if not rows:
        print("Error: No samples in metadata file")
        sys.exit(1)
    
    # Use first sample that hasn't been completed yet
    test_sample = None
    fastq_dir = Path("output/amalgkit/pogonomyrmex_barbatus/fastq/getfastq")
    for row in rows:
        sample_id = row.get('run', '').strip()
        if sample_id:
            sample_path = fastq_dir / sample_id
            # Check if already has FASTQ files
            if not sample_path.exists() or not any(sample_path.glob("*.fastq.gz")):
                test_sample = sample_id
                break
    
    if not test_sample:
        # Fall back to first sample
        test_sample = rows[0].get('run', '').strip()
        if not test_sample:
            print("Error: No valid sample ID found in metadata")
            sys.exit(1)

# Create filtered metadata with full row data
header = list(rows[0].keys())
test_row = rows[0] if test_sample == rows[0].get('run', '').strip() else next((r for r in rows if r.get('run', '').strip() == test_sample), rows[0])

metadata_path = Path(f"output/amalgkit/pogonomyrmex_barbatus/work/metadata/metadata.single.{test_sample}.tsv")
metadata_path.parent.mkdir(parents=True, exist_ok=True)

# Write full metadata row (not just run ID)
with open(metadata_path, 'w') as f:
    writer = csv.DictWriter(f, fieldnames=header, delimiter='\t')
    writer.writeheader()
    writer.writerow(test_row)

print(f"Testing heartbeat logging with sample: {test_sample}")
print("=" * 80)
print("Watch for heartbeat messages with size, rate, and current sample info")
print("=" * 80)
print()

# Run getfastq using step wrapper (handles pfd properly)
result = run_getfastq(
    {
        "metadata": str(metadata_path.absolute()),
        "out_dir": "output/amalgkit/pogonomyrmex_barbatus/fastq",
        "threads": 24,
        "aws": True,
        "gcp": True,
        "ncbi": True,
        "fastp": True,
        "pfd": False,  # Explicitly disable pfd
    },
    log_dir="output/amalgkit/pogonomyrmex_barbatus/logs",
)

print()
if hasattr(result, 'returncode'):
print(f"Test completed with return code: {result.returncode}")
else:
    print("Test completed")

