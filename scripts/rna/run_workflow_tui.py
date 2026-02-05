#!/usr/bin/env python3
"""Full workflow runner with TUI visualization.

This script runs the complete RNA-seq workflow (download → getfastq → quant)
with a real-time terminal interface showing per-sample progress through each stage.

Usage:
    python scripts/rna/run_workflow_tui.py --config config/amalgkit/amalgkit_pbarbatus_5sample.yaml
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

# Import setup utilities (must be before other imports)
sys.path.insert(0, str(Path(__file__).parent))
from _setup_utils import check_environment_or_exit, ensure_venv_activated

# Suppress optional dependency warnings until venv is ready
from metainformant.core.utils.optional_deps import enable_optional_warnings, suppress_optional_warnings

suppress_optional_warnings()

# Auto-setup and activate venv
ensure_venv_activated(auto_setup=True)
check_environment_or_exit(auto_setup=True)

# Now enable warnings since venv should be active
enable_optional_warnings()

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.core.engine.workflow_manager import WorkflowManager
from metainformant.core.utils.logging import get_logger

logger = get_logger("run_workflow_tui")


def load_samples_from_metadata(metadata_path: Path) -> list[dict]:
    """Load sample info from metadata TSV."""
    samples = []

    with open(metadata_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample_id = row.get("Run") or row.get("run") or row.get("SRR")
            if not sample_id:
                continue

            # Construct URL (prioritize AWS)
            url = None
            if row.get("AWS_Link"):
                url = row["AWS_Link"]
            elif row.get("sra_url"):
                url = row["sra_url"]
            else:
                # Fallback to constructing from SRA ID
                url = f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{sample_id}/{sample_id}"

            samples.append(
                {
                    "sample_id": sample_id,
                    "url": url,
                }
            )

    return samples


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run full RNA-seq workflow with TUI visualization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Path to workflow config YAML",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=5,
        help="Number of parallel download threads (default: 5)",
    )

    args = parser.parse_args()

    if not args.config.exists():
        logger.error(f"Config file not found: {args.config}")
        return 1

    # Load config to find metadata
    import yaml

    with open(args.config) as f:
        config = yaml.safe_load(f)

    work_dir = Path(config.get("work_dir", "output/amalgkit"))
    species = config.get("species", "unknown")

    # Find metadata file - work_dir from config already points to the work directory
    metadata_path = work_dir / "metadata" / "metadata_selected.tsv"

    if not metadata_path.exists():
        # Try with species subdirectory (alternate structure)
        metadata_path = work_dir / species / "work" / "metadata" / "metadata_selected.tsv"

    if not metadata_path.exists():
        logger.error(f"Metadata not found at {metadata_path}")
        logger.info("Run the metadata step first to select samples.")
        return 1

    # Load samples
    samples = load_samples_from_metadata(metadata_path)
    if not samples:
        logger.error("No samples found in metadata")
        return 1

    logger.info(f"Found {len(samples)} samples to process")

    # Initialize workflow manager
    manager = WorkflowManager(args.config, max_threads=args.threads)

    # Add samples - fastq is at same level as work (parent dir / fastq)
    base_dir = work_dir.parent  # Go up from 'work' to 'pbarbatus_test5'
    fastq_dir = base_dir / "fastq" / "getfastq"

    if not fastq_dir.exists():
        fastq_dir.mkdir(parents=True, exist_ok=True)

    for sample in samples:
        sample_id = sample["sample_id"]
        dest = fastq_dir / sample_id / f"{sample_id}.sra"
        dest.parent.mkdir(parents=True, exist_ok=True)

        manager.add_sample(sample_id, sample["url"], dest)

    # Run workflow
    print(f"\n🚀 Starting RNA-seq workflow for {len(samples)} samples...\n")
    results = manager.run()

    # Summary
    success = sum(1 for v in results.values() if v)
    failed = len(results) - success

    print(f"\n{'='*60}")
    print(f"Workflow Complete")
    print(f"{'='*60}")
    print(f"  ✓ Successful: {success}")
    print(f"  ✗ Failed: {failed}")
    print(f"  Total: {len(results)}")
    print()

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
