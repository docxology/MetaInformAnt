#!/usr/bin/env python3
"""Test script to quantify a single completed sample and clean up FASTQ files.

This script demonstrates the workflow:
1. Quantify sample using kallisto index
2. Verify quantification succeeded
3. Delete raw FASTQ files to free disk space

Usage:
    python3 scripts/rna/test_quantify_sample.py --sample SRR14740514 --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.core.io import read_delimited
from metainformant.core.logging import get_logger
from metainformant.rna.steps.getfastq import delete_sample_fastqs
from metainformant.rna.steps.quant import quantify_sample
from metainformant.rna.workflow import load_workflow_config

logger = get_logger("test_quantify_sample")


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Quantify a single sample and clean up FASTQ files",
    )
    parser.add_argument(
        "--sample",
        type=str,
        required=True,
        help="Sample ID (e.g., SRR14740514)",
    )
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Path to workflow config file",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without executing",
    )

    args = parser.parse_args()

    sample_id = args.sample
    config_path = args.config.resolve()

    if not config_path.exists():
        logger.error(f"Config file not found: {config_path}")
        return 1

    # Load config
    logger.info(f"Loading config from: {config_path}")
    cfg = load_workflow_config(config_path)

    # Get paths
    work_dir = cfg.work_dir
    fastq_dir_str = cfg.per_step.get("getfastq", {}).get("out_dir", str(work_dir / "fastq"))
    fastq_dir = Path(fastq_dir_str)
    if not fastq_dir.is_absolute():
        fastq_dir = fastq_dir.resolve()
    
    quant_dir_str = cfg.per_step.get("quant", {}).get("out_dir", str(work_dir / "quant"))
    quant_dir = Path(quant_dir_str)
    if not quant_dir.is_absolute():
        quant_dir = quant_dir.resolve()
    
    metadata_file = work_dir / "metadata" / "metadata.tsv"

    if not metadata_file.exists():
        logger.error(f"Metadata file not found: {metadata_file}")
        return 1

    # Read metadata
    logger.info(f"Reading metadata from: {metadata_file}")
    rows = list(read_delimited(metadata_file, delimiter="\t"))
    sample_rows = [row for row in rows if row.get("run") == sample_id]

    if not sample_rows:
        logger.error(f"Sample {sample_id} not found in metadata")
        return 1

    logger.info(f"Found metadata for {sample_id}")

    # Check if FASTQ files exist (check both getfastq subdirectory and direct structure)
    # Ensure paths are absolute
    fastq_dir = fastq_dir.resolve()
    sample_fastq_dir = fastq_dir / "getfastq" / sample_id
    if not sample_fastq_dir.exists():
        sample_fastq_dir = fastq_dir / sample_id

    fastq_files = list(sample_fastq_dir.glob("*.fastq*")) if sample_fastq_dir.exists() else []
    if not fastq_files:
        logger.warning(f"No FASTQ files found for {sample_id}")
        logger.warning(f"  Checked: {sample_fastq_dir}")
        logger.info("Checking if already quantified...")
        abundance_file = quant_dir / sample_id / "abundance.tsv"
        if abundance_file.exists():
            logger.info(f"Sample {sample_id} already quantified at {abundance_file}")
            return 0
        return 1

    logger.info(f"Found {len(fastq_files)} FASTQ files for {sample_id} in {sample_fastq_dir}")

    # Check if already quantified
    abundance_file = quant_dir / sample_id / "abundance.tsv"
    if abundance_file.exists():
        logger.info(f"Sample {sample_id} already quantified at {abundance_file}")
        if args.dry_run:
            logger.info("DRY RUN: Would delete FASTQ files")
            return 0
        logger.info("Deleting FASTQ files...")
        delete_sample_fastqs(sample_id, fastq_dir)
        logger.info("✓ FASTQ files deleted")
        return 0

    # Prepare quant params
    quant_params = dict(cfg.per_step.get("quant", {}))
    quant_params["out_dir"] = str(quant_dir.absolute())
    quant_params["threads"] = cfg.threads or 12

    # Inject index_dir if not already set
    if "index_dir" not in quant_params and "index-dir" not in quant_params:
        index_dir = work_dir / "index"
        if index_dir.exists():
            quant_params["index_dir"] = str(index_dir.absolute())
            logger.info(f"Using index: {quant_params['index_dir']}")
        else:
            logger.error(f"Index directory not found: {index_dir}")
            return 1

    # Note: amalgkit quant automatically finds FASTQ files based on metadata
    # and standard directory structure relative to work_dir
    # It looks in work_dir/fastq/getfastq/<SRR>/ or work_dir/fastq/<SRR>/
    # We need to pass work_dir to run_amalgkit so it can find the fastq directory

    log_dir = cfg.log_dir or (work_dir.parent / "logs")
    
    # Add work_dir to quant_params so quantify_sample can use it
    quant_params["work_dir"] = str(work_dir.absolute())

    if args.dry_run:
        logger.info("DRY RUN: Would quantify sample and delete FASTQ files")
        logger.info(f"  Sample: {sample_id}")
        logger.info(f"  Quant dir: {quant_dir}")
        logger.info(f"  Index dir: {quant_params.get('index_dir')}")
        logger.info(f"  FASTQ files: {len(fastq_files)}")
        return 0

    # Step 1: Quantify
    logger.info(f"Quantifying {sample_id}...")
    success, message, abundance_path = quantify_sample(
        sample_id=sample_id,
        metadata_rows=sample_rows,
        quant_params=quant_params,
        log_dir=log_dir,
        step_name=f"quant_{sample_id}",
    )

    if not success:
        logger.error(f"Quantification failed: {message}")
        return 1

    if not abundance_path or not abundance_path.exists():
        logger.error(f"Abundance file not found at expected location: {abundance_path}")
        return 1

    logger.info(f"✓ Quantification successful: {abundance_path}")

    # Step 2: Verify abundance file
    if abundance_file.exists():
        file_size = abundance_file.stat().st_size
        logger.info(f"✓ Abundance file verified: {file_size} bytes")
    else:
        logger.error(f"Abundance file not found at: {abundance_file}")
        return 1

    # Step 3: Delete FASTQ files
    logger.info(f"Deleting FASTQ files for {sample_id}...")
    delete_sample_fastqs(sample_id, fastq_dir)

    # Step 4: Verify deletion
    remaining_fastq = list(sample_fastq_dir.glob("*.fastq*")) if sample_fastq_dir.exists() else []
    if remaining_fastq:
        logger.warning(f"Some FASTQ files still exist: {remaining_fastq}")
        return 1

    logger.info(f"✓ FASTQ files deleted successfully")

    # Summary
    logger.info("=" * 80)
    logger.info("SUCCESS: Sample processed and cleaned up")
    logger.info(f"  Sample: {sample_id}")
    logger.info(f"  Abundance file: {abundance_file}")
    logger.info(f"  FASTQ files: Deleted")
    logger.info("=" * 80)

    return 0


if __name__ == "__main__":
    sys.exit(main())

