#!/usr/bin/env python3
"""Retry failed downloads for amalgkit workflow.

This script identifies failed samples from the workflow manifest and retries
their downloads with proper environment configuration.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.core.io import read_delimited
from metainformant.core.logging import get_logger
from metainformant.rna.workflow import load_workflow_config

logger = get_logger("retry_failed_downloads")


def find_repo_root(start_path: Path) -> Path:
    """Find repository root by looking for markers."""
    current = start_path.resolve()
    while current != current.parent:
        if any((current / marker).exists() for marker in [".git", "pyproject.toml", ".cursorrules"]):
            return current
        current = current.parent
    return start_path.resolve()


def setup_environment(repo_root: Path, work_dir: Path) -> None:
    """Set up environment variables for external drive usage."""
    sra_temp_dir = work_dir / "fastq" / "temp" / "sra"
    general_temp_dir = repo_root / "output" / ".tmp"
    
    sra_temp_dir.mkdir(parents=True, exist_ok=True)
    general_temp_dir.mkdir(parents=True, exist_ok=True)
    
    os.environ["TMPDIR"] = str(general_temp_dir)
    os.environ["TEMP"] = str(general_temp_dir)
    os.environ["TMP"] = str(general_temp_dir)
    os.environ["NCBI_SRA_REPOSITORY"] = str(sra_temp_dir)
    
    logger.info(f"Environment configured:")
    logger.info(f"  TMPDIR: {general_temp_dir}")
    logger.info(f"  NCBI_SRA_REPOSITORY: {sra_temp_dir}")


def get_failed_samples(manifest_path: Path) -> list[str]:
    """Extract failed sample IDs from manifest."""
    failed_samples = []
    
    if not manifest_path.exists():
        logger.warning(f"Manifest not found: {manifest_path}")
        return failed_samples
    
    try:
        with open(manifest_path) as f:
            for line in f:
                if not line.strip():
                    continue
                record = json.loads(line)
                if record.get("step") == "getfastq+quant (immediate)":
                    stats = record.get("params", {}).get("statistics", {})
                    failed_runs = stats.get("failed_runs", [])
                    failed_samples.extend(failed_runs)
                    break
    except Exception as e:
        logger.error(f"Error reading manifest: {e}")
    
    return failed_samples


def get_all_samples(metadata_path: Path) -> list[str]:
    """Get all sample IDs from metadata."""
    if not metadata_path.exists():
        logger.warning(f"Metadata not found: {metadata_path}")
        return []
    
    try:
        rows = list(read_delimited(metadata_path, delimiter="\t"))
        samples = [row.get("run") for row in rows if row.get("run")]
        return samples
    except Exception as e:
        logger.error(f"Error reading metadata: {e}")
        return []


def retry_downloads(
    config_path: Path,
    sample_ids: list[str] | None = None,
    *,
    max_retries: int = 3,
) -> int:
    """Retry downloads for specified samples."""
    repo_root = find_repo_root(config_path)
    config = load_workflow_config(config_path)
    work_dir = config.work_dir
    
    # Set up environment
    setup_environment(repo_root, work_dir)
    
    # Get samples to retry
    if sample_ids is None:
        # Get failed samples from manifest
        manifest_path = work_dir / "amalgkit.manifest.jsonl"
        sample_ids = get_failed_samples(manifest_path)
        
        if not sample_ids:
            # If no failed samples in manifest, check which samples don't have FASTQ files
            metadata_path = work_dir / "metadata" / "metadata.tsv"
            all_samples = get_all_samples(metadata_path)
            fastq_dir = Path(config.per_step.get("getfastq", {}).get("out_dir", work_dir / "fastq"))
            
            sample_ids = []
            for sample_id in all_samples:
                sample_dir = fastq_dir / "getfastq" / sample_id
                fastq_files = list(sample_dir.glob("*.fastq*")) if sample_dir.exists() else []
                if not fastq_files:
                    sample_ids.append(sample_id)
    
    if not sample_ids:
        logger.info("No samples to retry")
        return 0
    
    logger.info(f"Retrying downloads for {len(sample_ids)} samples")
    
    # Get getfastq parameters
    getfastq_params = config.per_step.get("getfastq", {})
    quant_params = config.per_step.get("quant", {})
    
    # Import step functions
    from metainformant.rna.steps.process_samples import run_download_quant_workflow
    
    # Create metadata file with just the samples to retry
    metadata_path = work_dir / "metadata" / "metadata.tsv"
    if not metadata_path.exists():
        logger.error(f"Metadata file not found: {metadata_path}")
        return 1
    
    all_rows = list(read_delimited(metadata_path, delimiter="\t"))
    retry_rows = [row for row in all_rows if row.get("run") in sample_ids]
    
    if not retry_rows:
        logger.warning("No matching rows found in metadata for retry samples")
        return 1
    
    # Create temporary metadata file
    temp_metadata = work_dir / "metadata" / "metadata.retry.tsv"
    from metainformant.core.io import write_delimited
    write_delimited(retry_rows, temp_metadata, delimiter="\t")
    
    logger.info(f"Created retry metadata with {len(retry_rows)} samples")
    
    # Run download+quant workflow
    try:
        stats = run_download_quant_workflow(
            temp_metadata,
            getfastq_params=getfastq_params,
            quant_params=quant_params,
            work_dir=work_dir,
            log_dir=config.log_dir or (work_dir.parent / "logs"),
            num_workers=1,  # Sequential for retries
            skip_completed=True,
        )
        
        logger.info("Retry results:")
        logger.info(f"  Processed: {stats['processed']}")
        logger.info(f"  Skipped: {stats['skipped']}")
        logger.info(f"  Failed: {stats['failed']}")
        
        if stats["failed"] > 0:
            logger.warning(f"  Failed samples: {stats['failed_runs'][:10]}{'...' if len(stats['failed_runs']) > 10 else ''}")
        
        return 0 if stats["failed"] == 0 else 1
    finally:
        # Clean up temp metadata
        if temp_metadata.exists():
            temp_metadata.unlink()


def main() -> int:
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Retry failed downloads for amalgkit workflow",
    )
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Path to workflow config file",
    )
    parser.add_argument(
        "--samples",
        nargs="+",
        help="Specific sample IDs to retry (default: all failed samples)",
    )
    parser.add_argument(
        "--max-retries",
        type=int,
        default=3,
        help="Maximum number of retry attempts (default: 3)",
    )
    
    args = parser.parse_args()
    
    config_path = Path(args.config).resolve()
    if not config_path.exists():
        logger.error(f"Config file not found: {config_path}")
        return 1
    
    return retry_downloads(
        config_path,
        sample_ids=args.samples,
        max_retries=args.max_retries,
    )


if __name__ == "__main__":
    sys.exit(main())

