"""Cleanup functions for RNA-seq workflows.

This module provides functions to clean up partial downloads, unquantified samples,
and fix file naming issues.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import Any

from ..core.logging import get_logger
from .workflow import load_workflow_config

logger = get_logger(__name__)


def find_partial_downloads(fastq_dir: Path, quant_dir: Path) -> list[tuple[str, Path, int]]:
    """Find samples with partial downloads that aren't quantified.

    Args:
        fastq_dir: Directory containing FASTQ files
        quant_dir: Directory containing quantification results

    Returns:
        List of (sample_id, sample_dir, size_mb) tuples
    """
    partial_samples = []

    if not fastq_dir.exists():
        return []

    # Check getfastq subdirectory
    getfastq_dir = fastq_dir / "getfastq"
    if getfastq_dir.exists():
        for sample_dir in getfastq_dir.glob("SRR*"):
            if not sample_dir.is_dir():
                continue

            sample_id = sample_dir.name

            # Check if quantified
            abundance_file = quant_dir / sample_id / "abundance.tsv"
            if abundance_file.exists():
                continue

            # Check for files
            has_files = False
            for pattern in ["*.fastq*", "*.sra"]:
                if list(sample_dir.glob(pattern)):
                    has_files = True
                    break

            if has_files:
                # Calculate size
                size_mb = sum(f.stat().st_size for f in sample_dir.rglob("*") if f.is_file()) / (1024 * 1024)
                partial_samples.append((sample_id, sample_dir, int(size_mb)))

    # Check direct structure
    for sample_dir in fastq_dir.glob("SRR*"):
        if not sample_dir.is_dir():
            continue

        sample_id = sample_dir.name

        # Check if quantified
        abundance_file = quant_dir / sample_id / "abundance.tsv"
        if abundance_file.exists():
            continue

        # Check for files
        has_files = False
        for pattern in ["*.fastq*", "*.sra"]:
            if list(sample_dir.glob(pattern)):
                has_files = True
                break

        if has_files:
            # Calculate size
            size_mb = sum(f.stat().st_size for f in sample_dir.rglob("*") if f.is_file()) / (1024 * 1024)
            partial_samples.append((sample_id, sample_dir, int(size_mb)))

    return partial_samples


def cleanup_partial_downloads(
    config_path: Path,
    *,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Clean up partial downloads for a species.

    Args:
        config_path: Path to species workflow config file
        dry_run: If True, only report what would be deleted

    Returns:
        Dictionary with 'deleted', 'freed_mb', 'errors' keys
    """
    cfg = load_workflow_config(config_path)
    fastq_dir = Path(cfg.per_step.get("getfastq", {}).get("out_dir", cfg.work_dir / "fastq"))
    quant_dir = Path(cfg.per_step.get("quant", {}).get("out_dir", cfg.work_dir / "quant"))

    partial_samples = find_partial_downloads(fastq_dir, quant_dir)

    if not partial_samples:
        return {"deleted": 0, "freed_mb": 0, "errors": 0}

    logger.info(f"Found {len(partial_samples)} partial downloads")

    deleted_count = 0
    freed_mb = 0
    errors = 0

    for sample_id, sample_dir, size_mb in partial_samples:
        try:
            if dry_run:
                logger.info(f"  [DRY RUN] Would delete {sample_id} ({size_mb}MB)")
            else:
                shutil.rmtree(sample_dir)
                logger.info(f"  🗑️  Deleted {sample_id} ({size_mb}MB)")
                deleted_count += 1
                freed_mb += size_mb
        except Exception as e:
            logger.warning(f"  ⚠️  Failed to delete {sample_dir}: {e}")
            errors += 1

    return {
        "deleted": deleted_count,
        "freed_mb": freed_mb,
        "errors": errors,
    }


def fix_abundance_naming(quant_dir: Path, sample_id: str) -> bool:
    """Create symlink from abundance.tsv to {SRR}_abundance.tsv for amalgkit merge.

    Args:
        quant_dir: Directory containing quantification results
        sample_id: Sample ID (e.g., 'SRR1234567')

    Returns:
        True if symlink was created or already exists
    """
    sample_dir = quant_dir / sample_id
    source = sample_dir / "abundance.tsv"
    target = sample_dir / f"{sample_id}_abundance.tsv"

    if not source.exists():
        return False

    if target.exists():
        return True

    try:
        target.symlink_to("abundance.tsv")  # Relative symlink
        return True
    except Exception as e:
        logger.warning(f"Failed to create symlink for {sample_id}: {e}")
        return False


def fix_abundance_naming_for_species(
    config_path: Path,
) -> tuple[int, int]:
    """Fix abundance naming for all samples in a species.

    Args:
        config_path: Path to species workflow config file

    Returns:
        Tuple of (created_count, already_exists_count)
    """
    cfg = load_workflow_config(config_path)
    quant_dir = Path(cfg.per_step.get("quant", {}).get("out_dir", cfg.work_dir / "quant"))

    if not quant_dir.exists():
        logger.error(f"Quant directory not found: {quant_dir}")
        return 0, 0

    created = 0
    already_exists = 0

    # Find all SRR* directories
    srr_dirs = sorted([d for d in quant_dir.iterdir() if d.is_dir() and d.name.startswith("SRR")])

    logger.info(f"Processing {len(srr_dirs)} samples")

    for srr_dir in srr_dirs:
        srr_id = srr_dir.name
        source = srr_dir / "abundance.tsv"
        target = srr_dir / f"{srr_id}_abundance.tsv"

        if not source.exists():
            continue

        if target.exists():
            if target.is_symlink():
                already_exists += 1
            continue

        if fix_abundance_naming(quant_dir, srr_id):
            created += 1

    logger.info(f"✅ Created {created} symlinks, {already_exists} already existed")
    return created, already_exists

