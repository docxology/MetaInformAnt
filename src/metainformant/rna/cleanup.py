"""RNA-seq workflow cleanup and maintenance utilities.

This module provides functions for cleaning up RNA-seq workflow artifacts,
managing disk space, and maintaining workflow state.
"""

from __future__ import annotations

from pathlib import Path
from typing import List

from metainformant.core import logging

logger = logging.get_logger(__name__)


def cleanup_partial_downloads(work_dir: Path) -> int:
    """Clean up partially downloaded FASTQ files.

    Args:
        work_dir: Workflow working directory

    Returns:
        Number of files cleaned up
    """
    cleanup_count = 0

    # Look for incomplete downloads (files ending with .tmp or incomplete FASTQ files)
    for tmp_file in work_dir.glob("**/*.tmp"):
        try:
            tmp_file.unlink()
            cleanup_count += 1
            logger.debug(f"Removed partial download: {tmp_file}")
        except Exception as e:
            logger.warning(f"Failed to remove {tmp_file}: {e}")

    # Look for incomplete FASTQ files (very small files that might be incomplete)
    for fastq_file in work_dir.glob("**/*.fastq"):
        try:
            if fastq_file.stat().st_size < 1000:  # Less than 1KB, likely incomplete
                fastq_file.unlink()
                cleanup_count += 1
                logger.debug(f"Removed incomplete FASTQ: {fastq_file}")
        except Exception as e:
            logger.warning(f"Failed to check {fastq_file}: {e}")

    logger.info(f"Cleaned up {cleanup_count} partial/incomplete files")
    return cleanup_count


def fix_abundance_naming_for_species(work_dir: Path, species: str) -> bool:
    """Fix abundance file naming inconsistencies for a species.

    Args:
        work_dir: Workflow working directory
        species: Species name

    Returns:
        True if fixes were applied, False otherwise
    """
    # This is a placeholder implementation
    # In practice, this would rename abundance files to consistent naming

    abundance_dir = work_dir / "abundance"
    if not abundance_dir.exists():
        return False

    logger.info(f"Checking abundance file naming for {species} in {abundance_dir}")

    # Placeholder - would implement actual renaming logic
    return False


def cleanup_unquantified_samples(work_dir: Path) -> List[str]:
    """Clean up samples that failed quantification.

    Args:
        work_dir: Workflow working directory

    Returns:
        List of cleaned up sample identifiers
    """
    cleaned_samples = []

    # Look for samples that have FASTQ files but no quantification results
    fastq_dir = work_dir / "fastq"
    quant_dir = work_dir / "quant"

    if not fastq_dir.exists() or not quant_dir.exists():
        return cleaned_samples

    # Get list of quantified samples
    quantified_samples = set()
    for quant_file in quant_dir.glob("*.sf"):  # Salmon quantification files
        sample_name = quant_file.stem
        quantified_samples.add(sample_name)

    # Check for unquantified FASTQ files
    for fastq_file in fastq_dir.glob("*_1.fastq"):
        sample_name = fastq_file.name.replace("_1.fastq", "")

        if sample_name not in quantified_samples:
            # Remove unquantified sample files
            try:
                fastq_file.unlink()
                paired_file = fastq_file.with_name(f"{sample_name}_2.fastq")
                if paired_file.exists():
                    paired_file.unlink()

                cleaned_samples.append(sample_name)
                logger.debug(f"Removed unquantified sample: {sample_name}")

            except Exception as e:
                logger.warning(f"Failed to remove files for {sample_name}: {e}")

    logger.info(f"Cleaned up {len(cleaned_samples)} unquantified samples")
    return cleaned_samples


def cleanup_workflow_artifacts(work_dir: Path, keep_logs: bool = True) -> Dict[str, int]:
    """Clean up workflow artifacts while optionally keeping logs.

    Args:
        work_dir: Workflow working directory
        keep_logs: Whether to keep log files

    Returns:
        Dictionary with cleanup statistics
    """
    stats = {
        'temp_files': 0,
        'cache_files': 0,
        'intermediate_files': 0,
        'log_files': 0
    }

    # Clean up temporary files
    for pattern in ['*.tmp', '*.temp', '*.swp']:
        for temp_file in work_dir.glob(f"**/{pattern}"):
            try:
                temp_file.unlink()
                stats['temp_files'] += 1
            except Exception:
                pass

    # Clean up cache files (but keep logs if requested)
    cache_dir = work_dir / ".cache"
    if cache_dir.exists():
        for cache_file in cache_dir.glob("*"):
            if not (keep_logs and cache_file.name.endswith('.log')):
                try:
                    if cache_file.is_file():
                        cache_file.unlink()
                        stats['cache_files'] += 1
                    elif cache_file.is_dir():
                        # Remove empty directories
                        try:
                            cache_file.rmdir()
                        except OSError:
                            pass  # Directory not empty
                except Exception:
                    pass

    # Clean up intermediate files that are no longer needed
    intermediate_patterns = [
        "**/intermediate_*",
        "**/*.intermediate",
        "**/temp_*"
    ]

    for pattern in intermediate_patterns:
        for intermediate_file in work_dir.glob(pattern):
            try:
                intermediate_file.unlink()
                stats['intermediate_files'] += 1
            except Exception:
                pass

    if not keep_logs:
        for log_file in work_dir.glob("**/*.log"):
            try:
                log_file.unlink()
                stats['log_files'] += 1
            except Exception:
                pass

    total_cleaned = sum(stats.values())
    logger.info(f"Cleaned up {total_cleaned} workflow artifacts")

    return stats


