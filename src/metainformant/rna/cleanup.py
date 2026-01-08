"""RNA-seq workflow cleanup and maintenance utilities.

This module provides functions for cleaning up RNA-seq workflow artifacts,
managing disk space, and maintaining workflow state.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

from metainformant.core import logging

logger = logging.get_logger(__name__)


def find_partial_downloads(fastq_dir: Path, quant_dir: Path) -> List[Tuple[str, Path, int]]:
    """Find samples with FASTQ files but no quantification results.

    Args:
        fastq_dir: Directory containing FASTQ downloads
        quant_dir: Directory containing quantification results

    Returns:
        List of tuples: (sample_id, sample_dir, size_mb)
    """
    partial = []

    if not fastq_dir.exists():
        return partial

    # Check for getfastq subdirectory (amalgkit structure)
    search_dirs = [fastq_dir]
    getfastq_subdir = fastq_dir / "getfastq"
    if getfastq_subdir.exists():
        search_dirs = [getfastq_subdir]

    for search_dir in search_dirs:
        for sample_dir in search_dir.iterdir():
            if not sample_dir.is_dir():
                continue

            sample_id = sample_dir.name

            # Check if sample has FASTQ files
            fastq_files = list(sample_dir.glob("*.fastq*"))
            if not fastq_files:
                continue

            # Check if sample is quantified
            quant_sample_dir = quant_dir / sample_id
            abundance_file = quant_sample_dir / "abundance.tsv"
            quant_sf = quant_sample_dir / "quant.sf"

            if abundance_file.exists() or quant_sf.exists():
                continue  # Already quantified, skip

            # Calculate size
            size_bytes = sum(f.stat().st_size for f in fastq_files)
            size_mb = size_bytes // (1024 * 1024)

            partial.append((sample_id, sample_dir, size_mb))

    return partial


def cleanup_partial_downloads(
    config_or_path: Union[Path, str],
    dry_run: bool = False
) -> Dict[str, Any]:
    """Clean up partially downloaded FASTQ files.

    Args:
        config_or_path: Path to config file or work directory
        dry_run: If True, only report what would be deleted

    Returns:
        Dictionary with cleanup results: {deleted, freed_mb, errors}
    """
    import shutil

    result = {"deleted": 0, "freed_mb": 0, "errors": 0}

    # Handle both config file path and direct work_dir
    config_path = Path(config_or_path)

    if config_path.is_file():
        # Load config to get directories
        from metainformant.rna.workflow import load_workflow_config
        config = load_workflow_config(config_path)
        work_dir = config.work_dir

        # Get step-specific directories
        steps = config.extra_config.get('per_step', {})
        fastq_dir = Path(steps.get('getfastq', {}).get('out_dir', work_dir / "fastq"))
        quant_dir = Path(steps.get('quant', {}).get('out_dir', work_dir / "quant"))
    else:
        # Treat as work_dir directly
        work_dir = config_path
        fastq_dir = work_dir / "fastq"
        quant_dir = work_dir / "quant"

    # Find partial downloads
    partial = find_partial_downloads(fastq_dir, quant_dir)

    for sample_id, sample_dir, size_mb in partial:
        if dry_run:
            logger.info(f"Would delete: {sample_dir} ({size_mb}MB)")
        else:
            try:
                shutil.rmtree(sample_dir)
                result["deleted"] += 1
                result["freed_mb"] += size_mb
                logger.debug(f"Deleted: {sample_dir}")
            except Exception as e:
                logger.warning(f"Failed to delete {sample_dir}: {e}")
                result["errors"] += 1

    logger.info(f"Cleanup: {result['deleted']} deleted, {result['freed_mb']}MB freed, {result['errors']} errors")
    return result


def fix_abundance_naming(quant_dir: Path, sample_id: str) -> bool:
    """Create symlink for abundance file with sample ID prefix.

    Args:
        quant_dir: Quantification output directory
        sample_id: Sample identifier

    Returns:
        True if symlink created or already exists, False if source missing
    """
    sample_dir = quant_dir / sample_id
    source = sample_dir / "abundance.tsv"
    target = sample_dir / f"{sample_id}_abundance.tsv"

    if target.exists():
        return True  # Already exists

    if not source.exists():
        return False  # Source missing

    try:
        target.symlink_to(source.name)
        return True
    except Exception as e:
        logger.warning(f"Failed to create symlink {target}: {e}")
        return False


def fix_abundance_naming_for_species(
    config_or_path: Union[Path, str],
    species: str = None
) -> Tuple[int, int]:
    """Fix abundance file naming for all samples in a species workflow.

    Args:
        config_or_path: Config file path or work directory
        species: Species name (optional, for logging)

    Returns:
        Tuple of (created_count, already_exists_count)
    """
    config_path = Path(config_or_path)

    if config_path.is_file():
        from metainformant.rna.workflow import load_workflow_config
        config = load_workflow_config(config_path)
        work_dir = config.work_dir
        steps = config.extra_config.get('per_step', {})
        quant_dir = Path(steps.get('quant', {}).get('out_dir', work_dir / "quant"))
    else:
        quant_dir = config_path / "quant"

    if not quant_dir.exists():
        return (0, 0)

    created = 0
    already_exists = 0

    for sample_dir in quant_dir.iterdir():
        if not sample_dir.is_dir():
            continue

        sample_id = sample_dir.name
        target = sample_dir / f"{sample_id}_abundance.tsv"

        if target.exists():
            already_exists += 1
        elif fix_abundance_naming(quant_dir, sample_id):
            created += 1

    logger.info(f"Abundance naming: {created} created, {already_exists} already exist")
    return (created, already_exists)


def cleanup_unquantified_samples(work_dir: Path) -> List[str]:
    """Clean up samples that failed quantification.

    Args:
        work_dir: Workflow working directory

    Returns:
        List of cleaned up sample identifiers
    """
    cleaned_samples = []

    fastq_dir = work_dir / "fastq"
    quant_dir = work_dir / "quant"

    if not fastq_dir.exists() or not quant_dir.exists():
        return cleaned_samples

    # Get quantified samples
    quantified_samples = set()
    for quant_file in quant_dir.glob("*/abundance.tsv"):
        quantified_samples.add(quant_file.parent.name)
    for quant_file in quant_dir.glob("*/quant.sf"):
        quantified_samples.add(quant_file.parent.name)

    # Find and remove unquantified FASTQ files
    for fastq_file in fastq_dir.glob("*_1.fastq"):
        sample_name = fastq_file.name.replace("_1.fastq", "")

        if sample_name not in quantified_samples:
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

    # Clean up cache files
    cache_dir = work_dir / ".cache"
    if cache_dir.exists():
        for cache_file in cache_dir.glob("*"):
            if not (keep_logs and cache_file.name.endswith('.log')):
                try:
                    if cache_file.is_file():
                        cache_file.unlink()
                        stats['cache_files'] += 1
                except Exception:
                    pass

    # Clean up intermediate files
    intermediate_patterns = ["**/intermediate_*", "**/*.intermediate", "**/temp_*"]
    for pattern in intermediate_patterns:
        for f in work_dir.glob(pattern):
            try:
                f.unlink()
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









