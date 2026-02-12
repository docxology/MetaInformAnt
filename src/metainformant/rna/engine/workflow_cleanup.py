"""Workflow cleanup and disk management utilities.

Extracted from workflow.py to reduce its complexity. These functions handle
FASTQ/SRA cleanup after quantification, disk space monitoring, and
temporary file management.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Set, Tuple

from metainformant.core.utils import logging

if TYPE_CHECKING:
    from metainformant.rna.engine.workflow import AmalgkitWorkflowConfig

logger = logging.get_logger(__name__)


def check_disk_space(path: Path, min_free_gb: float = 10.0) -> Tuple[bool, float]:
    """Check if there's sufficient disk space.

    Args:
        path: Path to check disk space for
        min_free_gb: Minimum free space required in GB

    Returns:
        Tuple of (is_sufficient, free_gb)
    """
    try:
        stat = os.statvfs(path)
        free_gb = (stat.f_bavail * stat.f_frsize) / (1024**3)

        if free_gb < min_free_gb:
            logger.warning(f"Low disk space: {free_gb:.2f} GB free (minimum: {min_free_gb} GB) at {path}")
            return False, free_gb
        else:
            logger.debug(f"Disk space check: {free_gb:.2f} GB free at {path}")
            return True, free_gb
    except Exception as e:
        logger.warning(f"Could not check disk space: {e}")
        return True, 100.0  # Assume OK if we can't check


def check_disk_space_or_fail(path: Path, min_free_gb: float = 5.0, step_name: str = "") -> float:
    """Check disk space and raise RuntimeError if insufficient.

    Args:
        path: Path to check
        min_free_gb: Minimum free GB required
        step_name: Step name for error message

    Returns:
        Free GB available

    Raises:
        RuntimeError: If disk space is below minimum
    """
    ok, free_gb = check_disk_space(path, min_free_gb)
    if not ok and free_gb >= 0:
        raise RuntimeError(
            f"CRITICAL: Disk space too low to continue {step_name}. "
            f"Only {free_gb:.2f} GB free (need {min_free_gb} GB minimum). "
            f"Free up disk space and re-run with --steps {step_name} to resume."
        )
    return free_gb


def cleanup_temp_files(tmp_dir: Path, max_size_gb: float = 50.0) -> None:
    """Clean up temporary files if directory gets too large.

    Args:
        tmp_dir: Temporary directory to clean
        max_size_gb: Maximum size in GB before cleanup
    """
    if not tmp_dir.exists():
        return

    try:
        total_size = sum(f.stat().st_size for f in tmp_dir.rglob("*") if f.is_file())
        size_gb = total_size / (1024**3)

        if size_gb > max_size_gb:
            logger.warning(f"Temporary directory {tmp_dir} is {size_gb:.2f} GB (max: {max_size_gb} GB), cleaning up...")
            for item in tmp_dir.iterdir():
                if item.is_file():
                    item.unlink()
                elif item.is_dir():
                    shutil.rmtree(item)
            logger.info(f"Cleaned up temporary directory {tmp_dir}")
    except Exception as e:
        logger.warning(f"Could not clean up temporary directory {tmp_dir}: {e}")


def cleanup_incorrectly_placed_sra_files(getfastq_dir: Path) -> None:
    """Find and move SRA files from wrong locations to correct location.

    SRA Toolkit's prefetch command sometimes downloads SRA files to default
    locations instead of the specified output directory.

    Args:
        getfastq_dir: Directory where amalgkit expects SRA files
    """
    default_locations = [
        Path.home() / "ncbi" / "public" / "sra",
        Path("/tmp") / "ncbi" / "public" / "sra",
    ]

    moved_count = 0
    for default_loc in default_locations:
        if not default_loc.exists():
            continue

        for sra_file in default_loc.rglob("*.sra"):
            try:
                sample_id = sra_file.stem
                target_dir = getfastq_dir / sample_id
                target_dir.mkdir(parents=True, exist_ok=True)
                target_file = target_dir / sra_file.name

                if not target_file.exists():
                    logger.info(f"Moving SRA file from wrong location: {sra_file} -> {target_file}")
                    shutil.move(str(sra_file), str(target_file))
                    moved_count += 1
                else:
                    logger.debug(f"Target SRA file already exists, removing duplicate: {sra_file}")
                    sra_file.unlink()
            except Exception as e:
                logger.warning(f"Could not move SRA file {sra_file}: {e}")

    if moved_count > 0:
        logger.info(f"Moved {moved_count} SRA files from wrong locations to correct location")


def cleanup_fastqs(config: AmalgkitWorkflowConfig, sample_ids: List[str]) -> None:
    """Delete FASTQ files for specific samples after quantification.

    Args:
        config: Workflow configuration containing work_dir and step settings
        sample_ids: List of sample identifiers (SRR accessions) to clean up
    """
    getfastq_conf_dir = None
    if config.per_step and "getfastq" in config.per_step:
        gf_out = config.per_step["getfastq"].get("out_dir")
        if gf_out:
            getfastq_conf_dir = Path(gf_out)

    for sample_id in sample_ids:
        if not sample_id:
            continue
        paths = [
            config.work_dir / "fastq" / "getfastq" / sample_id,
            config.work_dir / "fastq" / sample_id,
            config.work_dir / "getfastq" / sample_id,
        ]

        if getfastq_conf_dir:
            paths.append(getfastq_conf_dir / "getfastq" / sample_id)

        for p in paths:
            if p.exists() and p.is_dir():
                try:
                    shutil.rmtree(p, ignore_errors=True)
                except Exception:
                    pass  # Best effort cleanup


def get_quantified_samples(config: AmalgkitWorkflowConfig) -> Set[str]:
    """Get set of sample IDs that already have successful quantification results.

    Args:
        config: Workflow configuration

    Returns:
        Set of sample IDs that have abundance.tsv files
    """
    quantified: Set[str] = set()

    steps_config = config.extra_config.get("steps", {})
    quant_dir_raw = steps_config.get("quant", {}).get("out_dir", config.work_dir / "quant")
    quant_dir = Path(quant_dir_raw)

    if not quant_dir.exists():
        return quantified

    patterns = ["**/abundance.tsv", "**/*_abundance.tsv"]
    for pattern in patterns:
        for abundance_file in quant_dir.glob(pattern):
            sample_id = abundance_file.parent.name
            if sample_id.startswith(("SRR", "ERR", "DRR")):
                if abundance_file.stat().st_size > 100:
                    quantified.add(sample_id)

    if quantified:
        logger.info(f"Found {len(quantified)} samples already quantified")

    return quantified


def cleanup_after_quant(config: AmalgkitWorkflowConfig, dry_run: bool = False) -> Dict[str, Any]:
    """Delete FASTQ and SRA files for samples with successful quantification.

    Args:
        config: Workflow configuration
        dry_run: If True, only report what would be deleted without deleting

    Returns:
        Dictionary with cleanup statistics
    """
    result: Dict[str, Any] = {
        "samples_cleaned": 0,
        "fastq_files_deleted": 0,
        "sra_files_deleted": 0,
        "bytes_freed": 0,
        "errors": [],
    }

    quantified_samples = get_quantified_samples(config)
    if not quantified_samples:
        logger.info("No quantified samples found - skipping cleanup")
        return result

    steps_config = config.extra_config.get("steps", {})
    fastq_dir_raw = steps_config.get("getfastq", {}).get("out_dir", config.work_dir / "fastq")
    fastq_dir = Path(fastq_dir_raw)

    if fastq_dir.name != "getfastq":
        getfastq_subdir = fastq_dir / "getfastq"
        if getfastq_subdir.exists():
            fastq_dir = getfastq_subdir

    if not fastq_dir.exists():
        logger.info(f"FASTQ directory does not exist: {fastq_dir}")
        return result

    logger.info(f"Cleaning up FASTQ/SRA files for {len(quantified_samples)} quantified samples")
    if dry_run:
        logger.info("DRY RUN - no files will be deleted")

    for sample_id in quantified_samples:
        sample_dir = fastq_dir / sample_id
        if not sample_dir.exists():
            continue

        files_to_delete: List[Path] = []

        for pattern in ["*.fastq.gz", "*.fastq", "*.fq.gz", "*.fq", "*.amalgkit.fastq.gz"]:
            files_to_delete.extend(sample_dir.glob(pattern))

        for pattern in ["*.sra", "*.sra.part"]:
            files_to_delete.extend(sample_dir.glob(pattern))

        if not files_to_delete:
            continue

        sample_bytes_freed = 0
        for file_path in files_to_delete:
            try:
                file_size = file_path.stat().st_size
                if not dry_run:
                    file_path.unlink()
                sample_bytes_freed += file_size

                if file_path.suffix in (".sra", ".part"):
                    result["sra_files_deleted"] += 1
                else:
                    result["fastq_files_deleted"] += 1

                logger.debug(
                    f"{'Would delete' if dry_run else 'Deleted'}: {file_path.name} ({file_size / 1024**2:.1f} MB)"
                )
            except Exception as e:
                error_msg = f"Failed to delete {file_path}: {e}"
                result["errors"].append(error_msg)
                logger.warning(error_msg)

        if sample_bytes_freed > 0:
            result["samples_cleaned"] += 1
            result["bytes_freed"] += sample_bytes_freed
            logger.info(f"{'Would clean' if dry_run else 'Cleaned'} {sample_id}: {sample_bytes_freed / 1024**3:.2f} GB")

    total_gb = result["bytes_freed"] / 1024**3
    action = "Would free" if dry_run else "Freed"
    logger.info(
        f"Cleanup complete: {result['samples_cleaned']} samples, "
        f"{result['fastq_files_deleted']} FASTQ files, "
        f"{result['sra_files_deleted']} SRA files, "
        f"{action} {total_gb:.2f} GB"
    )

    if result["errors"]:
        logger.warning(f"{len(result['errors'])} errors during cleanup")

    return result


def filter_metadata_for_unquantified(
    config: AmalgkitWorkflowConfig, source_metadata: Path, output_metadata: Path
) -> int:
    """Filter metadata to include only samples not yet quantified.

    Args:
        config: Workflow configuration
        source_metadata: Path to source metadata TSV
        output_metadata: Path to write filtered metadata

    Returns:
        Number of samples remaining after filtering
    """
    import csv

    quantified = get_quantified_samples(config)

    if not quantified:
        if source_metadata.exists():
            shutil.copy2(source_metadata, output_metadata)
            with open(source_metadata, "r") as f:
                return sum(1 for _ in f) - 1
        return 0

    try:
        with open(source_metadata, "r", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            fieldnames = reader.fieldnames
            rows = list(reader)
    except Exception as e:
        logger.error(f"Could not read metadata: {e}")
        return 0

    filtered_rows = []
    for row in rows:
        run_id = row.get("run", "")
        if run_id not in quantified:
            filtered_rows.append(row)
        else:
            logger.debug(f"Skipping already-quantified sample: {run_id}")

    skipped_count = len(rows) - len(filtered_rows)
    logger.info(f"Filtered metadata: {len(filtered_rows)} remaining, {skipped_count} already quantified")

    if not filtered_rows:
        logger.info("All samples already quantified - nothing to process")
        return 0

    try:
        with open(output_metadata, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(filtered_rows)
        logger.info(f"Created filtered metadata with {len(filtered_rows)} unprocessed samples: {output_metadata}")
    except Exception as e:
        logger.error(f"Could not write filtered metadata: {e}")
        return 0

    return len(filtered_rows)
