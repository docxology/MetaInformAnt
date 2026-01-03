"""Disk space monitoring and file system utilities for METAINFORMANT.

This module provides utilities for monitoring disk usage, managing file cleanup,
and ensuring adequate disk space for bioinformatics operations.
"""

from __future__ import annotations

import os
import shutil
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core.logging import get_logger

logger = get_logger(__name__)


def get_disk_usage(path: str | Path) -> Dict[str, Any]:
    """Get disk usage statistics for a path.

    Args:
        path: Path to check disk usage for

    Returns:
        Dictionary with disk usage statistics:
        - total: Total disk space in bytes
        - used: Used disk space in bytes
        - free: Free disk space in bytes
        - percent: Percentage used
        - path: Path checked

    Raises:
        OSError: If disk usage cannot be determined
    """
    path = Path(path)

    try:
        stat = os.statvfs(path)

        # Calculate sizes in bytes
        total = stat.f_bsize * stat.f_blocks
        free = stat.f_bsize * stat.f_bavail
        used = total - free

        # Calculate percentage
        percent = (used / total * 100) if total > 0 else 0

        return {
            'total': total,
            'used': used,
            'free': free,
            'percent': round(percent, 1),
            'path': str(path),
        }

    except OSError as e:
        logger.error(f"Failed to get disk usage for {path}: {e}")
        raise


def get_free_space(path: str | Path) -> int:
    """Get available disk space in bytes.

    Args:
        path: Path to check

    Returns:
        Free space in bytes

    Raises:
        OSError: If free space cannot be determined
    """
    usage = get_disk_usage(path)
    return usage['free']


def get_recommended_temp_dir(repo_root: Path) -> Path:
    """Get recommended temporary directory, preferring locations with more space.

    Args:
        repo_root: Repository root directory

    Returns:
        Path to recommended temporary directory

    Notes:
        Tries in order:
        1. /tmp if it has >1GB free space
        2. repo_root/.tmp if repo has >5GB free space
        3. Falls back to /tmp regardless of space
    """
    tmp_dir = Path("/tmp")
    repo_tmp = repo_root / ".tmp"

    try:
        # Check /tmp first
        if tmp_dir.exists():
            free_space = get_free_space(tmp_dir)
            if free_space > 1_000_000_000:  # >1GB
                return tmp_dir

        # Check repo .tmp
        repo_tmp.mkdir(parents=True, exist_ok=True)
        free_space = get_free_space(repo_tmp)
        if free_space > 5_000_000_000:  # >5GB
            return repo_tmp

    except (OSError, Exception) as e:
        logger.warning(f"Error checking disk space: {e}")

    # Fall back to /tmp
    return tmp_dir


def ensure_disk_space(
    path: str | Path,
    required_bytes: int,
    cleanup_threshold: float = 0.9,
) -> bool:
    """Ensure adequate disk space is available, with optional cleanup.

    Args:
        path: Path to check
        required_bytes: Minimum required free space in bytes
        cleanup_threshold: Usage threshold (0-1) above which to attempt cleanup

    Returns:
        True if adequate space is available

    Raises:
        OSError: If disk space cannot be checked
        RuntimeError: If cleanup fails to free enough space
    """
    usage = get_disk_usage(path)

    if usage['free'] >= required_bytes:
        return True

    logger.warning(
        f"Insufficient disk space: {usage['free']} bytes free, "
        f"{required_bytes} bytes required"
    )

    # Attempt cleanup if usage is high
    if usage['percent'] / 100 >= cleanup_threshold:
        logger.info("Attempting automatic cleanup of temp files")
        cleaned = cleanup_temp_files(path, max_age_hours=24)

        if cleaned > 0:
            logger.info(f"Cleaned up {cleaned} temp files")

            # Recheck space after cleanup
            usage = get_disk_usage(path)
            if usage['free'] >= required_bytes:
                return True

    if usage['free'] < required_bytes:
        raise RuntimeError(
            f"Insufficient disk space: {usage['free']} bytes free, "
            f"{required_bytes} bytes required. "
            f"Usage: {usage['percent']}%"
        )

    return True


def cleanup_temp_files(
    directory: str | Path,
    max_age_hours: int = 24,
    patterns: Optional[List[str]] = None,
) -> int:
    """Clean up old temporary files.

    Args:
        directory: Directory to clean
        max_age_hours: Maximum age of files to keep (hours)
        patterns: File patterns to clean (default: common temp patterns)

    Returns:
        Number of files removed
    """
    if patterns is None:
        patterns = ["*.tmp", "*.temp", "*.log", "*.cache", "*.bak"]

    directory = Path(directory)
    if not directory.exists():
        return 0

    cutoff_time = time.time() - (max_age_hours * 3600)
    removed_count = 0

    try:
        for pattern in patterns:
            for file_path in directory.glob(pattern):
                if file_path.is_file():
                    try:
                        # Check file age
                        stat = file_path.stat()
                        if stat.st_mtime < cutoff_time:
                            file_path.unlink()
                            removed_count += 1

                    except OSError as e:
                        logger.warning(f"Failed to remove temp file {file_path}: {e}")

        # Also clean empty directories
        for dir_path in directory.rglob("*"):
            if dir_path.is_dir() and not any(dir_path.iterdir()):
                try:
                    dir_path.rmdir()
                    logger.debug(f"Removed empty directory {dir_path}")
                except OSError:
                    pass  # Directory might not be empty anymore

    except Exception as e:
        logger.error(f"Error during temp file cleanup: {e}")

    return removed_count


def cleanup_old_files(
    directory: str | Path,
    max_age_days: int = 30,
    exclude_patterns: Optional[List[str]] = None,
) -> int:
    """Clean up old files from directory.

    Args:
        directory: Directory to clean
        max_age_days: Maximum age of files to keep (days)
        exclude_patterns: Patterns to exclude from cleanup

    Returns:
        Number of files removed
    """
    directory = Path(directory)
    if not directory.exists():
        return 0

    if exclude_patterns is None:
        exclude_patterns = []

    cutoff_time = time.time() - (max_age_days * 24 * 3600)
    removed_count = 0

    try:
        for file_path in directory.rglob("*"):
            if file_path.is_file():
                # Check if file matches exclude patterns
                excluded = False
                for pattern in exclude_patterns:
                    if file_path.match(pattern):
                        excluded = True
                        break

                if not excluded:
                    try:
                        stat = file_path.stat()
                        if stat.st_mtime < cutoff_time:
                            file_path.unlink()
                            removed_count += 1
                            logger.debug(f"Removed old file: {file_path}")

                    except OSError as e:
                        logger.warning(f"Failed to remove old file {file_path}: {e}")

    except Exception as e:
        logger.error(f"Error during old file cleanup: {e}")

    return removed_count


def get_directory_size(path: str | Path) -> int:
    """Get total size of directory in bytes.

    Args:
        path: Directory path

    Returns:
        Total size in bytes

    Raises:
        OSError: If directory cannot be accessed
    """
    path = Path(path)
    if not path.exists():
        return 0

    total_size = 0

    try:
        for file_path in path.rglob("*"):
            if file_path.is_file():
                try:
                    total_size += file_path.stat().st_size
                except OSError:
                    pass  # Skip files we can't access

    except OSError as e:
        logger.error(f"Failed to calculate directory size for {path}: {e}")
        raise

    return total_size


def get_largest_files(
    directory: str | Path,
    n: int = 10,
) -> List[Tuple[Path, int]]:
    """Find largest files in directory.

    Args:
        directory: Directory to scan
        n: Number of largest files to return

    Returns:
        List of (path, size) tuples for largest files
    """
    directory = Path(directory)
    if not directory.exists():
        return []

    files_with_sizes = []

    try:
        for file_path in directory.rglob("*"):
            if file_path.is_file():
                try:
                    size = file_path.stat().st_size
                    files_with_sizes.append((file_path, size))
                except OSError:
                    pass

        # Sort by size descending and return top n
        files_with_sizes.sort(key=lambda x: x[1], reverse=True)
        return files_with_sizes[:n]

    except Exception as e:
        logger.error(f"Error finding largest files in {directory}: {e}")
        return []


def monitor_disk_space(
    path: str | Path,
    warning_threshold: float = 0.8,
    critical_threshold: float = 0.95,
) -> Dict[str, Any]:
    """Monitor disk space and return status.

    Args:
        path: Path to monitor
        warning_threshold: Usage threshold for warning (0-1)
        critical_threshold: Usage threshold for critical (0-1)

    Returns:
        Dictionary with monitoring results:
        - status: 'ok', 'warning', or 'critical'
        - usage: Disk usage statistics
        - message: Status message
    """
    try:
        usage = get_disk_usage(path)
        percent_used = usage['percent'] / 100

        if percent_used >= critical_threshold:
            status = 'critical'
            message = f"Critical disk usage: {usage['percent']}% used"
        elif percent_used >= warning_threshold:
            status = 'warning'
            message = f"High disk usage: {usage['percent']}% used"
        else:
            status = 'ok'
            message = f"Disk usage normal: {usage['percent']}% used"

        return {
            'status': status,
            'usage': usage,
            'message': message,
        }

    except Exception as e:
        logger.error(f"Disk monitoring failed for {path}: {e}")
        return {
            'status': 'error',
            'usage': None,
            'message': f"Monitoring failed: {e}",
        }


def safe_remove_directory(
    directory: str | Path,
    confirm_size_mb: Optional[int] = None,
) -> bool:
    """Safely remove directory with confirmation.

    Args:
        directory: Directory to remove
        confirm_size_mb: If provided, only remove if directory size is below this limit

    Returns:
        True if directory was removed
    """
    directory = Path(directory)

    if not directory.exists():
        return True

    try:
        # Check size if confirmation requested
        if confirm_size_mb is not None:
            size_mb = get_directory_size(directory) / (1024 * 1024)
            if size_mb > confirm_size_mb:
                logger.warning(
                    f"Directory {directory} is {size_mb:.1f} MB, "
                    f"larger than limit {confirm_size_mb} MB. Not removing."
                )
                return False

        shutil.rmtree(directory)
        logger.info(f"Removed directory: {directory}")
        return True

    except Exception as e:
        logger.error(f"Failed to remove directory {directory}: {e}")
        return False





