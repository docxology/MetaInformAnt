"""Disk space monitoring and management utilities."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Optional


def get_disk_usage(path: Path) -> tuple[float, float, float, str]:
    """Get disk usage statistics for a path.
    
    Args:
        path: Path to check (can be file or directory)
        
    Returns:
        Tuple of (total_gb, used_gb, free_gb, percent_used_str)
        Returns (0, 0, 0, "0%") if unable to determine
    """
    try:
        # Use shutil for cross-platform compatibility
        stat = shutil.disk_usage(path)
        total_gb = stat.total / (1024 ** 3)
        used_gb = stat.used / (1024 ** 3)
        free_gb = stat.free / (1024 ** 3)
        percent = (stat.used / stat.total) * 100
        percent_str = f"{percent:.1f}%"
        
        return total_gb, used_gb, free_gb, percent_str
    except Exception:
        # Fallback to df command on Unix systems
        try:
            result = subprocess.run(
                ["df", "-h", str(path)],
                capture_output=True,
                text=True,
                timeout=5,
            )
            if result.returncode == 0:
                lines = result.stdout.strip().split("\n")
                if len(lines) > 1:
                    parts = lines[1].split()
                    if len(parts) >= 5:
                        # Parse output: Filesystem Size Used Avail Use% Mounted
                        size_str = parts[1]
                        used_str = parts[2]
                        avail_str = parts[3]
                        percent_str = parts[4]
                        
                        # Convert to GB
                        def parse_size(s: str) -> float:
                            """Parse size string like '468G' or '1.5T' to GB."""
                            s = s.upper().replace(",", "")
                            if s.endswith("K"):
                                return float(s[:-1]) / (1024 * 1024)
                            elif s.endswith("M"):
                                return float(s[:-1]) / 1024
                            elif s.endswith("G"):
                                return float(s[:-1])
                            elif s.endswith("T"):
                                return float(s[:-1]) * 1024
                            else:
                                return float(s) / (1024 ** 3)
                        
                        total_gb = parse_size(size_str)
                        used_gb = parse_size(used_str)
                        free_gb = parse_size(avail_str)
                        
                        return total_gb, used_gb, free_gb, percent_str
        except Exception:
            pass
        
        return 0.0, 0.0, 0.0, "0%"


def check_disk_space(path: Path, min_free_gb: float = 10.0, min_free_percent: float = 5.0) -> tuple[bool, str]:
    """Check if sufficient disk space is available.
    
    Args:
        path: Path to check disk space for
        min_free_gb: Minimum free space required in GB
        min_free_percent: Minimum free space required as percentage
        
    Returns:
        Tuple of (is_sufficient: bool, message: str)
    """
    total_gb, used_gb, free_gb, percent_str = get_disk_usage(path)
    
    if total_gb == 0:
        return False, "Unable to determine disk space"
    
    free_percent = (free_gb / total_gb) * 100
    
    # Check both absolute and percentage thresholds
    if free_gb < min_free_gb:
        return False, (
            f"Insufficient disk space: {free_gb:.1f}GB free "
            f"(need at least {min_free_gb}GB). "
            f"Usage: {percent_str} ({used_gb:.1f}GB / {total_gb:.1f}GB)"
        )
    
    if free_percent < min_free_percent:
        return False, (
            f"Insufficient disk space: {free_percent:.1f}% free "
            f"(need at least {min_free_percent}%). "
            f"Usage: {percent_str} ({used_gb:.1f}GB / {total_gb:.1f}GB)"
        )
    
    return True, (
        f"Disk space OK: {free_gb:.1f}GB free ({free_percent:.1f}%). "
        f"Usage: {percent_str} ({used_gb:.1f}GB / {total_gb:.1f}GB)"
    )


def get_disk_space_info(path: Path) -> dict[str, float | str]:
    """Get comprehensive disk space information.
    
    Args:
        path: Path to check
        
    Returns:
        Dictionary with disk space information
    """
    total_gb, used_gb, free_gb, percent_str = get_disk_usage(path)
    
    return {
        "total_gb": total_gb,
        "used_gb": used_gb,
        "free_gb": free_gb,
        "percent_used": percent_str,
        "percent_free": f"{100 - float(percent_str.rstrip('%')):.1f}%",
    }


def detect_drive_size_category(path: Path) -> str:
    """Detect drive size category based on available free space.
    
    Args:
        path: Path to check (directory or file on the drive)
        
    Returns:
        "large" (> 1TB free), "medium" (500GB-1TB free), or "small" (< 500GB free)
    """
    total_gb, used_gb, free_gb, _ = get_disk_usage(path)
    
    if total_gb == 0:
        # Can't determine, assume small to be safe
        return "small"
    
    if free_gb > 1024:  # > 1TB free
        return "large"
    elif free_gb > 500:  # 500GB-1TB free
        return "medium"
    else:  # < 500GB free
        return "small"


def get_recommended_batch_size(path: Path, sample_size_gb: float = 1.5, safety_buffer: float = 0.3) -> int:
    """Calculate recommended batch size based on available disk space.
    
    Args:
        path: Path to check (directory or file on the drive)
        sample_size_gb: Estimated size per sample in GB (default: 1.5GB)
        safety_buffer: Safety buffer as fraction (default: 0.3 = 30% buffer)
        
    Returns:
        Recommended batch size (number of samples)
    """
    total_gb, used_gb, free_gb, _ = get_disk_usage(path)
    
    if total_gb == 0 or free_gb < 10:
        # Can't determine or very low space, return conservative default
        return 8
    
    # Calculate available space after buffer
    available_gb = free_gb * (1 - safety_buffer)
    
    # Calculate batch size (round down to be safe)
    batch_size = int(available_gb / sample_size_gb)
    
    # Apply limits based on drive category
    category = detect_drive_size_category(path)
    if category == "large":
        # Large drives: 50-100 samples
        return min(max(batch_size, 50), 100)
    elif category == "medium":
        # Medium drives: 20-30 samples
        return min(max(batch_size, 20), 30)
    else:
        # Small drives: 8-12 samples
        return min(max(batch_size, 8), 12)


def get_recommended_temp_dir(repo_root: Path) -> Path:
    """Get recommended temporary directory location based on drive size.
    
    Prefers external drive if it's large enough, otherwise falls back to system temp.
    
    Args:
        repo_root: Root directory of the repository
        
    Returns:
        Path to recommended temporary directory (directory is created if needed)
    """
    import os
    import tempfile
    
    # Check if TMPDIR is explicitly set
    tmpdir_env = os.environ.get("TMPDIR")
    if tmpdir_env:
        tmp_path = Path(tmpdir_env)
        tmp_path.mkdir(parents=True, exist_ok=True)
        return tmp_path
    
    # Check output directory location (usually on external drive)
    output_dir = repo_root / "output"
    if output_dir.exists():
        category = detect_drive_size_category(output_dir)
        if category in ("large", "medium"):
            # Use external drive for temp files
            temp_dir = output_dir / ".tmp"
            temp_dir.mkdir(parents=True, exist_ok=True)
            return temp_dir
    
    # Fall back to system temp
    temp_dir = Path(tempfile.gettempdir())
    temp_dir.mkdir(parents=True, exist_ok=True)
    return temp_dir

