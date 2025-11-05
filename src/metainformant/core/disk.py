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

