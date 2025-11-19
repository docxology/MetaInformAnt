"""Filesystem detection and compatibility utilities.

Provides functions to detect filesystem type, check symlink support,
and determine appropriate cache directories for tools like uv.
"""

from __future__ import annotations

import os
import platform
import subprocess
import tempfile
from pathlib import Path
from typing import Optional


def detect_filesystem_type(path: Path | str) -> str:
    """Detect the filesystem type for a given path.
    
    Uses platform-specific methods to determine filesystem type:
    - Linux: Uses `df -T` command
    - macOS: Uses `df -T` command
    - Windows: Uses `fsutil` command
    
    Args:
        path: Path to check (file or directory)
        
    Returns:
        Filesystem type string (e.g., 'ext4', 'exfat', 'ntfs', 'apfs', 'fat32')
        Returns 'unknown' if detection fails
    """
    path_str = str(path)
    
    # Resolve to absolute path
    try:
        path_obj = Path(path_str).resolve()
        if path_obj.is_file():
            path_str = str(path_obj.parent)
        else:
            path_str = str(path_obj)
    except Exception:
        pass
    
    system = platform.system().lower()
    
    if system == "linux":
        try:
            # Use df -T to get filesystem type
            result = subprocess.run(
                ["df", "-T", path_str],
                capture_output=True,
                text=True,
                timeout=5,
            )
            if result.returncode == 0:
                lines = result.stdout.strip().split("\n")
                if len(lines) > 1:
                    # Second line contains filesystem info
                    parts = lines[1].split()
                    if len(parts) >= 2:
                        fs_type = parts[1].lower()
                        return fs_type
        except Exception:
            pass
    
    elif system == "darwin":  # macOS
        try:
            # macOS df -T works similarly
            result = subprocess.run(
                ["df", "-T", path_str],
                capture_output=True,
                text=True,
                timeout=5,
            )
            if result.returncode == 0:
                lines = result.stdout.strip().split("\n")
                if len(lines) > 1:
                    parts = lines[1].split()
                    if len(parts) >= 2:
                        fs_type = parts[1].lower()
                        return fs_type
        except Exception:
            pass
    
    elif system == "windows":
        try:
            # Use fsutil to get filesystem type
            drive = path_str[0:2] if len(path_str) >= 2 else "C:"
            result = subprocess.run(
                ["fsutil", "fsinfo", "volumeinfo", drive],
                capture_output=True,
                text=True,
                timeout=5,
            )
            if result.returncode == 0:
                for line in result.stdout.split("\n"):
                    if "File System Name" in line or "FileSystem" in line:
                        parts = line.split(":")
                        if len(parts) > 1:
                            fs_type = parts[1].strip().lower()
                            return fs_type
        except Exception:
            pass
    
    return "unknown"


def supports_symlinks(path: Path | str) -> bool:
    """Check if the filesystem at the given path supports symlinks.
    
    Uses both filesystem type detection and a practical test to determine
    symlink support.
    
    Args:
        path: Path to check (file or directory)
        
    Returns:
        True if filesystem supports symlinks, False otherwise
    """
    # Known filesystems that don't support symlinks
    no_symlink_fs = {"exfat", "fat32", "fat", "vfat", "msdos"}
    
    fs_type = detect_filesystem_type(path)
    
    if fs_type in no_symlink_fs:
        return False
    
    # Known filesystems that support symlinks
    symlink_fs = {
        "ext4", "ext3", "ext2", "xfs", "btrfs", "zfs",
        "ntfs", "apfs", "hfs+", "hfs",
    }
    
    if fs_type in symlink_fs:
        return True
    
    # For unknown filesystems, try a practical test
    try:
        path_obj = Path(path)
        if path_obj.is_file():
            test_dir = path_obj.parent
        else:
            test_dir = path_obj
        
        # Create a temporary file
        test_file = test_dir / ".symlink_test_temp"
        test_file.write_text("test")
        
        # Try to create a symlink
        test_link = test_dir / ".symlink_test_link"
        try:
            test_link.symlink_to(test_file)
            # Success - symlinks are supported
            test_link.unlink()
            test_file.unlink()
            return True
        except (OSError, NotImplementedError):
            # Failed - symlinks not supported
            test_file.unlink()
            return False
    except Exception:
        # If test fails, assume no symlink support to be safe
        return False


def get_uv_cache_dir(repo_root: Path | str | None = None) -> Path:
    """Get the appropriate UV cache directory based on filesystem support.
    
    Returns `/tmp/uv-cache` for filesystems without symlink support (FAT/exFAT),
    otherwise returns `.uv-cache` relative to repo root or current directory.
    
    Args:
        repo_root: Root directory of the repository (optional)
        
    Returns:
        Path to UV cache directory
    """
    # Check if UV_CACHE_DIR is explicitly set
    uv_cache_env = os.environ.get("UV_CACHE_DIR")
    if uv_cache_env:
        cache_path = Path(uv_cache_env)
        cache_path.mkdir(parents=True, exist_ok=True)
        return cache_path
    
    # Determine base path for cache directory
    if repo_root:
        base_path = Path(repo_root)
    else:
        # Try to find repo root by looking for pyproject.toml
        current = Path.cwd()
        for parent in [current] + list(current.parents):
            if (parent / "pyproject.toml").exists():
                base_path = parent
                break
        else:
            base_path = current
    
    # Check if filesystem supports symlinks
    if not supports_symlinks(base_path):
        # Use /tmp/uv-cache for filesystems without symlink support
        cache_path = Path("/tmp/uv-cache")
        cache_path.mkdir(parents=True, exist_ok=True)
        return cache_path
    
    # Use .uv-cache in repo root for filesystems with symlink support
    cache_path = base_path / ".uv-cache"
    cache_path.mkdir(parents=True, exist_ok=True)
    return cache_path


def get_venv_location(repo_root: Path | str | None = None) -> Path:
    """Get the appropriate virtual environment location based on filesystem support.
    
    Returns `/tmp/metainformant_venv` for filesystems without symlink support,
    otherwise returns `.venv` relative to repo root.
    
    Args:
        repo_root: Root directory of the repository (optional)
        
    Returns:
        Path to virtual environment directory
    """
    # Check if UV_PROJECT_ENVIRONMENT is explicitly set (UV's preferred env var)
    uv_venv_env = os.environ.get("UV_PROJECT_ENVIRONMENT")
    if uv_venv_env:
        return Path(uv_venv_env)
    
    # Check if METAINFORMANT_VENV is explicitly set (for compatibility)
    venv_env = os.environ.get("METAINFORMANT_VENV")
    if venv_env:
        return Path(venv_env)
    
    # Determine base path
    if repo_root:
        base_path = Path(repo_root)
    else:
        # Try to find repo root
        current = Path.cwd()
        for parent in [current] + list(current.parents):
            if (parent / "pyproject.toml").exists():
                base_path = parent
                break
        else:
            base_path = current
    
    # Check if filesystem supports symlinks
    if not supports_symlinks(base_path):
        # Use /tmp/metainformant_venv for filesystems without symlink support
        return Path("/tmp/metainformant_venv")
    
    # Use .venv in repo root for filesystems with symlink support
    return base_path / ".venv"


def is_fat_filesystem(path: Path | str) -> bool:
    """Check if the path is on a FAT filesystem (FAT32, exFAT, etc.).
    
    Args:
        path: Path to check
        
    Returns:
        True if filesystem is FAT-based, False otherwise
    """
    fs_type = detect_filesystem_type(path)
    fat_types = {"exfat", "fat32", "fat", "vfat", "msdos"}
    return fs_type in fat_types

