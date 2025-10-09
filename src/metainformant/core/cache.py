from __future__ import annotations

import time
from pathlib import Path
from typing import Any

from . import io as core_io
from .text import safe_filename


def _key_to_path(cache_dir: Path, key: str) -> Path:
    # keep slashes as hierarchy but sanitize each segment
    parts = [safe_filename(seg) for seg in key.split("/") if seg]
    path = cache_dir.joinpath(*parts)
    if path.suffix not in {".json", ".gz"}:
        path = path.with_suffix(".json")
    return path


def cache_json(cache_dir: Path, key: str, obj: Any) -> Path:
    path = _key_to_path(cache_dir, key)
    core_io.ensure_directory(path.parent)
    core_io.dump_json(obj, path)
    return path


def load_cached_json(cache_dir: Path, key: str, *, ttl_seconds: int) -> Any | None:
    path = _key_to_path(cache_dir, key)
    if not path.exists():
        return None
    age = time.time() - path.stat().st_mtime
    if ttl_seconds <= 0 or age > ttl_seconds:
        return None
    return core_io.load_json(path)


def get_json_cache(cache_file: Path, default: Any = None, max_age_seconds: float | None = None) -> Any:
    """Get data from JSON cache file with optional TTL check.

    Args:
        cache_file: Path to cache file
        default: Default value if cache miss or expired
        max_age_seconds: Maximum age in seconds, None for no TTL check

    Returns:
        Cached data or default value
    """
    if not cache_file.exists():
        return default

    if max_age_seconds is not None:
        age = time.time() - cache_file.stat().st_mtime
        if age > max_age_seconds:
            return default

    try:
        return core_io.load_json(cache_file)
    except Exception:
        return default


def set_json_cache(cache_file: Path, data: Any) -> None:
    """Set data in JSON cache file.

    Args:
        cache_file: Path to cache file
        data: Data to cache
    """
    core_io.ensure_directory(cache_file.parent)
    core_io.dump_json(data, cache_file)


def clear_cache_dir(cache_dir: Path) -> None:
    """Clear all cache files in a directory."""
    if cache_dir.exists():
        for file in cache_dir.rglob("*"):
            if file.is_file() and file.suffix in {".json", ".gz"}:
                file.unlink()


def get_cache_info(cache_dir: Path) -> Dict[str, Any]:
    """Get information about cache directory contents.

    Args:
        cache_dir: Path to cache directory

    Returns:
        Dictionary with cache statistics
    """
    if not cache_dir.exists():
        return {"exists": False, "total_files": 0, "total_size": 0}

    total_files = 0
    total_size = 0

    for file in cache_dir.rglob("*"):
        if file.is_file() and file.suffix in {".json", ".gz"}:
            total_files += 1
            total_size += file.stat().st_size

    return {
        "exists": True,
        "total_files": total_files,
        "total_size": total_size,
        "directory": str(cache_dir)
    }
