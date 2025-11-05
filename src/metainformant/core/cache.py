from __future__ import annotations

import threading
import time
from collections import OrderedDict
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


# Thread-safe cache operations
_cache_locks: dict[str, threading.Lock] = {}
_cache_lock = threading.Lock()


def _get_cache_lock(key: str) -> threading.Lock:
    """Get or create a lock for a cache key."""
    with _cache_lock:
        if key not in _cache_locks:
            _cache_locks[key] = threading.Lock()
        return _cache_locks[key]


def cache_json(cache_dir: Path, key: str, obj: Any) -> Path:
    """Cache an object as JSON with thread-safe access.
    
    Args:
        cache_dir: Cache directory
        key: Cache key
        obj: Object to cache
        
    Returns:
        Path to cached file
    """
    path = _key_to_path(cache_dir, key)
    lock = _get_cache_lock(str(path))
    with lock:
        core_io.ensure_directory(path.parent)
        core_io.dump_json(obj, path, atomic=True)
    return path


def load_cached_json(cache_dir: Path, key: str, *, ttl_seconds: int) -> Any | None:
    """Load cached JSON with thread-safe access.
    
    Args:
        cache_dir: Cache directory
        key: Cache key
        ttl_seconds: Time to live in seconds
        
    Returns:
        Cached object or None if not found or expired
    """
    path = _key_to_path(cache_dir, key)
    lock = _get_cache_lock(str(path))
    with lock:
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


def get_cache_info(cache_dir: Path) -> dict[str, Any]:
    """Get information about cache directory contents.

    Args:
        cache_dir: Path to cache directory

    Returns:
        Dictionary with cache statistics including hit/miss ratios if available
    """
    if not cache_dir.exists():
        return {"exists": False, "total_files": 0, "total_size": 0, "hits": 0, "misses": 0}

    total_files = 0
    total_size = 0
    oldest_file_time = None
    newest_file_time = None

    for file in cache_dir.rglob("*"):
        if file.is_file() and file.suffix in {".json", ".gz"}:
            total_files += 1
            file_size = file.stat().st_size
            total_size += file_size
            mtime = file.stat().st_mtime
            if oldest_file_time is None or mtime < oldest_file_time:
                oldest_file_time = mtime
            if newest_file_time is None or mtime > newest_file_time:
                newest_file_time = mtime

    result = {
        "exists": True,
        "total_files": total_files,
        "total_size": total_size,
        "directory": str(cache_dir)
    }
    
    if oldest_file_time is not None:
        result["oldest_file_age"] = time.time() - oldest_file_time
    if newest_file_time is not None:
        result["newest_file_age"] = time.time() - newest_file_time

    return result


class LRUCache:
    """Thread-safe LRU cache with TTL support.
    
    Attributes:
        max_size: Maximum number of items in cache
        ttl_seconds: Time to live in seconds (0 for no expiration)
    """
    
    def __init__(self, max_size: int = 100, ttl_seconds: int = 0):
        """Initialize LRU cache.
        
        Args:
            max_size: Maximum number of items
            ttl_seconds: Time to live in seconds (0 for no expiration)
        """
        self.max_size = max_size
        self.ttl_seconds = ttl_seconds
        self._cache: OrderedDict[str, tuple[Any, float]] = OrderedDict()
        self._lock = threading.Lock()
        self._hits = 0
        self._misses = 0
    
    def get(self, key: str) -> Any | None:
        """Get item from cache.
        
        Args:
            key: Cache key
            
        Returns:
            Cached value or None if not found or expired
        """
        with self._lock:
            if key not in self._cache:
                self._misses += 1
                return None
            
            value, timestamp = self._cache[key]
            
            # Check TTL
            if self.ttl_seconds > 0 and (time.time() - timestamp) > self.ttl_seconds:
                del self._cache[key]
                self._misses += 1
                return None
            
            # Move to end (most recently used)
            self._cache.move_to_end(key)
            self._hits += 1
            return value
    
    def set(self, key: str, value: Any) -> None:
        """Set item in cache.
        
        Args:
            key: Cache key
            value: Value to cache
        """
        with self._lock:
            if key in self._cache:
                # Update existing item
                self._cache[key] = (value, time.time())
                self._cache.move_to_end(key)
            else:
                # Add new item
                if len(self._cache) >= self.max_size:
                    # Remove oldest item
                    self._cache.popitem(last=False)
                self._cache[key] = (value, time.time())
    
    def clear(self) -> None:
        """Clear all items from cache."""
        with self._lock:
            self._cache.clear()
            self._hits = 0
            self._misses = 0
    
    def get_stats(self) -> dict[str, Any]:
        """Get cache statistics.
        
        Returns:
            Dictionary with cache statistics
        """
        with self._lock:
            total = self._hits + self._misses
            hit_rate = (self._hits / total) if total > 0 else 0.0
            return {
                "size": len(self._cache),
                "max_size": self.max_size,
                "hits": self._hits,
                "misses": self._misses,
                "hit_rate": hit_rate,
                "ttl_seconds": self.ttl_seconds
            }
