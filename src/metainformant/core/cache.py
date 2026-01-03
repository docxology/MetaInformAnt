"""JSON-based caching with TTL support for METAINFORMANT.

This module provides thread-safe JSON caching with automatic cleanup of expired entries.
"""

from __future__ import annotations

import json
import os
import threading
import time
from pathlib import Path
from typing import Any, Dict, Optional, Union

from metainformant.core import io
from metainformant.core.logging import get_logger

logger = get_logger(__name__)


class CacheEntry:
    """Represents a single cache entry with TTL."""

    def __init__(self, value: Any, ttl_seconds: int | None = None):
        """Initialize cache entry.

        Args:
            value: Value to cache
            ttl_seconds: Time-to-live in seconds, None for no expiration
        """
        self.value = value
        self.created_at = time.time()
        self.ttl_seconds = ttl_seconds
        self.expires_at = self.created_at + ttl_seconds if ttl_seconds else None

    def is_expired(self) -> bool:
        """Check if entry is expired."""
        if self.expires_at is None:
            return False
        return time.time() > self.expires_at

    def time_to_expiry(self) -> float | None:
        """Get seconds until expiry, None if no expiration."""
        if self.expires_at is None:
            return None
        return max(0, self.expires_at - time.time())

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            'value': self.value,
            'created_at': self.created_at,
            'ttl_seconds': self.ttl_seconds,
            'expires_at': self.expires_at,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> CacheEntry:
        """Create from dictionary."""
        entry = cls.__new__(cls)
        entry.value = data['value']
        entry.created_at = data['created_at']
        entry.ttl_seconds = data['ttl_seconds']
        entry.expires_at = data['expires_at']
        return entry


class JsonCache:
    """Thread-safe JSON-based cache with TTL support.

    This cache stores data as JSON files on disk with automatic expiration
    and cleanup of stale entries. All operations are thread-safe.

    Example:
        cache = JsonCache("cache_dir", ttl_seconds=3600)

        # Store data
        cache.set("key", {"data": "value"})

        # Retrieve data
        data = cache.get("key")

        # Clear expired entries
        cache.cleanup_expired()
    """

    def __init__(self, cache_dir: str | Path, ttl_seconds: int = 3600):
        """Initialize JSON cache.

        Args:
            cache_dir: Directory to store cache files
            ttl_seconds: Default time-to-live for cache entries
        """
        self.cache_dir = Path(cache_dir)
        self.default_ttl = ttl_seconds
        self._lock = threading.RLock()
        self._ensure_cache_dir()

    def _ensure_cache_dir(self) -> None:
        """Ensure cache directory exists."""
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _get_cache_path(self, key: str) -> Path:
        """Get filesystem path for cache key."""
        # Sanitize key for filesystem
        safe_key = "".join(c for c in key if c.isalnum() or c in "._-")
        if not safe_key:
            safe_key = "default"
        return self.cache_dir / f"{safe_key}.json"

    def get(self, key: str, default: Any = None) -> Any:
        """Retrieve value from cache.

        Args:
            key: Cache key
            default: Default value if key not found or expired

        Returns:
            Cached value or default
        """
        with self._lock:
            cache_path = self._get_cache_path(key)

            if not cache_path.exists():
                return default

            try:
                data = io.load_json(cache_path)
                entry = CacheEntry.from_dict(data)

                if entry.is_expired():
                    # Remove expired entry
                    cache_path.unlink(missing_ok=True)
                    return default

                return entry.value

            except (json.JSONDecodeError, KeyError, ValueError) as e:
                logger.warning(f"Failed to load cache entry {key}: {e}")
                # Remove corrupted entry
                cache_path.unlink(missing_ok=True)
                return default

    def set(self, key: str, value: Any, ttl_seconds: int | None = None) -> None:
        """Store value in cache.

        Args:
            key: Cache key
            value: Value to cache (must be JSON serializable)
            ttl_seconds: Time-to-live override, uses default if None
        """
        with self._lock:
            ttl = ttl_seconds if ttl_seconds is not None else self.default_ttl
            entry = CacheEntry(value, ttl)

            cache_path = self._get_cache_path(key)

            try:
                io.dump_json(entry.to_dict(), cache_path)
            except Exception as e:
                logger.error(f"Failed to save cache entry {key}: {e}")

    def delete(self, key: str) -> bool:
        """Delete cache entry.

        Args:
            key: Cache key

        Returns:
            True if entry was deleted, False if not found
        """
        with self._lock:
            cache_path = self._get_cache_path(key)
            if cache_path.exists():
                cache_path.unlink()
                return True
            return False

    def clear(self) -> int:
        """Clear all cache entries.

        Returns:
            Number of entries cleared
        """
        with self._lock:
            count = 0
            for cache_file in self.cache_dir.glob("*.json"):
                try:
                    cache_file.unlink()
                    count += 1
                except OSError as e:
                    logger.warning(f"Failed to remove cache file {cache_file}: {e}")

            return count

    def cleanup_expired(self) -> int:
        """Remove expired cache entries.

        Returns:
            Number of expired entries removed
        """
        with self._lock:
            count = 0
            for cache_file in self.cache_dir.glob("*.json"):
                try:
                    data = io.load_json(cache_file)
                    entry = CacheEntry.from_dict(data)

                    if entry.is_expired():
                        cache_file.unlink()
                        count += 1

                except (json.JSONDecodeError, KeyError, ValueError, OSError) as e:
                    # Remove corrupted files
                    try:
                        cache_file.unlink()
                        count += 1
                    except OSError:
                        pass
                    logger.warning(f"Removed corrupted cache file {cache_file}: {e}")

            return count

    def size(self) -> int:
        """Get number of cache entries.

        Returns:
            Number of cache files
        """
        with self._lock:
            return len(list(self.cache_dir.glob("*.json")))

    def stats(self) -> Dict[str, Any]:
        """Get cache statistics.

        Returns:
            Dictionary with cache statistics
        """
        with self._lock:
            total_files = 0
            valid_entries = 0
            expired_entries = 0
            corrupted_files = 0

            for cache_file in self.cache_dir.glob("*.json"):
                total_files += 1
                try:
                    data = io.load_json(cache_file)
                    entry = CacheEntry.from_dict(data)

                    if entry.is_expired():
                        expired_entries += 1
                    else:
                        valid_entries += 1

                except (json.JSONDecodeError, KeyError, ValueError):
                    corrupted_files += 1

            return {
                'total_files': total_files,
                'valid_entries': valid_entries,
                'expired_entries': expired_entries,
                'corrupted_files': corrupted_files,
                'cache_dir': str(self.cache_dir),
                'default_ttl': self.default_ttl,
            }





