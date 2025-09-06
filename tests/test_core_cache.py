from __future__ import annotations

from pathlib import Path

from metainformant.core import cache as core_cache


def test_json_cache_roundtrip(tmp_path: Path) -> None:
    cache_dir = tmp_path / "cache"
    key = "example/key"
    data = {"v": 1}
    core_cache.cache_json(cache_dir, key, data)
    loaded = core_cache.load_cached_json(cache_dir, key, ttl_seconds=3600)
    assert loaded == data


def test_json_cache_ttl_expired(tmp_path: Path) -> None:
    cache_dir = tmp_path / "cache"
    key = "k"
    data = {"x": 2}
    core_cache.cache_json(cache_dir, key, data)
    # ttl 0 => always expired
    loaded = core_cache.load_cached_json(cache_dir, key, ttl_seconds=0)
    assert loaded is None
