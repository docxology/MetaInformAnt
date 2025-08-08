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


