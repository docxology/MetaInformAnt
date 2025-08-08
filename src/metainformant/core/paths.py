from __future__ import annotations

from pathlib import Path


def expand_and_resolve(path: str | Path) -> Path:
    """Expand user (~) and resolve to absolute path without requiring existence."""
    return Path(path).expanduser().resolve(strict=False)


def is_within(path: str | Path, parent: str | Path) -> bool:
    """Return True if path is inside parent directory (after resolving)."""
    p = expand_and_resolve(path)
    r = expand_and_resolve(parent)
    try:
        p.relative_to(r)
        return True
    except ValueError:
        return False


