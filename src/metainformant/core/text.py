from __future__ import annotations

import re
from pathlib import Path

_WHITESPACE_RE = re.compile(r"\s+")
_SLUG_INVALID_RE = re.compile(r"[^a-z0-9-]+")


def normalize_whitespace(s: str) -> str:
    """Collapse all whitespace to single spaces and strip ends."""
    return _WHITESPACE_RE.sub(" ", s).strip()


def slugify(s: str) -> str:
    """Lowercase, normalize spaces to dashes, and drop invalid URL chars."""
    s = normalize_whitespace(s).lower().replace(" ", "-")
    s = _SLUG_INVALID_RE.sub("", s)
    s = re.sub(r"-+", "-", s)
    return s.strip("-")


def safe_filename(name: str) -> str:
    """Make a safe filename, preserving extension when present."""
    p = Path(name)
    stem = slugify(p.stem)
    suffix = p.suffix
    return f"{stem}{suffix}"


