"""Explicit lossless FASTQ compression settings for bounded acquisition workers."""

from __future__ import annotations

import os
from pathlib import Path


def compression_level() -> int:
    """Read a strict pigz level; retain level 6 unless an operator selects another."""
    raw = os.environ.get("AMALGKIT_PIPELINE_COMPRESSION_LEVEL", "6")
    try:
        level = int(raw)
    except ValueError as exc:
        raise ValueError("AMALGKIT_PIPELINE_COMPRESSION_LEVEL must be an integer from 1 to 9") from exc
    if not 1 <= level <= 9:
        raise ValueError("AMALGKIT_PIPELINE_COMPRESSION_LEVEL must be an integer from 1 to 9")
    return level


def pigz_command(source: Path, *, threads: int, level: int) -> list[str]:
    """Build the native compression command with a validated level and thread budget."""
    if type(level) is not int or not 1 <= level <= 9:
        raise ValueError("compression level must be an integer from 1 to 9")
    if type(threads) is not int or threads <= 0:
        raise ValueError("compression threads must be a positive integer")
    return ["pigz", "-f", f"-{level}", "-p", str(threads), str(source)]
