from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class RNAPipelineConfig:
    work_dir: Path


def summarize_curate_tables(curate_dir: Path) -> dict[str, int]:
    """Return simple counts of expected TSV tables in a curate directory.

    This uses existing repository outputs if present.
    """
    counts: dict[str, int] = {}
    if not curate_dir.is_dir():
        return counts
    for tsv in curate_dir.glob("**/*.tsv"):
        key = tsv.name
        counts[key] = counts.get(key, 0) + 1
    return counts


