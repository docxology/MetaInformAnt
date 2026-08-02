"""RNA-seq pipeline utilities and high-level workflow orchestration.

This module provides high-level functions for RNA-seq analysis pipelines,
including result summarization and workflow coordination.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict

from metainformant.core.utils import logging
from metainformant.rna.core.configs import RNAPipelineConfig

logger = logging.get_logger(__name__)

__all__ = ["RNAPipelineConfig", "summarize_finalize_tables"]


def summarize_finalize_tables(finalize_dir: str | Path) -> Dict[str, int]:
    """Summarize current Amalgkit finalization outputs.

    Scans the finalization directory for TSV files and current metadata,
    exclusion, and statistics artifacts.

    Args:
        finalize_dir: Path to the finalization output directory

    Returns:
        Dictionary mapping filename to count of occurrences.
        Special keys:
            - "_error": Error message if directory doesn't exist or is inaccessible
            - "_total": Total file count

    Example:
        >>> counts = summarize_finalize_tables("output/amalgkit/finalize/")
        >>> print(counts)
        {'metadata.tsv': 1, 'tc.tsv': 1, 'stats.json': 1, '_total': 3}
    """
    finalize_path = Path(finalize_dir)
    if not finalize_path.exists():
        logger.warning(f"Finalize directory does not exist: {finalize_path}")
        return {"_error": f"Directory not found: {finalize_path}", "_total": 0}

    counts = {}

    # Count only TSV files by their basename
    for tsv_file in finalize_path.glob("**/*.tsv"):
        name = tsv_file.name
        counts[name] = counts.get(name, 0) + 1

    other_expected = ["stats.json"]

    for expected_file in other_expected:
        if (finalize_path / expected_file).exists():
            counts[expected_file] = counts.get(expected_file, 0) + 1

    total = sum(counts.values())
    counts["_total"] = total
    logger.info(f"Summarized {total} files in finalize directory")
    return counts
