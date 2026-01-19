"""RNA-seq pipeline utilities and high-level workflow orchestration.

This module provides high-level functions for RNA-seq analysis pipelines,
including result summarization and workflow coordination.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Any

from metainformant.core import logging
from metainformant.rna.core.configs import RNAPipelineConfig

logger = logging.get_logger(__name__)


def summarize_curate_tables(curate_dir: str | Path) -> Dict[str, int]:
    """Summarize the contents of amalgkit curate output directory.

    Args:
        curate_dir: Path to the curate output directory

    Returns:
        Dictionary mapping file patterns to counts

    Example:
        >>> counts = summarize_curate_tables("output/amalgkit/curate/")
        >>> print(counts)
        {'metadata.tsv': 1, 'tc.tsv': 1, 'stats.json': 1}
    """
    curate_path = Path(curate_dir)
    if not curate_path.exists():
        logger.warning(f"Curate directory does not exist: {curate_path}")
        return {}

    counts = {}

    # Count only TSV files by their basename
    for tsv_file in curate_path.glob("**/*.tsv"):
        name = tsv_file.name
        counts[name] = counts.get(name, 0) + 1

    # Specifically check for other expected amalgkit curate outputs (non-TSV)
    # like stats.json and outlier_samples.txt
    other_expected = [
        "stats.json",
        "outlier_samples.txt",
    ]

    for expected_file in other_expected:
        if (curate_path / expected_file).exists():
            counts[expected_file] = counts.get(expected_file, 0) + 1

    logger.info(f"Summarized {sum(counts.values())} files in curate directory")
    return counts
