"""Metadata filtering utilities for RNA-seq workflows.

This module provides functions to filter metadata tables based on selection criteria,
ensuring downstream steps only process selected samples.
"""

from __future__ import annotations

import pandas as pd
from pathlib import Path
from typing import Optional

from metainformant.core import logging

logger = logging.get_logger(__name__)


def filter_selected_metadata(
    metadata_path: str | Path,
    output_path: Optional[str | Path] = None,
    require_sampled: bool = True,
    require_qualified: bool = True,
    exclude_excluded: bool = True
) -> Path:
    """Filter metadata to only selected samples.

    Args:
        metadata_path: Path to full metadata.tsv file
        output_path: Path for filtered metadata file.
                    Default: metadata_path.parent / "metadata_selected.tsv"
        require_sampled: Only include rows where is_sampled=yes
        require_qualified: Only include rows where is_qualified=yes
        exclude_excluded: Exclude rows where exclusion column is not empty

    Returns:
        Path to filtered metadata file

    Raises:
        FileNotFoundError: If metadata_path doesn't exist
        ValueError: If no samples meet the filtering criteria
    """
    metadata_path = Path(metadata_path)

    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

    # Determine output path
    if output_path is None:
        output_path = metadata_path.parent / "metadata_selected.tsv"
    else:
        output_path = Path(output_path)

    logger.info(f"Filtering metadata from {metadata_path} to {output_path}")

    # Read metadata
    try:
        df = pd.read_csv(metadata_path, sep='\t', low_memory=False)
        logger.debug(f"Loaded {len(df)} samples from metadata")
    except Exception as e:
        raise ValueError(f"Failed to read metadata file: {e}")

    # Check required columns exist
    required_cols = []
    if require_sampled:
        required_cols.append('is_sampled')
    if require_qualified:
        required_cols.append('is_qualified')
    if exclude_excluded:
        required_cols.append('exclusion')

    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Required columns missing from metadata: {missing_cols}")

    # Apply filters
    original_count = len(df)

    if require_sampled:
        df = df[df['is_sampled'] == 'yes']
        logger.debug(f"After is_sampled=yes filter: {len(df)} samples")

    if require_qualified:
        df = df[df['is_qualified'] == 'yes']
        logger.debug(f"After is_qualified=yes filter: {len(df)} samples")

    if exclude_excluded:
        # Exclude rows where exclusion column has values (not empty/NaN and not "no")
        df = df[df['exclusion'].isna() | (df['exclusion'] == '') | (df['exclusion'] == 'no')]
        logger.debug(f"After exclusion filter: {len(df)} samples")

    filtered_count = len(df)

    if filtered_count == 0:
        raise ValueError("No samples meet the filtering criteria")

    if filtered_count == original_count:
        logger.warning("Filtering criteria had no effect - all samples were already selected")

    # Write filtered metadata
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path, sep='\t', index=False)
        logger.info(f"Created filtered metadata with {filtered_count} samples: {output_path}")
    except Exception as e:
        raise ValueError(f"Failed to write filtered metadata: {e}")

    return output_path


def get_sample_count(metadata_path: str | Path) -> int:
    """Get the number of samples in a metadata file.

    Args:
        metadata_path: Path to metadata file

    Returns:
        Number of samples (excluding header row)

    Raises:
        FileNotFoundError: If metadata_path doesn't exist
    """
    metadata_path = Path(metadata_path)

    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

    try:
        df = pd.read_csv(metadata_path, sep='\t', low_memory=False)
        return len(df)
    except Exception as e:
        raise ValueError(f"Failed to read metadata file: {e}")


def validate_filtered_metadata(
    original_path: str | Path,
    filtered_path: str | Path,
    expected_count: Optional[int] = None
) -> bool:
    """Validate that filtered metadata looks correct.

    Args:
        original_path: Path to original metadata file
        filtered_path: Path to filtered metadata file
        expected_count: Expected number of samples in filtered file

    Returns:
        True if validation passes

    Raises:
        ValueError: If validation fails
    """
    try:
        original_count = get_sample_count(original_path)
        filtered_count = get_sample_count(filtered_path)

        logger.info(f"Validation: original={original_count}, filtered={filtered_count}")

        if filtered_count > original_count:
            raise ValueError("Filtered metadata has more samples than original")

        if filtered_count == 0:
            raise ValueError("Filtered metadata is empty")

        if expected_count is not None and filtered_count != expected_count:
            logger.warning(f"Filtered count ({filtered_count}) differs from expected ({expected_count})")

        # Check that filtered metadata has required columns
        df = pd.read_csv(filtered_path, sep='\t', nrows=1)
        required_cols = ['run', 'experiment', 'scientific_name']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Filtered metadata missing required columns: {missing_cols}")

        return True

    except Exception as e:
        raise ValueError(f"Metadata validation failed: {e}")
