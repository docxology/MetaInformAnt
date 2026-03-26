"""Phenotype/trait data loading utilities for GWAS analysis.

Provides standardised loading of trait values from CSV phenotype files,
with automatic column detection and robust numeric coercion.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import List, Optional, Tuple

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


def load_traits(
    phenotype_path: str | Path,
    *,
    trait_column: Optional[str] = None,
    id_column: Optional[str] = None,
) -> List[float]:
    """Load numeric trait values from a CSV phenotype file.

    When *trait_column* is ``None``, the last column whose header is not
    a known ID column (``sample_id``, ``sample``, ``id``) is selected.

    Args:
        phenotype_path: Path to a CSV file with a header row.
        trait_column: Explicit column name to use as trait values.
        id_column: Column name for sample IDs (unused by this function
            but reserved for future ``load_traits_with_ids``).

    Returns:
        List of float trait values (rows with non-numeric values are skipped).
    """
    traits: List[float] = []
    id_columns = {"sample_id", "sample", "id"}

    with open(phenotype_path) as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            logger.warning("Empty phenotype file: %s", phenotype_path)
            return traits

        # Resolve trait column name
        col = (
            trait_column if trait_column and trait_column in reader.fieldnames else None
        )
        if not col:
            for c in reversed(reader.fieldnames):
                if c.lower() not in id_columns:
                    col = c
                    break
        if not col:
            logger.warning(
                "No suitable trait column found in %s (headers: %s)",
                phenotype_path,
                reader.fieldnames,
            )
            return traits

        logger.info("Loading trait column '%s' from %s", col, phenotype_path)
        for row in reader:
            try:
                traits.append(float(row[col]))
            except (ValueError, KeyError):
                continue

    logger.info("Loaded %d trait values from %s", len(traits), phenotype_path)
    return traits


def load_traits_with_ids(
    phenotype_path: str | Path,
    *,
    trait_column: Optional[str] = None,
    id_column: Optional[str] = None,
) -> Tuple[List[str], List[float]]:
    """Load sample IDs and trait values from a CSV phenotype file.

    Args:
        phenotype_path: Path to CSV with header.
        trait_column: Column for trait values (auto-detected if None).
        id_column: Column for sample IDs (auto-detected if None).

    Returns:
        Tuple of (sample_ids, trait_values).
    """
    id_columns = {"sample_id", "sample", "id"}
    sample_ids: List[str] = []
    trait_values: List[float] = []

    with open(phenotype_path) as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            return sample_ids, trait_values

        # Resolve ID column
        id_col = id_column
        if not id_col:
            for c in reader.fieldnames:
                if c.lower() in id_columns:
                    id_col = c
                    break

        # Resolve trait column
        trait_col = trait_column
        if not trait_col:
            for c in reversed(reader.fieldnames):
                if c.lower() not in id_columns:
                    trait_col = c
                    break

        if not id_col or not trait_col:
            return sample_ids, trait_values

        for row in reader:
            try:
                val = float(row[trait_col])
                sample_ids.append(row[id_col])
                trait_values.append(val)
            except (ValueError, KeyError):
                continue

    return sample_ids, trait_values
