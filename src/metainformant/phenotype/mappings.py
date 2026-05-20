"""Apis mellifera phenotype-to-GWAS mapping constants and helpers.

Centralises the domain-specific biological-group and phenotype-linkage
mappings used by the BeeWAS GWAS study.  Moving these out of the
orchestrator script (``04b_phenotype_analysis.py``) makes them
importable by any downstream module that needs consistent caste/strain
semantics.
"""

from __future__ import annotations

from typing import Dict

try:
    import pandas as pd

    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

# ── Biological group mapping ─────────────────────────────────────────────
# Maps the short caste codes in the raw data to human-readable
# biological group labels used in publications and visualisations.
BIOLOGICAL_GROUP_MAP: Dict[str, str] = {
    "WORK": "Worker",
    "ITW": "In vitro worker",
    "ITQ": "In vitro queen",
    "IV": "In vivo (Queen)",
    "G": "Grafted",
}

# ── Phenotype linkage mapping ────────────────────────────────────────────
# Maps caste codes to the binary Worker/Queen phenotype axis (W / Q)
# used for GWAS trait assignment.
PHENOTYPE_LINK_MAP: Dict[str, str] = {
    "WORK": "W",
    "ITW": "W",
    "ITQ": "Q",
    "IV": "Q",
    "G": "Q",
}

# ── Strain metadata ──────────────────────────────────────────────────────
STRAIN_FULL_NAMES: Dict[str, str] = {
    "C": "Carnica",
    "I": "Italian (Ligustica)",
    "M": "Mellifera mellifera",
    "R": "Russian",
}

STRAIN_PALETTE: Dict[str, str] = {
    "C": "#2196F3",
    "I": "#FF9800",
    "M": "#4CAF50",
    "R": "#F44336",
}

STRAIN_ORDER = ["C", "I", "M", "R"]


def map_biological_groups(df: "pd.DataFrame", caste_column: str = "caste") -> "pd.DataFrame":
    """Add a ``biological_group`` column to *df* based on BIOLOGICAL_GROUP_MAP.

    Args:
        df: DataFrame with a caste column.
        caste_column: Name of the column containing caste codes.

    Returns:
        The same DataFrame with the new column added (in-place).
    """
    if not HAS_PANDAS:
        raise ImportError("pandas is required for map_biological_groups")
    df["biological_group"] = df[caste_column].map(BIOLOGICAL_GROUP_MAP)
    return df


def map_phenotype_links(df: "pd.DataFrame", caste_column: str = "caste") -> "pd.DataFrame":
    """Add a ``phenotype_link`` column mapping castes to W/Q phenotype axis.

    Args:
        df: DataFrame with a caste column.
        caste_column: Name of the column containing caste codes.

    Returns:
        The same DataFrame with the new column added (in-place).
    """
    if not HAS_PANDAS:
        raise ImportError("pandas is required for map_phenotype_links")
    df["phenotype_link"] = df[caste_column].map(PHENOTYPE_LINK_MAP)
    return df
