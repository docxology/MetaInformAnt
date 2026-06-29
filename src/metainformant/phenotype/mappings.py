"""Apis mellifera phenotype-to-GWAS mapping constants and helpers.

Centralises the domain-specific biological-group and phenotype-linkage
mappings used by the BeeWAS GWAS study.  Keeping these in the package makes
them importable by downstream modules that need consistent biological-group and
strain semantics.
"""

from __future__ import annotations

from typing import Dict

try:
    import pandas as pd

    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

# ── Biological group mapping ─────────────────────────────────────────────
# Maps the short biological-group codes in the raw data to human-readable
# biological group labels used in publications and visualisations.
BIOLOGICAL_GROUP_MAP: Dict[str, str] = {
    "WORK": "Worker",
    "ITW": "In vitro worker",
    "ITQ": "In vitro queen",
    "IV": "In vivo (Queen)",
    "G": "Grafted",
}

# ── Phenotype linkage mapping ────────────────────────────────────────────
# Maps biological-group codes to the binary Worker/Queen phenotype axis (W / Q)
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


def map_biological_groups(
    df: "pd.DataFrame", biological_group_code_column: str = "biological_group_code"
) -> "pd.DataFrame":
    """Add a ``biological_group`` column to *df* based on BIOLOGICAL_GROUP_MAP.

    Args:
        df: DataFrame with a biological-group code column.
        biological_group_code_column: Name of the column containing biological-group codes.

    Returns:
        The same DataFrame with the new column added (in-place).
    """
    if not HAS_PANDAS:
        raise ImportError("pandas is required for map_biological_groups")
    df["biological_group"] = df[biological_group_code_column].map(BIOLOGICAL_GROUP_MAP)
    return df


def map_phenotype_links(
    df: "pd.DataFrame", biological_group_code_column: str = "biological_group_code"
) -> "pd.DataFrame":
    """Add a ``phenotype_link`` column mapping biological groups to W/Q axis.

    Args:
        df: DataFrame with a biological-group code column.
        biological_group_code_column: Name of the column containing biological-group codes.

    Returns:
        The same DataFrame with the new column added (in-place).
    """
    if not HAS_PANDAS:
        raise ImportError("pandas is required for map_phenotype_links")
    df["phenotype_link"] = df[biological_group_code_column].map(PHENOTYPE_LINK_MAP)
    return df
