"""Tissue metadata normalization for RNA-seq workflows.

This module provides functions to normalize tissue metadata values to canonical
forms while preserving original values. Supports YAML-based synonym mappings
and sample-specific patches.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Set

import pandas as pd
import yaml

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def load_tissue_mapping(mapping_path: str | Path) -> Dict[str, List[str]]:
    """Load tissue synonym mapping from YAML file.

    Args:
        mapping_path: Path to tissue_mapping.yaml

    Returns:
        Dict mapping canonical tissue names to list of synonyms
    """
    path = Path(mapping_path)
    if not path.exists():
        logger.warning(f"Tissue mapping file not found: {mapping_path}")
        return {}

    with open(path) as f:
        mapping = yaml.safe_load(f)

    # Filter out non-mapping keys (like comments)
    return {k: v for k, v in mapping.items() if isinstance(v, list)}


def load_tissue_patches(patches_path: str | Path) -> Dict[str, Dict[str, str]]:
    """Load tissue patches from YAML file.

    Args:
        patches_path: Path to tissue_patches.yaml

    Returns:
        Dict with 'samples', 'bioprojects', 'biosamples' patch mappings
    """
    path = Path(patches_path)
    if not path.exists():
        logger.warning(f"Tissue patches file not found: {patches_path}")
        return {"samples": {}, "bioprojects": {}, "biosamples": {}}

    with open(path) as f:
        patches = yaml.safe_load(f) or {}

    return {
        "samples": patches.get("samples") or {},
        "bioprojects": patches.get("bioprojects") or {},
        "biosamples": patches.get("biosamples") or {},
    }


def build_synonym_lookup(mapping: Dict[str, List[str]]) -> Dict[str, str]:
    """Build reverse lookup from synonym to canonical name.

    Args:
        mapping: Dict mapping canonical names to synonym lists

    Returns:
        Dict mapping each synonym (lowercase) to its canonical name
    """
    lookup = {}
    for canonical, synonyms in mapping.items():
        for synonym in synonyms:
            # Case-insensitive lookup
            lookup[synonym.lower().strip()] = canonical
    return lookup


def normalize_tissue(raw_value: str, synonym_lookup: Dict[str, str], default: str = "unknown") -> str:
    """Normalize a single tissue value to canonical form.

    Uses exact match first, then prefix-based fallback for long-form
    descriptions (e.g. "fat body of 1 queen; kept in a group..." → fat_body).

    Args:
        raw_value: Original tissue value from metadata
        synonym_lookup: Dict mapping synonyms to canonical names
        default: Value to return if no mapping found

    Returns:
        Canonical tissue name or default if unmapped
    """
    if not raw_value or pd.isna(raw_value):
        return ""

    # Clean and lowercase for lookup
    cleaned = str(raw_value).lower().strip()

    # Exact match
    result = synonym_lookup.get(cleaned, "")
    if result:
        return result

    # Prefix-based fallback: check if value starts with a known synonym
    # Sort by length (longest first) to prefer more specific matches
    for synonym, canonical in sorted(synonym_lookup.items(), key=lambda x: -len(x[0])):
        if cleaned.startswith(synonym) and len(synonym) >= 3:
            return canonical

    return ""


def apply_tissue_normalization(
    df: pd.DataFrame,
    mapping_path: str | Path,
    patches_path: Optional[str | Path] = None,
    tissue_column: str = "tissue",
    run_column: str = "run",
    bioproject_column: str = "bioproject",
    biosample_column: str = "biosample",
    output_column: str = "tissue_normalized",
) -> pd.DataFrame:
    """Apply tissue normalization to a metadata DataFrame.

    Adds a new column with normalized tissue values. Never modifies
    the original tissue column.

    Args:
        df: Metadata DataFrame
        mapping_path: Path to tissue_mapping.yaml
        patches_path: Optional path to tissue_patches.yaml
        tissue_column: Name of original tissue column
        run_column: Name of run accession column
        bioproject_column: Name of bioproject column
        biosample_column: Name of biosample column
        output_column: Name for normalized tissue column

    Returns:
        DataFrame with added tissue_normalized column
    """
    # Load configurations
    mapping = load_tissue_mapping(mapping_path)
    synonym_lookup = build_synonym_lookup(mapping)

    patches = {}
    if patches_path:
        patches = load_tissue_patches(patches_path)

    # Create copy to avoid modifying original
    result = df.copy()

    # Initialize normalized column
    result[output_column] = ""

    # Track statistics
    stats = {
        "total": len(df),
        "normalized": 0,
        "patched_sample": 0,
        "patched_bioproject": 0,
        "patched_biosample": 0,
        "unmapped": 0,
        "empty": 0,
    }

    for idx, row in result.iterrows():
        run_id = row.get(run_column, "")
        bioproject = row.get(bioproject_column, "")
        biosample = row.get(biosample_column, "")
        raw_tissue = row.get(tissue_column, "")

        normalized = ""

        # Priority 1: Sample-specific patch
        if run_id and run_id in patches.get("samples", {}):
            normalized = patches["samples"][run_id]
            stats["patched_sample"] += 1

        # Priority 2: Biosample patch
        elif biosample and biosample in patches.get("biosamples", {}):
            normalized = patches["biosamples"][biosample]
            stats["patched_biosample"] += 1

        # Priority 3: Bioproject patch
        elif bioproject and bioproject in patches.get("bioprojects", {}):
            normalized = patches["bioprojects"][bioproject]
            stats["patched_bioproject"] += 1

        # Priority 4: Synonym mapping
        else:
            normalized = normalize_tissue(raw_tissue, synonym_lookup)
            if normalized:
                stats["normalized"] += 1
            elif raw_tissue and not pd.isna(raw_tissue) and str(raw_tissue).strip():
                stats["unmapped"] += 1
            else:
                stats["empty"] += 1

        result.at[idx, output_column] = normalized

    # Log statistics
    logger.info(f"Tissue normalization complete:")
    logger.info(f"  Total samples: {stats['total']}")
    logger.info(f"  Normalized via mapping: {stats['normalized']}")
    logger.info(f"  Patched (sample): {stats['patched_sample']}")
    logger.info(f"  Patched (bioproject): {stats['patched_bioproject']}")
    logger.info(f"  Patched (biosample): {stats['patched_biosample']}")
    logger.info(f"  Unmapped (needs review): {stats['unmapped']}")
    logger.info(f"  Empty tissue field: {stats['empty']}")

    return result


def get_unmapped_tissues(df: pd.DataFrame, mapping_path: str | Path, tissue_column: str = "tissue") -> Dict[str, int]:
    """Get tissue values that don't have mappings.

    Args:
        df: Metadata DataFrame
        mapping_path: Path to tissue_mapping.yaml
        tissue_column: Name of tissue column

    Returns:
        Dict mapping unmapped tissue values to their counts
    """
    mapping = load_tissue_mapping(mapping_path)
    synonym_lookup = build_synonym_lookup(mapping)

    unmapped = {}

    for value in df[tissue_column]:
        if not value or pd.isna(value):
            continue

        cleaned = str(value).lower().strip()
        if cleaned not in synonym_lookup:
            raw = str(value).strip()
            unmapped[raw] = unmapped.get(raw, 0) + 1

    # Sort by count descending
    return dict(sorted(unmapped.items(), key=lambda x: -x[1]))


def get_missing_tissue_samples(
    df: pd.DataFrame, tissue_column: str = "tissue", run_column: str = "run", bioproject_column: str = "bioproject"
) -> pd.DataFrame:
    """Get samples with missing tissue information.

    Args:
        df: Metadata DataFrame
        tissue_column: Name of tissue column
        run_column: Name of run column
        bioproject_column: Name of bioproject column

    Returns:
        DataFrame with samples missing tissue, grouped by bioproject
    """
    missing = df[df[tissue_column].isna() | (df[tissue_column].astype(str).str.strip() == "")][
        [run_column, bioproject_column]
    ].copy()

    return missing


def get_canonical_tissues(mapping_path: str | Path) -> List[str]:
    """Get list of canonical tissue names from mapping.

    Args:
        mapping_path: Path to tissue_mapping.yaml

    Returns:
        Sorted list of canonical tissue names
    """
    mapping = load_tissue_mapping(mapping_path)
    return sorted(mapping.keys())
