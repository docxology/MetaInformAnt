"""Sample metadata loading, validation, and merging for GWAS analysis.

This module provides utilities for handling sample metadata (TSV files),
merging with phenotype data, validating completeness, and extracting
geographic and population information for downstream GWAS workflows.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def load_sample_metadata(metadata_file: Union[str, Path]) -> Dict[str, Any]:
    """Load sample metadata from a TSV file.

    Expects a tab-separated file with a header row. The first column must be
    ``sample_id``. All remaining columns are treated as metadata fields and
    stored per sample.

    Args:
        metadata_file: Path to the TSV metadata file.

    Returns:
        Dictionary with keys:
        - status: "success" or "error"
        - metadata: ``{sample_id: {col: value, ...}}``
        - columns: list of column names (excluding sample_id)
        - n_samples: number of samples loaded
        - message: (on error) human-readable description
    """
    metadata_path = Path(metadata_file)

    if not metadata_path.exists():
        logger.error(f"Metadata file not found: {metadata_path}")
        return {
            "status": "error",
            "message": f"Metadata file not found: {metadata_path}",
            "metadata": {},
            "columns": [],
            "n_samples": 0,
        }

    try:
        metadata: Dict[str, Dict[str, Any]] = {}
        columns: List[str] = []

        with open(metadata_path, "r", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")

            if reader.fieldnames is None:
                logger.warning(f"Empty or headerless metadata file: {metadata_path}")
                return {
                    "status": "error",
                    "message": "Metadata file has no header row",
                    "metadata": {},
                    "columns": [],
                    "n_samples": 0,
                }

            fieldnames = list(reader.fieldnames)

            if "sample_id" not in fieldnames:
                logger.error(f"Missing 'sample_id' column in {metadata_path}")
                return {
                    "status": "error",
                    "message": "Missing required 'sample_id' column",
                    "metadata": {},
                    "columns": [],
                    "n_samples": 0,
                }

            columns = [c for c in fieldnames if c != "sample_id"]

            for row in reader:
                sample_id = row.get("sample_id", "").strip()
                if not sample_id:
                    continue
                sample_data: Dict[str, Any] = {}
                for col in columns:
                    sample_data[col] = row.get(col, "")
                metadata[sample_id] = sample_data

        logger.info(f"Loaded metadata for {len(metadata)} samples " f"with {len(columns)} columns from {metadata_path}")

        return {
            "status": "success",
            "metadata": metadata,
            "columns": columns,
            "n_samples": len(metadata),
        }

    except Exception as exc:
        logger.error(f"Failed to load metadata from {metadata_path}: {exc}")
        return {
            "status": "error",
            "message": f"Failed to load metadata: {exc}",
            "metadata": {},
            "columns": [],
            "n_samples": 0,
        }


def merge_metadata_with_phenotypes(
    metadata: Dict[str, Dict[str, Any]],
    phenotypes: Dict[str, Any],
) -> Dict[str, Any]:
    """Merge sample metadata with phenotype data.

    For each sample present in *metadata*, any matching key in *phenotypes*
    is added under a ``"phenotype"`` sub-key. Samples without phenotype data
    receive ``"phenotype": None``.

    Args:
        metadata: ``{sample_id: {col: value, ...}}`` as returned by
            :func:`load_sample_metadata`.
        phenotypes: ``{sample_id: phenotype_value}`` mapping.

    Returns:
        Dictionary with keys:
        - status: "success"
        - merged: ``{sample_id: {**metadata_fields, "phenotype": value}}``
        - n_matched: number of samples with phenotype data
        - n_unmatched: number of samples missing phenotype data
    """
    merged: Dict[str, Dict[str, Any]] = {}
    n_matched = 0
    n_unmatched = 0

    for sample_id, sample_meta in metadata.items():
        entry = dict(sample_meta)
        if sample_id in phenotypes:
            entry["phenotype"] = phenotypes[sample_id]
            n_matched += 1
        else:
            entry["phenotype"] = None
            n_unmatched += 1
        merged[sample_id] = entry

    logger.info(f"Merged metadata with phenotypes: {n_matched} matched, " f"{n_unmatched} unmatched")

    return {
        "status": "success",
        "merged": merged,
        "n_matched": n_matched,
        "n_unmatched": n_unmatched,
    }


def validate_metadata(
    metadata: Dict[str, Dict[str, Any]],
    sample_ids: List[str],
) -> Dict[str, Any]:
    """Validate metadata completeness against a list of expected sample IDs.

    Reports samples that are present in *sample_ids* but missing from
    *metadata* (``missing_samples``), and samples in *metadata* that are
    not in *sample_ids* (``extra_samples``).

    Args:
        metadata: ``{sample_id: {col: value, ...}}``
        sample_ids: Expected sample IDs (e.g. from genotype data).

    Returns:
        Dictionary with keys:
        - status: "success"
        - valid: ``True`` if no missing samples
        - missing_samples: sample IDs expected but absent from metadata
        - extra_samples: sample IDs in metadata but not in expected list
        - completeness: fraction of expected samples present (0.0 -- 1.0)
    """
    expected = set(sample_ids)
    present = set(metadata.keys())

    missing_samples = sorted(expected - present)
    extra_samples = sorted(present - expected)

    completeness = 0.0
    if expected:
        completeness = len(expected & present) / len(expected)

    valid = len(missing_samples) == 0

    if missing_samples:
        logger.warning(f"Metadata validation: {len(missing_samples)} samples missing " f"(e.g. {missing_samples[:3]})")
    if extra_samples:
        logger.info(f"Metadata validation: {len(extra_samples)} extra samples in metadata " f"not in expected list")

    return {
        "status": "success",
        "valid": valid,
        "missing_samples": missing_samples,
        "extra_samples": extra_samples,
        "completeness": completeness,
    }


def get_population_labels(
    metadata: Dict[str, Dict[str, Any]],
    column: str = "population",
) -> Dict[str, str]:
    """Extract population or strain labels from metadata.

    Args:
        metadata: ``{sample_id: {col: value, ...}}``
        column: Metadata column to use for population label (default:
            ``"population"``).

    Returns:
        ``{sample_id: population_label}`` mapping.  Samples where the
        column is missing or empty are omitted.
    """
    labels: Dict[str, str] = {}

    for sample_id, sample_meta in metadata.items():
        value = sample_meta.get(column, "")
        if isinstance(value, str):
            value = value.strip()
        if value:
            labels[sample_id] = str(value)

    logger.info(f"Extracted population labels for {len(labels)}/{len(metadata)} samples " f"using column '{column}'")

    return labels


def get_geographic_coordinates(
    metadata: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Extract geographic coordinates from metadata.

    Looks for ``latitude`` and ``longitude`` columns. Non-numeric values
    are silently skipped with a debug-level log message.

    Args:
        metadata: ``{sample_id: {col: value, ...}}``

    Returns:
        List of dictionaries, each with keys ``sample_id``, ``latitude``,
        and ``longitude`` (floats).
    """
    coordinates: List[Dict[str, Any]] = []

    for sample_id, sample_meta in metadata.items():
        lat_raw = sample_meta.get("latitude", "")
        lon_raw = sample_meta.get("longitude", "")

        if not lat_raw or not lon_raw:
            continue

        try:
            lat = float(lat_raw)
            lon = float(lon_raw)
        except (ValueError, TypeError):
            logger.debug(f"Skipping non-numeric coordinates for {sample_id}: " f"lat={lat_raw!r}, lon={lon_raw!r}")
            continue

        coordinates.append(
            {
                "sample_id": sample_id,
                "latitude": lat,
                "longitude": lon,
            }
        )

    logger.info(f"Extracted geographic coordinates for {len(coordinates)}/{len(metadata)} samples")

    return coordinates
