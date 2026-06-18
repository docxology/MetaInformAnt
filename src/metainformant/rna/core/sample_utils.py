"""Shared RNA workflow helpers for sample IDs and quantification outputs."""

from __future__ import annotations

import csv
from collections.abc import Mapping
from pathlib import Path
from typing import Any

SAMPLE_ID_COLUMNS = ("run", "sra_run", "sample_id", "accession")


def extract_sample_id(row: Mapping[str, Any], fallback: str | None = None) -> str:
    """Return a sample/run identifier from a metadata row.

    The RNA workflow receives metadata from several sources whose accession
    column names vary by case and provider. This helper keeps the accepted
    column set consistent across validation, streaming, and internal download
    paths.
    """
    lower_to_key = {str(key).lower(): key for key in row.keys()}
    for candidate in SAMPLE_ID_COLUMNS:
        key = lower_to_key.get(candidate)
        if key is None:
            continue
        value = row.get(key)
        if value is not None and str(value).strip():
            return str(value).strip()
    return fallback or ""


def read_sample_ids_from_metadata(metadata_path: Path) -> list[str]:
    """Read sample IDs from an amalgkit-style TSV metadata file."""
    sample_ids: list[str] = []
    with open(metadata_path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample_id = extract_sample_id(row)
            if sample_id:
                sample_ids.append(sample_id)
    return sample_ids


def quantification_file_candidates(sample_quant_dir: Path, sample_id: str | None = None) -> list[Path]:
    """Return known per-sample quantification output paths in priority order."""
    candidates = [sample_quant_dir / "abundance.tsv"]
    if sample_id:
        candidates.append(sample_quant_dir / f"{sample_id}_abundance.tsv")
    candidates.append(sample_quant_dir / "quant.sf")
    return candidates


def find_quantification_file(
    sample_quant_dir: Path,
    sample_id: str | None = None,
    *,
    require_nonempty: bool = True,
) -> Path | None:
    """Find a recognized quantification output for one sample directory."""
    for candidate in quantification_file_candidates(sample_quant_dir, sample_id):
        if not candidate.exists():
            continue
        if require_nonempty and candidate.stat().st_size <= 0:
            continue
        return candidate

    for candidate in sorted(sample_quant_dir.glob("*_abundance.tsv")):
        if require_nonempty and candidate.stat().st_size <= 0:
            continue
        return candidate

    return None
