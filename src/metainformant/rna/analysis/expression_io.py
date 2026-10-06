"""Streaming finalized-matrix validation and descriptive expression profiles."""

from __future__ import annotations
import gzip
from pathlib import Path
import numpy as np
import pandas as pd


def _open_text(path: Path):
    """Open plain or gzip-compressed tabular text using one interface."""

    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def validate_expression_matrix(
    path: Path,
    *,
    min_samples: int = 1,
    chunk_size: int = 4096,
) -> tuple[int, int]:
    """Validate and dimension a finalized expression matrix.

    The cross-species contract is deliberately stricter than "the file exists":
    the header must contain unique, non-empty sample identifiers; feature IDs
    must be unique and non-empty; every value must be finite and non-negative;
    and every sample must contain at least one positive value.  The table is
    streamed in chunks so validation does not require a second full in-memory
    copy of a large finalized matrix.
    """

    if min_samples < 1:
        raise ValueError("min_samples must be at least 1")
    if chunk_size < 1:
        raise ValueError("chunk_size must be at least 1")
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Expression matrix is missing or empty: {path}")

    with _open_text(path) as handle:
        raw_header = handle.readline()
    if not raw_header:
        raise ValueError(f"Expression matrix has no header: {path}")
    header = raw_header.rstrip("\r\n").split("\t")
    if len(header) < min_samples + 1:
        requirement = (
            "At least two sample columns"
            if min_samples == 2
            else f"At least {min_samples} sample column(s)"
        )
        raise ValueError(f"{requirement} are required in expression matrix: {path}")
    if not header[0].strip():
        raise ValueError(f"Expression matrix has an empty feature-ID column: {path}")
    sample_ids = [value.strip() for value in header[1:]]
    if any(not value for value in sample_ids):
        raise ValueError(f"Expression matrix has an empty sample identifier: {path}")
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError(f"Expression matrix has duplicate sample identifiers: {path}")

    seen_features: set[str] = set()
    nonzero_samples = np.zeros(len(sample_ids), dtype=bool)
    feature_count = 0
    try:
        chunks = pd.read_csv(
            path,
            sep="\t",
            compression="infer",
            index_col=0,
            dtype=str,
            keep_default_na=False,
            chunksize=chunk_size,
        )
        for chunk in chunks:
            if chunk.shape[1] != len(sample_ids):
                raise ValueError(
                    f"Expression matrix row width does not match its header: {path}"
                )
            feature_ids = [str(value).strip() for value in chunk.index]
            if any(not value for value in feature_ids):
                raise ValueError(
                    f"Expression matrix has an empty feature identifier: {path}"
                )
            local_features = set(feature_ids)
            feature_index = pd.Index(feature_ids)
            local_duplicates = set(feature_index[feature_index.duplicated()].tolist())
            duplicates = (local_features & seen_features) | local_duplicates
            if duplicates:
                detail = sorted(duplicates)[0]
                raise ValueError(
                    f"Expression matrix has duplicate feature identifier ({detail}): {path}"
                )
            seen_features.update(local_features)

            numeric = chunk.apply(pd.to_numeric, errors="coerce")
            values = numeric.to_numpy(dtype=float)
            if not np.isfinite(values).all():
                raise ValueError(
                    f"Expression matrix contains non-numeric, missing, or non-finite values: {path}"
                )
            if (values < 0).any():
                raise ValueError(f"Expression matrix contains negative values: {path}")
            nonzero_samples |= np.any(values > 0, axis=0)
            feature_count += len(chunk)
    except (
        OSError,
        UnicodeError,
        pd.errors.ParserError,
        pd.errors.EmptyDataError,
    ) as exc:
        raise ValueError(f"Could not parse expression matrix {path}: {exc}") from exc

    if feature_count == 0:
        raise ValueError(f"Expression matrix has no feature rows: {path}")
    if not nonzero_samples.all():
        empty_samples = [
            sample_ids[i] for i, present in enumerate(nonzero_samples) if not present
        ]
        raise ValueError(
            "Expression matrix has all-zero sample column(s): "
            + ", ".join(empty_samples[:5])
            + (" ..." if len(empty_samples) > 5 else "")
        )
    return len(sample_ids), feature_count


def load_expression_profile(finalize_path: str) -> pd.Series:
    """Load a finalized matrix and compute a mean expression feature profile."""
    validate_expression_matrix(Path(finalize_path), min_samples=2)
    df = pd.read_csv(finalize_path, sep="\t", index_col=0, compression="infer")
    numeric = df.apply(pd.to_numeric, errors="raise")
    profile = numeric.mean(axis=1)
    if profile.empty or not np.isfinite(profile.to_numpy(dtype=float)).all():
        raise ValueError(f"No finite expression values found in {finalize_path}")
    if (profile.to_numpy(dtype=float) < 0).any():
        raise ValueError(f"Negative expression values found in {finalize_path}")
    return profile


def compute_profile_quality(species_profiles: dict[str, pd.Series]) -> pd.DataFrame:
    """Summarize profile validity without treating missing values as zeros."""

    rows: list[dict[str, float | int | str]] = []
    for species, profile in sorted(species_profiles.items()):
        values = pd.to_numeric(profile, errors="coerce").to_numpy(dtype=float)
        finite = np.isfinite(values)
        negative = finite & (values < 0)
        if negative.any():
            raise ValueError(f"Negative expression values found for {species}")
        positive = finite & (values > 0)
        zero = finite & (values == 0)
        rows.append(
            {
                "species": species,
                "total_features": int(values.size),
                "finite_features": int(finite.sum()),
                "positive_features": int(positive.sum()),
                "zero_features": int(zero.sum()),
                "nonfinite_features": int((~finite).sum()),
                "positive_fraction_finite": float(positive.sum() / finite.sum())
                if finite.any()
                else np.nan,
                "mean_positive_expression": float(values[positive].mean())
                if positive.any()
                else np.nan,
                "median_positive_expression": float(np.median(values[positive]))
                if positive.any()
                else np.nan,
            }
        )
    return pd.DataFrame(rows)
