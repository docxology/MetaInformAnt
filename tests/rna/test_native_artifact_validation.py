"""Real numerical controls for native cross-artifact evidence gates."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from metainformant.rna.analysis.expression_io import compute_profile_quality
from metainformant.rna.analysis.native_artifact_validation import (
    validate_divergence_stability_table,
    validate_profile_quality_table,
)


def quality_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    manifest = pd.DataFrame({"species_title": ["A", "B"], "features": [3, 3]})
    quality = compute_profile_quality({"A": pd.Series([1.0, 0.0, np.nan]), "B": pd.Series([2.0, 3.0, 4.0])})
    return manifest, quality


def stability_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    matrix = pd.DataFrame([[0.0, 1.0], [1.0, 0.0]], index=["A", "B"], columns=["A", "B"])
    # A point estimate need not fall inside a resampling-sensitivity interval.
    table = pd.DataFrame(
        [
            {
                "species_a": "A",
                "species_b": "B",
                "point_estimate": 1.0,
                "sensitivity_median": 0.35,
                "sensitivity_lower": 0.2,
                "sensitivity_upper": 0.7,
                "sensitivity_iqr": 0.1,
                "median_absolute_shift": 0.6,
                "replicate_count": 20,
            }
        ]
    )
    return matrix, table


def test_generated_quality_and_sensitivity_tables_are_valid() -> None:
    manifest, quality = quality_data()
    validate_profile_quality_table(manifest, quality)
    matrix, stability = stability_data()
    validate_divergence_stability_table(matrix, stability)


@pytest.mark.parametrize(
    "field,value",
    [
        ("total_features", 4),
        ("finite_features", 3),
        ("zero_features", 2),
        ("positive_features", -1),
        ("nonfinite_features", 0.5),
        ("positive_fraction_finite", 0.9),
        ("mean_positive_expression", float("inf")),
        ("median_positive_expression", 0),
    ],
)
def test_quality_rejects_unreconciled_or_invalid_evidence(field: str, value: object) -> None:
    manifest, quality = quality_data()
    quality[field] = quality[field].astype(object)
    quality.loc[0, field] = value
    with pytest.raises(ValueError):
        validate_profile_quality_table(manifest, quality)


def test_unavailable_quality_summaries_cannot_be_imputed_zero() -> None:
    manifest = pd.DataFrame({"species_title": ["A"], "features": [3]})
    quality = compute_profile_quality({"A": pd.Series([np.nan, np.inf, -np.inf])})
    validate_profile_quality_table(manifest, quality)
    quality["mean_positive_expression"] = 0.0
    with pytest.raises(ValueError):
        validate_profile_quality_table(manifest, quality)


@pytest.mark.parametrize("kind", ["empty", "duplicate", "unknown"])
def test_quality_requires_the_nonempty_declared_species_partition(kind: str) -> None:
    manifest, quality = quality_data()
    if kind == "empty":
        manifest, quality = manifest.iloc[:0], quality.iloc[:0]
    elif kind == "duplicate":
        quality.loc[1, "species"] = "A"
    else:
        quality.loc[1, "species"] = "C"
    with pytest.raises(ValueError):
        validate_profile_quality_table(manifest, quality)


@pytest.mark.parametrize(
    "field,value",
    [
        ("point_estimate", 0.9),
        ("point_estimate", True),
        ("sensitivity_lower", 0.8),
        ("sensitivity_upper", float("inf")),
        ("sensitivity_median", 0.1),
        ("sensitivity_iqr", -0.1),
        ("median_absolute_shift", 2.1),
        ("replicate_count", 19),
        ("replicate_count", 20.5),
    ],
)
def test_stability_rejects_invalid_bounds_counts_and_point_values(field: str, value: object) -> None:
    matrix, table = stability_data()
    table[field] = table[field].astype(object)
    table.loc[0, field] = value
    with pytest.raises(ValueError):
        validate_divergence_stability_table(matrix, table)


@pytest.mark.parametrize("kind", ["empty", "duplicate", "unknown", "diagonal", "asymmetric"])
def test_stability_requires_unique_complete_pairs_and_a_valid_matrix(kind: str) -> None:
    matrix, table = stability_data()
    if kind == "empty":
        table = table.iloc[:0]
    elif kind == "duplicate":
        table = pd.concat([table, table], ignore_index=True)
    elif kind == "unknown":
        table.loc[0, "species_b"] = "C"
    elif kind == "diagonal":
        table.loc[0, "species_b"] = "A"
    else:
        matrix.loc["A", "B"] = 0.5
    with pytest.raises(ValueError):
        validate_divergence_stability_table(matrix, table)
