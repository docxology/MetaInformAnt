"""Cross-artifact contracts for descriptive native-feature diagnostics."""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd


def _numeric(frame: pd.DataFrame) -> pd.DataFrame:
    if any(isinstance(value, (bool, np.bool_)) for value in frame.to_numpy(dtype=object).flat):
        raise ValueError("numeric artifacts cannot contain boolean values")
    return frame.apply(pd.to_numeric, errors="raise")


def validate_profile_quality_table(manifest: pd.DataFrame, quality: pd.DataFrame) -> None:
    """Require a complete species partition with reconciled feature counts."""
    count_fields = [
        "total_features",
        "finite_features",
        "positive_features",
        "zero_features",
        "nonfinite_features",
    ]
    required = {
        "species",
        *count_fields,
        "positive_fraction_finite",
        "mean_positive_expression",
        "median_positive_expression",
    }
    if not {"species_title", "features"} <= set(manifest) or not required <= set(quality):
        raise ValueError("profile-quality or manifest columns are incomplete")
    titles = manifest["species_title"]
    if (
        manifest.empty
        or not titles.map(lambda value: isinstance(value, str)).all()
        or titles.isna().any()
        or titles.duplicated().any()
        or titles.astype(str).str.strip().eq("").any()
    ):
        raise ValueError("manifest species must be nonempty and unique")
    if (
        quality["species"].isna().any()
        or quality["species"].duplicated().any()
        or set(quality["species"]) != set(titles)
    ):
        raise ValueError("profile-quality species differ from the manifest")
    counts = _numeric(quality[count_fields]).to_numpy(dtype=float)
    if not np.isfinite(counts).all() or (counts < 0).any() or (counts != np.floor(counts)).any():
        raise ValueError("profile-quality counts must be finite nonnegative integers")
    total, finite, positive, zero, nonfinite = counts.T
    expected = (
        _numeric(manifest.set_index("species_title")[["features"]])["features"]
        .reindex(quality["species"])
        .to_numpy(dtype=float)
    )
    if (
        (total < 1).any()
        or not np.array_equal(total, expected)
        or not np.array_equal(finite + nonfinite, total)
        or not np.array_equal(positive + zero, finite)
    ):
        raise ValueError("profile-quality feature counts do not reconcile")
    fraction = _numeric(quality[["positive_fraction_finite"]]).to_numpy(dtype=float).ravel()
    expected_fraction = np.divide(positive, finite, out=np.full_like(positive, np.nan), where=finite > 0)
    if not np.allclose(fraction, expected_fraction, rtol=0, atol=1e-12, equal_nan=True):
        raise ValueError("profile-quality fractions disagree with their denominator")
    summaries = _numeric(quality[["mean_positive_expression", "median_positive_expression"]]).to_numpy(dtype=float)
    if (
        not np.isfinite(summaries[positive > 0]).all()
        or (summaries[positive > 0] <= 0).any()
        or not np.isnan(summaries[positive == 0]).all()
    ):
        raise ValueError("profile-quality positive summaries must retain unavailable values")


def validate_divergence_stability_table(matrix: pd.DataFrame, stability: pd.DataFrame) -> None:
    """Require every unordered pair, matching point values and bounded sensitivity."""
    required = {
        "species_a",
        "species_b",
        "point_estimate",
        "sensitivity_median",
        "sensitivity_lower",
        "sensitivity_upper",
        "sensitivity_iqr",
        "median_absolute_shift",
        "replicate_count",
    }
    if (
        len(matrix) < 2
        or not all(isinstance(label, str) and label.strip() for label in matrix.index)
        or matrix.index.has_duplicates
        or matrix.columns.tolist() != matrix.index.tolist()
        or matrix.index.isna().any()
        or not required <= set(stability)
    ):
        raise ValueError("stability matrix or columns are incomplete")
    values = _numeric(matrix).to_numpy(dtype=float)
    if (
        not np.isfinite(values).all()
        or (values < 0).any()
        or (values > 2).any()
        or not np.allclose(values, values.T, rtol=0, atol=1e-12)
        or not np.allclose(np.diag(values), 0, rtol=0, atol=1e-12)
    ):
        raise ValueError("stability requires a valid symmetric divergence matrix")
    pairs = []
    for row in stability.itertuples(index=False):
        if pd.isna(row.species_a) or pd.isna(row.species_b) or row.species_a == row.species_b:
            raise ValueError("stability pairs must contain two declared species")
        pairs.append(frozenset((row.species_a, row.species_b)))
    expected_pairs = {frozenset(pair) for pair in combinations(matrix.index, 2)}
    if len(pairs) != len(expected_pairs) or len(set(pairs)) != len(pairs) or set(pairs) != expected_pairs:
        raise ValueError("stability rows do not cover the unique species pairs")
    numeric_fields = sorted(required - {"species_a", "species_b"})
    numeric = _numeric(stability[numeric_fields])
    array = numeric.to_numpy(dtype=float)
    if not np.isfinite(array).all() or (array < 0).any():
        raise ValueError("stability values must be finite and nonnegative")
    bounded = numeric.drop(columns="replicate_count")
    if (
        (bounded.to_numpy(dtype=float) > 2).any()
        or (numeric["sensitivity_lower"] > numeric["sensitivity_upper"]).any()
        or (numeric["sensitivity_median"] < numeric["sensitivity_lower"]).any()
        or (numeric["sensitivity_median"] > numeric["sensitivity_upper"]).any()
    ):
        raise ValueError("stability sensitivity values are outside their bounds")
    replicates = numeric["replicate_count"].to_numpy(dtype=float)
    if (replicates < 20).any() or (replicates != np.floor(replicates)).any() or len(set(replicates)) != 1:
        raise ValueError("stability requires a common integer replicate count of at least 20")
    for row, point in zip(stability.itertuples(index=False), numeric["point_estimate"]):
        if not np.isclose(point, matrix.loc[row.species_a, row.species_b], rtol=0, atol=1e-12):
            raise ValueError("stability point estimate differs from the divergence matrix")
