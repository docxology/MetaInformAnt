"""Multivariate phenotype statistics.

Reusable helpers for building phenotype trait matrices, distance matrices,
ordination coordinates, multivariate group tests, dispersion summaries, and
axis loadings. Domain-specific code should prepare biologically meaningful
traits and grouping columns before calling these functions.
"""

from __future__ import annotations

from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

from metainformant.ecology.analysis.indicators import permanova
from metainformant.ecology.analysis.ordination import pcoa

try:
    from statsmodels.multivariate.manova import MANOVA

    HAS_STATSMODELS_MANOVA = True
except Exception:  # pragma: no cover - depends on optional statsmodels import path
    MANOVA = None
    HAS_STATSMODELS_MANOVA = False


def _as_list(values: Iterable[str] | None) -> list[str]:
    return [str(value) for value in values or []]


def build_standardized_trait_matrix(
    data: pd.DataFrame,
    traits: Sequence[str],
    *,
    id_columns: Sequence[str] = (),
    min_nonmissing: int = 3,
    min_unique: int = 2,
    drop_incomplete: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build a z-scored numeric phenotype matrix plus trait metadata.

    Traits are retained only when they are present, numeric enough, variable,
    and have enough non-missing observations. Standardization uses population
    standard deviation (``ddof=0``) for deterministic matrix construction.
    """
    id_columns = _as_list(id_columns)
    rows: list[dict[str, Any]] = []
    matrix = data[[column for column in id_columns if column in data.columns]].copy()
    included: list[str] = []
    for trait in traits:
        trait = str(trait)
        if trait not in data.columns:
            rows.append(
                {
                    "trait": trait,
                    "included": False,
                    "exclusion_reason": "missing_column",
                    "n_nonmissing": 0,
                    "n_unique": 0,
                    "mean": np.nan,
                    "sd": np.nan,
                }
            )
            continue
        numeric = pd.to_numeric(data[trait], errors="coerce")
        nonmissing = int(numeric.notna().sum())
        unique = int(numeric.nunique(dropna=True))
        mean = float(numeric.mean()) if nonmissing else np.nan
        sd = float(numeric.std(ddof=0)) if nonmissing else np.nan
        reason = ""
        if nonmissing < min_nonmissing:
            reason = "insufficient_nonmissing"
        elif unique < min_unique:
            reason = "constant_or_singleton"
        elif not np.isfinite(sd) or sd <= 0:
            reason = "zero_or_invalid_sd"
        included_flag = reason == ""
        if included_flag:
            matrix[trait] = (numeric - mean) / sd
            included.append(trait)
        rows.append(
            {
                "trait": trait,
                "included": included_flag,
                "exclusion_reason": reason,
                "n_nonmissing": nonmissing,
                "n_unique": unique,
                "mean": mean,
                "sd": sd,
            }
        )
    if drop_incomplete and included:
        matrix = matrix.dropna(subset=included).reset_index(drop=True)
    return matrix, pd.DataFrame(rows)


def phenotype_distance_matrix(
    matrix: pd.DataFrame,
    trait_columns: Sequence[str],
    *,
    id_column: str | None = None,
) -> pd.DataFrame:
    """Compute a deterministic standardized Euclidean phenotype distance matrix."""
    trait_columns = [trait for trait in _as_list(trait_columns) if trait in matrix.columns]
    if not trait_columns or len(matrix) < 2:
        labels = matrix[id_column].astype(str).tolist() if id_column and id_column in matrix.columns else matrix.index.astype(str).tolist()
        return pd.DataFrame(np.zeros((len(matrix), len(matrix))), index=labels, columns=labels)
    values = matrix[trait_columns].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(values).all(axis=1)
    values = values[valid]
    labels_source = matrix.loc[valid, id_column].astype(str) if id_column and id_column in matrix.columns else matrix.index[valid].astype(str)
    labels = labels_source.tolist()
    if len(values) == 0:
        return pd.DataFrame()
    diff = values[:, None, :] - values[None, :, :]
    distances = np.sqrt(np.sum(diff**2, axis=2) / len(trait_columns))
    return pd.DataFrame(distances, index=labels, columns=labels)


def run_pcoa(distance_matrix: pd.DataFrame, *, n_components: int = 2) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run PCoA and return coordinate and variance tables."""
    if distance_matrix.empty:
        return pd.DataFrame(), pd.DataFrame()
    actual_components = min(max(1, int(n_components)), max(1, len(distance_matrix)))
    result = pcoa(distance_matrix.to_numpy(dtype=float).tolist(), n_components=actual_components)
    labels = distance_matrix.index.astype(str).tolist()
    coordinates = pd.DataFrame(result["coordinates"], columns=[f"PCoA{i + 1}" for i in range(len(result["coordinates"][0]))])
    coordinates.insert(0, "entity_id", labels)
    variance = pd.DataFrame(
        {
            "axis": [f"PCoA{i + 1}" for i in range(len(result.get("variance_explained", [])))],
            "eigenvalue": result.get("eigenvalues", []),
            "variance_explained": result.get("variance_explained", []),
        }
    )
    return coordinates, variance


def _aligned_labels(
    distance_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    term: str,
    *,
    id_column: str,
    min_group_size: int,
) -> tuple[pd.DataFrame, list[str], str]:
    if term not in metadata.columns:
        return pd.DataFrame(), [], "missing_term"
    if id_column not in metadata.columns:
        return pd.DataFrame(), [], "missing_id_column"
    labels = distance_matrix.index.astype(str).tolist()
    meta = metadata.copy()
    meta[id_column] = meta[id_column].astype(str)
    meta = meta.set_index(id_column).reindex(labels)
    term_values = meta[term].astype("string")
    valid = term_values.notna() & term_values.ne("")
    counts = term_values[valid].value_counts()
    allowed = set(counts[counts >= min_group_size].index.astype(str))
    keep = valid & term_values.isin(allowed)
    kept_labels = [label for label, ok in zip(labels, keep.tolist()) if ok]
    if len(kept_labels) < 3:
        return pd.DataFrame(), [], "insufficient_observations"
    groups = term_values.loc[kept_labels].astype(str).tolist()
    if len(set(groups)) < 2:
        return pd.DataFrame(), [], "insufficient_groups"
    return distance_matrix.loc[kept_labels, kept_labels], groups, "tested"


def run_permanova_terms(
    distance_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    terms: Sequence[str],
    *,
    id_column: str = "entity_id",
    n_permutations: int = 499,
    seed: int = 20260618,
    min_group_size: int = 2,
) -> pd.DataFrame:
    """Run one PERMANOVA row per grouping term."""
    rows: list[dict[str, Any]] = []
    for offset, term in enumerate(_as_list(terms)):
        subset, groups, status = _aligned_labels(distance_matrix, metadata, term, id_column=id_column, min_group_size=min_group_size)
        row: dict[str, Any] = {
            "test_family": "permanova",
            "term": term,
            "n": int(len(groups)),
            "n_groups": int(len(set(groups))) if groups else 0,
            "groups_tested": ",".join(sorted(set(groups))) if groups else "",
            "f_statistic": np.nan,
            "p_value": np.nan,
            "r_squared": np.nan,
            "n_permutations": int(n_permutations),
            "test_status": status,
        }
        if status == "tested":
            try:
                result = permanova(
                    subset.to_numpy(dtype=float).tolist(),
                    groups,
                    n_permutations=n_permutations,
                    seed=seed + offset,
                )
                row.update(
                    {
                        "f_statistic": float(result["f_statistic"]),
                        "p_value": float(result["p_value"]),
                        "r_squared": float(result["r_squared"]),
                        "test_status": "tested",
                    }
                )
            except Exception as exc:
                row["test_status"] = f"failed:{type(exc).__name__}"
        rows.append(row)
    return pd.DataFrame(rows)


def axis_trait_loadings(
    trait_matrix: pd.DataFrame,
    ordination: pd.DataFrame,
    trait_columns: Sequence[str],
    *,
    id_column: str = "entity_id",
) -> pd.DataFrame:
    """Compute Pearson correlations between traits and ordination axes."""
    trait_columns = [trait for trait in _as_list(trait_columns) if trait in trait_matrix.columns]
    axis_columns = [column for column in ordination.columns if column.startswith("PCoA")]
    if not trait_columns or not axis_columns or id_column not in trait_matrix or id_column not in ordination:
        return pd.DataFrame()
    merged = trait_matrix[[id_column, *trait_columns]].merge(ordination[[id_column, *axis_columns]], on=id_column, how="inner")
    rows: list[dict[str, Any]] = []
    for trait in trait_columns:
        values = pd.to_numeric(merged[trait], errors="coerce")
        for axis in axis_columns:
            axis_values = pd.to_numeric(merged[axis], errors="coerce")
            valid = values.notna() & axis_values.notna()
            if valid.sum() >= 3 and float(values[valid].std(ddof=0)) > 0.0 and float(axis_values[valid].std(ddof=0)) > 0.0:
                loading = float(values[valid].corr(axis_values[valid]))
            else:
                loading = np.nan
            rows.append({"trait": trait, "axis": axis, "loading": loading, "abs_loading": abs(loading) if np.isfinite(loading) else np.nan})
    return pd.DataFrame(rows)


def multivariate_dispersion(
    ordination: pd.DataFrame,
    metadata: pd.DataFrame,
    term: str,
    *,
    id_column: str = "entity_id",
) -> pd.DataFrame:
    """Summarize group centroids and distances to centroid in ordination space."""
    axis_columns = [column for column in ordination.columns if column.startswith("PCoA")]
    if not axis_columns or term not in metadata.columns or id_column not in ordination or id_column not in metadata:
        return pd.DataFrame()
    merged = ordination[[id_column, *axis_columns]].merge(metadata[[id_column, term]], on=id_column, how="inner").dropna(subset=[term])
    rows: list[dict[str, Any]] = []
    for group, sub in merged.groupby(term, sort=True):
        coords = sub[axis_columns].to_numpy(dtype=float)
        if len(coords) == 0:
            continue
        centroid = coords.mean(axis=0)
        distances = np.sqrt(((coords - centroid) ** 2).sum(axis=1))
        row: dict[str, Any] = {
            "term": term,
            "group": str(group),
            "n": int(len(sub)),
            "mean_distance_to_centroid": float(distances.mean()),
            "median_distance_to_centroid": float(np.median(distances)),
        }
        for axis, value in zip(axis_columns, centroid):
            row[f"{axis}_centroid"] = float(value)
        rows.append(row)
    return pd.DataFrame(rows)


def fit_manova_terms(
    data: pd.DataFrame,
    response_columns: Sequence[str],
    terms: Sequence[str],
    *,
    min_observations: int = 8,
    min_group_size: int = 2,
) -> pd.DataFrame:
    """Fit simple MANOVA models for each categorical term when possible."""
    responses = [column for column in _as_list(response_columns) if column in data.columns]
    rows: list[dict[str, Any]] = []
    for term in _as_list(terms):
        row: dict[str, Any] = {
            "test_family": "manova",
            "term": term,
            "responses": ",".join(responses),
            "n": 0,
            "n_groups": 0,
            "wilks_lambda": np.nan,
            "pillai_trace": np.nan,
            "f_statistic": np.nan,
            "p_value": np.nan,
            "test_status": "not_tested",
        }
        if term not in data.columns or not responses:
            row["test_status"] = "missing_term_or_responses"
            rows.append(row)
            continue
        sub = data[[term, *responses]].copy()
        for response in responses:
            sub[response] = pd.to_numeric(sub[response], errors="coerce")
        sub = sub.dropna()
        counts = sub[term].astype(str).value_counts()
        allowed = set(counts[counts >= min_group_size].index)
        sub = sub[sub[term].astype(str).isin(allowed)].copy()
        row["n"] = int(len(sub))
        row["n_groups"] = int(sub[term].nunique())
        if len(sub) < min_observations or sub[term].nunique() < 2:
            row["test_status"] = "insufficient_data"
            rows.append(row)
            continue
        if not HAS_STATSMODELS_MANOVA or MANOVA is None:
            row["test_status"] = "statsmodels_unavailable"
            rows.append(row)
            continue
        try:
            formula = " + ".join(responses) + f" ~ C({term})"
            result = MANOVA.from_formula(formula, data=sub).mv_test()
            stat = result.results[f"C({term})"]["stat"]
            row["wilks_lambda"] = float(stat.loc["Wilks' lambda", "Value"])
            row["pillai_trace"] = float(stat.loc["Pillai's trace", "Value"])
            row["f_statistic"] = float(stat.loc["Pillai's trace", "F Value"])
            row["p_value"] = float(stat.loc["Pillai's trace", "Pr > F"])
            row["test_status"] = "tested"
        except Exception as exc:
            row["test_status"] = f"failed:{type(exc).__name__}"
        rows.append(row)
    return pd.DataFrame(rows)


__all__ = [
    "axis_trait_loadings",
    "build_standardized_trait_matrix",
    "fit_manova_terms",
    "multivariate_dispersion",
    "phenotype_distance_matrix",
    "run_pcoa",
    "run_permanova_terms",
]
