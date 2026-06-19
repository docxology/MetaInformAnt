from __future__ import annotations

import numpy as np
import pandas as pd

from metainformant.phenotype.analysis.multivariate import (
    axis_trait_loadings,
    build_standardized_trait_matrix,
    fit_manova_terms,
    phenotype_distance_matrix,
    run_pcoa,
    run_permanova_terms,
)


def test_standardized_trait_matrix_records_exclusions_and_zscores() -> None:
    data = pd.DataFrame(
        {
            "entity_id": ["s1", "s2", "s3", "s4"],
            "group": ["A", "A", "B", "B"],
            "trait_a": [1.0, 2.0, 4.0, 5.0],
            "trait_b": [2.0, 4.0, 8.0, 10.0],
            "constant": [7.0, 7.0, 7.0, 7.0],
        }
    )

    matrix, metadata = build_standardized_trait_matrix(
        data,
        ["trait_a", "trait_b", "constant", "missing"],
        id_columns=["entity_id", "group"],
        min_nonmissing=3,
    )

    assert set(matrix.columns) >= {"entity_id", "group", "trait_a", "trait_b"}
    assert np.isclose(matrix["trait_a"].mean(), 0.0)
    assert np.isclose(matrix["trait_a"].std(ddof=0), 1.0)
    by_trait = metadata.set_index("trait")
    assert bool(by_trait.loc["trait_a", "included"]) is True
    assert by_trait.loc["constant", "exclusion_reason"] == "constant_or_singleton"
    assert by_trait.loc["missing", "exclusion_reason"] == "missing_column"


def test_distance_permanova_pcoa_and_loadings_are_deterministic() -> None:
    data = pd.DataFrame(
        {
            "entity_id": ["s1", "s2", "s3", "s4", "s5", "s6"],
            "group": ["A", "A", "A", "B", "B", "B"],
            "trait_a": [1.0, 1.4, 1.8, 6.0, 6.4, 6.8],
            "trait_b": [0.2, 0.3, 0.4, 2.0, 2.1, 2.2],
        }
    )
    matrix, _ = build_standardized_trait_matrix(data, ["trait_a", "trait_b"], id_columns=["entity_id", "group"])
    distance = phenotype_distance_matrix(matrix, ["trait_a", "trait_b"], id_column="entity_id")

    assert distance.shape == (6, 6)
    assert np.allclose(distance.to_numpy(), distance.to_numpy().T)
    assert np.allclose(np.diag(distance.to_numpy()), 0.0)

    first = run_permanova_terms(distance, matrix[["entity_id", "group"]], ["group"], n_permutations=49, seed=7)
    second = run_permanova_terms(distance, matrix[["entity_id", "group"]], ["group"], n_permutations=49, seed=7)
    assert first[["f_statistic", "p_value", "r_squared", "test_status"]].to_dict("records") == second[
        ["f_statistic", "p_value", "r_squared", "test_status"]
    ].to_dict("records")
    assert first.loc[0, "test_status"] == "tested"
    assert 0.0 <= first.loc[0, "p_value"] <= 1.0

    ordination, variance = run_pcoa(distance, n_components=2)
    assert len(ordination) == 6
    assert {"PCoA1", "PCoA2"}.issubset(ordination.columns)
    assert len(variance) == 2

    loadings = axis_trait_loadings(matrix, ordination, ["trait_a", "trait_b"], id_column="entity_id")
    assert set(loadings["axis"]) == {"PCoA1", "PCoA2"}
    assert loadings["loading"].notna().any()
    assert loadings["abs_loading"].max() > 0.5


def test_manova_underpowered_status_is_explicit() -> None:
    data = pd.DataFrame(
        {
            "group": ["A", "A", "B", "B"],
            "trait_a": [1.0, 2.0, 3.0, 4.0],
            "trait_b": [2.0, 3.0, 4.0, 5.0],
        }
    )

    rows = fit_manova_terms(data, ["trait_a", "trait_b"], ["group"], min_observations=8, min_group_size=2)

    assert len(rows) == 1
    assert rows.loc[0, "test_family"] == "manova"
    assert rows.loc[0, "test_status"] == "insufficient_data"
