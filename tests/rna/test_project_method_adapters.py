"""Independent numerical and real-file controls for shared project methods."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from statsmodels.stats.proportion import proportion_confint

from metainformant.rna.analysis.counting_statistics import wilson_interval
from metainformant.rna.analysis.expression_io import (
    compute_profile_quality,
    load_expression_profile,
    validate_expression_matrix,
)
from metainformant.rna.analysis.ortholog_diagnostics import (
    audit_gene_drop_reasons,
    classify_orthogroup_cardinality,
)
from metainformant.rna.analysis.ortholog_profiles import (
    compute_ortholog_profile_divergence,
)
from metainformant.rna.analysis.statistics_io import read_analysis_provenance


def test_compressed_matrix_profile_and_chunked_dimensions(tmp_path: Path) -> None:
    path = tmp_path / "expression.tsv.gz"
    with gzip.open(path, "wt") as f:
        f.write("feature\tA\tB\ng1\t1\t3\ng2\t4\t0\n")
    assert validate_expression_matrix(path, min_samples=2, chunk_size=1) == (2, 2)
    assert load_expression_profile(str(path)).to_dict() == {"g1": 2.0, "g2": 2.0}


@pytest.mark.parametrize("rows", ["g1\t1\t1\ng1\t2\t2\n", "g1\t1\t0\ng2\t1\t0\n"])
def test_matrix_rejects_cross_chunk_duplicate_or_empty_sample(tmp_path: Path, rows: str) -> None:
    path = tmp_path / "matrix.tsv"
    path.write_text("feature\tA\tB\n" + rows)
    with pytest.raises(ValueError):
        validate_expression_matrix(path, chunk_size=1)


def test_quality_preserves_nonfinite_counts_instead_of_imputing_zero() -> None:
    table = compute_profile_quality({"species": pd.Series([1, 0, np.nan, np.inf])})
    assert table.loc[0, ["finite_features", "zero_features", "nonfinite_features"]].tolist() == [2, 1, 2]
    assert table.loc[0, "positive_fraction_finite"] == 0.5


def test_orthology_audit_reconciles_each_mapping_failure(tmp_path: Path) -> None:
    path = tmp_path / "orthogroups.tsv"
    path.write_text("orthogroup\t1\t2\nog1\tg1,g2,g3,g4\th1\n")
    taxonomy = {"1": "A", "2": "B"}
    classes = classify_orthogroup_cardinality(path, taxonomy)
    assert classes.loc[0, "cardinality_class"] == "one_to_many"
    audit = audit_gene_drop_reasons(
        path,
        {"g2": "p2", "g3": "p3", "g4": "p4", "h1": "q1"},
        {"p3": "r3", "p4": "r4", "q1": "s1"},
        {"A": {"r4": "t4"}, "B": {"s1": "u1"}},
        taxonomy,
    ).set_index("species")
    assert audit.loc[
        "A",
        [
            "genes_seen",
            "genes_retained",
            "dropped_no_protein",
            "dropped_no_rna_map",
            "dropped_no_transcript",
        ],
    ].tolist() == [4, 1, 1, 1, 1]


def test_ortholog_mean_profile_distance_and_explicit_first_transcript() -> None:
    mapping = pd.DataFrame(
        {
            "A": ["a1,unused", "a2", "a3"],
            "B": ["b1", "b2", "b3"],
            "missing": ["c1", "c2", "c3"],
        },
        index=["g1", "g2", "g3"],
    )
    profiles = {
        "A": pd.Series({"a1": 1.0, "unused": 999.0, "a2": 2.0, "a3": 3.0}),
        "B": pd.Series({"b1": 3.0, "b2": 2.0, "b3": 1.0}),
    }
    distance, overlap = compute_ortholog_profile_divergence(mapping, profiles, min_shared_orthologs=3)
    assert distance.loc["A", "B"] == pytest.approx(2.0)
    assert overlap.loc["A", "B"] == 3
    assert np.isnan(distance.loc["A", "missing"])
    assert overlap.loc["A", "missing"] == 0


def test_constant_or_insufficient_ortholog_profiles_are_unavailable() -> None:
    mapping = pd.DataFrame({"A": ["a", "b"], "B": ["a", "b"]}, index=["g1", "g2"])
    profiles = {
        "A": pd.Series({"a": 1.0, "b": 1.0}),
        "B": pd.Series({"a": 2.0, "b": 3.0}),
    }
    distances, _ = compute_ortholog_profile_divergence(mapping, profiles, min_shared_orthologs=2)
    assert np.isnan(distances.loc["A", "B"])
    profiles["A"]["a"] = np.inf
    with pytest.raises(ValueError):
        compute_ortholog_profile_divergence(mapping, profiles, min_shared_orthologs=2)


def test_wilson_interval_matches_independent_implementation() -> None:
    assert wilson_interval(7, 20) == pytest.approx(proportion_confint(7, 20, method="wilson"))
    assert wilson_interval(np.int64(7), np.int64(20)) == pytest.approx(proportion_confint(7, 20, method="wilson"))
    with pytest.raises(ValueError):
        wilson_interval(True, 20)
    assert wilson_interval(1, 0) is None


@pytest.mark.parametrize("z", [float("nan"), float("inf"), 0, -1])
def test_wilson_interval_rejects_invalid_uncertainty_parameter(z: float) -> None:
    with pytest.raises(ValueError):
        wilson_interval(7, 20, z)


@pytest.mark.parametrize("entries", [False, 0, ""])
def test_contract_reader_rejects_falsey_nonlist_sensitivities(tmp_path: Path, entries: object) -> None:
    path = tmp_path / "contract.json"
    path.write_text(
        json.dumps(
            {
                "analysis_id": "control",
                "estimand": "profile",
                "replicate_unit": "species",
                "random_seed": 42,
                "resampling_count": 100,
                "null_model": "permutation",
                "sensitivity_analyses": entries,
            }
        )
    )
    with pytest.raises(ValueError, match="JSON list"):
        read_analysis_provenance(path)
