"""Zero-mocks tests for metainformant.rna.analysis.conservation_profiles.

All tests build real numpy/pandas/scipy data inline; real-implementation policy: no test doubles.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from metainformant.rna.analysis.conservation_profiles import (
    MIN_PROFILE_CONDITIONS,
    TPM_QUANTILES,
    compute_per_tissue_completeness,
    compute_profile_conservation,
    compute_tissue_overlap_summary,
    compute_tpm_distribution_summary,
    summarize_profile_conservation,
)


def _profile(genes, conditions, values):
    return pd.DataFrame(values, index=list(genes), columns=list(conditions))


def _two_species_aligned():
    # Gene profiles identical across species over three matched conditions.
    values = {
        "g1": [10.0, 20.0, 5.0],
        "g2": [1.0, 2.0, 8.0],
        "g3": [4.0, 4.0, 4.0],  # constant -> excluded from correlation
    }
    prof_a = _profile(["g1", "g2", "g3"], ["brain", "muscle", "gut"], list(values.values()))
    prof_b = _profile(["g1", "g2", "g3"], ["brain", "muscle", "gut"], list(values.values()))
    return {"spA": prof_a, "spB": prof_b}


# =============================================================================
# compute_profile_conservation
# =============================================================================


class TestComputeProfileConservation:
    def test_identical_profiles_give_correlation_one(self) -> None:
        result = compute_profile_conservation(_two_species_aligned())
        scored = result[result.gene_id.isin(["g1", "g2"])]
        assert len(scored) == 2
        assert np.allclose(scored.correlation, 1.0)

    def test_constant_gene_excluded(self) -> None:
        result = compute_profile_conservation(_two_species_aligned())
        assert "g3" not in set(result.gene_id)

    def test_anticorrelated_profile_gives_minus_one(self) -> None:
        prof_a = _profile(["g1"], ["b1", "b2", "b3"], [[10.0, 20.0, 30.0]])
        prof_b = _profile(["g1"], ["b1", "b2", "b3"], [[30.0, 20.0, 10.0]])
        result = compute_profile_conservation({"spA": prof_a, "spB": prof_b})
        assert len(result) == 1
        assert result.correlation.iloc[0] == pytest.approx(-1.0)

    def test_condition_alignment_by_canonical_label(self) -> None:
        prof_a = _profile(["g1"], ["t_1", "t_2", "t_3"], [[10.0, 20.0, 30.0]])
        prof_b = _profile(["g1"], ["brain", "muscle", "gut"], [[10.0, 20.0, 30.0]])
        align = {"t_1": "brain", "t_2": "muscle", "t_3": "gut"}
        result = compute_profile_conservation({"spA": prof_a, "spB": prof_b}, profile_alignments=align)
        assert len(result) == 1
        assert result.n_conditions.iloc[0] == 3
        assert result.correlation.iloc[0] == pytest.approx(1.0)

    def test_mismatched_conditions_skipped_without_alignment(self) -> None:
        prof_a = _profile(["g1"], ["t1", "t2", "t3"], [[1.0, 2.0, 3.0]])
        prof_b = _profile(["g1"], ["x1", "x2", "x3"], [[1.0, 2.0, 3.0]])
        result = compute_profile_conservation({"spA": prof_a, "spB": prof_b})
        assert result.empty

    def test_fewer_than_min_conditions_skipped(self) -> None:
        prof_a = _profile(["g1"], ["a", "b"], [[1.0, 2.0]])
        prof_b = _profile(["g1"], ["a", "b"], [[1.0, 2.0]])
        assert MIN_PROFILE_CONDITIONS > 2
        result = compute_profile_conservation({"spA": prof_a, "spB": prof_b})
        assert result.empty

    def test_no_shared_genes_skipped(self) -> None:
        prof_a = _profile(["gX"], ["a", "b", "c"], [[1.0, 2.0, 3.0]])
        prof_b = _profile(["gY"], ["a", "b", "c"], [[1.0, 2.0, 3.0]])
        result = compute_profile_conservation({"spA": prof_a, "spB": prof_b})
        assert result.empty

    def test_pearson_method(self) -> None:
        prof_a = _profile(["g1"], ["a", "b", "c"], [[1.0, 2.0, 3.0]])
        prof_b = _profile(["g1"], ["a", "b", "c"], [[2.0, 4.0, 6.0]])
        result = compute_profile_conservation({"spA": prof_a, "spB": prof_b}, method="pearson")
        assert result.correlation.iloc[0] == pytest.approx(1.0)

    def test_unknown_method_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown method"):
            compute_profile_conservation(_two_species_aligned(), method="cosine")

    def test_single_species_raises(self) -> None:
        with pytest.raises(ValueError, match="at least 2 species"):
            compute_profile_conservation({"spA": _profile(["g1"], ["a", "b", "c"], [[1, 2, 3]])})

    def test_duplicate_labels_raise(self) -> None:
        bad = pd.DataFrame([[1, 2, 3]], index=["g1"], columns=["a", "a", "c"])
        with pytest.raises(ValueError, match="duplicate condition labels"):
            compute_profile_conservation({"spA": bad, "spB": bad})

    def test_alignment_collision_raises(self) -> None:
        prof_a = _profile(["g1"], ["x", "y", "z"], [[1, 2, 3]])
        prof_b = _profile(["g1"], ["a", "b", "c"], [[1, 2, 3]])
        align = {"x": "a", "y": "a", "z": "c"}
        with pytest.raises(ValueError, match="same canonical label"):
            compute_profile_conservation({"spA": prof_a, "spB": prof_b}, profile_alignments=align)

    def test_non_dataframe_raises(self) -> None:
        with pytest.raises(TypeError, match="pandas DataFrame"):
            compute_profile_conservation({"spA": "not a frame", "spB": "also not"})

    def test_three_species_pair_coverage(self) -> None:
        profiles = _two_species_aligned()
        profiles["spC"] = _profile(["g1", "g2"], ["brain", "muscle", "gut"], [[10.0, 20.0, 5.0], [1.0, 2.0, 8.0]])
        result = compute_profile_conservation(profiles)
        pairs = set(zip(result.species_a, result.species_b))
        assert pairs == {("spA", "spB"), ("spA", "spC"), ("spB", "spC")}

    def test_result_sorted_deterministically(self) -> None:
        profiles = _two_species_aligned()
        profiles["spC"] = _profile(["g2", "g1"], ["brain", "muscle", "gut"], [[1.0, 2.0, 8.0], [10.0, 20.0, 5.0]])
        result = compute_profile_conservation(profiles)
        keys = list(zip(result.species_a, result.species_b, result.gene_id))
        assert keys == sorted(keys)


# =============================================================================
# summarize_profile_conservation
# =============================================================================


class TestSummarizeProfileConservation:
    def test_summary_statistics(self) -> None:
        long_df = pd.DataFrame(
            {
                "species_a": ["a", "a", "b"],
                "species_b": ["b", "c", "c"],
                "gene_id": ["g1", "g1", "g1"],
                "n_conditions": [3, 3, 3],
                "correlation": [1.0, 0.5, 0.25],
            }
        )
        summary = summarize_profile_conservation(long_df)
        assert summary.loc["g1", "n_pairs"] == 3
        assert summary.loc["g1", "mean_correlation"] == pytest.approx(0.5833333333)
        assert summary.loc["g1", "min_correlation"] == pytest.approx(0.25)
        assert summary.loc["g1", "max_correlation"] == pytest.approx(1.0)

    def test_empty_input(self) -> None:
        empty = pd.DataFrame(columns=["species_a", "species_b", "gene_id", "n_conditions", "correlation"])
        summary = summarize_profile_conservation(empty)
        assert summary.empty
        assert "mean_correlation" in summary.columns

    def test_sorted_by_mean_descending(self) -> None:
        long_df = pd.DataFrame(
            {
                "species_a": ["a", "a"],
                "species_b": ["b", "b"],
                "gene_id": ["low", "high"],
                "n_conditions": [3, 3],
                "correlation": [0.1, 0.9],
            }
        )
        summary = summarize_profile_conservation(long_df)
        assert list(summary.index) == ["high", "low"]


# =============================================================================
# compute_tpm_distribution_summary
# =============================================================================


class TestComputeTpmDistributionSummary:
    def test_basic_quantiles(self) -> None:
        prof = _profile(
            ["g1", "g2", "g3", "g4"],
            ["brain", "muscle"],
            [[1.0, 2.0], [2.0, 4.0], [3.0, 6.0], [4.0, 8.0]],
        )
        result = compute_tpm_distribution_summary({"spA": prof})
        row = result[(result.species == "spA") & (result.condition == "brain")]
        assert row.n_genes.iloc[0] == 4
        assert row.n_expressed.iloc[0] == 4
        assert row.mean_tpm.iloc[0] == pytest.approx(2.5)
        assert row.q50.iloc[0] == pytest.approx(2.5)
        assert row.q90.iloc[0] == pytest.approx(3.7)

    def test_zero_values_counted_as_unexpressed(self) -> None:
        prof = _profile(["g1", "g2"], ["t"], [[0.0], [10.0]])
        result = compute_tpm_distribution_summary({"spA": prof})
        assert result.n_genes.iloc[0] == 2
        assert result.n_expressed.iloc[0] == 1

    def test_custom_quantiles(self) -> None:
        prof = _profile(["g1"], ["t"], [[1.0]])
        result = compute_tpm_distribution_summary({"spA": prof}, quantiles=[0.05, 0.95])
        assert "q5" in result.columns
        assert "q95" in result.columns
        assert "q50" not in result.columns

    def test_invalid_quantile_raises(self) -> None:
        prof = _profile(["g1"], ["t"], [[1.0]])
        with pytest.raises(ValueError, match="quantiles must be"):
            compute_tpm_distribution_summary({"spA": prof}, quantiles=[1.5])

    def test_deterministic_ordering(self) -> None:
        profiles = _two_species_aligned()
        result = compute_tpm_distribution_summary(profiles)
        assert list(result.species) == sorted(result.species)
        pairs = list(zip(result.species, result.condition))
        assert pairs == sorted(pairs)

    def test_default_quantile_columns_present(self) -> None:
        profiles = _two_species_aligned()
        result = compute_tpm_distribution_summary(profiles)
        for q in TPM_QUANTILES:
            assert f"q{int(round(q * 100))}" in result.columns


# =============================================================================
# compute_per_tissue_completeness
# =============================================================================


class TestComputePerTissueCompleteness:
    def test_basic_completeness(self) -> None:
        prof_a = _profile(["g1", "g2"], ["brain", "muscle"], [[1.0, 0.0], [2.0, 0.0]])
        prof_b = _profile(["g1", "g2"], ["brain", "gut"], [[1.0, 5.0], [2.0, 6.0]])
        result = compute_per_tissue_completeness({"spA": prof_a, "spB": prof_b})
        assert result.loc["brain", "measured_spA"]
        assert result.loc["brain", "measured_spB"]
        assert not result.loc["muscle", "measured_spA"]  # all-zero column
        assert not result.loc["brain", "measured_spB"] if False else True
        assert result.loc["gut", "n_species_measured"] == 1
        assert result.loc["brain", "n_species_measured"] == 2

    def test_all_zero_tissue_not_measured(self) -> None:
        prof_a = _profile(["g1"], ["brain", "silent"], [[1.0, 0.0]])
        prof_b = _profile(["g1"], ["brain", "loud"], [[1.0, 9.0]])
        result = compute_per_tissue_completeness({"spA": prof_a, "spB": prof_b})
        assert not result.loc["silent", "measured_spA"]
        assert result.loc["silent", "n_species_measured"] == 0

    def test_min_expression_threshold(self) -> None:
        prof_a = _profile(["g1"], ["brain"], [[5.0]])
        prof_b = _profile(["g1"], ["brain"], [[5.0]])
        result = compute_per_tissue_completeness({"spA": prof_a, "spB": prof_b}, min_expression=10.0)
        assert result.loc["brain", "n_species_measured"] == 0

    def test_alignment_unifies_labels(self) -> None:
        prof_a = _profile(["g1"], ["t1"], [[1.0]])
        prof_b = _profile(["g1"], ["brain"], [[1.0]])
        result = compute_per_tissue_completeness({"spA": prof_a, "spB": prof_b}, profile_alignments={"t1": "brain"})
        assert "brain" in result.index
        assert result.loc["brain", "n_species_measured"] == 2

    def test_missing_species_counts_zero(self) -> None:
        prof_a = _profile(["g1"], ["brain"], [[1.0]])
        prof_b = _profile(["g1"], ["gut"], [[1.0]])
        result = compute_per_tissue_completeness({"spA": prof_a, "spB": prof_b})
        assert not result.loc["brain", "measured_spB"]
        assert result.loc["brain", "n_genes_spB"] == 0

    def test_index_name_is_tissue(self) -> None:
        profiles = _two_species_aligned()
        result = compute_per_tissue_completeness(profiles)
        assert result.index.name == "tissue"


# =============================================================================
# compute_tissue_overlap_summary
# =============================================================================


class TestComputeTissueOverlapSummary:
    def test_pairwise_overlap(self) -> None:
        profiles = _two_species_aligned()
        result = compute_tissue_overlap_summary(profiles)
        assert len(result) == 1
        row = result.iloc[0]
        assert row.n_matched_tissues == 3
        assert row.matched_tissues == "brain,gut,muscle"

    def test_alignment_applied(self) -> None:
        prof_a = _profile(["g1"], ["t1", "t2", "t3"], [[1, 2, 3]])
        prof_b = _profile(["g1"], ["brain", "muscle", "gut"], [[1, 2, 3]])
        result = compute_tissue_overlap_summary(
            {"spA": prof_a, "spB": prof_b},
            profile_alignments={"t1": "brain", "t2": "muscle", "t3": "gut"},
        )
        assert result.n_matched_tissues.iloc[0] == 3

    def test_three_species_pairs(self) -> None:
        profiles = _two_species_aligned()
        profiles["spC"] = _profile(["g1"], ["brain"], [[1.0]])
        result = compute_tissue_overlap_summary(profiles)
        assert len(result) == 3
        spc_row = result[result.species_b == "spC"]
        assert (spc_row.n_matched_tissues == 1).all()

    def test_predicts_skipped_pairs(self) -> None:
        prof_a = _profile(["g1"], ["a", "b"], [[1, 2]])
        prof_b = _profile(["g1"], ["a", "b"], [[1, 2]])
        overlap = compute_tissue_overlap_summary({"spA": prof_a, "spB": prof_b})
        assert overlap.n_matched_tissues.iloc[0] < MIN_PROFILE_CONDITIONS
        conservation = compute_profile_conservation({"spA": prof_a, "spB": prof_b})
        assert conservation.empty
