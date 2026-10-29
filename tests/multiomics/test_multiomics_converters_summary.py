"""Converters choose their transform from the data's scale, not from its zeros; summaries keep real feature names.

Regression tests for ``metainformant.multiomics.analysis``:

* ``from_rna_expression(normalize=True)`` log-transformed (with ``log1p``) only when *every* value was strictly
  positive. ``log1p`` is defined at zero, so a single zero count (the usual case in RNA-seq) or a single missing
  value switched the transform off and the matrix was z-scored on the raw scale instead, without a word;
* ``from_metabolomics(normalize=True)`` skipped its log the same way when one value was missing, and needed
  scikit-learn although a pandas z-score does the same;
* ``from_epigenome_data`` promised "processed" data and had ``pass`` branches per data type: it returns a copy;
* ``MultiOmicsDescriptiveSummary.to_dict`` listed the first features' means under invented names
  (``feature_0``, ``feature_1``, ...) instead of the block's columns.

Real computation; no test doubles (real-implementation policy).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from metainformant.multiomics.analysis.descriptive_summary import (
    MultiOmicsDescriptiveSummary,
    OmicsBlockSummary,
    multi_omics_descriptive_summary,
    summarize_omics_block,
)
from metainformant.multiomics.analysis.integration import from_epigenome_data, from_metabolomics, from_rna_expression


def _zscore(frame: pd.DataFrame) -> pd.DataFrame:
    """Column z-scores with the population standard deviation, ignoring NaN (the converters' convention)."""
    return (frame - frame.mean(axis=0)) / frame.std(axis=0, ddof=0)


def _counts(seed: int = 0, shape: tuple[int, int] = (12, 5)) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        rng.poisson(6.0, size=shape).astype(float),
        index=[f"S{i}" for i in range(shape[0])],
        columns=[f"g{j}" for j in range(shape[1])],
    )


class TestRnaExpressionTransform:
    def test_counts_with_zeros_are_log_transformed(self) -> None:
        counts = _counts()
        counts.iloc[0, 0] = 0.0
        counts.iloc[3, 2] = 0.0

        result = from_rna_expression(counts, normalize=True)

        pd.testing.assert_frame_equal(result, _zscore(np.log1p(counts)), atol=1e-12, check_exact=False)

    def test_a_missing_value_does_not_switch_the_log_off(self) -> None:
        counts = _counts(seed=1)
        counts.iloc[2, 1] = np.nan

        result = from_rna_expression(counts, normalize=True)

        pd.testing.assert_frame_equal(result, _zscore(np.log1p(counts)), atol=1e-12, check_exact=False)
        assert result.isna().sum().sum() == 1  # the gap stays a gap

    def test_the_same_matrix_is_treated_alike_with_and_without_a_zero(self) -> None:
        # the first and the second matrix differ in one count only; their transforms must not differ in kind
        positive = _counts(seed=2) + 1.0
        with_zero = positive.copy()
        with_zero.iloc[0, 0] = 0.0

        difference = (from_rna_expression(positive) - from_rna_expression(with_zero)).abs().max().max()

        assert difference < 5.0  # on the raw scale the z-scores of counts up to ~15 would have been another matter
        pd.testing.assert_frame_equal(
            from_rna_expression(with_zero), _zscore(np.log1p(with_zero)), atol=1e-12, check_exact=False
        )

    def test_data_with_a_negative_value_is_only_z_scored(self) -> None:
        centred = _counts(seed=3) - 6.0

        result = from_rna_expression(centred, normalize=True)

        pd.testing.assert_frame_equal(result, _zscore(centred), atol=1e-12, check_exact=False)

    def test_normalize_false_returns_the_numbers_untouched(self) -> None:
        counts = _counts(seed=4)

        pd.testing.assert_frame_equal(from_rna_expression(counts, normalize=False), counts)

    def test_the_docstring_states_the_rule(self) -> None:
        doc = " ".join((from_rna_expression.__doc__ or "").split())

        assert "log1p" in doc and "non-negative" in doc


class TestMetabolomicsTransform:
    @staticmethod
    def _intensities(seed: int = 5) -> pd.DataFrame:
        rng = np.random.default_rng(seed)
        return pd.DataFrame(rng.lognormal(mean=4.0, sigma=1.0, size=(14, 4)), columns=[f"m{j}" for j in range(4)])

    def test_positive_data_is_logged_then_z_scored(self) -> None:
        data = self._intensities()

        result = from_metabolomics(data, normalize=True)

        pd.testing.assert_frame_equal(result, _zscore(np.log(data)), atol=1e-12, check_exact=False)

    def test_a_missing_value_does_not_switch_the_log_off(self) -> None:
        data = self._intensities(seed=6)
        data.iloc[4, 2] = np.nan

        result = from_metabolomics(data, normalize=True)

        pd.testing.assert_frame_equal(result, _zscore(np.log(data)), atol=1e-12, check_exact=False)
        assert np.isnan(result.iloc[4, 2]) and result.isna().sum().sum() == 1

    def test_data_with_a_zero_cannot_be_logged_and_is_only_z_scored(self) -> None:
        data = self._intensities(seed=7)
        data.iloc[0, 0] = 0.0

        result = from_metabolomics(data, normalize=True)

        pd.testing.assert_frame_equal(result, _zscore(data), atol=1e-12, check_exact=False)


class TestEpigenomeConverter:
    @pytest.mark.parametrize("data_type", ["methylation", "chipseq", "atacseq", "anything else"])
    def test_it_returns_an_unchanged_copy_for_every_data_type(self, data_type: str) -> None:
        data = pd.DataFrame({"cg1": [0.1, 0.9], "cg2": [0.4, np.nan]}, index=["S1", "S2"])

        result = from_epigenome_data(data, data_type=data_type)

        pd.testing.assert_frame_equal(result, data)
        assert result is not data

    def test_the_docstring_says_nothing_is_transformed(self) -> None:
        doc = " ".join((from_epigenome_data.__doc__ or "").split())

        assert "No transformation is applied" in doc


class TestSummaryKeepsFeatureNames:
    @staticmethod
    def _blocks() -> dict[str, pd.DataFrame]:
        rng = np.random.default_rng(8)
        samples = [f"S{i}" for i in range(10)]
        return {
            "rna": pd.DataFrame(rng.normal(size=(10, 3)), index=samples, columns=["TP53", "BRCA1", "MYC"]),
            "protein": pd.DataFrame(rng.normal(size=(10, 2)), index=samples, columns=[101, 102]),
        }

    def test_the_block_summary_carries_the_column_labels_as_text(self) -> None:
        block = summarize_omics_block(self._blocks()["rna"], "rna")
        numeric = summarize_omics_block(self._blocks()["protein"], "protein")

        assert block.feature_names == ["TP53", "BRCA1", "MYC"]
        assert numeric.feature_names == ["101", "102"]

    def test_to_dict_reports_the_means_under_the_real_feature_names(self) -> None:
        blocks = self._blocks()
        summary = multi_omics_descriptive_summary(blocks)

        payload = summary.to_dict()

        rna = next(b for b in payload["blocks"] if b["name"] == "rna")
        assert list(rna["mean_head"]) == ["TP53", "BRCA1", "MYC"]  # not feature_0, feature_1, feature_2
        for name, value in rna["mean_head"].items():
            assert value == pytest.approx(round(float(blocks["rna"][name].mean()), 6))
        protein = next(b for b in payload["blocks"] if b["name"] == "protein")
        assert list(protein["mean_head"]) == ["101", "102"]

    def test_only_the_first_ten_features_are_listed_and_they_are_the_first_ten_columns(self) -> None:
        rng = np.random.default_rng(9)
        wide = pd.DataFrame(rng.normal(size=(6, 14)), columns=[f"gene{j:02d}" for j in range(14)])

        payload = multi_omics_descriptive_summary({"rna": wide}).to_dict()

        assert list(payload["blocks"][0]["mean_head"]) == [f"gene{j:02d}" for j in range(10)]

    def test_a_missing_mean_is_reported_as_none_under_its_name(self) -> None:
        data = pd.DataFrame({"a": [1.0, 2.0, 3.0], "b": [np.nan, np.nan, np.nan]})

        payload = multi_omics_descriptive_summary({"x": data}).to_dict()

        assert payload["blocks"][0]["mean_head"] == {"a": 2.0, "b": None}

    def test_a_summary_built_by_hand_without_names_falls_back_to_positional_names(self) -> None:
        block = OmicsBlockSummary(
            name="manual",
            n_samples=3,
            n_features=2,
            mean=np.array([1.0, 2.0]),
            std=np.array([0.5, 0.5]),
            missing_fraction=0.0,
            quantiles={},
        )

        payload = MultiOmicsDescriptiveSummary(blocks=[block], block_correlations=[]).to_dict()

        assert payload["blocks"][0]["mean_head"] == {"feature_0": 1.0, "feature_1": 2.0}
