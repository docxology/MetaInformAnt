"""Tests for metabolomics identification and quantification."""
from __future__ import annotations

import numpy as np
import pytest

from metainformant.metabolomics.analysis.identification import (
    AdductMatch,
    COMMON_ADDUCTS,
    MetaboliteMatch,
    cosine_spectral_similarity,
    differential_abundance,
    fold_change,
    identify_metabolites,
    identify_with_adducts,
    missing_value_imputation,
    normalize_intensities,
)


class TestIdentifyMetabolites:
    """Tests for m/z-based metabolite identification."""

    def test_exact_match(self) -> None:
        """Exact m/z match should yield score ~1.0."""
        db = {"glucose": 180.0634, "fructose": 180.0634}
        observed = np.array([180.0634])
        results = identify_metabolites(observed, db, ppm_tolerance=10.0)
        assert len(results) == 1
        assert len(results[0]) == 2  # both glucose and fructose match
        assert results[0][0].score == pytest.approx(1.0, abs=0.01)

    def test_no_match_outside_tolerance(self) -> None:
        """m/z outside tolerance should yield no matches."""
        db = {"glucose": 180.0634}
        observed = np.array([200.0])
        results = identify_metabolites(observed, db, ppm_tolerance=5.0)
        assert len(results) == 1
        assert len(results[0]) == 0

    def test_multiple_observations(self) -> None:
        """Multiple observed m/z values return one result list each."""
        db = {"A": 100.0, "B": 200.0, "C": 300.0}
        observed = np.array([100.0, 200.0, 999.0])
        results = identify_metabolites(observed, db, ppm_tolerance=20.0)
        assert len(results) == 3
        assert len(results[0]) >= 1  # A matches
        assert len(results[1]) >= 1  # B matches
        assert len(results[2]) == 0  # no match for 999

    def test_match_dataclass_fields(self) -> None:
        """MetaboliteMatch should have correct field values."""
        db = {"alanine": 89.0935}
        observed = np.array([89.0935])
        results = identify_metabolites(observed, db, ppm_tolerance=10.0)
        match = results[0][0]
        assert isinstance(match, MetaboliteMatch)
        assert match.matched_name == "alanine"
        assert match.delta_ppm < 1.0
        assert 0.0 <= match.score <= 1.0

    def test_ppm_tolerance_sensitivity(self) -> None:
        """Tighter tolerance should find fewer matches."""
        db = {"A": 100.000, "B": 100.001}
        observed = np.array([100.0005])
        loose = identify_metabolites(observed, db, ppm_tolerance=50.0)
        tight = identify_metabolites(observed, db, ppm_tolerance=1.0)
        assert len(loose[0]) >= len(tight[0])


class TestNormalizeIntensities:
    """Tests for intensity normalization methods."""

    def setup_method(self) -> None:
        self.data = np.array([[100.0, 200.0, 300.0], [50.0, 100.0, 150.0]])

    def test_total_ion_count(self) -> None:
        result = normalize_intensities(self.data, method="total_ion_count")
        assert result.shape == self.data.shape
        # Column sums should be approximately equal after TIC normalization
        col_sums = result.sum(axis=0)
        assert np.std(col_sums) / np.mean(col_sums) < 0.01

    def test_median_normalization(self) -> None:
        result = normalize_intensities(self.data, method="median")
        assert result.shape == self.data.shape
        assert np.all(result > 0)

    def test_log2_transform(self) -> None:
        result = normalize_intensities(self.data, method="log2")
        expected = np.log2(self.data + 1)
        np.testing.assert_allclose(result, expected)

    def test_pareto_scaling(self) -> None:
        result = normalize_intensities(self.data, method="pareto")
        assert result.shape == self.data.shape
        # Row means should be ~0 after centering
        assert np.abs(result.mean(axis=1)).max() < 1e-10

    def test_unknown_method_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown normalization"):
            normalize_intensities(self.data, method="invalid")

    def test_zero_column_handling(self) -> None:
        """All-zero column should not cause division by zero."""
        data = np.array([[0.0, 100.0], [0.0, 50.0]])
        result = normalize_intensities(data, method="total_ion_count")
        assert not np.any(np.isnan(result))
        assert not np.any(np.isinf(result))


class TestFoldChange:
    """Tests for fold change computation."""

    def test_equal_groups(self) -> None:
        data = np.array([[100.0, 100.0, 200.0, 200.0]])
        fc = fold_change(data, group_a=[0, 1], group_b=[2, 3])
        assert fc[0] == pytest.approx(-1.0)  # log2(100/200) = -1

    def test_positive_fold_change(self) -> None:
        data = np.array([[200.0, 200.0, 100.0, 100.0]])
        fc = fold_change(data, group_a=[0, 1], group_b=[2, 3])
        assert fc[0] == pytest.approx(1.0)  # log2(200/100) = 1

    def test_zero_handling(self) -> None:
        """Zero means should not produce NaN/Inf."""
        data = np.array([[0.0, 0.0, 100.0, 100.0]])
        fc = fold_change(data, group_a=[0, 1], group_b=[2, 3])
        assert np.isfinite(fc[0])


class TestDifferentialAbundance:
    """Tests for t-test differential abundance."""

    def test_significant_difference(self) -> None:
        """Clearly different groups should have small p-value."""
        rng = np.random.default_rng(42)
        n_met = 5
        data = np.zeros((n_met, 10))
        data[:, :5] = rng.normal(10, 1, (n_met, 5))
        data[:, 5:] = rng.normal(20, 1, (n_met, 5))
        t_stats, p_vals = differential_abundance(data, group_a=[0, 1, 2, 3, 4], group_b=[5, 6, 7, 8, 9])
        assert all(p < 0.05 for p in p_vals)
        assert all(t < 0 for t in t_stats)  # group A < group B

    def test_no_difference(self) -> None:
        """Identical values in both groups should give t-stat near zero."""
        # Use truly identical groups to avoid random variation
        base = np.array([[10.0, 20.0, 30.0]])
        data = np.hstack([base, base])  # 1 metabolite, 6 samples (3+3 identical)
        t_stats, p_vals = differential_abundance(data, group_a=[0, 1, 2], group_b=[3, 4, 5])
        assert abs(t_stats[0]) < 1e-10

    def test_output_shapes(self) -> None:
        data = np.random.default_rng(0).random((10, 8))
        t_stats, p_vals = differential_abundance(data, group_a=[0, 1, 2, 3], group_b=[4, 5, 6, 7])
        assert t_stats.shape == (10,)
        assert p_vals.shape == (10,)


class TestCosineSpectralSimilarity:
    """Tests for spectral similarity scoring."""

    def test_identical_spectra(self) -> None:
        spec = np.array([100.0, 50.0, 25.0, 10.0])
        assert cosine_spectral_similarity(spec, spec) == pytest.approx(1.0)

    def test_orthogonal_spectra(self) -> None:
        a = np.array([1.0, 0.0, 0.0])
        b = np.array([0.0, 1.0, 0.0])
        assert cosine_spectral_similarity(a, b) == pytest.approx(0.0)

    def test_zero_spectrum(self) -> None:
        a = np.array([1.0, 2.0])
        b = np.zeros(2)
        assert cosine_spectral_similarity(a, b) == 0.0

    def test_similarity_range(self) -> None:
        rng = np.random.default_rng(42)
        a = rng.random(100)
        b = rng.random(100)
        sim = cosine_spectral_similarity(a, b)
        assert 0.0 <= sim <= 1.0


class TestIdentifyWithAdducts:
    """Tests for adduct-aware identification."""

    def test_protonated_adduct(self) -> None:
        """[M+H]+ should find neutral mass = m/z - 1.007276."""
        neutral_mass = 180.0634  # glucose
        db = {"glucose": neutral_mass}
        observed_mz = np.array([neutral_mass + COMMON_ADDUCTS["[M+H]+"]])
        results = identify_with_adducts(observed_mz, db, ppm_tolerance=10.0, ion_mode="positive")
        assert len(results) == 1
        matches = results[0]
        assert any(m.matched_name == "glucose" and m.adduct_type == "[M+H]+" for m in matches)

    def test_negative_mode(self) -> None:
        """Negative mode should use negative adducts."""
        db = {"glucose": 180.0634}
        mz = 180.0634 + COMMON_ADDUCTS["[M-H]-"]
        results = identify_with_adducts(np.array([mz]), db, ppm_tolerance=10.0, ion_mode="negative")
        assert len(results) == 1
        if results[0]:
            assert all("-" in m.adduct_type for m in results[0])

    def test_custom_adducts(self) -> None:
        db = {"test_met": 100.0}
        custom = {"[M+2H]2+": 1.00728}
        results = identify_with_adducts(np.array([101.00728]), db, adducts=custom, ppm_tolerance=10.0)
        assert len(results) == 1


class TestMissingValueImputation:
    """Tests for missing value imputation."""

    def test_min_half_imputation(self) -> None:
        data = np.array([[10.0, 0.0, 30.0], [0.0, 20.0, 40.0]])
        result = missing_value_imputation(data, method="min_half")
        assert result[0, 1] == pytest.approx(5.0)  # min(10,30)/2
        assert result[1, 0] == pytest.approx(10.0)  # min(20,40)/2
        assert result[0, 0] == 10.0  # unchanged
        assert result[1, 2] == 40.0  # unchanged

    def test_median_imputation(self) -> None:
        data = np.array([[10.0, 0.0, 30.0]])
        result = missing_value_imputation(data, method="median")
        assert result[0, 1] == pytest.approx(20.0)  # median(10,30)

    def test_knn_imputation(self) -> None:
        """KNN should fill missing values using correlated rows."""
        data = np.array([
            [10.0, 20.0, 0.0, 40.0],
            [10.0, 20.0, 30.0, 40.0],
            [5.0, 10.0, 15.0, 20.0],
        ])
        result = missing_value_imputation(data, method="knn")
        # Row 0 missing value at col 2 should be imputed from neighbors
        assert result[0, 2] > 0  # no longer zero
        assert np.isfinite(result[0, 2])

    def test_unknown_method_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown imputation"):
            missing_value_imputation(np.ones((2, 3)), method="unknown")

    def test_no_missing_values(self) -> None:
        """Data without zeros/NaN should be unchanged."""
        data = np.array([[10.0, 20.0], [30.0, 40.0]])
        result = missing_value_imputation(data, method="min_half")
        np.testing.assert_array_equal(result, data)
