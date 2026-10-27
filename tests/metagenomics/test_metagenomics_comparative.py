"""Tests for metagenomics comparative submodule.

Tests differential abundance, CLR transform, indicator species, biomarkers.
Uses real implementations -- real-implementation policy.
"""

from __future__ import annotations

import json
import math
import os
import random
import subprocess
import sys
import textwrap

import numpy as np
import pytest
from scipy import stats as scipy_stats

from metainformant.metagenomics.comparative.differential_abundance import (
    biomarker_discovery,
    clr_transform,
    differential_abundance,
    effect_size_analysis,
    indicator_species,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def count_matrix_two_groups() -> tuple[list[list[int]], list[int], list[str]]:
    """Count matrix with clear differential abundance between two groups.

    Group 0 (samples 0-3): taxon_A dominant.
    Group 1 (samples 4-7): taxon_B dominant.
    """
    counts = [
        # Group 0: high taxon_A, low taxon_B
        [100, 5, 20, 10],
        [90, 8, 18, 12],
        [110, 3, 22, 9],
        [95, 6, 19, 11],
        # Group 1: low taxon_A, high taxon_B
        [5, 100, 15, 12],
        [8, 95, 17, 10],
        [3, 110, 14, 13],
        [6, 105, 16, 11],
    ]
    groups = [0, 0, 0, 0, 1, 1, 1, 1]
    taxa_names = ["taxon_A", "taxon_B", "taxon_C", "taxon_D"]
    return counts, groups, taxa_names


@pytest.fixture
def simple_count_matrix() -> tuple[list[list[int]], list[int], list[str]]:
    """Minimal count matrix for basic function checks."""
    counts = [
        [10, 20, 30],
        [15, 25, 35],
        [5, 50, 10],
        [8, 45, 12],
    ]
    groups = [0, 0, 1, 1]
    taxa_names = ["sp_X", "sp_Y", "sp_Z"]
    return counts, groups, taxa_names


# ---------------------------------------------------------------------------
# Tests: differential_abundance
# ---------------------------------------------------------------------------


class TestDifferentialAbundance:
    """Tests for differential abundance analysis."""

    def test_aldex2_like_basic(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        results = differential_abundance(counts, groups, taxa, method="aldex2_like")
        assert isinstance(results, list)
        assert len(results) == len(taxa)
        for r in results:
            assert "taxon" in r
            assert "log2fc" in r
            assert "p_value" in r
            assert "adjusted_p" in r
            assert "effect_size" in r

    def test_ancom_like_basic(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        results = differential_abundance(counts, groups, taxa, method="ancom_like")
        assert len(results) == len(taxa)

    def test_simple_deseq_basic(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        results = differential_abundance(counts, groups, taxa, method="simple_deseq")
        assert len(results) == len(taxa)

    def test_sorted_by_adjusted_p(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        results = differential_abundance(counts, groups, taxa)
        adjusted_ps = [r["adjusted_p"] for r in results]
        assert adjusted_ps == sorted(adjusted_ps)

    def test_invalid_method_raises(
        self,
        simple_count_matrix: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = simple_count_matrix
        with pytest.raises(ValueError, match="Invalid method"):
            differential_abundance(counts, groups, taxa, method="bogus")

    def test_mismatched_groups_raises(self) -> None:
        with pytest.raises(ValueError, match="groups length"):
            differential_abundance([[1, 2]], [0, 0], ["a", "b"])

    def test_three_groups_raises(self) -> None:
        with pytest.raises(ValueError, match="exactly 2 unique"):
            differential_abundance(
                [[1, 2], [3, 4], [5, 6]],
                [0, 1, 2],
                ["a", "b"],
            )


# ---------------------------------------------------------------------------
# Tests: clr_transform
# ---------------------------------------------------------------------------


class TestClrTransform:
    """Tests for centered log-ratio transformation."""

    def test_rows_are_centered(self) -> None:
        counts = [[10, 20, 30], [5, 50, 10]]
        result = clr_transform(counts)
        assert len(result) == 2
        for row in result:
            assert len(row) == 3
            # CLR rows should sum to approximately 0
            assert abs(sum(row)) < 1e-6

    def test_pseudocount_avoids_log_zero(self) -> None:
        counts = [[0, 0, 100]]
        result = clr_transform(counts, pseudocount=0.5)
        assert len(result) == 1
        # No NaN or Inf values
        for v in result[0]:
            assert math.isfinite(v)

    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError, match="not be empty"):
            clr_transform([])


# ---------------------------------------------------------------------------
# Tests: indicator_species
# ---------------------------------------------------------------------------


class TestIndicatorSpecies:
    """Tests for indicator species analysis."""

    def test_returns_species_group_associations(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        results = indicator_species(counts, groups, taxa, n_permutations=99, seed=42)
        assert isinstance(results, list)
        assert len(results) == len(taxa)
        for r in results:
            assert "taxon" in r
            assert "indicator_value" in r
            assert "p_value" in r
            assert "associated_group" in r
            assert 0.0 <= r["indicator_value"] <= 1.0
            assert 0.0 <= r["p_value"] <= 1.0

    def test_sorted_by_indicator_value(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        results = indicator_species(counts, groups, taxa, n_permutations=49, seed=42)
        indvals = [r["indicator_value"] for r in results]
        assert indvals == sorted(indvals, reverse=True)


# ---------------------------------------------------------------------------
# Tests: effect_size_analysis
# ---------------------------------------------------------------------------


class TestEffectSizeAnalysis:
    """Tests for LEfSe-style effect size ranking."""

    def test_returns_lda_scores(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        results = effect_size_analysis(counts, groups, taxa)
        assert isinstance(results, list)
        assert len(results) == len(taxa)
        for r in results:
            assert "taxon" in r
            assert "lda_score" in r
            assert "p_value" in r
            assert "direction" in r

    def test_sorted_by_abs_lda(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        results = effect_size_analysis(counts, groups, taxa)
        lda_scores = [abs(r["lda_score"]) for r in results]
        assert lda_scores == sorted(lda_scores, reverse=True)


# ---------------------------------------------------------------------------
# Tests: biomarker_discovery
# ---------------------------------------------------------------------------


class TestBiomarkerDiscovery:
    """Tests for ML-based biomarker discovery."""

    def test_returns_expected_keys(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        result = biomarker_discovery(counts, groups, taxa, seed=42)
        assert isinstance(result, dict)
        assert "selected_taxa" in result
        assert "importances" in result
        assert "cv_accuracy" in result
        assert "model_summary" in result
        assert isinstance(result["selected_taxa"], list)
        assert isinstance(result["importances"], dict)
        assert all(t in taxa for t in result["importances"])

    def test_importances_are_nonneg(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        result = biomarker_discovery(counts, groups, taxa, seed=42)
        for imp in result["importances"].values():
            assert imp >= 0.0

    def test_invalid_method_raises(
        self,
        simple_count_matrix: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = simple_count_matrix
        with pytest.raises(ValueError, match="Invalid method"):
            biomarker_discovery(counts, groups, taxa, method="xgboost")


class TestBiomarkerDiscoveryValidation:
    """Input validation for biomarker discovery."""

    def test_mismatched_groups_length_raises(
        self,
        simple_count_matrix: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, _, taxa = simple_count_matrix
        with pytest.raises(ValueError, match="groups length"):
            biomarker_discovery(counts, [0, 0], taxa)


class TestDifferentialAbundanceDirection:
    """Fold-change direction must match group dominance."""

    def test_log2fc_direction(
        self,
        count_matrix_two_groups: tuple[list[list[int]], list[int], list[str]],
    ) -> None:
        counts, groups, taxa = count_matrix_two_groups
        results = differential_abundance(counts, groups, taxa, method="aldex2_like")
        by_taxon = {r["taxon"]: r for r in results}
        # log2fc = log2((mean_group0 + 1) / (mean_group1 + 1)):
        # taxon_A dominates group 0 -> positive; taxon_B dominates group 1 -> negative
        assert by_taxon["taxon_A"]["log2fc"] > 0
        assert by_taxon["taxon_B"]["log2fc"] < 0


# ---------------------------------------------------------------------------
# Helpers for the regression tests below
# ---------------------------------------------------------------------------


def _subprocess_env() -> dict[str, str]:
    """Environment for a child interpreter that imports the same ``metainformant`` as this process.

    ``tests/conftest.py`` prepends the checkout's ``src`` to ``PYTHONPATH``, so inheriting the
    variable could silently select a different copy of the package than the one under test.
    """
    import metainformant

    package_root = os.path.dirname(os.path.dirname(os.path.abspath(metainformant.__file__)))
    return {**os.environ, "PYTHONPATH": package_root}


def _two_group_counts(n_per_group: int, n_taxa: int = 6, seed: int = 0) -> tuple[list[list[int]], list[int], list[str]]:
    """Seeded count matrix with a shift in the first two taxa between the groups."""
    rng = random.Random(seed)
    counts = []
    groups = []
    for group in (0, 1):
        for _ in range(n_per_group):
            row = [rng.randint(20, 60) for _ in range(n_taxa)]
            if group == 1:
                row[0] += 80
                row[1] = max(1, row[1] - 15)
            counts.append(row)
            groups.append(group)
    return counts, groups, [f"taxon_{i}" for i in range(n_taxa)]


# ---------------------------------------------------------------------------
# Regression tests: cross-validation folds are bounded by class size, not group count
# ---------------------------------------------------------------------------


class TestBiomarkerCrossValidationFolds:
    """cv_folds used to be capped at the number of groups (2 for a two-group comparison)."""

    @pytest.fixture(autouse=True)
    def _require_sklearn(self) -> None:
        pytest.importorskip("sklearn")

    def test_requested_folds_are_used(self) -> None:
        counts, groups, taxa = _two_group_counts(10)
        result = biomarker_discovery(counts, groups, taxa, cv_folds=5, seed=1)
        assert result["model_summary"]["cv_folds"] == 5
        assert len(result["model_summary"]["cv_scores"]) == 5

    def test_folds_are_bounded_by_the_smallest_class(self) -> None:
        rows_a, _, taxa = _two_group_counts(6, seed=1)
        rows_b, _, _ = _two_group_counts(4, seed=2)
        counts = rows_a[:6] + rows_b[4:]  # 6 samples of group 0, 4 of group 1
        groups = [0] * 6 + [1] * 4
        result = biomarker_discovery(counts, groups, taxa, cv_folds=10, seed=1)
        assert result["model_summary"]["cv_folds"] == 4
        assert len(result["model_summary"]["cv_scores"]) == 4

    def test_default_five_folds(self) -> None:
        counts, groups, taxa = _two_group_counts(8)
        assert biomarker_discovery(counts, groups, taxa, seed=1)["model_summary"]["cv_folds"] == 5


# ---------------------------------------------------------------------------
# Regression tests: Kruskal-Wallis without SciPy covers every group
# ---------------------------------------------------------------------------


class TestKruskalWallisFallback:
    """The SciPy-free fallback compared only the first two groups, with the wrong tail."""

    GROUPS = [[1.0, 2.0, 3.0, 4.0], [3.5, 4.5, 5.5, 6.5], [10.0, 11.0, 12.0, 13.0]]

    def test_three_groups_use_all_groups(self) -> None:
        from metainformant.metagenomics.comparative.differential_abundance import _kruskal_wallis_pure

        expected = scipy_stats.kruskal(*self.GROUPS).pvalue
        assert expected == pytest.approx(0.0097098, rel=1e-4)
        assert _kruskal_wallis_pure(self.GROUPS) == pytest.approx(expected, rel=1e-9)

    def test_two_groups_match_scipy(self) -> None:
        from metainformant.metagenomics.comparative.differential_abundance import _kruskal_wallis_pure

        a, b = self.GROUPS[0], self.GROUPS[1]
        assert _kruskal_wallis_pure([a, b]) == pytest.approx(scipy_stats.kruskal(a, b).pvalue, rel=1e-9)

    def test_random_groups_with_ties_match_scipy(self) -> None:
        from metainformant.metagenomics.comparative.differential_abundance import _kruskal_wallis_pure

        rng = random.Random(4)
        for _ in range(60):
            n_groups = rng.randint(2, 5)
            groups = [[float(rng.randint(0, 12)) for _ in range(rng.randint(3, 9))] for _ in range(n_groups)]
            if len({v for g in groups for v in g}) < 2:
                continue
            expected = scipy_stats.kruskal(*groups).pvalue
            assert _kruskal_wallis_pure(groups) == pytest.approx(expected, rel=1e-8, abs=1e-14)

    def test_undefined_cases_return_one(self) -> None:
        from metainformant.metagenomics.comparative.differential_abundance import _kruskal_wallis_pure

        assert _kruskal_wallis_pure([[1.0, 1.0], [1.0, 1.0]]) == 1.0  # all values identical
        assert _kruskal_wallis_pure([[1.0, 2.0]]) == 1.0  # one group
        assert _kruskal_wallis_pure([[1.0, 2.0], []]) == 1.0  # empty group

    def test_effect_size_analysis_p_values_without_scipy(self) -> None:
        """Run the real fallback path in a child interpreter in which SciPy cannot be imported."""
        rng = random.Random(8)
        counts = []
        groups = []
        for group, shift in enumerate((0, 25, 90)):
            for _ in range(5):
                counts.append([rng.randint(5, 40) + shift, rng.randint(5, 40), rng.randint(5, 40) + shift // 2])
                groups.append(group)
        taxa = ["t0", "t1", "t2"]

        code = textwrap.dedent(f"""
            import json, sys
            sys.modules["scipy"] = None  # makes `import scipy` raise ImportError
            from metainformant.metagenomics.comparative import differential_abundance as da
            assert da.HAS_SCIPY is False
            rows = da.effect_size_analysis({counts!r}, {groups!r}, {taxa!r})
            print(json.dumps({{r["taxon"]: r["p_value"] for r in rows}}))
            """)
        run = subprocess.run(
            [sys.executable, "-c", code],
            env=_subprocess_env(),
            capture_output=True,
            text=True,
            check=True,
            timeout=180,
        )
        child = json.loads(run.stdout.strip().splitlines()[-1])

        clr = clr_transform(counts)
        for j, taxon in enumerate(taxa):
            samples = [[clr[i][j] for i in range(len(groups)) if groups[i] == g] for g in (0, 1, 2)]
            assert child[taxon] == pytest.approx(scipy_stats.kruskal(*samples).pvalue, rel=1e-8)


class TestChiSquareSurvival:
    """The pure-Python chi-square tail matches SciPy, including tiny tails."""

    @pytest.mark.parametrize("df", [1, 2, 3, 4, 7, 12, 40])
    def test_matches_scipy(self, df: int) -> None:
        from metainformant.metagenomics._special import chi2_sf

        for x in [0.0, 0.01, 0.5, 1.0, 3.84, 9.0, 25.0, 80.0, 200.0]:
            expected = float(scipy_stats.chi2.sf(x, df))
            assert chi2_sf(x, df) == pytest.approx(expected, rel=1e-9, abs=1e-300)

    def test_gamma_tails_match_scipy(self) -> None:
        from scipy import special

        from metainformant.metagenomics._special import regularized_gamma_p, regularized_gamma_q

        for a in (0.5, 1.0, 2.5, 10.0, 60.0, 400.0):
            for x in (0.01, 0.5, 3.0, 20.0, 80.0, 500.0):
                assert regularized_gamma_p(a, x) == pytest.approx(float(special.gammainc(a, x)), rel=1e-9, abs=1e-300)
                assert regularized_gamma_q(a, x) == pytest.approx(float(special.gammaincc(a, x)), rel=1e-9, abs=1e-300)

    def test_invalid_shape_raises(self) -> None:
        from metainformant.metagenomics._special import regularized_gamma_p, regularized_gamma_q

        with pytest.raises(ValueError, match="positive"):
            regularized_gamma_p(0.0, 1.0)
        with pytest.raises(ValueError, match="positive"):
            regularized_gamma_q(-1.0, 1.0)


# ---------------------------------------------------------------------------
# Regression tests: numpy arrays are accepted (no truthiness tests on inputs)
# ---------------------------------------------------------------------------


class TestNumpyArrayInputs:
    """`if not counts` raised "truth value of an array is ambiguous" for numpy input."""

    @pytest.fixture
    def data(self) -> tuple[list[list[int]], list[int], list[str]]:
        return _two_group_counts(6)

    def test_clr_transform(self, data: tuple[list[list[int]], list[int], list[str]]) -> None:
        counts, _, _ = data
        np.testing.assert_allclose(clr_transform(np.array(counts)), clr_transform(counts))

    @pytest.mark.parametrize("method", ["aldex2_like", "ancom_like", "simple_deseq"])
    def test_differential_abundance(self, data: tuple[list[list[int]], list[int], list[str]], method: str) -> None:
        counts, groups, taxa = data
        from_lists = differential_abundance(counts, groups, taxa, method=method, n_monte_carlo=8, seed=3)
        from_arrays = differential_abundance(
            np.array(counts), np.array(groups), np.array(taxa), method=method, n_monte_carlo=8, seed=3
        )
        assert from_arrays == from_lists

    def test_indicator_species(self, data: tuple[list[list[int]], list[int], list[str]]) -> None:
        counts, groups, taxa = data
        from_lists = indicator_species(counts, groups, taxa, n_permutations=49, seed=5)
        from_arrays = indicator_species(np.array(counts), np.array(groups), np.array(taxa), n_permutations=49, seed=5)
        assert from_arrays == from_lists

    def test_effect_size_analysis(self, data: tuple[list[list[int]], list[int], list[str]]) -> None:
        counts, groups, taxa = data
        assert effect_size_analysis(np.array(counts), np.array(groups), np.array(taxa)) == effect_size_analysis(
            counts, groups, taxa
        )

    def test_biomarker_discovery(self, data: tuple[list[list[int]], list[int], list[str]]) -> None:
        counts, groups, taxa = data
        from_lists = biomarker_discovery(counts, groups, taxa, seed=2)
        from_arrays = biomarker_discovery(np.array(counts), np.array(groups), np.array(taxa), seed=2)
        assert from_arrays["importances"] == from_lists["importances"]
        assert from_arrays["selected_taxa"] == from_lists["selected_taxa"]

    def test_empty_array_is_reported_as_empty(self) -> None:
        with pytest.raises(ValueError, match="not be empty"):
            clr_transform(np.zeros((0, 3)))
        with pytest.raises(ValueError, match="not be empty"):
            differential_abundance(np.zeros((0, 3)), [], ["a", "b", "c"])

    def test_one_dimensional_counts_are_rejected_clearly(self) -> None:
        with pytest.raises(ValueError, match="2-D"):
            clr_transform(np.array([1, 2, 3]))
