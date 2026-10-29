"""Exact Student's t tail and the SciPy-free Welch test of the metagenomics package.

The Welch test used a normal approximation when SciPy was not installed, which is badly
anticonservative for the group sizes of microbiome studies (a handful of samples). The tail is
now the regularized incomplete beta function. Uses real implementations; SciPy is only an
oracle and each use is guarded with ``pytest.importorskip``.
"""

from __future__ import annotations

import json
import math
import os
import random
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from metainformant.metagenomics.comparative.differential_abundance import _welch_t_test


def _special():
    """The module under test, imported lazily so a missing function fails one test, not the collection."""
    from metainformant.metagenomics import _special as special

    return special


def _package_root() -> str:
    """Directory that holds the ``metainformant`` package this process imported."""
    import metainformant

    return os.path.dirname(os.path.dirname(os.path.abspath(metainformant.__file__)))


def _child_env() -> dict[str, str]:
    """Child environment that imports exactly the package under test (the conftest path is not inherited)."""
    return {**os.environ, "PYTHONPATH": _package_root()}


def _cauchy_two_sided(t: float) -> float:
    """Closed form of the df = 1 tail with full relative precision: (2 / pi) * atan(1 / |t|)."""
    if t == 0:
        return 1.0
    return 2.0 / math.pi * math.atan(1.0 / abs(t))


def _df2_two_sided(t: float) -> float:
    """Closed form of the df = 2 tail, ``1 - |t| / sqrt(2 + t**2)``, rationalised so large |t| stays exact."""
    root = math.sqrt(2.0 + t * t)
    return 2.0 / (root * (root + abs(t)))


class TestStudentTAgainstClosedForms:
    """df = 1 and df = 2 have elementary closed forms, independent of SciPy."""

    @pytest.mark.parametrize("t", [0.0, 1e-6, 0.3, 1.0, 2.5, 10.0, 80.0, 1e3, 1e6, 1e9])
    def test_cauchy_tail(self, t: float) -> None:
        assert _special().student_t_two_sided_p(t, 1.0) == pytest.approx(_cauchy_two_sided(t), rel=1e-12)

    @pytest.mark.parametrize("t", [0.0, 1e-6, 0.3, 1.0, 2.5, 10.0, 80.0, 1e3, 1e6, 1e9])
    def test_df2_tail(self, t: float) -> None:
        assert _special().student_t_two_sided_p(t, 2.0) == pytest.approx(_df2_two_sided(t), rel=1e-12)

    def test_survival_function_is_one_sided_and_symmetric(self) -> None:
        special = _special()
        for df in (1.0, 2.0, 7.5, 40.0):
            for t in (0.0, 0.4, 1.7, 6.0):
                upper = special.student_t_sf(t, df)
                assert upper + special.student_t_sf(-t, df) == pytest.approx(1.0, abs=1e-14)
                assert upper == pytest.approx(0.5 * special.student_t_two_sided_p(t, df), rel=1e-12)

    def test_extreme_tails_do_not_lose_precision(self) -> None:
        """``1 - cdf`` rounds to zero long before the true tail does."""
        special = _special()
        assert special.student_t_two_sided_p(1e6, 5.0) > 0.0
        p = special.student_t_two_sided_p(1e3, 2.0)
        assert p == pytest.approx(_df2_two_sided(1e3), rel=1e-12)
        assert 0.0 < p < 1e-5

    def test_huge_degrees_of_freedom_approach_the_normal_tail(self) -> None:
        special = _special()
        for df in (1e10, 1e14, math.inf):
            for t in (0.5, 2.0, 5.0):
                assert special.student_t_two_sided_p(t, df) == pytest.approx(math.erfc(t / math.sqrt(2.0)), rel=1e-8)

    def test_limits_and_invalid_input(self) -> None:
        special = _special()
        assert special.student_t_two_sided_p(math.inf, 3.0) == 0.0
        assert special.student_t_sf(math.inf, 3.0) == 0.0
        assert special.student_t_sf(-math.inf, 3.0) == 1.0
        assert math.isnan(special.student_t_two_sided_p(math.nan, 3.0))
        for df in (0.0, -1.0, math.nan):
            with pytest.raises(ValueError, match="df"):
                special.student_t_two_sided_p(1.0, df)


class TestStudentTAgainstScipy:
    """The tail matches SciPy wherever SciPy itself is accurate."""

    DFS = [0.5, 1.0, 2.0, 3.7, 10.0, 30.5, 120.0, 5000.0, 1e6]
    TS = [0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 15.0, 30.0]

    def test_survival_function_and_two_sided_p(self) -> None:
        stats = pytest.importorskip("scipy.stats")
        special = _special()
        for df in self.DFS:
            for t in self.TS:
                for signed in (t, -t):
                    assert special.student_t_sf(signed, df) == pytest.approx(float(stats.t.sf(signed, df)), rel=1e-9)
                assert special.student_t_two_sided_p(t, df) == pytest.approx(float(2.0 * stats.t.sf(t, df)), rel=1e-9)

    def test_incomplete_beta_matches_scipy(self) -> None:
        scipy_special = pytest.importorskip("scipy.special")
        special = _special()
        rng = random.Random(12)
        for _ in range(300):
            a, b, x = 10 ** rng.uniform(-1.5, 3.0), 10 ** rng.uniform(-1.5, 3.0), rng.random()
            lower = float(scipy_special.betainc(a, b, x))
            upper = float(scipy_special.betainc(b, a, 1.0 - x))
            assert special.regularized_beta_p(a, b, x) == pytest.approx(lower, rel=1e-8, abs=1e-250)
            assert special.regularized_beta_q(a, b, x) == pytest.approx(upper, rel=1e-8, abs=1e-250)

    def test_incomplete_beta_identities(self) -> None:
        special = _special()
        for a, b, x in [(0.5, 0.5, 0.3), (2.0, 5.0, 0.9), (30.0, 1.5, 0.999), (1.0, 1.0, 0.25)]:
            total = special.regularized_beta_p(a, b, x) + special.regularized_beta_q(a, b, x)
            assert total == pytest.approx(1.0, abs=1e-14)
        assert special.regularized_beta_p(1.0, 1.0, 0.25) == pytest.approx(0.25, rel=1e-13)  # uniform distribution
        assert special.regularized_beta_p(2.0, 3.0, 0.0) == 0.0
        assert special.regularized_beta_p(2.0, 3.0, 1.0) == 1.0

    def test_incomplete_beta_rejects_bad_arguments(self) -> None:
        special = _special()
        with pytest.raises(ValueError):
            special.regularized_beta_p(0.0, 1.0, 0.5)
        with pytest.raises(ValueError):
            special.regularized_beta_q(1.0, 1.0, 1.5)


class TestWelchTTest:
    """``_welch_t_test`` is exact, matches SciPy, and reports NaN when undefined."""

    def test_matches_scipy_for_random_samples(self) -> None:
        stats = pytest.importorskip("scipy.stats")
        rng = random.Random(0)
        for _ in range(300):
            a = [rng.gauss(rng.uniform(-2, 2), rng.uniform(0.05, 4)) for _ in range(rng.randint(2, 30))]
            b = [rng.gauss(rng.uniform(-2, 2), rng.uniform(0.05, 4)) for _ in range(rng.randint(2, 30))]
            expected = float(stats.ttest_ind(a, b, equal_var=False).pvalue)
            assert _welch_t_test(a, b) == pytest.approx(expected, rel=1e-9, abs=1e-300)

    def test_hand_computed_small_samples(self) -> None:
        # means 2 and 7/3, variances 1 and 7/3: t = -1/sqrt(10) = -0.31623, Welch df = 100/29 = 3.44828,
        # and P(|T_3.44828| >= 0.31623) = 0.77003 (the normal tail would say 0.75183)
        a, b = [1.0, 2.0, 3.0], [1.0, 2.0, 4.0]
        assert _welch_t_test(a, b) == pytest.approx(0.7700256383911763, rel=1e-9)
        assert _special().student_t_two_sided_p(-1.0 / math.sqrt(10.0), 100.0 / 29.0) == pytest.approx(
            0.7700256383911763, rel=1e-9
        )

    def test_undefined_tests_are_nan_not_one(self) -> None:
        assert math.isnan(_welch_t_test([1.0], [1.0, 2.0, 3.0]))  # one sample has no variance
        assert math.isnan(_welch_t_test([], [1.0, 2.0]))
        assert math.isnan(_welch_t_test([2.0, 2.0, 2.0], [2.0, 2.0]))  # no variance anywhere: 0 / 0
        assert math.isnan(_welch_t_test([1.0, 1.0], [5.0, 5.0]))  # no variance anywhere: x / 0

    def test_one_constant_group_is_still_testable(self) -> None:
        """Without variance in one group Welch's df collapses to ``n_b - 1`` and t to a one-sample t."""
        a, b = [3.0, 3.0, 3.0, 3.0], [1.0, 2.0, 2.5, 4.0, 6.5]
        mean_b = sum(b) / len(b)
        var_b = sum((x - mean_b) ** 2 for x in b) / (len(b) - 1)
        t_stat = (3.0 - mean_b) / math.sqrt(var_b / len(b))
        assert _welch_t_test(a, b) == pytest.approx(_special().student_t_two_sided_p(t_stat, len(b) - 1), rel=1e-12)

    def test_result_does_not_depend_on_the_unit_of_the_data(self) -> None:
        """Squared variances overflow or underflow for data of magnitude 1e200 or 1e-200 unless it is rescaled."""
        a, b = [1.0, 2.0, 3.5, 2.2], [3.0, 4.5, 2.0, 5.1, 3.9]
        expected = _welch_t_test(a, b)
        assert 0.0 < expected < 1.0
        for unit in (1e-200, 1e-12, 1e12, 1e200):
            assert _welch_t_test([x * unit for x in a], [x * unit for x in b]) == pytest.approx(expected, rel=1e-9)

    def test_non_finite_values_are_undefined(self) -> None:
        assert math.isnan(_welch_t_test([1.0, math.nan, 3.0], [1.0, 2.0, 4.0]))
        assert math.isnan(_welch_t_test([1.0, 2.0, 3.0], [1.0, math.inf, 4.0]))

    def test_small_groups_are_not_judged_by_a_normal_tail(self) -> None:
        """3 vs 3 samples with |t| = 3.41 and Welch df = 3.81: exact p = 0.0291, normal tail 0.00064."""
        a, b = [10.0, 11.0, 12.0], [13.0, 14.0, 15.5]
        assert _welch_t_test(a, b) == pytest.approx(0.029145305862049366, rel=1e-9)


class TestWelchWithoutScipy:
    """Run the real SciPy-free path in a child interpreter in which SciPy cannot be imported."""

    def test_child_p_values_match_scipy(self, tmp_path: Path) -> None:
        stats = pytest.importorskip("scipy.stats")
        rng = random.Random(5)
        cases = []
        for _ in range(25):
            n_a, n_b = rng.randint(2, 6), rng.randint(2, 6)  # microbiome-sized groups
            cases.append(
                (
                    [rng.gauss(0.0, 1.0) for _ in range(n_a)],
                    [rng.gauss(rng.uniform(-1.5, 1.5), rng.uniform(0.5, 2.0)) for _ in range(n_b)],
                )
            )
        code = textwrap.dedent(f"""
            import json, os, sys
            sys.modules["scipy"] = None  # makes `import scipy` raise ImportError
            import metainformant
            assert os.path.dirname(os.path.dirname(os.path.abspath(metainformant.__file__))) == {_package_root()!r}
            from metainformant.metagenomics.comparative import differential_abundance as da
            assert da.HAS_SCIPY is False
            print(json.dumps([da._welch_t_test(a, b) for a, b in {cases!r}]))
            """)
        run = subprocess.run(
            [sys.executable, "-c", code],
            env=_child_env(),
            cwd=tmp_path,
            capture_output=True,
            text=True,
            check=True,
            timeout=180,
        )
        child = json.loads(run.stdout.strip().splitlines()[-1])
        expected = [float(stats.ttest_ind(a, b, equal_var=False).pvalue) for a, b in cases]
        assert child == pytest.approx(expected, rel=1e-9, abs=1e-300)

    def test_differential_abundance_p_values_do_not_depend_on_scipy(self, tmp_path: Path) -> None:
        rng = random.Random(21)
        counts = [[rng.randint(5, 60) + (30 if g and j < 2 else 0) for j in range(7)] for g in (0, 0, 0, 1, 1, 1)]
        groups = [0, 0, 0, 1, 1, 1]
        taxa = [f"t{j}" for j in range(7)]
        code = textwrap.dedent(f"""
            import json, os, sys
            sys.modules["scipy"] = None
            import metainformant
            assert os.path.dirname(os.path.dirname(os.path.abspath(metainformant.__file__))) == {_package_root()!r}
            from metainformant.metagenomics.comparative.differential_abundance import differential_abundance
            rows = differential_abundance({counts!r}, {groups!r}, {taxa!r}, method="simple_deseq")
            print(json.dumps({{r["taxon"]: r["p_value"] for r in rows}}))
            """)
        run = subprocess.run(
            [sys.executable, "-c", code],
            env=_child_env(),
            cwd=tmp_path,
            capture_output=True,
            text=True,
            check=True,
            timeout=180,
        )
        without_scipy = json.loads(run.stdout.strip().splitlines()[-1])

        from metainformant.metagenomics.comparative.differential_abundance import differential_abundance

        with_scipy = {r["taxon"]: r["p_value"] for r in differential_abundance(counts, groups, taxa, "simple_deseq")}
        assert without_scipy == pytest.approx(with_scipy, rel=1e-10, abs=1e-300)
