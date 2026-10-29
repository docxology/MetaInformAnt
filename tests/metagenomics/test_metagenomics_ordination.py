"""PCoA with honest negative eigenvalues, Kruskal's non-metric MDS, and the ordination plot.

NMDS used to be metric MDS: gradient descent on raw stress between the inter-point distances and
the dissimilarities themselves (no isotonic step), one random start, with ``1 - stress`` reported as
"variance explained". PCoA in the plot used power iteration, which returns the eigenvalue of
largest magnitude, so a large negative eigenvalue of a non-Euclidean distance matrix hijacked an
axis. Uses real implementations; NumPy, SciPy and scikit-learn are oracles, each guarded with
``pytest.importorskip``.
"""

from __future__ import annotations

import inspect
import json
import math
import os
import random
import subprocess
import sys
import textwrap
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import pytest
from matplotlib.text import Annotation

from metainformant.metagenomics.diversity.metrics import ordination, permanova
from metainformant.metagenomics.visualization.plots import plot_ordination

np = pytest.importorskip("numpy")


# ---------------------------------------------------------------------------
# Helpers and oracles (independent of the code under test)
# ---------------------------------------------------------------------------


def _package_root() -> str:
    import metainformant

    return os.path.dirname(os.path.dirname(os.path.abspath(metainformant.__file__)))


def _run_child(code: str, cwd: Path) -> str:
    guarded = (
        "import os, metainformant\n"
        f"assert os.path.dirname(os.path.dirname(os.path.abspath(metainformant.__file__))) == {_package_root()!r}\n"
    ) + code
    run = subprocess.run(
        [sys.executable, "-c", guarded],
        env={**os.environ, "PYTHONPATH": _package_root()},
        cwd=cwd,
        capture_output=True,
        text=True,
        check=True,
        timeout=300,
    )
    return run.stdout.strip().splitlines()[-1]


def _pairwise(points: object) -> object:
    arr = np.asarray(points, dtype=float)
    return np.linalg.norm(arr[:, None, :] - arr[None, :, :], axis=2)


def _double_centred(distances: object) -> object:
    d = np.asarray(distances, dtype=float)
    n = len(d)
    j = np.eye(n) - np.ones((n, n)) / n
    return -0.5 * j @ (d**2) @ j


def _procrustes_error(planted: object, found: object) -> float:
    """Residual after centring, unit-scaling and the best rotation/reflection (0 is a perfect match)."""
    x = np.asarray(planted, dtype=float)
    y = np.asarray(found, dtype=float)
    x = x - x.mean(axis=0)
    y = y - y.mean(axis=0)
    x = x / np.linalg.norm(x)
    y = y / np.linalg.norm(y)
    u, _, vt = np.linalg.svd(y.T @ x)
    return float(np.linalg.norm(x - y @ (u @ vt)))


def _kruskal_stress1(dissimilarities: object, config: object) -> float:
    """Kruskal's stress-1 of a configuration, with an independent isotonic fit (ties by distance)."""
    isotonic = pytest.importorskip("sklearn.isotonic")
    delta = np.asarray(dissimilarities, dtype=float)
    iu = np.triu_indices(len(delta), 1)
    d = _pairwise(config)[iu]
    order = np.lexsort((d, delta[iu]))
    fitted = np.empty_like(d)
    fitted[order] = isotonic.IsotonicRegression(increasing=True).fit_transform(np.arange(len(d)), d[order])
    return float(np.sqrt(((d - fitted) ** 2).sum() / (d**2).sum()))


def _complete_bipartite(m: int) -> list[list[float]]:
    """Graph metric of K(m,m): distance 1 across the two sides and 2 within a side (not Euclidean)."""
    n = 2 * m
    d = [[0.0 if i == j else (1.0 if (i < m) != (j < m) else 2.0) for j in range(n)] for i in range(n)]
    return d


def _cycle_four() -> list[list[float]]:
    return [[0, 1, 2, 1], [1, 0, 1, 2], [2, 1, 0, 1], [1, 2, 1, 0]]


def _gradient_bray_curtis(seed: int, n: int = 22, taxa: int = 16) -> object:
    """Bray-Curtis dissimilarities of samples along an environmental gradient (continuous, no ties)."""
    rng = np.random.default_rng(seed)
    x = np.sort(rng.uniform(0, 10, n))
    mu = np.linspace(-1, 11, taxa)
    mean = 60 * np.exp(-((x[:, None] - mu[None, :]) ** 2) / (2 * 1.6**2))
    counts = rng.poisson(mean + 0.5).astype(float)
    numerator = np.abs(counts[:, None, :] - counts[None, :, :]).sum(axis=2)
    denominator = (counts[:, None, :] + counts[None, :, :]).sum(axis=2)
    return numerator / denominator


def _planted(n: int = 14, seed: int = 3) -> object:
    rng = np.random.default_rng(seed)
    return rng.normal(size=(n, 2)) * np.array([3.0, 1.5])


# ---------------------------------------------------------------------------
# PCoA
# ---------------------------------------------------------------------------


class TestPcoaEigenvalues:
    """PCoA keeps the whole spectrum, so negative eigenvalues are reported, not hidden."""

    def test_non_euclidean_matrix_reports_its_negative_eigenvalue(self) -> None:
        # K(4,4): B has eigenvalues 2 (six times), 0 and -4 (hand-derived from the circulant block structure)
        result = ordination(_complete_bipartite(4), method="pcoa", n_components=2)
        assert result["all_eigenvalues"] == pytest.approx([2.0] * 6 + [0.0, -4.0], abs=1e-9)
        assert result["n_negative_eigenvalues"] == 1
        assert result["negative_eigenvalue_fraction"] == pytest.approx(4.0 / 16.0)
        assert result["eigenvalues"] == pytest.approx([2.0, 2.0], abs=1e-9)

    def test_variance_explained_is_relative_to_all_positive_eigenvalues(self) -> None:
        """It was divided by the sum of the *retained* eigenvalues, so two axes always summed to 1."""
        result = ordination(_complete_bipartite(4), method="pcoa", n_components=2)
        assert result["variance_explained"] == pytest.approx([2.0 / 12.0, 2.0 / 12.0], abs=1e-9)
        one_axis = ordination(_cycle_four(), method="pcoa", n_components=1)
        assert one_axis["all_eigenvalues"] == pytest.approx([2.0, 2.0, 0.0, -1.0], abs=1e-9)
        assert one_axis["variance_explained"] == pytest.approx([0.5], abs=1e-9)  # 2 of 2 + 2 positive

    def test_euclidean_matrix_has_no_negative_part(self) -> None:
        points = _planted(10, seed=1)
        result = ordination(_pairwise(points).tolist(), method="pcoa", n_components=2)
        assert result["n_negative_eigenvalues"] == 0
        assert result["negative_eigenvalue_fraction"] == 0.0
        assert sum(result["variance_explained"]) == pytest.approx(1.0, abs=1e-9)  # exactly two positive axes
        assert _procrustes_error(points, result["coordinates"]) < 1e-6

    def test_axes_have_the_eigenvalue_as_their_sum_of_squares(self) -> None:
        result = ordination(_complete_bipartite(3), method="pcoa", n_components=2)
        coords = np.asarray(result["coordinates"])
        assert (coords**2).sum(axis=0) == pytest.approx(result["eigenvalues"], abs=1e-9)
        assert float(coords[:, 0] @ coords[:, 1]) == pytest.approx(0.0, abs=1e-9)

    def test_eigenvalues_match_an_independent_decomposition(self) -> None:
        rng = np.random.default_rng(8)
        a = rng.uniform(1, 2, size=(9, 9))
        d = np.triu(a, 1)
        d = d + d.T
        expected = np.sort(np.linalg.eigvalsh(_double_centred(d)))[::-1]
        result = ordination(d.tolist(), method="pcoa", n_components=3)
        assert result["all_eigenvalues"] == pytest.approx(expected.tolist(), abs=1e-9)

    def test_pure_python_eigensolver_matches_lapack(self) -> None:
        from metainformant.metagenomics.diversity.metrics import _jacobi_eigh

        rng = np.random.default_rng(4)
        for n in (2, 5, 9):
            a = rng.normal(size=(n, n))
            a = a + a.T
            values, vectors = _jacobi_eigh(a.tolist())
            assert values == pytest.approx(np.linalg.eigvalsh(a).tolist(), abs=1e-10)
            v = np.asarray(vectors)
            assert v @ np.diag(values) @ v.T == pytest.approx(a, abs=1e-9)
            assert v.T @ v == pytest.approx(np.eye(n), abs=1e-10)

    def test_without_numpy_the_spectrum_is_still_exact(self, tmp_path: Path) -> None:
        """The fallback used power iteration, which picks the eigenvalue of largest magnitude: here -4."""
        rng = np.random.default_rng(6)
        a = rng.uniform(1, 2, size=(7, 7))
        generic = np.triu(a, 1)
        generic = (generic + generic.T).tolist()  # distinct eigenvalues, so the axes are determined
        code = textwrap.dedent(f"""
            import json, sys
            sys.modules["numpy"] = None  # makes `import numpy` raise ImportError
            from metainformant.metagenomics.diversity import metrics
            assert metrics.HAS_NUMPY is False
            print(json.dumps({{
                "bipartite": metrics.ordination({_complete_bipartite(4)!r}, method="pcoa", n_components=2),
                "generic": metrics.ordination({generic!r}, method="pcoa", n_components=2),
            }}))
            """)
        child = json.loads(_run_child(code, tmp_path))

        # K(4,4): the plane inside the six-fold eigenvalue 2 is arbitrary, so compare what is determined
        bipartite = child["bipartite"]
        assert bipartite["all_eigenvalues"] == pytest.approx([2.0] * 6 + [0.0, -4.0], abs=1e-9)
        assert bipartite["variance_explained"] == pytest.approx([2.0 / 12.0] * 2, abs=1e-9)
        axes = np.asarray(bipartite["coordinates"])
        assert (axes**2).sum(axis=0) == pytest.approx([2.0, 2.0], abs=1e-9)  # no axis was zeroed by the -4

        parent = ordination(generic, method="pcoa", n_components=2)
        assert child["generic"]["all_eigenvalues"] == pytest.approx(parent["all_eigenvalues"], abs=1e-9)
        assert _pairwise(child["generic"]["coordinates"]) == pytest.approx(_pairwise(parent["coordinates"]), abs=1e-6)


class TestDistanceMatrixValidation:
    """A matrix that is not a distance matrix is rejected instead of being silently ordinated."""

    GOOD = [[0.0, 1.0, 2.0], [1.0, 0.0, 1.5], [2.0, 1.5, 0.0]]

    def _with(self, i: int, j: int, value: float, symmetric: bool = True) -> list[list[float]]:
        matrix = [list(row) for row in self.GOOD]
        matrix[i][j] = value
        if symmetric:
            matrix[j][i] = value
        return matrix

    @pytest.mark.parametrize("method", ["pcoa", "nmds"])
    def test_rejects_invalid_matrices(self, method: str) -> None:
        with pytest.raises(ValueError, match="symmetric"):
            ordination(self._with(0, 1, 5.0, symmetric=False), method=method)
        with pytest.raises(ValueError, match="negative"):
            ordination(self._with(0, 2, -1.0), method=method)
        with pytest.raises(ValueError, match="not finite"):
            ordination(self._with(1, 2, math.nan), method=method)
        with pytest.raises(ValueError, match="not finite"):
            ordination(self._with(1, 2, math.inf), method=method)
        with pytest.raises(ValueError, match="zero diagonal"):
            ordination(self._with(1, 1, 0.5), method=method)
        with pytest.raises(ValueError, match="square"):
            ordination([[0.0, 1.0, 2.0], [1.0, 0.0]], method=method)
        with pytest.raises(ValueError, match="empty"):
            ordination([], method=method)

    def test_numpy_input_and_parameter_checks(self) -> None:
        assert ordination(np.asarray(self.GOOD), method="pcoa")["coordinates"] == ordination(self.GOOD)["coordinates"]
        with pytest.raises(ValueError, match="n_components"):
            ordination(self.GOOD, n_components=0)
        with pytest.raises(ValueError, match="n_init"):
            ordination(self.GOOD, method="nmds", n_init=-1)


# ---------------------------------------------------------------------------
# NMDS
# ---------------------------------------------------------------------------


class TestNmdsIsKruskalNonMetric:
    """Isotonic regression of distances on dissimilarities, minimising stress-1."""

    def test_reports_stress_and_no_variance_or_eigenvalues(self) -> None:
        result = ordination(_gradient_bray_curtis(0).tolist(), method="nmds", n_components=2, seed=3)
        assert "variance_explained" not in result  # was 1 - stress, which is not a share of anything
        assert "eigenvalues" not in result
        assert 0.0 <= result["stress"] < 0.1
        assert result["stress_by_start"][0] >= result["stress"] - 1e-15
        assert min(result["stress_by_start"]) == pytest.approx(result["stress"], rel=1e-9, abs=1e-12)

    def test_reported_stress_is_the_kruskal_stress_of_the_returned_configuration(self) -> None:
        dissimilarities = _gradient_bray_curtis(1)
        result = ordination(dissimilarities.tolist(), method="nmds", n_components=2, seed=5)
        independent = _kruskal_stress1(dissimilarities, result["coordinates"])
        assert result["stress"] == pytest.approx(independent, rel=1e-6, abs=1e-10)

    def test_recovers_a_planted_configuration_up_to_a_monotone_distortion(self) -> None:
        planted = _planted()
        distances = _pairwise(planted)
        for transform in (lambda d: d, lambda d: d**3, lambda d: np.exp(d / 3.0) - 1.0, np.sqrt):
            dissimilarities = transform(distances)
            result = ordination(dissimilarities.tolist(), method="nmds", n_components=2, seed=1)
            assert result["stress"] < 1e-3
            assert _procrustes_error(planted, result["coordinates"]) < 0.06

    def test_metric_mds_cannot_recover_a_distorted_configuration(self) -> None:
        """What the old implementation did (a fit to the dissimilarities themselves) fails on this input."""
        planted = _planted()
        dissimilarities = _pairwise(planted) ** 3
        classical = ordination(dissimilarities.tolist(), method="pcoa", n_components=2)
        assert _procrustes_error(planted, classical["coordinates"]) > 0.2

    def test_random_starts_depend_only_on_the_rank_order_of_the_dissimilarities(self) -> None:
        """A monotone rescaling changes nothing for a rank-based method (the classical start is not rank-based)."""
        distances = _gradient_bray_curtis(1)
        base = ordination(distances.tolist(), method="nmds", n_components=2, seed=2)
        for transform in (lambda d: d**3, lambda d: np.exp(d) - 1.0, lambda d: np.where(d > 0, 5.0 * d + 0.25, 0.0)):
            other = ordination(transform(distances).tolist(), method="nmds", n_components=2, seed=2)
            assert other["stress_by_start"][1:] == pytest.approx(base["stress_by_start"][1:], abs=1e-9)

    @pytest.mark.parametrize("data_seed", [0, 1, 2, 3])
    def test_stress_is_not_worse_than_scikit_learn(self, data_seed: int) -> None:
        manifold = pytest.importorskip("sklearn.manifold")
        dissimilarities = _gradient_bray_curtis(data_seed)
        parameters = inspect.signature(manifold.MDS.__init__).parameters
        options: dict[str, object] = {"n_components": 2, "n_init": 4, "max_iter": 2000, "eps": 1e-9}
        if "metric_mds" in parameters:  # scikit-learn >= 1.8
            options.update(metric_mds=False, metric="precomputed", init="random")
        else:
            options.update(metric=False, dissimilarity="precomputed")
        reference = []
        for seed in (0, 1, 2):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                embedding = manifold.MDS(random_state=seed, **options).fit_transform(dissimilarities)
            reference.append(_kruskal_stress1(dissimilarities, embedding))
        ours = ordination(dissimilarities.tolist(), method="nmds", n_components=2, seed=7)
        assert ours["stress"] <= min(reference) + 1e-3
        assert _kruskal_stress1(dissimilarities, ours["coordinates"]) <= min(reference) + 1e-3

    def test_several_starts_beat_a_single_bad_one(self) -> None:
        dissimilarities = _gradient_bray_curtis(3).tolist()
        result = ordination(dissimilarities, method="nmds", n_components=2, seed=7)
        assert len(result["stress_by_start"]) == 1 + 5  # a classical start and five random ones
        assert max(result["stress_by_start"]) > 1.2 * result["stress"]  # some start ends in a worse local minimum
        classical_only = ordination(dissimilarities, method="nmds", n_components=2, n_init=0)
        assert len(classical_only["stress_by_start"]) == 1
        more = ordination(dissimilarities, method="nmds", n_components=2, seed=7, n_init=9)
        assert len(more["stress_by_start"]) == 10
        assert more["stress"] <= result["stress"] + 1e-12  # the first five random starts are shared

    def test_ties_in_the_dissimilarities_do_not_force_equal_distances(self) -> None:
        """Bray-Curtis is exactly 1 for samples that share nothing; those pairs may sit at different distances."""
        counts = np.zeros((9, 9))
        for i in range(9):
            counts[i, i] = 40
            counts[i, (i + 1) % 9] = 15 * (i % 3)
        numerator = np.abs(counts[:, None, :] - counts[None, :, :]).sum(axis=2)
        denominator = (counts[:, None, :] + counts[None, :, :]).sum(axis=2)
        dissimilarities = numerator / denominator
        assert (dissimilarities == 1.0).sum() > 20  # many tied pairs
        result = ordination(dissimilarities.tolist(), method="nmds", n_components=2, seed=1)
        assert 0.0 <= result["stress"] < 0.15
        assert _kruskal_stress1(dissimilarities, result["coordinates"]) == pytest.approx(
            result["stress"], rel=1e-6, abs=1e-10
        )

    def test_configuration_is_centred_and_on_principal_axes(self) -> None:
        coords = np.asarray(ordination(_gradient_bray_curtis(2).tolist(), method="nmds", n_components=2)["coordinates"])
        assert coords.mean(axis=0) == pytest.approx([0.0, 0.0], abs=1e-9)
        assert float(coords[:, 0] @ coords[:, 1]) == pytest.approx(0.0, abs=1e-8)
        assert (coords[:, 0] ** 2).sum() >= (coords[:, 1] ** 2).sum()

    def test_rms_distance_matches_the_dissimilarities(self) -> None:
        dissimilarities = _gradient_bray_curtis(0)
        coords = ordination(dissimilarities.tolist(), method="nmds", n_components=2)["coordinates"]
        iu = np.triu_indices(len(dissimilarities), 1)
        rms_config = math.sqrt(float((_pairwise(coords)[iu] ** 2).mean()))
        rms_input = math.sqrt(float((dissimilarities[iu] ** 2).mean()))
        assert rms_config == pytest.approx(rms_input, rel=1e-9)

    def test_degenerate_input_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="at least 3"):
            ordination([[0.0, 1.0], [1.0, 0.0]], method="nmds")
        with pytest.raises(ValueError, match="zero"):
            ordination([[0.0] * 4 for _ in range(4)], method="nmds")


class TestNmdsSeeding:
    """The seed picks the random starts; the best of them is the same configuration up to rotation."""

    DISSIMILARITIES = _gradient_bray_curtis(0).tolist()

    def test_same_seed_is_bit_identical(self) -> None:
        assert ordination(self.DISSIMILARITIES, method="nmds", seed=4) == ordination(
            self.DISSIMILARITIES, method="nmds", seed=4
        )

    def test_other_seed_changes_the_random_starts_not_the_solution(self) -> None:
        first = ordination(self.DISSIMILARITIES, method="nmds", seed=1)
        other = ordination(self.DISSIMILARITIES, method="nmds", seed=2)
        assert first["stress_by_start"][0] == pytest.approx(other["stress_by_start"][0], rel=1e-3)  # classical start
        assert first["stress_by_start"][1:] != other["stress_by_start"][1:]  # different random starts
        assert other["stress"] == pytest.approx(first["stress"], rel=1e-3)
        assert _procrustes_error(first["coordinates"], other["coordinates"]) < 0.02

    def test_seed_none_draws_fresh_entropy_but_stays_valid(self) -> None:
        result = ordination(self.DISSIMILARITIES, method="nmds", seed=None)
        assert 0.0 <= result["stress"] < 0.1


# ---------------------------------------------------------------------------
# plot_ordination
# ---------------------------------------------------------------------------


def _as_nested(matrix: list[list[float]]) -> dict[str, dict[str, float]]:
    names = [f"s{i}" for i in range(len(matrix))]
    return {a: {b: matrix[i][j] for j, b in enumerate(names)} for i, a in enumerate(names)}


def _notes(fig: object) -> list[str]:
    """Text on the axes other than the sample labels (those are annotations)."""
    return [t.get_text() for t in fig.axes[0].texts if not isinstance(t, Annotation)]  # type: ignore[attr-defined]


def _plotted_points(fig: object, n_samples: int) -> object:
    ax = fig.axes[0]  # type: ignore[attr-defined]
    points = np.vstack([np.asarray(c.get_offsets()).reshape(-1, 2) for c in ax.collections])
    assert points.shape == (n_samples, 2)
    return points


class TestPlotOrdinationPcoa:
    """The plotted PCoA is the eigh decomposition, not power iteration."""

    def test_large_negative_eigenvalue_does_not_hijack_an_axis(self) -> None:
        # K(4,4) has eigenvalues 2 (x6), 0 and -4: power iteration converged to the -4 axis and zeroed it
        fig = plot_ordination(_as_nested(_complete_bipartite(4)), method="PCoA")
        try:
            points = _plotted_points(fig, 8)
            assert (points**2).sum(axis=0) == pytest.approx([2.0, 2.0], abs=1e-6)  # eigenvalue = sum of squares
            assert float(points[:, 0] @ points[:, 1]) == pytest.approx(0.0, abs=1e-6)
        finally:
            plt.close(fig)

    def test_axes_are_labelled_with_the_share_of_positive_eigenvalue_variance(self) -> None:
        fig = plot_ordination(_as_nested(_complete_bipartite(4)), method="PCoA")
        try:
            ax = fig.axes[0]
            assert ax.get_xlabel() == "PC1 (16.7%)"
            assert ax.get_ylabel() == "PC2 (16.7%)"
            assert len(_notes(fig)) == 1
            assert "1 negative eigenvalues" in _notes(fig)[0] and "25%" in _notes(fig)[0]
        finally:
            plt.close(fig)

    def test_euclidean_input_has_no_negative_eigenvalue_note(self) -> None:
        points = _planted(8, seed=2)
        fig = plot_ordination(_as_nested(_pairwise(points).tolist()), method="PCoA")
        try:
            assert _notes(fig) == []
            assert _procrustes_error(points, _plotted_points(fig, 8)) < 1e-6
        finally:
            plt.close(fig)

    def test_pcoa_does_not_depend_on_the_seed(self) -> None:
        nested = _as_nested(_complete_bipartite(3))
        figs = [plot_ordination(nested, method="PCoA", seed=s) for s in (1, 2)]
        try:
            np.testing.assert_allclose(_plotted_points(figs[0], 6), _plotted_points(figs[1], 6))
        finally:
            for f in figs:
                plt.close(f)


class TestPlotOrdinationNmds:
    def test_stress_is_shown_on_the_plot(self) -> None:
        dissimilarities = _gradient_bray_curtis(0)
        fig = plot_ordination(_as_nested(dissimilarities.tolist()), method="NMDS", seed=3)
        try:
            ax = fig.axes[0]
            expected = ordination(dissimilarities.tolist(), method="nmds", n_components=2, seed=3)["stress"]
            assert _notes(fig) == [f"stress = {expected:.3f}"]
            assert ax.get_xlabel() == "NMDS1"
        finally:
            plt.close(fig)

    def test_seed_changes_the_starts_and_is_reproducible(self) -> None:
        nested = _as_nested(_gradient_bray_curtis(1).tolist())
        figs = [plot_ordination(nested, method="NMDS", seed=s) for s in (1, 1, 2)]
        try:
            first, repeat, other = (_plotted_points(f, 22) for f in figs)
            np.testing.assert_allclose(first, repeat)
            assert _procrustes_error(first, other) < 0.02  # the best of the starts is the same solution
        finally:
            for f in figs:
                plt.close(f)


class TestPlotOrdinationDistances:
    """A distance that is not given is an error, never a distance of zero."""

    NESTED = _as_nested(_complete_bipartite(2))

    def test_missing_distance_is_an_error(self) -> None:
        incomplete = {a: dict(row) for a, row in self.NESTED.items()}
        del incomplete["s0"]["s3"]
        del incomplete["s3"]["s0"]
        for method in ("PCoA", "NMDS"):
            with pytest.raises(ValueError, match="no distance between 's0' and 's3'"):
                plot_ordination(incomplete, method=method)

    def test_one_triangle_is_enough(self) -> None:
        upper = {a: {b: d for b, d in row.items() if b > a} for a, row in self.NESTED.items()}
        fig_full = plot_ordination(self.NESTED, method="PCoA")
        fig_upper = plot_ordination(upper, method="PCoA")
        try:
            np.testing.assert_allclose(_plotted_points(fig_full, 4), _plotted_points(fig_upper, 4), atol=1e-9)
        finally:
            plt.close(fig_full)
            plt.close(fig_upper)

    def test_inconsistent_triangles_are_an_error(self) -> None:
        skewed = {a: dict(row) for a, row in self.NESTED.items()}
        skewed["s0"]["s1"] = 9.0
        with pytest.raises(ValueError, match="symmetric"):
            plot_ordination(skewed, method="PCoA")

    def test_samples_named_only_in_inner_keys_are_not_dropped(self) -> None:
        upper = {"s0": {"s1": 1.0, "s2": 2.0}, "s1": {"s2": 1.5}}  # s2 never appears as an outer key
        fig = plot_ordination(upper, method="PCoA")
        try:
            assert _plotted_points(fig, 3).shape == (3, 2)
        finally:
            plt.close(fig)


# ---------------------------------------------------------------------------
# PERMANOVA on perfectly separated groups
# ---------------------------------------------------------------------------


class TestPermanovaPerfectSeparation:
    """No within-group variation gave pseudo-F = 0 and so p = 1 for the most separated design possible."""

    @staticmethod
    def _two_tight_groups(size: int = 5) -> tuple[list[list[float]], list[str]]:
        n = 2 * size
        matrix = [[0.0 if (i < size) == (j < size) else 1.0 for j in range(n)] for i in range(n)]
        return matrix, ["A"] * size + ["B"] * size

    def test_infinite_pseudo_f_and_small_p(self) -> None:
        matrix, groups = self._two_tight_groups()
        result = permanova(matrix, groups, n_permutations=199, seed=1)
        assert result["pseudo_f"] == math.inf
        assert result["r_squared"] == pytest.approx(1.0)
        # only the observed partition and its mirror image reproduce F = inf: 2 of 252 labelings
        assert result["p_value"] < 0.05

    def test_all_zero_distances_are_undefined(self) -> None:
        result = permanova([[0.0] * 4 for _ in range(4)], ["A", "A", "B", "B"], n_permutations=9, seed=1)
        assert math.isnan(result["pseudo_f"])
        assert math.isnan(result["p_value"])
        assert math.isnan(result["r_squared"])

    def test_no_residual_degrees_of_freedom_is_an_error(self) -> None:
        matrix = [[0.0, 1.0, 2.0], [1.0, 0.0, 1.5], [2.0, 1.5, 0.0]]
        with pytest.raises(ValueError, match="more samples"):
            permanova(matrix, ["A", "B", "C"])

    def test_ordinary_designs_are_unchanged(self) -> None:
        matrix = [
            [0.0, 0.1, 0.8, 0.9],
            [0.1, 0.0, 0.7, 0.85],
            [0.8, 0.7, 0.0, 0.1],
            [0.9, 0.85, 0.1, 0.0],
        ]
        result = permanova(matrix, ["A", "A", "B", "B"], n_permutations=99, seed=42)
        ss_total = sum(matrix[i][j] ** 2 for i in range(4) for j in range(i + 1, 4)) / 4
        ss_within = (0.1**2) / 2 + (0.1**2) / 2
        assert result["pseudo_f"] == pytest.approx(((ss_total - ss_within) / 1) / (ss_within / 2), rel=1e-12)
        assert result["r_squared"] == pytest.approx((ss_total - ss_within) / ss_total, rel=1e-12)


def test_random_generator_state_is_not_shared_between_ordinations() -> None:
    """The module draws from its own seeded generator, never the global one."""
    random.seed(123)
    before = random.getstate()
    ordination(_gradient_bray_curtis(0).tolist(), method="nmds", seed=5)
    assert random.getstate() == before
