"""Multi-omic methods keep their randomness private, run on every supported NumPy, and do not advertise a no-op.

Regression tests for ``metainformant.multiomics.methods``:

* the pure-Python path of ``factorization.joint_nmf`` called ``random.seed(42)``, which re-seeds the *global*
  generator of the whole process (every other user of ``random`` in the program was silently reset); it uses a
  private ``random.Random(42)`` now, which draws the same numbers;
* ``clustering.consensus_clustering`` integrated the consensus CDF with ``np.trapezoid``, which exists only from
  NumPy 2.0 while ``pyproject.toml`` allows ``numpy>=1.24`` (on NumPy 1.x every call raised ``AttributeError``);
  its docstring also promised hierarchical clustering of the consensus matrix, the code runs spectral clustering;
* ``factorization.similarity_network_fusion`` documented ``alpha`` as the kernel bandwidth although it never used it;
* ``consensus_clustering`` accepted ``n_resamples=0`` (no resampling, so every pair of samples had consensus 0 and
  every k a perfect PAC of 0.0) and ``k=1`` (one cluster: consensus 1 everywhere, PAC 0.0, always "optimal"), and
  ``multi_omic_clustering`` accepted more clusters than samples (the spectral route returned one cluster per sample
  with a silhouette of 0.0, the k-means routes an ``IndexError``).

The source guards exist because the NumPy 1.x failure cannot be reproduced on the NumPy that the test environment has
(they fail on the old modules whatever NumPy is installed) and because no module may seed the global generator.
Real computation; no test doubles (real-implementation policy).
"""

from __future__ import annotations

import ast
import random
import warnings
from pathlib import Path

import numpy as np
import pytest

from metainformant.multiomics import methods as methods_package
from metainformant.multiomics.methods import clustering
from metainformant.multiomics.methods.clustering import consensus_clustering, multi_omic_clustering
from metainformant.multiomics.methods.factorization import (
    _py_random_matrix,
    similarity_network_fusion,
)

PACKAGE_ROOTS = [
    Path(methods_package.__file__).resolve().parents[1],  # metainformant/multiomics
    Path(methods_package.__file__).resolve().parents[2] / "metabolomics",
    Path(methods_package.__file__).resolve().parents[2] / "quality" / "batch",
]
# attribute names that exist in only some of the NumPy releases the package supports (numpy>=1.24)
NUMPY_VERSION_SPECIFIC_NAMES = {"trapezoid", "trapz"}


def _python_sources() -> list[Path]:
    return sorted(source for root in PACKAGE_ROOTS for source in root.rglob("*.py"))


class TestRandomnessStaysPrivate:
    def test_drawing_a_random_matrix_leaves_the_global_generator_alone(self) -> None:
        random.seed(2024)
        before = random.getstate()

        _py_random_matrix(3, 4, random.Random(42))

        assert random.getstate() == before

    def test_the_private_generator_draws_what_the_global_seed_42_used_to_draw(self) -> None:
        random.seed(42)
        legacy = [[random.random() * 0.5 + 0.01 for _ in range(4)] for _ in range(3)]  # the former implementation

        assert _py_random_matrix(3, 4, random.Random(42)) == legacy

    def test_the_values_are_positive_and_below_the_documented_ceiling(self) -> None:
        matrix = _py_random_matrix(20, 20, random.Random(1))

        assert all(0.01 <= value < 0.51 for row in matrix for value in row)

    def test_no_module_of_the_three_packages_seeds_a_global_generator(self) -> None:
        offences: list[str] = []
        for source in _python_sources():
            for node in ast.walk(ast.parse(source.read_text(), filename=str(source))):
                if not (
                    isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == "seed"
                ):
                    continue
                target = node.func.value
                is_stdlib_random = isinstance(target, ast.Name) and target.id == "random"
                is_numpy_global = (
                    isinstance(target, ast.Attribute)
                    and target.attr == "random"
                    and isinstance(target.value, ast.Name)
                    and target.value.id in {"np", "numpy"}
                )
                if is_stdlib_random or is_numpy_global:
                    offences.append(f"{source.name}:{node.lineno}")

        assert offences == [], f"global random state is re-seeded at {offences}; use a private generator"


class TestNumpyVersionCompatibility:
    def test_no_module_of_the_three_packages_calls_a_numpy_function_missing_from_supported_versions(self) -> None:
        offences: list[str] = []
        for source in _python_sources():
            for node in ast.walk(ast.parse(source.read_text(), filename=str(source))):
                if (
                    isinstance(node, ast.Attribute)
                    and node.attr in NUMPY_VERSION_SPECIFIC_NAMES
                    and isinstance(node.value, ast.Name)
                    and node.value.id in {"np", "numpy"}
                ):
                    offences.append(f"{source.name}:{node.lineno} np.{node.attr}")

        assert offences == [], f"NumPy-version-specific calls while numpy>=1.24 is supported: {offences}"

    def test_the_trapezoid_helper_is_numpys_own_rule(self) -> None:
        rng = np.random.default_rng(5)
        x = np.sort(rng.uniform(0, 1, size=40))
        y = np.arange(1, 41) / 40.0
        reference = getattr(np, "trapezoid", None) or getattr(np, "trapz")

        assert clustering._trapezoid_area(y, x) == float(reference(y, x))
        assert clustering._trapezoid_area(np.array([1.0, 3.0, 2.0]), np.array([0.0, 1.0, 3.0])) == 2.0 + 5.0

    def test_the_consensus_cdf_area_is_the_trapezoid_area_of_its_own_consensus_matrix(self) -> None:
        rng = np.random.default_rng(6)
        data = np.vstack([rng.normal(0, 1, (10, 4)), rng.normal(5, 1, (10, 4))])

        result = consensus_clustering(data, k_range=[2], n_resamples=15)

        consensus = np.array(result["consensus_matrix"])
        values = np.sort(consensus[np.triu_indices(20, k=1)])
        cdf = np.arange(1, values.size + 1) / values.size
        by_hand = float(np.sum(np.diff(values) * (cdf[1:] + cdf[:-1]) / 2.0))
        assert result["optimal_k"] == 2
        assert result["cdf_area"][2] == pytest.approx(by_hand, rel=1e-12)


class TestSimilarityNetworkFusionAlpha:
    @staticmethod
    def _networks() -> list[np.ndarray]:
        rng = np.random.default_rng(7)
        networks = []
        for _ in range(2):
            points = rng.normal(size=(12, 3))
            distances = np.linalg.norm(points[:, None] - points[None], axis=2)
            networks.append(np.exp(-distances))
        return networks

    def test_the_default_alpha_is_silent(self) -> None:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            similarity_network_fusion(self._networks(), k_neighbors=4, n_iter=3)

        assert [str(w.message) for w in caught if "alpha" in str(w.message)] == []

    def test_another_alpha_says_that_it_does_nothing(self) -> None:
        networks = self._networks()
        baseline = similarity_network_fusion(networks, k_neighbors=4, n_iter=3)

        with pytest.warns(UserWarning, match=r"alpha=0\.9 has no effect"):
            changed = similarity_network_fusion(networks, k_neighbors=4, n_iter=3, alpha=0.9)

        assert changed["fused_network"] == baseline["fused_network"]  # the warning is true


class TestClusteringArguments:
    """Arguments that used to produce a perfect-looking answer from no computation."""

    @staticmethod
    def _two_blobs(n: int = 20) -> np.ndarray:
        rng = np.random.default_rng(3)
        return np.vstack([rng.normal(0.0, 1.0, (n // 2, 4)), rng.normal(6.0, 1.0, (n // 2, 4))])

    @pytest.mark.parametrize("n_resamples", [0, -3])
    def test_consensus_clustering_needs_at_least_one_resample(self, n_resamples: int) -> None:
        # with none, every pair of samples has consensus 0, PAC = 0.0 ("perfect") for every k
        with pytest.raises(ValueError, match=r"n_resamples must be at least 1"):
            consensus_clustering(self._two_blobs(), k_range=[2, 3], n_resamples=n_resamples)

    @pytest.mark.parametrize("k_range", [[1, 2, 3], [0], range(1, 4)])
    def test_a_single_cluster_is_not_a_candidate_k(self, k_range: object) -> None:
        # k = 1 puts everyone together: consensus 1 everywhere, PAC = 0.0, so it always "won"
        with pytest.raises(ValueError, match=r"every k in k_range must be at least 2"):
            consensus_clustering(self._two_blobs(), k_range=k_range, n_resamples=5)  # type: ignore[arg-type]

    def test_the_default_range_still_picks_the_two_blobs(self) -> None:
        result = consensus_clustering(self._two_blobs(), n_resamples=20)

        assert result["optimal_k"] == 2

    @pytest.mark.parametrize("method", ["snf", "concatenation", "late_integration"])
    def test_more_clusters_than_samples_raise_for_every_method(self, method: str) -> None:
        views = {
            "rna": np.random.default_rng(1).normal(size=(4, 5)),
            "protein": np.random.default_rng(2).normal(size=(4, 3)),
        }

        with pytest.raises(ValueError, match=r"Cannot form n_clusters=5 clusters from 4 samples"):
            multi_omic_clustering(views, n_clusters=5, method=method)

    def test_as_many_clusters_as_samples_is_still_allowed_for_snf(self) -> None:
        views = {
            "rna": np.random.default_rng(1).normal(size=(4, 5)),
            "protein": np.random.default_rng(2).normal(size=(4, 3)),
        }

        result = multi_omic_clustering(views, n_clusters=4, method="snf")

        assert len(result["labels"]) == 4
