"""Round-4 test-depth tests for singlecell data (preprocessing, integration) and doublet detection.

Covers previously untested public APIs:
- data.preprocessing: identify_highly_variable_genes, remove_batch_effects
- data.integration: combat_integration, integrate_multiple_batches (mnn/scanorama),
  evaluate_integration_quality
- doublet.detection: detect_doublets, DoubletResult

All tests use real computation on synthetic data (zero-mocks policy).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from metainformant.singlecell.data.integration import (
    combat_integration,
    evaluate_integration_quality,
    integrate_multiple_batches,
)
from metainformant.singlecell.data.preprocessing import (
    SingleCellData,
    identify_highly_variable_genes,
    log_transform,
    normalize_counts,
    remove_batch_effects,
)
from metainformant.singlecell.doublet.detection import DoubletResult, detect_doublets

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_batched_counts(n_per_batch: int = 30, n_genes: int = 50, seed: int = 0) -> SingleCellData:
    """Poisson counts with a batch offset between two groups."""
    rng = np.random.RandomState(seed)
    X = np.vstack([rng.poisson(2.0, (n_per_batch, n_genes)), rng.poisson(8.0, (n_per_batch, n_genes))]).astype(float)
    obs = pd.DataFrame(
        {"batch": ["b1"] * n_per_batch + ["b2"] * n_per_batch},
        index=[f"c{i}" for i in range(2 * n_per_batch)],
    )
    var = pd.DataFrame(index=[f"g{i}" for i in range(n_genes)])
    return SingleCellData(X, obs=obs, var=var)


def _make_batch_object(batch_id: int, n_cells: int = 20, n_genes: int = 15) -> SingleCellData:
    rng = np.random.RandomState(batch_id)
    X = rng.normal(size=(n_cells, n_genes)) + batch_id
    obs = pd.DataFrame({"batch": [f"b{batch_id}"] * n_cells}, index=[f"b{batch_id}_c{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=[f"g{i}" for i in range(n_genes)])
    return SingleCellData(X, obs=obs, var=var)


# ---------------------------------------------------------------------------
# identify_highly_variable_genes
# ---------------------------------------------------------------------------


class TestIdentifyHighlyVariableGenes:
    def _normalized(self) -> SingleCellData:
        return log_transform(normalize_counts(_make_batched_counts()))

    def test_marks_highly_variable_genes(self) -> None:
        data = identify_highly_variable_genes(self._normalized(), n_top_genes=10)
        assert "highly_variable" in data.var.columns
        assert data.var["highly_variable"].sum() == 10
        for col in ("gene_mean", "gene_variance"):
            assert col in data.var.columns

    def test_invalid_flavor_raises(self) -> None:
        with pytest.raises(Exception):
            identify_highly_variable_genes(self._normalized(), flavor="bogus")

    def test_rejects_non_singlecelldata(self) -> None:
        with pytest.raises(Exception):
            identify_highly_variable_genes("not data")  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# remove_batch_effects
# ---------------------------------------------------------------------------


class TestRemoveBatchEffects:
    def test_reduces_batch_mean_shift(self) -> None:
        data = _make_batched_counts()
        before = data.X[:30].mean(0) - data.X[30:].mean(0)
        corrected = remove_batch_effects(data.copy(), batch_key="batch")
        after = corrected.X[:30].mean(0) - corrected.X[30:].mean(0)
        assert corrected.X.shape == data.X.shape
        assert np.abs(after).mean() < np.abs(before).mean()

    def test_missing_batch_key_raises(self) -> None:
        with pytest.raises(Exception):
            remove_batch_effects(_make_batched_counts(), batch_key="missing")


# ---------------------------------------------------------------------------
# combat_integration
# ---------------------------------------------------------------------------


class TestCombatIntegration:
    def test_preserves_shape(self) -> None:
        data = _make_batched_counts()
        result = combat_integration(data.copy(), batch_key="batch")
        assert result.X.shape == data.X.shape

    def test_with_covariates(self) -> None:
        data = _make_batched_counts()
        data.obs["condition"] = ["ctrl"] * 30 + ["treated"] * 30
        result = combat_integration(data.copy(), batch_key="batch", covariates=["condition"])
        assert result.X.shape == data.X.shape


# ---------------------------------------------------------------------------
# integrate_multiple_batches (regression: mnn/scanorama branches previously
# dropped batch_key, raising TypeError)
# ---------------------------------------------------------------------------


class TestIntegrateMultipleBatches:
    @pytest.mark.parametrize("method", ["mnn", "scanorama"])
    def test_mnn_and_scanorama_concatenate(self, method: str) -> None:
        batches = [_make_batch_object(0), _make_batch_object(1)]
        integrated = integrate_multiple_batches(batches, integration_method=method)
        assert integrated.n_obs == 40
        assert integrated.n_vars == 15
        assert "batch" in integrated.obs.columns
        # scanorama/mnn rewrite batch labels to positional batch_0/batch_1
        assert set(integrated.obs["batch"].unique()) == {"batch_0", "batch_1"}

    def test_unsupported_method_raises(self) -> None:
        with pytest.raises(Exception):
            integrate_multiple_batches([_make_batch_object(0)], integration_method="bogus")

    def test_empty_list_raises(self) -> None:
        with pytest.raises(Exception):
            integrate_multiple_batches([], integration_method="mnn")


# ---------------------------------------------------------------------------
# evaluate_integration_quality
# ---------------------------------------------------------------------------


class TestEvaluateIntegrationQuality:
    def test_quality_metrics_structure(self) -> None:
        batches = [_make_batch_object(0), _make_batch_object(1)]
        integrated = integrate_multiple_batches(batches, integration_method="mnn")
        quality = evaluate_integration_quality(integrated, batch_key="batch")
        for key in (
            "n_batches",
            "batch_sizes",
            "mean_batch_mixing_score",
            "mean_within_batch_distance",
            "mean_between_batch_distance",
            "batch_distance_ratio",
        ):
            assert key in quality
        assert quality["n_batches"] == 2
        assert sum(quality["batch_sizes"]) == 40


# ---------------------------------------------------------------------------
# detect_doublets
# ---------------------------------------------------------------------------


class TestDetectDoublets:
    def _counts(self, seed: int = 0) -> np.ndarray:
        rng = np.random.RandomState(seed)
        return rng.poisson(np.random.RandomState(seed + 1).uniform(0.1, 5.0, (60, 30)))

    def test_result_structure(self) -> None:
        result = detect_doublets(self._counts(), random_state=0)
        assert isinstance(result, DoubletResult)
        assert len(result.scores) == 60
        assert 0.0 <= result.threshold <= 1.0
        assert result.n_doublets == int(np.sum(result.predicted_doublets))
        assert 0.0 <= result.doublet_rate <= 1.0
        assert len(result.synthetic_scores) > 0

    def test_deterministic_with_seed(self) -> None:
        r1 = detect_doublets(self._counts(), random_state=42)
        r2 = detect_doublets(self._counts(), random_state=42)
        assert list(r1.scores) == list(r2.scores)
        assert r1.n_doublets == r2.n_doublets

    def test_explicit_synthetic_count(self) -> None:
        result = detect_doublets(self._counts(), n_synthetic=50, random_state=0)
        assert len(result.synthetic_scores) == 50
