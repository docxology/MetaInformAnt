"""Tests for methylation simulation.

Tests cover simulate_methylation dataset structure, ground-truth DMR
injection, value ranges, determinism with seeds, and calculate_dmr_statistics.

All tests use real implementations (real-implementation policy) and small
configurations to keep runtime fast.
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.simulation.methylation.simulator import (
    MethylationSimulationConfig,
    calculate_dmr_statistics,
    simulate_methylation,
)


def _small_config(**overrides) -> MethylationSimulationConfig:
    """Return a small fast configuration for testing."""
    defaults = dict(
        n_samples=10,
        n_cpg_islands=3,
        n_gene_body_regions=2,
        n_promoters=2,
        dmr_fraction=0.3,
        random_seed=42,
    )
    defaults.update(overrides)
    return MethylationSimulationConfig(**defaults)


class TestSimulateMethylation:
    """Test methylation dataset simulation."""

    def test_dataset_shapes(self):
        """Test that beta values, labels, and masks have consistent shapes."""
        config = _small_config()
        dataset = simulate_methylation(config)
        n_sites = config.n_cpg_islands * 50 + config.n_gene_body_regions * 20 + config.n_promoters * 30

        assert dataset.beta_values.shape == (n_sites, config.n_samples)
        assert len(dataset.site_types) == n_sites
        assert len(dataset.dmr_mask) == n_sites
        assert len(dataset.group_labels) == config.n_samples
        assert dataset.config is config

    def test_group_labels_split_in_half(self):
        """Test that group labels are 0 for the first half and 1 for the rest."""
        config = _small_config(n_samples=7)
        dataset = simulate_methylation(config)
        labels = dataset.group_labels.tolist()
        assert labels == [0, 0, 0, 1, 1, 1, 1]

    def test_site_types_composition(self):
        """Test that region types are recorded in generation order."""
        config = _small_config()
        dataset = simulate_methylation(config)
        types = dataset.site_types
        assert types[: config.n_cpg_islands * 50] == ["island"] * (config.n_cpg_islands * 50)
        assert set(types) == {"island", "body", "promoter"}

    def test_beta_values_within_unit_interval(self):
        """Test that all beta values are clipped to [0, 1]."""
        dataset = simulate_methylation(_small_config(noise_std=0.5))
        assert dataset.beta_values.min() >= 0.0
        assert dataset.beta_values.max() <= 1.0

    def test_dmr_mask_matches_fraction(self):
        """Test that DMRs cover whole regions and match dmr_fraction."""
        config = _small_config(dmr_fraction=0.5)
        dataset = simulate_methylation(config)
        region_sizes = [50] * config.n_cpg_islands + [20] * config.n_gene_body_regions + [30] * config.n_promoters
        assert sum(region_sizes) == len(dataset.site_types)

        # The mask must be constant within each region
        n_dmr_regions = 0
        offset = 0
        for size in region_sizes:
            region_mask = dataset.dmr_mask[offset : offset + size]
            assert region_mask.all() or not region_mask.any()
            if region_mask.any():
                n_dmr_regions += 1
            offset += size
        assert offset == len(dataset.dmr_mask)
        assert n_dmr_regions == int(len(region_sizes) * config.dmr_fraction)

    def test_dmr_sites_shift_group_means(self):
        """Test that DMR sites show a mean shift between groups."""
        config = _small_config(dmr_fraction=1.0, dmr_effect_size=0.5)
        dataset = simulate_methylation(config)
        g0 = dataset.beta_values[:, dataset.group_labels == 0].mean(axis=1)
        g1 = dataset.beta_values[:, dataset.group_labels == 1].mean(axis=1)
        # Every region is a DMR, so group means should differ somewhere
        assert np.any(np.abs(g0 - g1) > 0.1)

    def test_deterministic_with_seed(self):
        """Test that a fixed seed reproduces the dataset exactly."""
        config = _small_config()
        dataset_a = simulate_methylation(config)
        dataset_b = simulate_methylation(config)
        assert np.array_equal(dataset_a.beta_values, dataset_b.beta_values)
        assert np.array_equal(dataset_a.dmr_mask, dataset_b.dmr_mask)

    def test_default_config_used_when_none(self):
        """Test that the default configuration is applied when config is None."""
        dataset = simulate_methylation()
        assert isinstance(dataset.config, MethylationSimulationConfig)
        assert dataset.config.n_samples == 100

    def test_dmrs_do_not_exceed_available_regions(self):
        """Test that dmr_fraction of 1.0 marks every region as a DMR."""
        config = _small_config(dmr_fraction=1.0)
        dataset = simulate_methylation(config)
        assert dataset.dmr_mask.all()


class TestCalculateDmrStatistics:
    """Test DMR summary statistics."""

    def test_statistics_keys_and_values(self):
        """Test that statistics match the underlying dataset."""
        dataset = simulate_methylation(_small_config())
        stats = calculate_dmr_statistics(dataset)

        expected_keys = {
            "n_sites",
            "n_samples",
            "n_dmr_sites",
            "mean_beta_overall",
            "mean_beta_group0",
            "mean_beta_group1",
        }
        assert set(stats) == expected_keys
        assert stats["n_sites"] == dataset.beta_values.shape[0]
        assert stats["n_samples"] == dataset.beta_values.shape[1]
        assert stats["n_dmr_sites"] == int(dataset.dmr_mask.sum())
        assert stats["mean_beta_overall"] == pytest.approx(float(dataset.beta_values.mean()))
        g0 = dataset.group_labels == 0
        g1 = dataset.group_labels == 1
        assert stats["mean_beta_group0"] == pytest.approx(float(dataset.beta_values[:, g0].mean()))
        assert stats["mean_beta_group1"] == pytest.approx(float(dataset.beta_values[:, g1].mean()))
