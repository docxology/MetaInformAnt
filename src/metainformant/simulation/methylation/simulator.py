"""Simulate DNA methylation patterns.

Generates synthetic methylation data for CpG islands, gene bodies, and
promoter regions using biologically motivated stochastic models.
Useful for benchmarking differential methylation analysis tools.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class CpGIslandConfig:
    """Configuration for a simulated CpG island.

    Attributes:
        n_cpgs: Number of CpG sites in the island.
        mean_methylation: Average beta value.
        spatial_correlation: Correlation between adjacent CpG sites.
    """

    n_cpgs: int = 50
    mean_methylation: float = 0.1
    spatial_correlation: float = 0.8


@dataclass
class MethylationSimulationConfig:
    """Configuration for methylation simulation.

    Attributes:
        n_samples: Number of samples to simulate.
        n_cpg_islands: Number of CpG islands.
        n_gene_body_regions: Number of gene body regions.
        n_promoters: Number of promoter regions.
        noise_std: Standard deviation of beta-value noise.
        dmr_fraction: Fraction of regions that are differentially methylated.
        dmr_effect_size: Mean effect size (delta beta) for DMRs.
        random_seed: Seed for reproducibility.
    """

    n_samples: int = 100
    n_cpg_islands: int = 200
    n_gene_body_regions: int = 300
    n_promoters: int = 200
    noise_std: float = 0.05
    dmr_fraction: float = 0.1
    dmr_effect_size: float = 0.3
    random_seed: int | None = None


@dataclass
class MethylationDataset:
    """A simulated methylation dataset.

    Attributes:
        beta_values: 2D array (sites × samples) of beta values in [0, 1].
        group_labels: 1D array of group assignments (0 or 1).
        site_types: List of site type labels ('island', 'body', 'promoter').
        dmr_mask: Boolean mask of truly differentially methylated sites.
        config: Configuration used to generate this dataset.
    """

    beta_values: np.ndarray
    group_labels: np.ndarray
    site_types: list[str]
    dmr_mask: np.ndarray
    config: MethylationSimulationConfig


def _simulate_correlated_beta(
    n_sites: int,
    n_samples: int,
    mean: float,
    correlation: float,
    noise_std: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Generate spatially correlated beta values for a region."""
    base = rng.beta(
        a=max(mean * 10, 0.1),
        b=max((1 - mean) * 10, 0.1),
        size=(1, n_samples),
    )
    offsets = np.zeros((n_sites, n_samples))
    offsets[0] = rng.normal(0, noise_std, n_samples)
    for i in range(1, n_sites):
        offsets[i] = correlation * offsets[i - 1] + rng.normal(
            0, noise_std * np.sqrt(1 - correlation**2), n_samples
        )
    result = base + offsets
    return np.clip(result, 0.0, 1.0)


def simulate_methylation(
    config: MethylationSimulationConfig | None = None,
) -> MethylationDataset:
    """Simulate a complete methylation dataset with ground truth DMRs.

    Generates beta values for CpG islands (typically hypomethylated),
    gene bodies (moderately methylated), and promoters (variable), then
    injects differential methylation between two groups.

    Args:
        config: Simulation configuration. Defaults to MethylationSimulationConfig().

    Returns:
        MethylationDataset with simulated data and ground truth.
    """
    if config is None:
        config = MethylationSimulationConfig()

    rng = np.random.default_rng(config.random_seed)

    n_group_a = config.n_samples // 2
    group_labels = np.array([0] * n_group_a + [1] * (config.n_samples - n_group_a))

    site_configs: list[tuple[str, int, float, float]] = []
    for _ in range(config.n_cpg_islands):
        site_configs.append(("island", 50, 0.1, 0.8))
    for _ in range(config.n_gene_body_regions):
        site_configs.append(("body", 20, 0.6, 0.5))
    for _ in range(config.n_promoters):
        site_configs.append(("promoter", 30, 0.3, 0.7))

    total_regions = len(site_configs)
    n_dmrs = int(total_regions * config.dmr_fraction)
    dmr_indices = set(rng.choice(total_regions, size=n_dmrs, replace=False))

    all_betas: list[np.ndarray] = []
    all_types: list[str] = []
    all_dmr: list[bool] = []

    for region_idx, (site_type, n_sites, mean_meth, corr) in enumerate(site_configs):
        region_beta = _simulate_correlated_beta(
            n_sites, config.n_samples, mean_meth, corr, config.noise_std, rng
        )

        is_dmr = region_idx in dmr_indices
        if is_dmr:
            effect = config.dmr_effect_size * rng.choice([-1, 1])
            region_beta[:, group_labels == 1] += effect
            region_beta = np.clip(region_beta, 0.0, 1.0)

        all_betas.append(region_beta)
        all_types.extend([site_type] * n_sites)
        all_dmr.extend([is_dmr] * n_sites)

    return MethylationDataset(
        beta_values=np.vstack(all_betas),
        group_labels=group_labels,
        site_types=all_types,
        dmr_mask=np.array(all_dmr),
        config=config,
    )


def calculate_dmr_statistics(
    dataset: MethylationDataset,
) -> dict[str, float]:
    """Compute summary statistics on the simulated methylation data.

    Args:
        dataset: A MethylationDataset.

    Returns:
        Dict with keys: n_sites, n_samples, n_dmr_sites,
        mean_beta_overall, mean_beta_group0, mean_beta_group1.
    """
    g0_mask = dataset.group_labels == 0
    g1_mask = dataset.group_labels == 1
    return {
        "n_sites": dataset.beta_values.shape[0],
        "n_samples": dataset.beta_values.shape[1],
        "n_dmr_sites": int(dataset.dmr_mask.sum()),
        "mean_beta_overall": float(dataset.beta_values.mean()),
        "mean_beta_group0": float(dataset.beta_values[:, g0_mask].mean()),
        "mean_beta_group1": float(dataset.beta_values[:, g1_mask].mean()),
    }
