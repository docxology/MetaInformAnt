"""RNA expression simulation utilities for generating synthetic transcriptomic data.

This module provides functions for simulating RNA-seq count data using
negative binomial distributions and modeling differential expression patterns.
All functions support reproducible results through random seed control.
"""

from __future__ import annotations

import random
from typing import Dict, List, Optional, Tuple, Any
import numpy as np

from metainformant.core import logging, validation, errors

logger = logging.get_logger(__name__)


def simulate_counts_negative_binomial(
    n_samples: int, n_features: int, means: np.ndarray, dispersions: np.ndarray, rng: random.Random | None = None
) -> np.ndarray:
    """Simulate RNA-seq counts using negative binomial distribution.

    Args:
        n_samples: Number of samples/cells to simulate
        n_features: Number of genes/features
        means: Mean expression values for each feature (shape: n_features,)
        dispersions: Dispersion parameters for each feature (shape: n_features,)
        rng: Random number generator for reproducibility

    Returns:
        Simulated count matrix (shape: n_samples x n_features)

    Raises:
        ValueError: If array shapes don't match or values are invalid
    """
    validation.validate_range(n_samples, min_val=1, name="n_samples")
    validation.validate_range(n_features, min_val=1, name="n_features")
    validation.validate_type(means, np.ndarray, "means")
    validation.validate_type(dispersions, np.ndarray, "dispersions")

    if means.shape != (n_features,):
        raise errors.ValidationError(f"means must have shape ({n_features},), got {means.shape}")

    if dispersions.shape != (n_features,):
        raise errors.ValidationError(f"dispersion must have shape ({n_features},), got {dispersions.shape}")

    if np.any(means < 0):
        raise errors.ValidationError("All mean values must be non-negative")

    if np.any(dispersions <= 0):
        raise errors.ValidationError("All dispersion values must be positive")

    if rng is None:
        rng = random.Random()

    # Set numpy random seed for reproducibility
    np.random.seed(rng.randint(0, 2**32))

    # Negative binomial parameters
    # mean = n * p / (1-p), variance = n * p / (1-p)^2
    # dispersion = 1/n, where n is number of successes
    # For RNA-seq, we use the parameterization: mean = mu, variance = mu + mu^2 * phi
    # where phi is the dispersion parameter

    counts = np.zeros((n_samples, n_features), dtype=int)

    for i in range(n_features):
        mu = means[i]
        phi = dispersions[i]

        if mu == 0:
            # No expression
            counts[:, i] = 0
            continue

        # Convert to scipy negative binomial parameterization
        # NB(n, p) where mean = n * (1-p)/p, variance = n * (1-p)/p^2
        # For RNA-seq: n = 1/phi, p = 1/(1 + mu * phi)
        n_trials = 1.0 / phi
        p_success = 1.0 / (1.0 + mu * phi)

        # Generate negative binomial samples
        # Use numpy's negative binomial (which uses n, p parameterization)
        counts[:, i] = np.random.negative_binomial(n_trials, p_success, size=n_samples)

    return counts


def simulate_differential_expression(
    n_samples: int, n_features: int, fold_changes: np.ndarray, rng: random.Random | None = None
) -> Tuple[np.ndarray, np.ndarray]:
    """Simulate differential expression with specified fold changes.

    Args:
        n_samples: Total number of samples
        n_features: Number of genes/features
        fold_changes: Fold changes for differential genes (0 = no change, positive/negative = upregulation/downregulation)
        rng: Random number generator for reproducibility

    Returns:
        Tuple of (expression_matrix, group_labels) where group_labels indicates condition (0 or 1)

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_range(n_samples, min_val=2, name="n_samples")
    validation.validate_range(n_features, min_val=1, name="n_features")
    validation.validate_type(fold_changes, np.ndarray, "fold_changes")

    n_diff = len(fold_changes)
    if n_diff >= n_features:
        raise errors.ValidationError(f"Cannot have {n_diff} differentially expressed genes out of {n_features} total")

    if rng is None:
        rng = random.Random()

    # Set numpy random seed
    np.random.seed(rng.randint(0, 2**32))

    # Split samples into two groups
    group_sizes = [n_samples // 2, n_samples - n_samples // 2]
    group_labels = np.concatenate([np.zeros(group_sizes[0]), np.ones(group_sizes[1])])

    # Generate baseline expression levels
    baseline_means = np.random.exponential(10, n_features)  # Exponential distribution for gene means

    # Create dispersion parameters (higher for low-expressed genes)
    dispersions = 1.0 / (baseline_means + 1.0)  # Inverse relationship with mean
    dispersions = np.clip(dispersions, 0.1, 10.0)  # Reasonable bounds

    # Modify means for differentially expressed genes in group 1
    modified_means = baseline_means.copy()
    for i, fc in enumerate(fold_changes):
        gene_idx = i  # First n_diff genes are differentially expressed
        if fc > 0:
            # Upregulation in group 1
            modified_means[gene_idx] *= fc
        elif fc < 0:
            # Downregulation in group 1
            modified_means[gene_idx] *= abs(fc)

    # Generate expression for group 0 (baseline)
    group0_means = baseline_means
    group0_counts = simulate_counts_negative_binomial(group_sizes[0], n_features, group0_means, dispersions, rng=rng)

    # Generate expression for group 1 (modified)
    group1_means = modified_means
    group1_counts = simulate_counts_negative_binomial(group_sizes[1], n_features, group1_means, dispersions, rng=rng)

    # Combine into single matrix
    expression_matrix = np.vstack([group0_counts, group1_counts])

    return expression_matrix, group_labels.astype(int)


def simulate_bulk_rnaseq(
    n_samples: int,
    n_genes: int,
    *,
    library_sizes: Optional[np.ndarray] = None,
    gene_means: Optional[np.ndarray] = None,
    rng: random.Random | None = None,
) -> np.ndarray:
    """Simulate bulk RNA-seq count data.

    Args:
        n_samples: Number of samples to simulate
        n_genes: Number of genes
        library_sizes: Library sizes for each sample (if None, randomly generated)
        gene_means: Mean expression for each gene (if None, randomly generated)
        rng: Random number generator for reproducibility

    Returns:
        Count matrix (samples x genes)

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_range(n_samples, min_val=1, name="n_samples")
    validation.validate_range(n_genes, min_val=1, name="n_genes")

    if rng is None:
        rng = random.Random()

    np.random.seed(rng.randint(0, 2**32))

    # Generate library sizes if not provided
    if library_sizes is None:
        # Library sizes follow a log-normal distribution
        library_sizes = np.random.lognormal(15, 0.5, n_samples)  # Mean ~1M reads

    if len(library_sizes) != n_samples:
        raise errors.ValidationError(f"library_sizes must have length {n_samples}")

    # Generate gene means if not provided
    if gene_means is None:
        # Gene means follow a power-law distribution (many low, few high)
        gene_means = np.random.power(0.5, n_genes) * 1000

    if len(gene_means) != n_genes:
        raise errors.ValidationError(f"gene_means must have length {n_genes}")

    # Calculate proportions for each gene in each sample
    # This simulates different expression patterns across samples
    proportions = np.zeros((n_samples, n_genes))

    for i in range(n_samples):
        # Add sample-specific variation
        sample_factors = np.random.lognormal(0, 0.5, n_genes)
        proportions[i, :] = gene_means * sample_factors

        # Normalize to sum to 1
        proportions[i, :] /= proportions[i, :].sum()

    # Generate multinomial counts
    counts = np.zeros((n_samples, n_genes), dtype=int)

    for i in range(n_samples):
        lib_size = int(library_sizes[i])
        probs = proportions[i, :]
        counts[i, :] = np.random.multinomial(lib_size, probs)

    return counts


def simulate_single_cell_rnaseq(
    n_cells: int, n_genes: int, *, n_cell_types: int = 5, dropout_rate: float = 0.3, rng: random.Random | None = None
) -> Tuple[np.ndarray, np.ndarray]:
    """Simulate single-cell RNA-seq data with cell types and dropout.

    Args:
        n_cells: Number of cells to simulate
        n_genes: Number of genes
        n_cell_types: Number of distinct cell types
        dropout_rate: Probability of dropout for each gene expression
        rng: Random number generator for reproducibility

    Returns:
        Tuple of (expression_matrix, cell_type_labels)

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_range(n_cells, min_val=1, name="n_cells")
    validation.validate_range(n_genes, min_val=1, name="n_genes")
    validation.validate_range(n_cell_types, min_val=1, name="n_cell_types")
    validation.validate_range(dropout_rate, min_val=0.0, max_val=1.0, name="dropout_rate")

    if rng is None:
        rng = random.Random()

    np.random.seed(rng.randint(0, 2**32))

    # Assign cells to cell types
    cells_per_type = n_cells // n_cell_types
    cell_type_labels = np.repeat(range(n_cell_types), cells_per_type)
    remaining = n_cells - len(cell_type_labels)
    if remaining > 0:
        cell_type_labels = np.concatenate([cell_type_labels, np.arange(remaining)])

    # Generate cell type-specific expression patterns
    cell_type_means = np.zeros((n_cell_types, n_genes))

    for ct in range(n_cell_types):
        # Each cell type has different expression patterns
        # Some genes are highly expressed in this cell type
        marker_genes = np.random.choice(n_genes, size=n_genes // 10, replace=False)
        cell_type_means[ct, marker_genes] = np.random.exponential(50, len(marker_genes))

        # Background expression for all genes
        background = np.random.exponential(5, n_genes)
        cell_type_means[ct, :] += background

    # Generate expression for each cell
    expression_matrix = np.zeros((n_cells, n_genes))

    for i in range(n_cells):
        cell_type = cell_type_labels[i]
        means = cell_type_means[cell_type, :]

        # Generate negative binomial counts
        dispersions = 1.0 / (means + 1.0)
        dispersions = np.clip(dispersions, 0.1, 10.0)

        # Simulate one cell (n_samples=1)
        counts = simulate_counts_negative_binomial(1, n_genes, means, dispersions, rng=rng)
        expression_matrix[i, :] = counts[0, :]

        # Apply dropout
        dropout_mask = np.random.random(n_genes) < dropout_rate
        expression_matrix[i, dropout_mask] = 0

    return expression_matrix.astype(int), cell_type_labels


def simulate_time_series_expression(
    n_timepoints: int, n_genes: int, *, oscillation_freq: Optional[np.ndarray] = None, rng: random.Random | None = None
) -> np.ndarray:
    """Simulate time-series gene expression with oscillatory patterns.

    Args:
        n_timepoints: Number of time points
        n_genes: Number of genes
        oscillation_freq: Oscillation frequencies for each gene (if None, randomly generated)
        rng: Random number generator for reproducibility

    Returns:
        Expression matrix (timepoints x genes)

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_range(n_timepoints, min_val=2, name="n_timepoints")
    validation.validate_range(n_genes, min_val=1, name="n_genes")

    if rng is None:
        rng = random.Random()

    np.random.seed(rng.randint(0, 2**32))

    time_points = np.linspace(0, 4 * np.pi, n_timepoints)  # Two full cycles

    # Generate oscillation frequencies
    if oscillation_freq is None:
        oscillation_freq = np.random.uniform(0.5, 2.0, n_genes)  # Random frequencies
    elif len(oscillation_freq) != n_genes:
        raise errors.ValidationError(f"oscillation_freq must have length {n_genes}")

    # Generate oscillatory expression
    expression = np.zeros((n_timepoints, n_genes))

    for t in range(n_timepoints):
        for g in range(n_genes):
            # Sine wave with random phase and amplitude
            phase = rng.uniform(0, 2 * np.pi)
            amplitude = rng.uniform(10, 100)
            baseline = rng.uniform(5, 20)

            expression[t, g] = baseline + amplitude * np.sin(oscillation_freq[g] * time_points[t] + phase)

            # Ensure non-negative
            expression[t, g] = max(0, expression[t, g])

    # Convert to counts (round and add noise)
    counts = np.random.poisson(expression)

    return counts.astype(int)


def simulate_spatial_expression(
    n_spots: int, n_genes: int, *, spatial_patterns: str = "random", rng: random.Random | None = None
) -> Tuple[np.ndarray, np.ndarray]:
    """Simulate spatially resolved RNA expression data.

    Args:
        n_spots: Number of spatial spots/locations
        n_genes: Number of genes
        spatial_patterns: Type of spatial pattern ("random", "gradient", "clusters")
        rng: Random number generator for reproducibility

    Returns:
        Tuple of (expression_matrix, coordinates) where coordinates are (x, y) positions

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_range(n_spots, min_val=1, name="n_spots")
    validation.validate_range(n_genes, min_val=1, name="n_genes")

    if spatial_patterns not in ["random", "gradient", "clusters"]:
        raise errors.ValidationError(f"Invalid spatial_patterns: {spatial_patterns}")

    if rng is None:
        rng = random.Random()

    np.random.seed(rng.randint(0, 2**32))

    # Generate spatial coordinates
    if spatial_patterns == "random":
        coordinates = np.random.uniform(0, 100, (n_spots, 2))
    elif spatial_patterns == "gradient":
        # Create a gradient across space
        x_coords = np.linspace(0, 100, int(np.sqrt(n_spots)))
        y_coords = np.linspace(0, 100, int(np.sqrt(n_spots)))
        xx, yy = np.meshgrid(x_coords, y_coords)
        coordinates = np.column_stack([xx.ravel()[:n_spots], yy.ravel()[:n_spots]])
    else:  # clusters
        # Create clustered spatial patterns
        n_clusters = 3
        coordinates = np.zeros((n_spots, 2))
        spots_per_cluster = n_spots // n_clusters

        for i in range(n_clusters):
            center_x, center_y = rng.uniform(20, 80), rng.uniform(20, 80)
            start_idx = i * spots_per_cluster
            end_idx = start_idx + spots_per_cluster if i < n_clusters - 1 else n_spots

            cluster_coords = np.random.normal([center_x, center_y], 10, (end_idx - start_idx, 2))
            coordinates[start_idx:end_idx] = cluster_coords

    # Generate expression based on spatial patterns
    expression_matrix = np.zeros((n_spots, n_genes))

    for g in range(n_genes):
        if spatial_patterns == "random":
            # Random spatial expression
            means = np.random.exponential(20, n_spots)
        elif spatial_patterns == "gradient":
            # Expression increases with x-coordinate
            means = 10 + coordinates[:, 0] * 0.5 + np.random.normal(0, 5, n_spots)
            means = np.maximum(means, 0)
        else:  # clusters
            # Cluster-specific expression
            means = np.zeros(n_spots)
            for i in range(n_clusters):
                cluster_mask = np.arange(n_spots) // spots_per_cluster == i
                if rng.random() < 0.7:  # 70% chance this gene is expressed in this cluster
                    means[cluster_mask] = rng.uniform(50, 200)
                else:
                    means[cluster_mask] = rng.uniform(1, 10)

        # Add spatial autocorrelation
        for i in range(n_spots):
            # Nearby spots have similar expression (simplified)
            nearby_indices = np.argsort(np.sum((coordinates - coordinates[i]) ** 2, axis=1))[:5]
            nearby_mean = np.mean(means[nearby_indices])
            means[i] = 0.7 * means[i] + 0.3 * nearby_mean

        # Generate counts
        dispersions = 1.0 / (means + 1.0)
        dispersions = np.clip(dispersions, 0.1, 10.0)

        counts = simulate_counts_negative_binomial(1, n_spots, means, dispersions, rng=rng)
        expression_matrix[:, g] = counts[0, :]

    return expression_matrix.astype(int), coordinates


def add_technical_noise(
    expression_matrix: np.ndarray,
    *,
    amplification_bias: float = 0.1,
    sequencing_depth: Optional[float] = None,
    rng: random.Random | None = None,
) -> np.ndarray:
    """Add technical noise to expression matrix.

    Args:
        expression_matrix: Input expression matrix
        amplification_bias: Strength of amplification bias (0-1)
        sequencing_depth: Target sequencing depth (if None, use current)
        rng: Random number generator for reproducibility

    Returns:
        Expression matrix with technical noise

    Raises:
        ValueError: If parameters are invalid
    """
    validation.validate_type(expression_matrix, np.ndarray, "expression_matrix")
    validation.validate_range(amplification_bias, min_val=0.0, max_val=1.0, name="amplification_bias")

    if rng is None:
        rng = random.Random()

    np.random.seed(rng.randint(0, 2**32))

    noisy_matrix = expression_matrix.copy().astype(float)

    # Add amplification bias (PCR amplification preferentially amplifies some transcripts)
    if amplification_bias > 0:
        for i in range(expression_matrix.shape[0]):
            bias_factors = np.random.lognormal(0, amplification_bias, expression_matrix.shape[1])
            noisy_matrix[i, :] *= bias_factors

    # Adjust sequencing depth
    if sequencing_depth is not None:
        current_depths = np.sum(noisy_matrix, axis=1)
        depth_factors = sequencing_depth / current_depths
        noisy_matrix *= depth_factors[:, np.newaxis]

    # Add Poisson noise (shot noise)
    noisy_matrix = np.random.poisson(noisy_matrix)

    return noisy_matrix.astype(int)
