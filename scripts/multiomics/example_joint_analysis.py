#!/usr/bin/env python3
"""Joint analysis methods example.

This script demonstrates joint PCA and joint NMF analysis on multi-omics data.
All outputs are written to output/multiomics/.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from metainformant.core import io as core_io
from metainformant.core import paths as core_paths
from metainformant.multiomics import MultiOmicsData, joint_nmf, joint_pca

# Set output directory
output_dir = Path("output/multiomics")
core_paths.ensure_directory(output_dir)

print("=" * 60)
print("Multi-Omics Joint Analysis Example")
print("=" * 60)

# Create synthetic data with some structure
np.random.seed(42)
n_samples = 30
n_genes = 150
n_proteins = 40

samples = [f"sample_{i:03d}" for i in range(n_samples)]

# Create correlated transcriptomics and proteomics data
print("\n1. Creating synthetic multi-omics data with structure...")

# Create latent factors that will be shared
latent_factors = np.random.randn(n_samples, 5)

# Transcriptomics: influenced by latent factors
transcriptomics_data = pd.DataFrame(
    latent_factors @ np.random.randn(5, n_genes) + np.random.randn(n_samples, n_genes) * 0.5,
    index=samples,
    columns=[f"GENE_{i:04d}" for i in range(n_genes)],
)

# Proteomics: also influenced by same latent factors (correlated)
proteomics_data = pd.DataFrame(
    latent_factors @ np.random.randn(5, n_proteins) + np.random.randn(n_samples, n_proteins) * 0.5,
    index=samples,
    columns=[f"PROTEIN_{i:04d}" for i in range(n_proteins)],
)

# Make data non-negative for NMF
transcriptomics_data = transcriptomics_data - transcriptomics_data.min().min() + 1
proteomics_data = proteomics_data - proteomics_data.min().min() + 1

print(f"   Transcriptomics shape: {transcriptomics_data.shape}")
print(f"   Proteomics shape: {proteomics_data.shape}")

# Create MultiOmicsData
omics_data = MultiOmicsData(
    transcriptomics=transcriptomics_data,
    proteomics=proteomics_data,
)
print(f"   Integrated {omics_data.n_samples} samples")

# Joint PCA
print("\n2. Performing Joint PCA...")
print("   - With equal layer weights")
embeddings_equal, loadings_equal, variance_equal = joint_pca(
    omics_data,
    n_components=10,
    layer_weights=None,  # Equal weights
    standardize=True,
)
print(f"   Embeddings shape: {embeddings_equal.shape}")
print(f"   Explained variance (first 5): {variance_equal[:5]}")

print("\n   - With different layer weights")
embeddings_weighted, loadings_weighted, variance_weighted = joint_pca(
    omics_data,
    n_components=10,
    layer_weights={"transcriptomics": 1.0, "proteomics": 2.0},  # Weight proteomics more
    standardize=True,
)
print(f"   Embeddings shape: {embeddings_weighted.shape}")
print(f"   Explained variance (first 5): {variance_weighted[:5]}")

# Save PCA results
pca_results = {
    "equal_weights": {
        "embeddings_shape": list(embeddings_equal.shape),
        "explained_variance": variance_equal.tolist(),
    },
    "weighted": {
        "embeddings_shape": list(embeddings_weighted.shape),
        "explained_variance": variance_weighted.tolist(),
    },
}
core_io.dump_json(pca_results, output_dir / "joint_pca_results.json", indent=2)
print(f"   Saved PCA results to {output_dir / 'joint_pca_results.json'}")

# Joint NMF
print("\n3. Performing Joint NMF...")
W, H = joint_nmf(
    omics_data,
    n_components=5,
    max_iter=200,
    regularization=0.01,
    random_state=42,
    tolerance=1e-6,
)
print(f"   Sample factors (W) shape: {W.shape}")
print(f"   Transcriptomics feature factors shape: {H['transcriptomics'].shape}")
print(f"   Proteomics feature factors shape: {H['proteomics'].shape}")

# Compute reconstruction error
print("\n4. Evaluating NMF reconstruction...")
reconstruction_error_transcript = np.sum(
    (omics_data.get_layer("transcriptomics").values - W @ H["transcriptomics"]) ** 2
)
reconstruction_error_protein = np.sum(
    (omics_data.get_layer("proteomics").values - W @ H["proteomics"]) ** 2
)
total_error = reconstruction_error_transcript + reconstruction_error_protein
print(f"   Total reconstruction error: {total_error:.6f}")
print(f"   Transcriptomics error: {reconstruction_error_transcript:.6f}")
print(f"   Proteomics error: {reconstruction_error_protein:.6f}")

# Save NMF results
nmf_results = {
    "sample_factors_shape": list(W.shape),
    "feature_factors_shapes": {k: list(v.shape) for k, v in H.items()},
    "reconstruction_errors": {
        "total": total_error,
        "transcriptomics": reconstruction_error_transcript,
        "proteomics": reconstruction_error_protein,
    },
}
core_io.dump_json(nmf_results, output_dir / "joint_nmf_results.json", indent=2)
print(f"   Saved NMF results to {output_dir / 'joint_nmf_results.json'}")

# Compare explained variance
print("\n5. Comparing analysis methods...")
print(f"   Joint PCA first component variance: {variance_equal[0]:.3f}")
print(f"   Joint PCA cumulative variance (first 5): {variance_equal[:5].sum():.3f}")
print(f"   NMF reconstruction error per sample: {total_error / omics_data.n_samples:.6f}")

print("\n" + "=" * 60)
print("Joint analysis example completed!")
print(f"All outputs saved to: {output_dir}")
print("=" * 60)


