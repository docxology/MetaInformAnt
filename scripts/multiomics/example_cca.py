#!/usr/bin/env python3
"""Canonical Correlation Analysis example.

This script demonstrates CCA between two omics layers to find correlated patterns.
All outputs are written to output/multiomics/.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from metainformant.core import io as core_io
from metainformant.core import paths as core_paths
from metainformant.multiomics import MultiOmicsData, canonical_correlation

# Set output directory
output_dir = Path("output/multiomics")
core_paths.ensure_directory(output_dir)

print("=" * 60)
print("Canonical Correlation Analysis Example")
print("=" * 60)

# Create synthetic data with known correlation structure
np.random.seed(42)
n_samples = 25
n_genes = 100
n_proteins = 50

samples = [f"sample_{i:03d}" for i in range(n_samples)]

print("\n1. Creating synthetic data with correlated structure...")

# Create shared latent factors
shared_factor = np.random.randn(n_samples, 1)

# Transcriptomics: strongly influenced by shared factor
transcriptomics_data = pd.DataFrame(
    shared_factor @ np.random.randn(1, n_genes) + np.random.randn(n_samples, n_genes) * 0.3,
    index=samples,
    columns=[f"GENE_{i:04d}" for i in range(n_genes)],
)

# Proteomics: also influenced by shared factor (correlated with transcriptomics)
proteomics_data = pd.DataFrame(
    shared_factor @ np.random.randn(1, n_proteins) + np.random.randn(n_samples, n_proteins) * 0.3,
    index=samples,
    columns=[f"PROTEIN_{i:04d}" for i in range(n_proteins)],
)

print(f"   Transcriptomics shape: {transcriptomics_data.shape}")
print(f"   Proteomics shape: {proteomics_data.shape}")

# Create MultiOmicsData
omics_data = MultiOmicsData(
    transcriptomics=transcriptomics_data,
    proteomics=proteomics_data,
)
print(f"   Integrated {omics_data.n_samples} samples")

# Perform CCA
print("\n2. Performing Canonical Correlation Analysis...")
X_c, Y_c, X_w, Y_w, correlations = canonical_correlation(
    omics_data,
    layer_pair=("transcriptomics", "proteomics"),
    n_components=10,
    regularization=0.01,
)

print(f"   Number of components: {len(correlations)}")
print(f"   Canonical correlations: {correlations[:5]}")
print(f"   First canonical correlation: {correlations[0]:.3f}")
print(f"   Transcriptomics weights shape: {X_w.shape}")
print(f"   Proteomics weights shape: {Y_w.shape}")
print(f"   Canonical variates shape (X): {X_c.shape}")
print(f"   Canonical variates shape (Y): {Y_c.shape}")

# Verify correlation
actual_correlation = np.corrcoef(X_c[:, 0], Y_c[:, 0])[0, 1]
print(f"\n   Verified correlation between first variates: {actual_correlation:.6f}")
print(f"   Expected correlation (from CCA): {correlations[0]:.6f}")

# Analyze which features contribute most
print("\n3. Analyzing feature contributions...")
# Top features for transcriptomics (first component)
top_genes_idx = np.argsort(np.abs(X_w[:, 0]))[::-1][:10]
top_genes = [omics_data.get_layer("transcriptomics").columns[i] for i in top_genes_idx]
print(f"   Top 10 genes contributing to first canonical variate:")
for i, gene in enumerate(top_genes, 1):
    weight = X_w[top_genes_idx[i - 1], 0]
    print(f"     {i}. {gene}: {weight:.4f}")

# Top features for proteomics (first component)
top_proteins_idx = np.argsort(np.abs(Y_w[:, 0]))[::-1][:10]
top_proteins = [omics_data.get_layer("proteomics").columns[i] for i in top_proteins_idx]
print(f"\n   Top 10 proteins contributing to first canonical variate:")
for i, protein in enumerate(top_proteins, 1):
    weight = Y_w[top_proteins_idx[i - 1], 0]
    print(f"     {i}. {protein}: {weight:.4f}")

# Save results
print("\n4. Saving CCA results...")
cca_results = {
    "canonical_correlations": correlations.tolist(),
    "n_components": len(correlations),
    "transcriptomics_weights_shape": list(X_w.shape),
    "proteomics_weights_shape": list(Y_w.shape),
    "canonical_variates_shape": {
        "transcriptomics": list(X_c.shape),
        "proteomics": list(Y_c.shape),
    },
    "top_features": {
        "transcriptomics": {
            "genes": top_genes[:5],
            "weights": X_w[top_genes_idx[:5], 0].tolist(),
        },
        "proteomics": {
            "proteins": top_proteins[:5],
            "weights": Y_w[top_proteins_idx[:5], 0].tolist(),
        },
    },
}
core_io.dump_json(cca_results, output_dir / "cca_results.json", indent=2)
print(f"   Saved CCA results to {output_dir / 'cca_results.json'}")

# Test with different regularization
print("\n5. Testing different regularization values...")
for reg in [0.001, 0.01, 0.1]:
    _, _, _, _, corr = canonical_correlation(
        omics_data,
        layer_pair=("transcriptomics", "proteomics"),
        n_components=5,
        regularization=reg,
    )
    print(f"   Regularization {reg}: first correlation = {corr[0]:.3f}")

print("\n" + "=" * 60)
print("CCA example completed!")
print(f"All outputs saved to: {output_dir}")
print("=" * 60)


