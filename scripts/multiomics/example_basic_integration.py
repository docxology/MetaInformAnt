#!/usr/bin/env python3
"""Basic multi-omics integration example.

This script demonstrates how to create and work with MultiOmicsData objects
using synthetic data. All outputs are written to output/multiomics/.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from metainformant.core import io as core_io
from metainformant.core import paths as core_paths
from metainformant.multiomics import MultiOmicsData, integrate_omics_data, joint_pca

# Set output directory
output_dir = Path("output/multiomics")
core_paths.ensure_directory(output_dir)

print("=" * 60)
print("Multi-Omics Basic Integration Example")
print("=" * 60)

# Create synthetic data
np.random.seed(42)
n_samples = 20
n_genes = 100
n_variants = 50

# Generate sample IDs
samples = [f"sample_{i:03d}" for i in range(n_samples)]

# Create genomics data (variant genotypes: 0, 1, 2)
print("\n1. Creating synthetic genomics data...")
genomics_data = pd.DataFrame(
    np.random.randint(0, 3, size=(n_samples, n_variants)),
    index=samples,
    columns=[f"SNP_{i:04d}" for i in range(n_variants)],
)
print(f"   Genomics shape: {genomics_data.shape}")

# Create transcriptomics data (gene expression)
print("\n2. Creating synthetic transcriptomics data...")
transcriptomics_data = pd.DataFrame(
    np.random.lognormal(mean=2, sigma=1, size=(n_samples, n_genes)),
    index=samples,
    columns=[f"GENE_{i:04d}" for i in range(n_genes)],
)
print(f"   Transcriptomics shape: {transcriptomics_data.shape}")

# Create proteomics data (protein abundance)
print("\n3. Creating synthetic proteomics data...")
proteomics_data = pd.DataFrame(
    np.random.lognormal(mean=1, sigma=0.5, size=(n_samples, 30)),
    index=samples,
    columns=[f"PROTEIN_{i:04d}" for i in range(30)],
)
print(f"   Proteomics shape: {proteomics_data.shape}")

# Save synthetic data to files
print("\n4. Saving synthetic data to files...")
genomics_path = output_dir / "genomics.csv"
transcriptomics_path = output_dir / "transcriptomics.csv"
proteomics_path = output_dir / "proteomics.csv"

core_io.write_csv(genomics_data, genomics_path)
core_io.write_csv(transcriptomics_data, transcriptomics_path)
core_io.write_csv(proteomics_data, proteomics_path)
print(f"   Saved data files to {output_dir}")

# Method 1: Create MultiOmicsData from DataFrames
print("\n5. Creating MultiOmicsData from DataFrames...")
omics_data = MultiOmicsData(
    genomics=genomics_data,
    transcriptomics=transcriptomics_data,
    proteomics=proteomics_data,
)
print(f"   Integrated {omics_data.n_samples} samples")
print(f"   Layers: {omics_data.layer_names}")

# Method 2: Load from files
print("\n6. Loading MultiOmicsData from files...")
data_dict = {
    "genomics": str(genomics_path),
    "transcriptomics": str(transcriptomics_path),
    "proteomics": str(proteomics_path),
}
omics_data_from_files = integrate_omics_data(data_dict)
print(f"   Loaded {omics_data_from_files.n_samples} samples")

# Save and load MultiOmicsData
print("\n7. Testing save/load functionality...")
save_path = output_dir / "integrated_data"
omics_data.save(save_path)
print(f"   Saved to {save_path}")

loaded_data = MultiOmicsData.load(save_path)
print(f"   Loaded from {save_path}: {loaded_data.n_samples} samples")

# Perform joint PCA
print("\n8. Performing joint PCA...")
embeddings, loadings, variance = joint_pca(omics_data, n_components=10)
print(f"   Embeddings shape: {embeddings.shape}")
print(f"   Explained variance (first 5 components): {variance[:5]}")
print(f"   Loadings keys: {list(loadings.keys())}")

# Save PCA results
print("\n9. Saving PCA results...")
pca_results = {
    "embeddings": embeddings.tolist(),
    "explained_variance": variance.tolist(),
    "loadings_shapes": {k: list(v.shape) for k, v in loadings.items()},
}
core_io.dump_json(pca_results, output_dir / "pca_results.json", indent=2)
print(f"   Saved PCA results to {output_dir / 'pca_results.json'}")

# Subset samples
print("\n10. Testing sample subsetting...")
subset_samples = ["sample_001", "sample_002", "sample_003"]
subset_data = omics_data.subset_samples(subset_samples)
print(f"   Original samples: {omics_data.n_samples}")
print(f"   Subset samples: {subset_data.n_samples}")

# Subset features
print("\n11. Testing feature subsetting...")
subset_features = {
    "transcriptomics": ["GENE_0000", "GENE_0001", "GENE_0002"],
    "genomics": ["SNP_0000", "SNP_0001"],
}
subset_features_data = omics_data.subset_features(subset_features)
print(f"   Transcriptomics features: {subset_features_data.get_layer('transcriptomics').shape[1]}")
print(f"   Genomics features: {subset_features_data.get_layer('genomics').shape[1]}")

print("\n" + "=" * 60)
print("Example completed successfully!")
print(f"All outputs saved to: {output_dir}")
print("=" * 60)

