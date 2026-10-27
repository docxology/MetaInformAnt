# Single-Cell Data Integration

The integration module provides tools for combining multiple single-cell datasets and correcting for batch effects.

## Dataset Integration

### integrate_multiple_batches()

Combine multiple datasets along the cell axis and correct for batch effects in one call. The expression matrices are stacked directly, so the datasets must contain the same genes in the same order (there are no gene intersection or union helpers; align the genes first). Cells are labelled `batch_0`, `batch_1`, ... in `obs['batch']` according to the position of their dataset in the list:

```python-snippet
from metainformant.singlecell.data.integration import integrate_multiple_batches

# Combine datasets from different samples/experiments
datasets = [data1, data2, data3]

integrated_data = integrate_multiple_batches(
    datasets,
    integration_method='scanorama'   # 'scanorama', 'mnn', 'harmony' or 'bbknn'
)

print(f"Integrated: {integrated_data.n_obs} cells × {integrated_data.n_vars} features")
```

The integrated object stores principal-component coordinates in `X` (at most 50 columns), so `n_vars` counts components rather than genes.

## Batch Correction

### remove_batch_effects()

Simple per-batch correction of the expression matrix:

```python-snippet
from metainformant.singlecell.data.preprocessing import remove_batch_effects

# Per-batch centering or standardization
corrected_data = remove_batch_effects(
    data,
    batch_key='batch',      # Batch column in obs
    method='regress_out'    # 'regress_out' (subtract batch means) or 'combat' (standardize within batch)
)
```

### combat_integration()

ComBat-style batch correction (a simplified per-gene adjustment of each batch to a reference batch; no additional dependencies are required):

```python-snippet
from metainformant.singlecell.data.integration import combat_integration

corrected_data = combat_integration(
    data,
    batch_key='batch',
    covariates=None        # Optional covariate names (recorded in data.uns['combat'])
)
```

### harmony_integration()

Harmony-style integration (a simplified implementation that centers each batch in PCA space; requires scikit-learn):

```python-snippet
from metainformant.singlecell.data.integration import harmony_integration

corrected_data = harmony_integration(
    data,
    batch_key='batch',
    n_components=50        # Number of components to compute
)
# Corrected coordinates are stored in corrected_data.obs as harmony_1 ... harmony_10
```

## Integration Workflows

### Standard Integration Pipeline

```python-snippet
from metainformant.singlecell.data.preprocessing import (
    calculate_qc_metrics,
    filter_cells,
    log_transform,
    normalize_counts,
)
from metainformant.singlecell.data.integration import integrate_multiple_batches
from metainformant.singlecell.analysis.nonlinear_methods import compute_neighbors, compute_umap

def integrate_datasets(datasets, batch_labels, method='scanorama'):
    """
    Standard integration workflow.
    
    Args:
        datasets: List of SingleCellData objects (same genes in the same order)
        batch_labels: List of batch identifiers  
        method: Integration method ('scanorama', 'mnn', 'harmony' or 'bbknn')
    """
    
    # 1. Standard preprocessing of each dataset
    #    (gene filtering must keep the same genes in every dataset)
    print("Preprocessing datasets...")
    processed = []
    for data in datasets:
        data = calculate_qc_metrics(data)
        data = filter_cells(data, min_genes=200, max_pct_mt=20)
        data = normalize_counts(data, method='total_count')
        data = log_transform(data)
        processed.append(data)
    
    # 2. Concatenate datasets and correct batch effects
    print(f"Applying {method} integration...")
    integrated = integrate_multiple_batches(processed, integration_method=method)
    
    # Replace the automatic batch_0, batch_1, ... labels
    label_map = {f"batch_{i}": label for i, label in enumerate(batch_labels)}
    integrated.obs['batch'] = integrated.obs['batch'].map(label_map)
    
    # 3. Neighbors and embedding (X now holds principal components)
    integrated = compute_neighbors(integrated)
    integrated = compute_umap(integrated)
    
    print("Integration complete!")
    return integrated

# Example usage
datasets = [sample1_data, sample2_data, sample3_data]
batch_labels = ['Ctrl', 'Treated_1', 'Treated_2']

integrated_data = integrate_datasets(datasets, batch_labels, method='scanorama')
```

## Quality Assessment

### evaluate_integration_quality()

Assess integration quality:

```python-snippet
from metainformant.singlecell.data.integration import evaluate_integration_quality

# Calculate integration metrics
metrics = evaluate_integration_quality(
    integrated_data,
    batch_key='batch',
    cluster_key='cluster'    # Optional cluster labels column
)

print("Integration Quality Metrics:")
print(f"  Mean batch mixing score: {metrics['mean_batch_mixing_score']:.3f}")  # share of neighbors from other batches
print(f"  Between/within batch distance ratio: {metrics['batch_distance_ratio']:.3f}")
print(f"  Clusters independent of batch: {metrics['cluster_batch_independence']}")  # needs cluster_key
```

## Advanced Integration Methods

### Multi-Modal Integration

```python
import numpy as np
from metainformant.singlecell.data.preprocessing import log_transform, normalize_counts
from metainformant.singlecell.analysis.pca_methods import compute_pca, select_hvgs

def integrate_multimodal_data(rna_data, adt_data, batch_key='batch'):
    """
    Integrate RNA and protein (ADT) data.
    
    Args:
        rna_data: Gene expression data
        adt_data: Antibody-derived tag data
        batch_key: Batch column
    """
    
    # 1. Separate processing for each modality
    # RNA processing
    rna_processed = normalize_counts(rna_data, method='total_count')
    rna_processed = log_transform(rna_processed)
    rna_processed = select_hvgs(rna_processed, n_top_genes=2000)
    
    # ADT processing (different normalization)
    adt_processed = normalize_counts(adt_data, method='centered_log_ratio')
    
    # 2. Dimensionality reduction per modality
    rna_processed = compute_pca(rna_processed, n_components=50)
    adt_processed = compute_pca(adt_processed, n_components=20)
    
    # 3. Joint embedding (simplified - would use more sophisticated methods)
    # Concatenate PCA coordinates
    joint_embedding = np.hstack([
        rna_processed.obsm['X_pca'],
        adt_processed.obsm['X_pca']
    ])
    
    # 4. Create integrated object
    integrated = rna_processed.copy()
    integrated.obsm['X_joint'] = joint_embedding
    
    return integrated
```

### Cross-Species Integration

```python
from metainformant.singlecell.data.integration import integrate_multiple_batches

def integrate_cross_species(human_data, mouse_data):
    """
    Integrate data across species using orthologous genes.
    
    Args:
        human_data: Human single-cell data
        mouse_data: Mouse single-cell data
    """
    
    # 1. Load orthology mapping (simplified)
    ortholog_mapping = {
        'human_gene': 'mouse_gene',
        # ... more mappings
    }
    
    # 2. Map mouse genes to human
    mouse_genes_mapped = [ortholog_mapping.get(g, g) for g in mouse_data.var.index]
    mouse_data.var.index = mouse_genes_mapped
    
    # 3. Find common genes
    common_genes = list(set(human_data.var.index) & set(mouse_data.var.index))
    
    # 4. Subset to common genes
    human_subset = human_data[:, common_genes].copy()
    mouse_subset = mouse_data[:, common_genes].copy()
    
    # 5. Add species labels
    human_subset.obs['species'] = 'Human'
    mouse_subset.obs['species'] = 'Mouse'
    
    # 6. Integrate (obs['species'] is kept; obs['batch'] holds batch_0 = human, batch_1 = mouse)
    integrated = integrate_multiple_batches(
        [human_subset, mouse_subset],
        integration_method='scanorama'
    )
    
    return integrated
```

## Batch Effect Visualization

### Before/After Comparison

```python
import matplotlib.pyplot as plt

def plot_batch_effect_comparison(data_before, data_after, batch_key='batch'):
    """Plot batch effects before and after correction."""
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Before correction - UMAP
    ax = axes[0, 0]
    for batch in data_before.obs[batch_key].unique():
        mask = data_before.obs[batch_key] == batch
        ax.scatter(
            data_before.obsm['X_umap'][mask, 0],
            data_before.obsm['X_umap'][mask, 1], 
            label=batch,
            alpha=0.6
        )
    ax.set_title('Before Correction (UMAP)')
    ax.legend()
    
    # After correction - UMAP
    ax = axes[0, 1]
    for batch in data_after.obs[batch_key].unique():
        mask = data_after.obs[batch_key] == batch
        ax.scatter(
            data_after.obsm['X_umap'][mask, 0],
            data_after.obsm['X_umap'][mask, 1],
            label=batch,
            alpha=0.6
        )
    ax.set_title('After Correction (UMAP)')
    ax.legend()
    
    # Before correction - PCA
    ax = axes[1, 0]
    for batch in data_before.obs[batch_key].unique():
        mask = data_before.obs[batch_key] == batch
        ax.scatter(
            data_before.obsm['X_pca'][mask, 0],
            data_before.obsm['X_pca'][mask, 1],
            label=batch,
            alpha=0.6
        )
    ax.set_title('Before Correction (PCA)')
    ax.legend()
    
    # After correction - PCA
    ax = axes[1, 1] 
    for batch in data_after.obs[batch_key].unique():
        mask = data_after.obs[batch_key] == batch
        ax.scatter(
            data_after.obsm['X_pca'][mask, 0],
            data_after.obsm['X_pca'][mask, 1],
            label=batch,
            alpha=0.6
        )
    ax.set_title('After Correction (PCA)')
    ax.legend()
    
    plt.tight_layout()
    return fig

# Usage
fig = plot_batch_effect_comparison(uncorrected_data, corrected_data)
plt.savefig('output/batch_correction_comparison.pdf', dpi=300)
plt.show()
```

## Method Selection Guidelines

### Choosing Integration Methods

```python
def recommend_integration_method(datasets, batch_info):
    """
    Recommend integration method based on data characteristics.
    
    Args:
        datasets: List of datasets
        batch_info: Dictionary with batch information
    """
    
    n_batches = len(datasets)
    total_cells = sum(d.n_obs for d in datasets)
    
    recommendations = []
    
    if n_batches <= 3 and total_cells < 50000:
        recommendations.append("scaling - Simple and fast")
    
    if n_batches > 3 or total_cells > 50000:
        recommendations.append("harmony - Good for large datasets") 
    
    if batch_info.get('technical_replicates', False):
        recommendations.append("combat - Good for technical batches")
    
    if batch_info.get('different_protocols', False):
        recommendations.append("harmony or combat - Handle protocol differences")
    
    return recommendations

# Example usage
batch_info = {
    'technical_replicates': False,
    'different_protocols': True,
    'batch_confounded': False
}

methods = recommend_integration_method(datasets, batch_info)
print("Recommended methods:")
for method in methods:
    print(f"  - {method}")
```

## Common Integration Issues

### Overcorrection Detection

```python
def detect_overcorrection(data_before, data_after, cell_type_key='cell_type'):
    """Detect if batch correction was too aggressive."""
    
    if cell_type_key not in data_before.obs.columns:
        print("Warning: No cell type information for overcorrection assessment")
        return None
    
    # Calculate within-cell-type batch mixing before and after
    from sklearn.metrics import silhouette_score
    
    # Before correction
    sil_before = silhouette_score(
        data_before.obsm['X_pca'],
        data_before.obs['batch']
    )
    
    # After correction  
    sil_after = silhouette_score(
        data_after.obsm['X_pca'],
        data_after.obs['batch']
    )
    
    # Cell type preservation
    ct_sil_before = silhouette_score(
        data_before.obsm['X_pca'],
        data_before.obs[cell_type_key]
    )
    
    ct_sil_after = silhouette_score(
        data_after.obsm['X_pca'], 
        data_after.obs[cell_type_key]
    )
    
    print(f"Batch separation - Before: {sil_before:.3f}, After: {sil_after:.3f}")
    print(f"Cell type separation - Before: {ct_sil_before:.3f}, After: {ct_sil_after:.3f}")
    
    if sil_after < 0.05 and (ct_sil_before - ct_sil_after) > 0.1:
        print("Warning: Possible overcorrection detected")
        return True
    
    return False
```

## Best Practices

### Data Preparation

1. **Quality control first**: Filter low-quality cells before integration
2. **Consistent preprocessing**: Use same normalization across datasets
3. **Feature selection**: Select HVGs on combined dataset
4. **Documentation**: Track all processing steps and batch information

### Method Selection

1. **Start simple**: Try scaling-based correction first
2. **Validate results**: Always assess integration quality
3. **Preserve biology**: Ensure cell types remain separable
4. **Multiple methods**: Compare different approaches

### Validation

1. **Visual inspection**: Plot embeddings colored by batch and cell type
2. **Quantitative metrics**: Calculate mixing entropy and silhouette scores
3. **Biological validation**: Check marker gene expression
4. **Downstream analysis**: Ensure clustering and trajectory results make sense

## Performance Tips

### Large Dataset Integration

```python
# For very large datasets (>500k cells)
# 1. Subsample for method optimization
sample_size = 50000
indices = np.random.choice(total_cells, sample_size, replace=False)

# 2. Optimize on subset
subset_integrated = integrate_datasets(
    [d[indices[i]:indices[i+1]] for i, d in enumerate(datasets)],
    batch_labels,
    method='harmony'
)

# 3. Apply to full dataset with optimized parameters
full_integrated = integrate_datasets(datasets, batch_labels, method='harmony')
```

## Testing

Integration functionality is tested as part of the single-cell test suite:

```bash
# Run single-cell tests (includes integration examples)
uv run pytest tests/singlecell/test_singlecell_*.py -v
```

## Related Documentation

- [Preprocessing](./preprocessing.md): For data preparation before integration
- [Dimensionality Reduction](./dimensionality.md): For post-integration analysis
- [Clustering](./clustering.md): For clustering integrated data
- [Visualization](./visualization.md): For plotting integration results
