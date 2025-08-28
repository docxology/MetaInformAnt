# Single-Cell Data Integration

The integration module provides tools for combining multiple single-cell datasets and correcting for batch effects.

## Dataset Integration

### concatenate_datasets()

Combine multiple datasets along the cell axis:

```python
from metainformant.singlecell.integration import concatenate_datasets

# Combine datasets from different samples/experiments
datasets = [data1, data2, data3]
batch_labels = ['Sample_1', 'Sample_2', 'Sample_3']

integrated_data = concatenate_datasets(
    datasets,
    batch_key='batch',      # Column name for batch labels
    batch_categories=batch_labels,
    join='outer'            # 'inner' for intersection, 'outer' for union
)

print(f"Integrated: {integrated_data.n_obs} cells × {integrated_data.n_vars} genes")
```

### intersect_datasets()

Keep only genes present in all datasets:

```python
from metainformant.singlecell.integration import intersect_datasets

# Find common genes across datasets
common_genes = intersect_datasets(datasets, return_indices=False)
print(f"Common genes: {len(common_genes)}")

# Apply intersection
integrated_data = concatenate_datasets(
    datasets, 
    batch_categories=batch_labels,
    join='inner'  # Equivalent to intersection
)
```

### union_datasets()

Keep all genes, filling missing values with zeros:

```python
from metainformant.singlecell.integration import union_datasets

# Union of all genes
all_genes = union_datasets(datasets, return_indices=False)
print(f"Total genes: {len(all_genes)}")

integrated_data = concatenate_datasets(
    datasets,
    batch_categories=batch_labels, 
    join='outer'  # Equivalent to union
)
```

## Batch Correction

### batch_correction_scaling()

Simple scaling-based batch correction:

```python
from metainformant.singlecell.integration import batch_correction_scaling

# Scale-based batch correction
corrected_data = batch_correction_scaling(
    data,
    batch_key='batch',      # Batch column in obs
    method='standardize'    # 'standardize', 'center', 'scale'
)
```

### batch_correction_combat()

ComBat batch correction (requires additional dependencies):

```python
from metainformant.singlecell.integration import batch_correction_combat

try:
    # ComBat batch correction
    corrected_data = batch_correction_combat(
        data,
        batch_key='batch',
        covariates=None,       # Optional covariates to preserve
        parametric=True        # Parametric vs non-parametric
    )
except ImportError:
    print("ComBat not available. Install with: pip install combat")
```

### batch_correction_harmony()

Harmony integration (requires harmony-pytorch):

```python
from metainformant.singlecell.integration import batch_correction_harmony

try:
    # Harmony batch correction
    corrected_data = batch_correction_harmony(
        data,
        batch_key='batch',
        n_components=50,       # Number of PCs to use
        theta=2.0,             # Diversity clustering penalty
        sigma=0.1              # Ridge regression penalty
    )
except ImportError:
    print("Harmony not available. Install with: pip install harmony-pytorch")
```

## Integration Workflows

### Standard Integration Pipeline

```python
from metainformant.singlecell.preprocessing import *
from metainformant.singlecell.dimensionality import *
from metainformant.singlecell.integration import *

def integrate_datasets(datasets, batch_labels, method='scaling'):
    """
    Standard integration workflow.
    
    Args:
        datasets: List of SingleCellData objects
        batch_labels: List of batch identifiers  
        method: Batch correction method
    """
    
    # 1. Concatenate datasets
    print("Concatenating datasets...")
    integrated = concatenate_datasets(
        datasets,
        batch_key='batch',
        batch_categories=batch_labels
    )
    
    # 2. Standard preprocessing
    print("Preprocessing integrated dataset...")
    integrated = calculate_qc_metrics(integrated)
    integrated = filter_cells(integrated, min_genes=200, max_pct_mt=20)
    integrated = filter_genes(integrated, min_cells=3)
    integrated = normalize_counts(integrated, method='total_count')
    integrated = log_transform(integrated)
    
    # 3. Feature selection (on full dataset)
    integrated = select_hvgs(integrated, n_top_genes=2000)
    
    # 4. Batch correction
    print(f"Applying {method} batch correction...")
    if method == 'scaling':
        integrated = batch_correction_scaling(integrated, batch_key='batch')
    elif method == 'combat':
        integrated = batch_correction_combat(integrated, batch_key='batch')
    elif method == 'harmony':
        integrated = batch_correction_harmony(integrated, batch_key='batch')
    
    # 5. Dimensionality reduction
    integrated = compute_pca(integrated, n_components=50)
    integrated = compute_neighbors(integrated)
    integrated = compute_umap(integrated)
    
    print("Integration complete!")
    return integrated

# Example usage
datasets = [sample1_data, sample2_data, sample3_data]
batch_labels = ['Ctrl', 'Treated_1', 'Treated_2']

integrated_data = integrate_datasets(datasets, batch_labels, method='scaling')
```

## Quality Assessment

### integration_metrics()

Assess integration quality:

```python
from metainformant.singlecell.integration import integration_metrics

# Calculate integration metrics
metrics = integration_metrics(
    integrated_data,
    batch_key='batch',
    label_key='cell_type',   # Cell type labels (if available)
    embedding='X_umap'       # Embedding to evaluate
)

print("Integration Quality Metrics:")
print(f"  Batch mixing entropy: {metrics['mixing_entropy']:.3f}")
print(f"  Silhouette score: {metrics['silhouette_score']:.3f}")
print(f"  ARI (batch vs clusters): {metrics['batch_ari']:.3f}")
```

### plot_integration_assessment()

Visualize integration quality:

```python
from metainformant.singlecell.integration import plot_integration_assessment

# Before and after batch correction comparison
fig = plot_integration_assessment(
    data_before=uncorrected_data,
    data_after=corrected_data,
    batch_key='batch',
    embedding='umap'
)
plt.show()
```

## Advanced Integration Methods

### Multi-Modal Integration

```python
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
    
    # 6. Integrate
    integrated = concatenate_datasets(
        [human_subset, mouse_subset],
        batch_key='species',
        batch_categories=['Human', 'Mouse']
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

Integration tests are available in the test suite:

```bash
# Run integration tests
uv run pytest tests/test_singlecell_integration.py -v
```

## Related Documentation

- [Preprocessing](./preprocessing.md): For data preparation before integration
- [Dimensionality Reduction](./dimensionality.md): For post-integration analysis
- [Clustering](./clustering.md): For clustering integrated data
- [Visualization](./visualization.md): For plotting integration results
