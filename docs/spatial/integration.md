# Multi-Sample Integration

Batch-correct and spatially align multiple tissue sections.

## Concatenation + BBKNN

```python
import scanpy as sc
from metainformant.spatial.integration import spatial_batch_correct

adatas = [load_visium(p) for p in dirs]
adata_joint = spatial_batch_correct(
    adatas,
    batch_key='sample_id',
    method='bbknn',        # also 'harmony' or 'scanorama'
    n_pcs=30,
    neighbors=15,
)
```

After integration, re-run spatial clustering on the corrected embedding:

```python
sc.pp.neighbors(adata_joint, use_rep='X_pca')
adata_joint.obs['integrated_domain'] = spatial_cluster(adata_joint, resolution=0.8)
```

## Transfer of labels

Project labels from a reference section to a new query:

```python
from metainformant.spatial.integration import transfer_labels
query = load_visium('new_sample/')
pred_labels = transfer_labels(
    reference=adata_joint,
    query=query,
    label_col='integrated_domain',
    method='knn',
)
query.obs['predicted_domain'] = pred_labels
```

## Spatial alignment (registration)

Rigid / non-rigid alignment of tissue shapes:

```python
from metainformant.spatial.integration import register_coordinates
aligned = register_coordinates(
    adatas,
    reference_idx=0,
    method=' affine',    # 'affine' or 'dense'
    landmark_key='landmarks',   # optional user-supplied points
)
```

Used when comparing histological landmarks across patients.

## Harmonization across platforms

If you mix Visium and Xenium data, down-sample Xenium to Visium-like spot spacing:

```python
from metainformant.spatial.integration import downsample_to_visium
visium_like = downsample_to_visium(xenium_adata, target_diameter=55.0)
```
