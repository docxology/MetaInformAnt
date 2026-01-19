# SPEC: SingleCell Scripts

Processing pipelines for scRNA-seq and single-cell multi-omics.

## Workflows

- `run_sc_pipeline.py`: End-to-end preprocessing, clustering, and UMAP generation.
- `infer_trajectories.py`: Orchestrates pseudotime analysis and lineage tracing.

## Standards

- **Performance**: High-memory tasks should utilize the `ParallelProcessor`.
- **Validation**: All input `h5ad` files are validated against the `anndata` schema.
