# Single-Cell Genomics Scripts

Single-cell RNA sequencing analysis and visualization workflow orchestrators.

## Directory Structure

```
scripts/singlecell/
├── run_singlecell_analysis.py     # Single-cell analysis workflow orchestrator
└── README.md                      # This file
```

## Single-Cell Analysis Workflow (`run_singlecell_analysis.py`)

Comprehensive single-cell analysis workflow orchestrator for scRNA-seq data processing and analysis.

**Features:**
- Quality control and filtering
- Normalization and scaling
- Dimensionality reduction
- Clustering and cell type identification
- Differential expression analysis
- Trajectory inference

**Usage:**
```bash
# Basic single-cell analysis
python3 scripts/singlecell/run_singlecell_analysis.py --input counts.h5ad --output output/singlecell/basic

# Full analysis pipeline
python3 scripts/singlecell/run_singlecell_analysis.py --input counts.h5ad --quality-control --normalize --cluster --find-markers
```

**Options:**
- `--input`: Input single-cell data (h5ad, CSV, or 10x format)
- `--output`: Output directory (defaults to output/singlecell/)
- `--quality-control`: Perform quality control filtering
- `--normalize`: Apply normalization and scaling
- `--cluster`: Perform cell clustering
- `--find-markers`: Identify cluster marker genes
- `--trajectory`: Perform trajectory inference
- `--threads`: Number of threads to use
- `--verbose`: Enable verbose logging

**Output Structure:**
```
output/singlecell/
├── quality_control/               # QC and filtering results
│   ├── qc_metrics.json
│   ├── filtered_cells.h5ad
│   └── qc_plots/
├── preprocessing/                 # Normalization results
│   ├── normalized_data.h5ad
│   ├── scaled_data.h5ad
│   └── preprocessing_summary.json
├── dimensionality_reduction/      # Dimension reduction results
│   ├── pca_components.h5ad
│   ├── umap_coordinates.csv
│   └── tsne_coordinates.csv
├── clustering/                    # Clustering results
│   ├── cluster_assignments.json
│   ├── cluster_markers.json
│   └── clustering_summary.json
├── differential_expression/       # DE analysis results
│   ├── de_results.json
│   ├── volcano_plots/
│   └── heatmap_plots/
├── trajectory_analysis/           # Trajectory inference results
│   ├── trajectory_graph.json
│   ├── pseudotime_values.json
│   └── trajectory_plots/
└── analysis_report.json           # Comprehensive analysis report
```

## Integration

Integrates with:
- **metainformant.singlecell**: Core single-cell analysis functionality
- **Scanpy**: Single-cell analysis framework
- **Core utilities**: I/O, logging, path management

## Dependencies

- **metainformant.singlecell**: Single-cell analysis module
- **scanpy**: Single-cell analysis
- **anndata**: Annotated data structures
- **matplotlib/seaborn**: Visualization support

## Related Documentation

- [Single-Cell Analysis Documentation](../../docs/singlecell/README.md)
- [METAINFORMANT CLI](../../docs/cli.md)

