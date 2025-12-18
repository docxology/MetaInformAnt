# Multi-Omics Integration Scripts

Multi-omics data integration and joint analysis workflow orchestrators.

## Directory Structure

```
scripts/multiomics/
├── run_multiomics_integration.py  # Main multi-omics integration orchestrator ⭐
├── example_basic_integration.py   # Basic integration examples
├── example_cca.py                 # Canonical correlation analysis examples
├── example_joint_analysis.py      # Joint analysis examples
└── README.md                      # This file
```

## Multi-Omics Integration (`run_multiomics_integration.py`)

Comprehensive multi-omics integration workflow orchestrator for joint analysis across genomics, transcriptomics, proteomics, and other omics layers.

**Features:**
- Multi-omics data loading and preprocessing
- Joint dimensionality reduction (PCA, NMF)
- Canonical correlation analysis
- Cross-omics correlation analysis
- Integration quality assessment

**Usage:**
```bash
# Basic integration with two omics layers
python3 scripts/multiomics/run_multiomics_integration.py --genomics genomics.csv --transcriptomics expression.tsv --output output/multiomics/basic

# Full analysis with joint PCA
python3 scripts/multiomics/run_multiomics_integration.py --genomics g.csv --transcriptomics t.tsv --proteomics p.csv --joint-pca --n-components 50

# Canonical correlation analysis
python3 scripts/multiomics/run_multiomics_integration.py --genomics g.csv --transcriptomics t.tsv --canonical-correlation --n-components 10
```

**Options:**
- `--genomics`: Genomics data file (CSV/TSV)
- `--transcriptomics`: Transcriptomics data file
- `--proteomics`: Proteomics data file
- `--methylomics`: Methylomics/epigenomics data file
- `--joint-pca`: Perform joint PCA analysis
- `--canonical-correlation`: Perform CCA analysis
- `--correlation-analysis`: Cross-omics correlation analysis
- `--output`: Output directory (defaults to output/multiomics/)
- `--threads`: Number of threads to use
- `--verbose`: Enable verbose logging

## Example Scripts

### Basic Integration (`example_basic_integration.py`)

Basic multi-omics integration examples demonstrating core functionality.

**Usage:**
```bash
# Run basic integration examples
python3 scripts/multiomics/example_basic_integration.py
```

### CCA Examples (`example_cca.py`)

Canonical correlation analysis examples and tutorials.

**Usage:**
```bash
# Run CCA examples
python3 scripts/multiomics/example_cca.py
```

### Joint Analysis (`example_joint_analysis.py`)

Joint analysis examples across multiple omics layers.

**Usage:**
```bash
# Run joint analysis examples
python3 scripts/multiomics/example_joint_analysis.py
```

**Output Structure:**
```
output/multiomics/
├── integrated_data/               # Preprocessed and integrated datasets
│   ├── normalized_genomics.csv
│   ├── normalized_transcriptomics.csv
│   └── integrated_metadata.json
├── dimensionality_reduction/      # Joint dimensionality reduction results
│   ├── joint_pca_components.csv
│   ├── joint_pca_loadings.json
│   └── explained_variance.json
├── canonical_correlation/         # CCA analysis results
│   ├── canonical_weights.json
│   ├── canonical_scores.csv
│   └── correlation_coefficients.json
├── correlation_analysis/          # Cross-omics correlations
│   ├── correlation_matrix.json
│   └── correlation_network.json
├── integration_plots/             # Generated visualizations
│   ├── joint_pca_plot.png
│   ├── correlation_heatmap.png
│   ├── cca_scatterplot.png
│   └── integration_quality.png
└── analysis_report.json           # Comprehensive analysis report
```

## Key Features

✅ **Multi-Omics Support**: Genomics, transcriptomics, proteomics, epigenomics
✅ **Advanced Integration**: Joint PCA, CCA, correlation analysis
✅ **Scalable Processing**: Efficient algorithms for large datasets
✅ **Quality Assessment**: Integration quality metrics and diagnostics
✅ **Comprehensive Visualization**: Rich plotting for integrated results
✅ **Modular Design**: Individual analysis components

## Integration

Integrates with:
- **metainformant.multiomics**: Core multi-omics integration functionality
- **Scientific computing**: NumPy, SciPy, scikit-learn
- **Statistical analysis**: pandas, statsmodels
- **Visualization**: matplotlib, seaborn, plotly
- **Core utilities**: I/O, logging, configuration management

## Dependencies

- **metainformant.multiomics**: Multi-omics integration module
- **scikit-learn**: Machine learning and dimensionality reduction
- **pandas/numpy**: Data manipulation and statistics
- **matplotlib/seaborn**: Visualization support

## Related Documentation

- [Multi-Omics Documentation](../../docs/multiomics/README.md)
- [Integration Methods](../../docs/multiomics/integration_methods.md)
- [Data Preparation](../../docs/multiomics/data_preparation.md)
- [METAINFORMANT CLI](../../docs/cli.md)

