# Stage 2 — Results Analysis

## Purpose

Reads the processed dataset from Stage 1 and computes statistical summaries, optional correlation analysis, and optional PCA. Outputs tables to `results/tables/`.

## Script

`scripts/02_analyze_results.py`

## Inputs

| Source | Path | Description |
| :--- | :--- | :--- |
| Processed data | `data/processed/processed_data.csv` | Output of Stage 1 |
| Config | `config/default.yaml` | Analysis method, number of components, significance threshold |

## Outputs

| Path | Description |
| :--- | :--- |
| `results/tables/summary_statistics.csv` | Descriptive stats (mean, std, min, max, CV) per feature |
| `results/tables/correlation_matrix.csv` | Pairwise Pearson correlations (`method=summary` or `correlation`) |
| `results/tables/pca_summary.json` | Explained variance ratio per component (`method=pca`) |
| `results/tables/analysis_metadata.json` | Run provenance (method, row/column counts) |
| `logs/02_analyze_results.log` | Structured execution log |

## Config Keys

| Key | Default | Effect |
| :--- | :--- | :--- |
| `analysis.method` | `"summary"` | `summary` (stats + correlation), `correlation` (correlation only), `pca` |
| `analysis.n_components` | `10` | PCA components to retain (PCA mode only) |
| `analysis.significance_threshold` | `0.05` | Alpha (reserved for domain fork extension) |
| `analysis.correction` | `"bonferroni"` | Multiple-testing method (reserved for domain fork extension) |

## Usage

```bash
# Standard run (idempotent)
uv run scripts/02_analyze_results.py --config config/default.yaml

# Force rerun
uv run scripts/02_analyze_results.py --config config/default.yaml --force

# PCA mode
# Set analysis.method: "pca" in your config, then:
uv run scripts/02_analyze_results.py --config config/pca_experiment.yaml
```

## Analysis Methods

### `summary` (default)
Runs both descriptive statistics **and** Pearson correlation matrix.

### `correlation`
Runs Pearson correlation matrix only (faster for large feature sets).

### `pca`
Runs SVD-based PCA (no scikit-learn dependency). Outputs:
- `n_components` principal component explained variance ratios
- `total_variance_explained`: sum of retained variance
- Written to `results/tables/pca_summary.json`

## Idempotency

Checks for `results/tables/summary_statistics.csv`. Skips if present unless `--force` is passed.

## Extension Hooks

Replace or extend the analysis with MetaInformAnt modules:

```python
from metainformant.dna.population import nucleotide_diversity, tajimas_d
from metainformant.ml.classification import train_classifier
```

## Troubleshooting

| Symptom | Cause | Fix |
| :--- | :--- | :--- |
| Exits non-zero | Stage 1 not run | Run `01_process_data.py` first |
| `Insufficient numeric columns for PCA` | All columns non-numeric or missing | Check raw data schema |
| All correlations NaN | Data not normalised | Enable `processing.normalize: true` in config |
