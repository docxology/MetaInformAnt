# Stage 1 — Data Processing

## Purpose

Ingests all raw CSV files from `data/raw/`, applies configured quality filters, and writes a clean unified dataset to `data/processed/processed_data.csv`.

## Script

`scripts/01_process_data.py`

## Inputs

| Source | Path | Description |
| :--- | :--- | :--- |
| Raw CSVs | `data/raw/*.csv` | One or more measurement files in standard format |
| Config | `config/default.yaml` | Processing thresholds and output paths |
| Metadata | `data/raw/metadata.yaml` | Provenance record (informational only) |

### Expected CSV Schema

| Column | Type | Required | Description |
| :--- | :--- | :--- | :--- |
| `sample_id` | str | [DONE] | Unique observation identifier |
| `group` | str | [DONE] | Experimental group label |
| `feature_*` | float | [DONE] | Numeric measurements |

## Outputs

| Path | Description |
| :--- | :--- |
| `data/processed/processed_data.csv` | Filtered, optionally normalised dataset |
| `logs/01_process_data.log` | Structured execution log with timing |

## Config Keys

| Key | Default | Effect |
| :--- | :--- | :--- |
| `processing.missing_fraction_max` | `0.2` | Drop columns exceeding this missing-value fraction |
| `processing.min_sample_count` | `5` | Drop rows with fewer non-null values than this |
| `processing.normalize` | `true` | Apply per-column z-score normalisation |
| `paths.data_raw` | `data/raw/` | Source directory |
| `paths.data_processed` | `data/processed/` | Output directory |

## Usage

```bash
# Standard run (idempotent — skips if output exists)
uv run scripts/01_process_data.py --config config/default.yaml

# Force reprocessing
uv run scripts/01_process_data.py --config config/default.yaml --force

# With custom config
uv run scripts/01_process_data.py --config config/my_experiment.yaml
```

## Idempotency

The script checks for `data/processed/processed_data.csv` at startup. If it exists, processing is skipped unless `--force` is passed. This makes pipeline reruns cheap.

## Processing Steps

1. **Discover** all `*.csv` files in `paths.data_raw`
2. **Load** each file into a DataFrame; attach `_source_file` column for traceability
3. **Concatenate** all frames into one unified table
4. **Filter columns**: drop numeric columns where missing fraction > `missing_fraction_max`
5. **Filter rows**: drop rows with fewer non-null values than `min_sample_count`
6. **Normalise**: if `normalize=true`, apply z-score (zero mean, unit variance) per numeric column
7. **Write** to `paths.data_processed/processed_data.csv`

## Extension Hooks

To add domain-specific processing (e.g., read trimming, sequence quality filtering):

```python
# In 01_process_data.py → process_data()
# After the normalisation block, add:
from metainformant.dna.sequences import validate_dna_sequence
# ... your domain logic here
```

## Troubleshooting

| Symptom | Cause | Fix |
| :--- | :--- | :--- |
| `No CSV files found in data/raw/` | Raw data missing | Run `uv run scripts/99_create_synthetic_data.py` first |
| Output empty / zero rows | Filters too strict | Increase `min_sample_count` or `missing_fraction_max` |
| `UnicodeDecodeError` reading CSV | Non-UTF-8 encoding | Pre-convert files: `iconv -f latin1 -t utf-8 input.csv > output.csv` |
