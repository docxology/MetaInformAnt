# Data Management Guide

This document describes the data layout, file format conventions, immutability rules, and provenance practices for the template bioinformatics project.

## Directory Conventions

```text
data/
├── raw/             ← ⚠ IMMUTABLE — treat as read-only after acquisition
│   ├── samples_A.csv
│   ├── samples_B.csv
│   └── metadata.yaml
└── processed/       ← Written by Stage 1; inputs to Stages 2–3
    └── processed_data.csv
```

> [!CAUTION]
> **`data/raw/` is immutable.** No script should ever write to, modify, or delete files in `data/raw/`. If raw data must be corrected, make a documented copy in `data/raw/corrected/` and update `metadata.yaml`.

## File Formats

### Raw Input: CSV

Raw datasets in `data/raw/*.csv` must have:

| Column | Type | Description |
| :--- | :--- | :--- |
| `sample_id` | str | Unique identifier per observation (e.g., `S0001`) |
| `group` | str | Experimental group label (e.g., `control`, `treatment_A`) |
| `feature_NN` | float | Numeric measurements (any number of columns) |

Missing values should be represented as empty cells (not `NA`, `null`, or `NaN` strings).

### Raw Input: `metadata.yaml`

Every `data/raw/` directory should contain a `metadata.yaml` provenance file:

```yaml
generated_at: "2026-03-19T16:25:00Z"
generator: "scripts/99_create_synthetic_data.py"
purpose: "End-to-end pipeline testing"
datasets:
  - file: "samples_A.csv"
    n_rows: 200
    n_features: 8
  - file: "samples_B.csv"
    n_rows: 100
    n_features: 8
notes: "Replace this with a description of your real data acquisition."
```

For real projects, record the data source URL, download date, and checksum.

### Processed Output: `processed_data.csv`

Produced by Stage 1. Same schema as raw CSVs but:
- Missing-value-filtered columns removed (above `processing.missing_fraction_max`)
- Low-count rows dropped (below `processing.min_sample_count`)
- Optionally z-score normalised per numeric column
- `_source_file` column added for traceability

## Acquisition Workflow (Real Data)

For real projects, replace synthetic data generation with a documented acquisition process:

```bash
# Option A: direct download
curl -O https://example.org/datasets/my_data.csv
mv my_data.csv data/raw/

# Option B: SRA download via metainformant
uv run -c "from metainformant.dna.ncbi import download_sra_run; ..."

# Always record provenance
echo "source: https://..." >> data/raw/metadata.yaml
```

## Immutability Enforcement

The `.gitignore` excludes all file content from `data/raw/`, `data/processed/`, `results/`, and `logs/` to prevent large data commits. Only `.gitkeep` placeholders are tracked:

```text
data/raw/*
!data/raw/.gitkeep
```

This means:
- Data is **never committed to the repository**.
- Data must be re-acquired or re-generated on each new machine.
- The `metadata.yaml` file (committed separately if needed) documents how to reproduce the data.

## Provenance and Reproducibility

Each pipeline run appends to log files in `logs/`. These logs record:
- Exact timestamp of execution
- Config file path and key parameter values
- Input file counts and shapes
- Output file paths and sizes
- Elapsed time

To reconstruct a previous analysis:
1. Check the git commit hash at analysis time: `git log --oneline -1`
2. Restore the config file at that commit.
3. Re-acquire raw data using `metadata.yaml`.
4. Re-run `./run.sh`.

## Adding New Data Modalities

To extend the pipeline with new data types (e.g., FASTQ files, VCFs):

1. Add new subdirectories under `data/raw/` (e.g., `data/raw/sequences/`).
2. Register the new path in `config/default.yaml` under `paths`.
3. Create a new numbered script to handle the new modality.
4. Update `metadata.yaml` to document the new data source.
5. Update this document.
