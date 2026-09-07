# Stage 99 — Synthetic Data Generator

## Purpose

Generates realistic synthetic CSV datasets and a provenance `metadata.yaml` in `data/raw/` to enable complete end-to-end pipeline testing **without depending on external data sources**.

This script is a test utility — it is not part of the scientific pipeline. Run it on a new machine or in CI to generate test data before running Stages 1–3.

## Script

`scripts/99_create_synthetic_data.py`

## Outputs

| Path | Description |
| :--- | :--- |
| `data/raw/samples_A.csv` | Primary synthetic dataset (`--n-samples` rows) |
| `data/raw/samples_B.csv` | Secondary dataset (50% of primary size) |
| `data/raw/metadata.yaml` | Provenance record for generated files |

### Output CSV Schema

| Column | Description |
| :--- | :--- |
| `sample_id` | Sequential identifier (`S0001`, `S0002`, …) |
| `group` | Randomly assigned from `{control, treatment_A, treatment_B}` |
| `feature_01` … `feature_NN` | Normally distributed floats; ~5% cells randomly set to NaN |

## Usage

```bash
# Default: 200-row dataset, 8 features, seed=2026
uv run scripts/99_create_synthetic_data.py

# Custom size and seed
uv run scripts/99_create_synthetic_data.py --n-samples 1000 --n-features 20 --seed 42

# Regenerate even if data already exists
uv run scripts/99_create_synthetic_data.py --force

# With custom config (to resolve data/raw/ path)
uv run scripts/99_create_synthetic_data.py --config config/my_experiment.yaml
```

## CLI Arguments

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--config` | `config/default.yaml` | Config YAML to resolve `paths.data_raw` |
| `--n-samples` | `200` | Number of rows in `samples_A.csv` |
| `--n-features` | `8` | Number of numeric feature columns |
| `--seed` | `2026` | NumPy random seed for reproducibility |
| `--force` | — | Overwrite existing files |

## Idempotency

If any `samples_*.csv` files already exist in `data/raw/`, the script exits cleanly without overwriting. Use `--force` to regenerate.

## Use in Tests

The test suite (`tests/test_pipeline.py`) invokes this script via subprocess in the `TestEndToEnd` class:

```python
run_script(SCRIPTS / "99_create_synthetic_data.py", config_path,
           "--n-samples", "60", "--seed", "99")
```

## Notes

> [!WARNING]
> Synthetic data is entirely fictional. **Do not use for scientific conclusions.** Replace with real acquisition (see [Data Management Guide](../data.md)) before production use.

The `metadata.yaml` produced by this script clearly records `purpose: "Synthetic data for end-to-end pipeline testing"`.
