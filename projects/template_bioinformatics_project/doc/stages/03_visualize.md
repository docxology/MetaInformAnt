# Stage 3 — Visualisation

## Purpose

Reads processed data and summary tables produced by Stages 1–2 and renders publication-quality figures to `results/figures/`.

## Script

`scripts/03_visualize.py`

## Inputs

| Source | Path | Description |
| :--- | :--- | :--- |
| Processed data | `data/processed/processed_data.csv` | Required |
| Correlation matrix | `results/tables/correlation_matrix.csv` | Optional; heatmap only rendered if present |
| Config | `config/default.yaml` | Visual parameters |

## Outputs

| Path | Description |
| :--- | :--- |
| `results/figures/distribution_grid.{fmt}` | Histogram grid — one panel per numeric feature |
| `results/figures/correlation_heatmap.{fmt}` | Lower-triangle correlation heatmap (if correlation data exists) |
| `logs/03_visualize.log` | Structured execution log |

`{fmt}` is controlled by `visualization.file_format` (default: `png`).

## Config Keys

| Key | Default | Effect |
| :--- | :--- | :--- |
| `visualization.dpi` | `150` | Figure resolution (dots per inch) |
| `visualization.figure_width` | `10` | Base width in inches |
| `visualization.figure_height` | `6` | Base height in inches |
| `visualization.color_palette` | `"viridis"` | Any valid matplotlib/seaborn palette name |
| `visualization.file_format` | `"png"` | Output format: `png`, `pdf`, `svg` |

## Usage

```bash
# Standard run
uv run scripts/03_visualize.py --config config/default.yaml

# Force regeneration
uv run scripts/03_visualize.py --config config/default.yaml --force

# PDF output for publication
# Set visualization.file_format: "pdf" in your config, then rerun
```

## Figures Produced

### `distribution_grid.png`
Grid of histograms, one per numeric column in `processed_data.csv`.
- Layout: up to 4 columns, auto-expanded rows.
- Colour: per-column from the configured palette.
- NaN values excluded from each histogram.

### `correlation_heatmap.png`
Lower-triangle Pearson correlation heatmap.
- Only rendered if `results/tables/correlation_matrix.csv` exists.
- Annotated with correlation values when ≤ 12 features (otherwise unlabeled for readability).
- Coloured on `coolwarm` diverging palette centred at 0.

## Idempotency

Checks for all expected output figures. Skips if all are present unless `--force` is passed.

## Extension Hooks

Add new plots using the same config-driven pattern:

```python
def plot_my_figure(data: pd.DataFrame, output_path: Path, vcfg: dict) -> None:
    fig, ax = plt.subplots(figsize=(vcfg["figure_width"], vcfg["figure_height"]))
    # ... plotting logic
    fig.savefig(output_path, dpi=vcfg["dpi"], bbox_inches="tight")
    plt.close(fig)
```

Then call it inside `run_visualize()` and register the expected path in the idempotency check.

## Troubleshooting

| Symptom | Cause | Fix |
| :--- | :--- | :--- |
| Script exits non-zero | Stage 1 not run | Run `01_process_data.py` first |
| Heatmap not produced | Correlation matrix missing | Run `02_analyze_results.py` with `method: summary` first |
| Blank / empty figure file | No numeric columns in data | Check raw data schema and filtering thresholds |
| Matplotlib backend error | Headless environment | Script sets `Agg` backend automatically — should not occur |
