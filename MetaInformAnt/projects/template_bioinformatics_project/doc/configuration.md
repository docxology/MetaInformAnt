# Configuration Reference

Complete reference for every key in `config/default.yaml`.

---

## `metadata`

Human-readable project metadata. These values appear in log headers and provenance files. They have no effect on pipeline logic.

| Key | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `project_name` | str | `"template_bioinformatics_project"` | Identifier for this project instance |
| `version` | str | `"0.1.0"` | Project version string (SemVer recommended) |
| `description` | str | `"..."` | Free-text description included in logs |
| `author` | str | `""` | Analyst or team name |
| `date` | str | `""` | ISO 8601 date of the experiment |
| `organism` | str | `""` | Organism name (e.g. `"Apis mellifera"`) |
| `genome_accession` | str | `""` | NCBI genome accession (if applicable) |

---

## `paths`

All file I/O is routed through these keys. Paths are relative to the project root. Scripts create missing directories automatically.

> [!IMPORTANT]
> Never write output to paths outside these configured directories. Stray files in the project root indicate a misconfiguration.

| Key | Default | Contents |
| :--- | :--- | :--- |
| `data_raw` | `data/raw/` | Immutable raw input files (CSV, YAML metadata) |
| `data_processed` | `data/processed/` | Cleaned, normalised data from Stage 1 |
| `results` | `results/` | Parent for all result subtrees |
| `results_figures` | `results/figures/` | Matplotlib figures from Stage 3 |
| `results_tables` | `results/tables/` | Summary tables and JSON from Stage 2 |
| `logs` | `logs/` | Per-stage `.log` files |
| `config` | `config/default.yaml` | Self-reference (used for provenance logging) |

---

## `processing`

Controls Stage 1 (`01_process_data.py`) behaviour.

| Key | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `threads` | int | `4` | Number of worker threads (reserved for parallelism in domain forks) |
| `filtering_threshold` | float | `10.0` | Generic threshold value (domain-specific meaning; log-scaled by default) |
| `normalize` | bool | `true` | If `true`, apply z-score normalisation per numeric column |
| `min_sample_count` | int | `5` | Minimum non-null values per row; rows below this are dropped |
| `missing_fraction_max` | float | `0.2` | Drop columns where the fraction of missing values exceeds this |

---

## `analysis`

Controls Stage 2 (`02_analyze_results.py`) behaviour.

| Key | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `method` | str | `"summary"` | Analysis type: `summary` (descriptive stats + correlation), `correlation` (correlation only), `pca` (PCA via SVD) |
| `n_components` | int | `10` | Number of principal components to retain (PCA mode only) |
| `significance_threshold` | float | `0.05` | Alpha level for statistical tests (reserved for domain forks) |
| `correction` | str | `"bonferroni"` | Multiple-testing correction method: `bonferroni`, `fdr` (reserved for domain forks) |

---

## `visualization`

Controls Stage 3 (`03_visualize.py`) output format and appearance.

| Key | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `dpi` | int | `150` | Figure resolution in dots per inch |
| `figure_width` | float | `10` | Base figure width in inches |
| `figure_height` | float | `6` | Base figure height in inches |
| `color_palette` | str | `"viridis"` | Matplotlib/seaborn colour palette name |
| `file_format` | str | `"png"` | Output format: `png`, `pdf`, or `svg` |

---

## `logging`

Controls log format and verbosity across all scripts.

| Key | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `level` | str | `"INFO"` | Log verbosity: `DEBUG`, `INFO`, `WARNING`, `ERROR` |
| `format` | str | `"%(asctime)s [%(levelname)s] %(name)s — %(message)s"` | Python `logging` format string |

---

## Example: Custom Experiment Config

Create `config/my_experiment.yaml` with only the keys you wish to override:

```yaml
metadata:
  project_name: "my_species_analysis"
  organism: "Drosophila melanogaster"
  date: "2026-03-19"

processing:
  normalize: false
  missing_fraction_max: 0.1

analysis:
  method: "pca"
  n_components: 5

visualization:
  dpi: 300
  file_format: "pdf"
```

Then run:

```bash
./run.sh --config config/my_experiment.yaml
```
