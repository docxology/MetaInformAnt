# Comprehensive RNA-seq Pipeline Guide

This guide details the METAINFORMANT RNA-seq pipeline, covering single-species execution, multi-species orchestration, configuration, and troubleshooting.

## 1. Pipeline Overview

The pipeline automates the processing of RNA-seq data using `amalgkit`, providing a robust, restartable, and parallelized workflow.

### Key Features

* **Automated Data Retrieval**: Direct download from ENA/SRA with fallback mechanisms.
* **Quality Control**: Integrated `fastp` support.
* **Quantification**: `kallisto`-based pseudo-alignment.
* **Analysis**: Differential expression, tissue specificity, and cross-species comparison.
* **Resource Management**: Automatic cleanup of intermediate FASTQ files to save disk space.

## 2. Orchestration Strategies

We provide two primary orchestration scripts to handle different workflow needs.

### A. Single-Species / High-Volume: `run_workflow.py`

**Best for**: Production runs, large datasets (e.g., *Apis mellifera*), and deep troubleshooting.

* **Script**: `scripts/rna/run_workflow.py`
* **Behavior**: Processes a single species configuration end-to-end.
* **Parallelism**: Controlled via `num_download_workers` in config.

**Usage**:

```bash
# Run standard workflow
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_apis_mellifera_all.yaml

# Check status
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_apis_mellifera_all.yaml --status

# Cleanup unquantified samples
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_apis_mellifera_all.yaml --cleanup-unquantified
```

### B. Multi-Species / Batch: `orchestrate_species.py`

**Best for**: Processing multiple species sequentially, automated nightly runs.

* **Script**: `scripts/rna/orchestrate_species.py`
* **Behavior**: Iterates through a list of species in a master config, running `amalgkit` steps for each.
* **Isolation**: Ensures failure in one species does not crash the entire batch (soft failure).

**Usage**:

```bash
# Run multi-species orchestration
uv run python scripts/rna/orchestrate_species.py --config config/amalgkit/cross_species.yaml
```

**Configuration (`cross_species.yaml`)**:

```yaml
work_dir: output/amalgkit
species:
  - Ooceraea_biroi
  - Linepithema_humile
  - Solenopsis_invicta
steps:
  - metadata
  - getfastq
  - quant
  - merge
```

## 3. Configuration

Configuration is managed via YAML files in `config/amalgkit/`.

### Example Species Config

```yaml
work_dir: blue/amalgkit/species_name/work
threads: 12

steps:
  getfastq:
    num_download_workers: 8  # Parallel downloads
    threads: 4               # Threads per download/process
    aws: yes                 # Use AWS SRA mirror
  quant:
    index_dir: path/to/index
    fasta_dir: path/to/fasta
  merge:
    out_dir: blue/amalgkit/species_name/work  # IMPORTANT: Must match quant output root
```

### Path Resolution

* **`quant`**: Outputs to `{work_dir}/quant/{sample_id}/`.
* **`merge`**: Looks for quantification results in `{out_dir}/quant/`.
* **Critical**: Ensure `merge:out_dir` points to the *parent* directory containing the `quant` folder (usually `work_dir`), NOT `work_dir/merged`.

## 4. Troubleshooting

### `IsADirectoryError` in `merge`

* **Cause**: `merge` step looking for input files in the wrong directory structure.
* **Fix**: Set `merge:out_dir` to the same path as `quant:out_dir` (the `work` directory).
  * *Incorrect*: `out_dir: blue/amalgkit/species/merged`
  * *Correct*: `out_dir: blue/amalgkit/species/work`

### `curate` Failures

* **Cause**: R script errors (e.g., missing columns, NA values).
* **Fix**: Check `blue/amalgkit/species/logs/curate.stderr.log`.
  * Use `skip_curation: yes` in config to bypass if curation is non-critical for initial analysis.

### Download Failures

* **Cause**: Network timeouts or SRA issues.
* **Fix**: `run_workflow.py` automatically retries. For persistent failures, try reducing `num_download_workers`.

## 5. Development & Testing

* Tests: Located in `tests/`.
  * `test_rna_amalgkit_steps.py`: Unit tests for individual steps.
  * `test_rna_orchestration_new.py`: Tests for `PipelineOrchestrator`.
  * `test_rna_curate.py`: Integration tests for `curate` wrapper.
* **Mocking Policy**: We use **Zero-Mock** policy where possible, but for CLI wrappers we mock `subprocess.run` to verify command construction.

---
*Generated February 2026 for METAINFORMANT RNA-seq V2*
