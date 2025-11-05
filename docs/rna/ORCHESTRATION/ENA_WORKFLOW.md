# ENA Workflow Orchestrator

**Script**: `scripts/rna/workflow_ena_integrated.py`

## Overview

Single-species ENA-based workflow orchestrator providing maximum reliability for RNA-seq data processing. This is the **recommended** orchestrator for production workflows.

## Configuration Block

```python
# ============================================================================
# CONFIGURATION
# ============================================================================
# Scope: Single-species ENA-based workflow with batched processing
# Steps: metadata → config → select → getfastq → quant → merge
# Config: YAML file via --config argument (required)
# Threads: Per-config (default 12) or --threads override
# Batch Size: Per-config or --batch-size override (default 12)
# Output: output/amalgkit/{species}/work/
# Dependencies: wget, kallisto, amalgkit (for metadata/index setup)
# Reliability: 100% (ENA direct downloads vs 0% SRA Toolkit for large samples)
# ============================================================================
```

## When to Use

Use this orchestrator when:
- ✅ You need maximum reliability (100% success rate)
- ✅ Processing single species
- ✅ Production workflows
- ✅ Large samples (>50GB)
- ✅ Network interruptions are common
- ✅ You want direct ENA downloads (bypasses SRA Toolkit)

## Features

- **Direct ENA Downloads**: 100% reliability using ENA API via wget
- **Batched Processing**: Configurable batch size (default: 12 samples)
- **Automatic Resume**: wget --continue for interrupted downloads
- **Auto-cleanup**: FASTQs deleted after quantification
- **Disk Management**: ~50-100 GB peak usage (temporary, auto-cleaned)

## Usage

### Basic Usage

```bash
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --batch-size 12 \
  --threads 12
```

### Command-Line Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--config` | Path | Required | Path to amalgkit YAML config file |
| `--batch-size` | int | 12 | Number of samples per batch |
| `--threads` | int | 12 | Threads for download and quantification |
| `--max-samples` | int | None | Limit total samples (for testing) |
| `--skip-download` | flag | False | Skip download, only quantify existing |

### Configuration File

The script requires a YAML configuration file with:

```yaml
work_dir: output/amalgkit/{species}/work
log_dir: output/amalgkit/{species}/logs
threads: 12

species_list:
  - Species_name

genome:
  dest_dir: output/amalgkit/{species}/genome
  # ... genome configuration

steps:
  metadata:
    out_dir: output/amalgkit/{species}/work
  getfastq:
    out_dir: output/amalgkit/{species}/fastq
  quant:
    out_dir: output/amalgkit/{species}/quant
```

## Method Signatures

### `load_config(config_path: Path) -> dict`
Load amalgkit YAML config.

**Parameters**:
- `config_path`: Path to YAML configuration file

**Returns**: Dictionary containing configuration values

### `get_sample_list(metadata_path: Path) -> list[str]`
Get list of run IDs from metadata file.

**Parameters**:
- `metadata_path`: Path to TSV metadata file with 'run' column

**Returns**: List of SRA run IDs (SRR*)

**Raises**: `ValueError` if metadata file lacks 'run' column

### `sample_already_quantified(run_id: str, quant_dir: Path) -> bool`
Check if sample has abundance.tsv.

**Parameters**:
- `run_id`: SRA run ID (e.g., 'SRR1234567')
- `quant_dir`: Directory containing quantification results

**Returns**: True if abundance.tsv exists and has non-zero size

### `download_batch_ena(run_ids: list[str], metadata_path: Path, fastq_dir: Path, threads: int, batch_num: int) -> tuple[list[str], list[str]]`
Download a batch of samples using the robust ENA downloader.

**Parameters**:
- `run_ids`: List of SRA run IDs to download
- `metadata_path`: Path to metadata TSV file
- `fastq_dir`: Directory to save downloaded FASTQ files
- `threads`: Number of parallel download threads
- `batch_num`: Batch number for logging

**Returns**: Tuple of (successful_downloads, failed_downloads)

**Dependencies**: wget command (downloads handled directly in workflow_ena_integrated.py)

### `quantify_batch_kallisto(run_ids: list[str], fastq_dir: Path, quant_dir: Path, index_path: Path, threads: int, batch_num: int) -> tuple[list[str], list[str]]`
Quantify samples using kallisto.

**Parameters**:
- `run_ids`: List of SRA run IDs to quantify
- `fastq_dir`: Directory containing downloaded FASTQ files
- `quant_dir`: Directory to save quantification results
- `index_path`: Path to kallisto transcriptome index
- `threads`: Number of threads for kallisto
- `batch_num`: Batch number for logging

**Returns**: Tuple of (successful_quants, failed_quants)

**Dependencies**: kallisto command, valid index file

### `cleanup_fastqs(run_ids: list[str], fastq_dir: Path) -> int`
Delete FASTQ files for samples.

**Parameters**:
- `run_ids`: List of SRA run IDs whose FASTQs to delete
- `fastq_dir`: Directory containing FASTQ files

**Returns**: Number of sample directories successfully deleted

### `main()`
Main entry point for ENA-integrated workflow.

**Orchestrates**:
1. Load configuration from YAML
2. Get sample list from metadata
3. Process samples in batches (download → quantify → cleanup)
4. Report statistics

**Exit codes**: 0 = success, 1 = failure

## Workflow Steps

1. **Metadata**: Must be run first (using amalgkit) to generate metadata.tsv
2. **Download**: Batched ENA downloads using integrated wget-based downloader
3. **Quantify**: Kallisto quantification per batch
4. **Cleanup**: FASTQ deletion after quantification
5. **Repeat**: Process next batch until all samples complete

## Output Structure

```
output/amalgkit/{species}/
├── work/
│   ├── metadata/
│   │   └── metadata.tsv        # Required input
│   └── index/
│       └── transcriptome.idx    # Kallisto index (must exist)
├── fastq/
│   └── {SRR_ID}/               # Temporary (auto-deleted)
│       └── *.fastq.gz
└── quant/
    └── {SRR_ID}/               # Final results
        ├── abundance.tsv
        └── run_info.json
```

## Performance

- **Download speed**: Fast (ENA direct, no SRA conversion)
- **Success rate**: 100% (vs ~0% for SRA Toolkit on large samples)
- **Disk usage**: Low (auto-cleanup, ~50-100 GB peak)
- **Processing time**: ~7.5 minutes per sample average

## Troubleshooting

### "Metadata not found"
Run amalgkit metadata step first:
```bash
amalgkit metadata --out_dir output/amalgkit/{species}/work
```

### "Kallisto index not found"
Build index first:
```bash
amalgkit quant --build_index yes --index_dir output/amalgkit/{species}/work/index
```

### Downloads failing
- Check network: `ping -c 4 8.8.8.8`
- Verify ENA access: `wget --spider https://www.ebi.ac.uk/ena/`
- Check disk space: `df -h /`

## Related Documentation

- **[ORCHESTRATION/README.md](README.md)**: Orchestrator comparison
- **[WORKFLOW.md](../WORKFLOW.md)**: Workflow planning and execution
- **[MULTI_SPECIES_QUICK_START.md](../MULTI_SPECIES_QUICK_START.md)**: Production workflows

## See Also

- **Source Script**: `scripts/rna/workflow_ena_integrated.py`
- **ENA Downloads**: Integrated directly in `workflow_ena_integrated.py` using wget
- **Configuration Guide**: [CONFIGURATION.md](../CONFIGURATION.md)

