# amalgkit getfastq: FASTQ File Generation from SRA

## Purpose

Downloads SRA files and converts them to FASTQ format for RNA-seq quantification. This step handles the **data acquisition** phase, retrieving raw sequencing reads from NCBI SRA or cloud repositories.

## Overview

The `getfastq` step:
- Downloads SRA files from NCBI, AWS, or GCP
- Extracts FASTQ files using parallel-fastq-dump or fasterq-dump
- Performs quality filtering with fastp
- Handles both single-end and paired-end libraries
- Supports batch processing for HPC environments
- Automatically removes SRA files after extraction (configurable)

## Usage

### Basic Usage

```bash
amalgkit getfastq \
  --out_dir output/amalgkit/work/fastq \
  --metadata output/amalgkit/work/metadata/metadata.tsv \
  --threads 8 \
  --pfd yes \
  --fastp yes

# Note: FASTQ files will be in output/amalgkit/work/fastq/getfastq/{SRR_ID}/
```

**Important**: The metadata file must be in row-per-sample format (with a `run` column), NOT a pivot table. The step automatically looks for `metadata.filtered.tissue.tsv` or `metadata.tsv` in the metadata directory.

### Python API

```python
from metainformant.rna.steps.getfastq import run

result = run(
    params={
        "out_dir": "output/amalgkit/work/fastq",
        "metadata": "output/amalgkit/work/metadata/metadata.tsv",
        "threads": 8,
        "pfd": "yes",
        "fastp": "yes",
    },
    work_dir="output/amalgkit/work",
    log_dir="output/amalgkit/work/logs",
)
```

### Configuration File

```yaml
steps:
  getfastq:
    out_dir: output/amalgkit/amellifera/fastq
    # Note: FASTQ files will be in {out_dir}/getfastq/{SRR_ID}/
    # Example: output/amalgkit/amellifera/fastq/getfastq/SRR12345678/SRR12345678_1.fastq.gz
    metadata: output/amalgkit/amellifera/work/metadata/metadata.tsv  # Row-per-sample format
    threads: 8
    pfd: yes                # Use parallel-fastq-dump (auto-detected if available)
    fastp: yes              # Quality filtering (default: no, set to yes to enable)
    accelerate: true        # Enable cloud acceleration (METAINFORMANT-specific)
    ncbi: yes
    aws: yes
    gcp: yes
    min_read_length: 25
    max_bp: 10000000000  # Optional: Limit download size per sample (10GB)
```

## Parameters

### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--out_dir` | PATH | `./` | Directory for intermediate and output files. |
| `--metadata` | PATH | `inferred` | Path to metadata.tsv. Default: `out_dir/metadata/metadata.tsv` |
| `--threads` | INT | `1` | Number of threads for parallel processing. |
| `--redo` | yes/no | `no` | Re-download behavior. See "When to Use `redo`" section below for details. |
| `--batch` | INT | `None` | Process only one SRA record (1-based index). For HPC array jobs. |
| `--entrez_email` | email | `` | Email for NCBI Entrez API. |
| `--id` | STR | `None` | Single BioProject/BioSample/SRR ID to download directly. |
| `--id_list` | PATH | `None` | File containing list of SRA IDs (one per line). |
| `--layout` | STR | `auto` | Library layout: `single`, `paired`, or `auto` (prefers paired). |
| `--max_bp` | INT | `999999999999999` | Target sequence size (bp) to download per sample. Use to limit download size for testing or disk space constraints. After quantification, only small abundance files (~few MB) remain, so large downloads are temporary. |
| `--min_read_length` | INT | `25` | Minimum read length after quality filtering. |
| `--pfd` | yes/no | `yes` | Use parallel-fastq-dump for SRA extraction. |
| `--pfd_exe` | PATH | `parallel-fastq-dump` | Path to parallel-fastq-dump executable. |
| `--prefetch_exe` | PATH | `prefetch` | Path to prefetch executable. |
| `--fastp` | yes/no | `yes` | Run fastp for quality filtering and adapter trimming. |
| `--fastp_exe` | PATH | `fastp` | Path to fastp executable. |
| `--fastp_option` | STR | `-j /dev/null -h /dev/null` | Custom fastp options. |
| `--remove_sra` | yes/no | `yes` | Delete SRA files after FASTQ extraction. |
| `--remove_tmp` | yes/no | `yes` | Remove temporary files. |
| `--pfd_print` | yes/no | `yes` | Show parallel-fastq-dump stdout/stderr. |
| `--fastp_print` | yes/no | `yes` | Show fastp stdout/stderr. |
| `--sci_name` | STR | `None` | Species name filter (if BioProject spans multiple species). |
| `--ncbi` | yes/no | `yes` | Download from NCBI cloud. |
| `--aws` | yes/no | `yes` | Download from Amazon Web Services (AWS). |
| `--gcp` | yes/no | `yes` | Download from Google Cloud Platform (GCP). |
| `--read_name` | STR | `default` | Read name formatting: `default` or `trinity` (for Trinity assembler). |
| `--entrez_additional_search_term` | STR | `None` | Additional Entrez search terms to restrict SRA entries. |
| `--tol` | FLOAT | `1` | Acceptable percentage loss of reads relative to --max_bp. |

### METAINFORMANT-Specific Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `show_progress` | bool | `true` | Enable live progress tracking with file size monitoring. Shows per-thread download progress, rates (MB/s), and elapsed time. |
| `progress_update_interval` | float | `2.0` | How often to update progress display in seconds. Lower values provide more frequent updates but use more CPU. |
| `progress_style` | str | `"bar"` | Progress display style: `"bar"` for tqdm progress bars (requires tqdm), or `"text"` for simple text updates. |

## When to Use `redo: yes` vs `redo: no`

The `redo` parameter controls whether amalgkit re-downloads SRA files that already exist locally.

### Default Behavior: `redo: no` (Recommended)

**Use `redo: no` for normal operations.** This is the default and recommended setting.

With `redo: no`:
- ✅ **Skips already-downloaded SRA files** - Amalgkit automatically detects existing `.sra` files and skips re-downloading them
- ✅ **Saves bandwidth and time** - Only downloads missing files
- ✅ **Still processes existing files** - Even if SRA files exist, amalgkit will still:
  - Extract FASTQ files if they don't exist
  - Run quality filtering with fastp if enabled
  - Generate output files as needed

**Example output with `redo: no`:**
```
Processing SRA ID: SRR34065661
Previously-downloaded sra file was detected at: .../SRR34065661.sra
Time elapsed for 1st-round sequence extraction: SRR34065661, 0.0 sec
```

### When to Use `redo: yes`

**Use `redo: yes` only in specific situations:**

1. **Corrupted files** - If SRA files are corrupted and need to be re-downloaded
2. **Data refresh** - If you want to ensure you have the latest version from SRA (rarely needed, as SRA data is typically stable)
3. **Testing** - When testing download functionality or troubleshooting download issues
4. **Configuration changes** - If you changed download settings (e.g., switching from NCBI to AWS) and want fresh downloads

**Warning**: Using `redo: yes` will:
- ⚠️ Re-download ALL files, even if they already exist
- ⚠️ Use significant bandwidth and time
- ⚠️ May overwrite existing files

### File Detection Behavior

Amalgkit detects existing files by checking for `.sra` files in the expected location:
- Pattern: `{out_dir}/getfastq/{SRR_ID}/{SRR_ID}.sra`
- Detection message: `"Previously-downloaded sra file was detected at: {path}"`
- Skip message: `"Previously-downloaded sra file was not detected. New sra file will be downloaded."`

### Summary Logging

After the `getfastq` step completes, METAINFORMANT logs a summary showing:
- Number of files skipped (already exists)
- Number of files downloaded
- Total samples processed

Example summary:
```
Step getfastq summary: 1 file(s) skipped (already exists), 14 file(s) downloaded (total: 15 samples)
```

## Input Requirements

### Prerequisites

- **Metadata Table**: Row-per-sample format from `amalgkit select` (typically `metadata.filtered.tissue.tsv` or `metadata.tsv`). **Note**: Pivot tables (like `pivot_qualified.tsv`) are NOT supported - they lack the required `run` column.
- **SRA Toolkit**: `prefetch`, `fasterq-dump`, or `parallel-fastq-dump`
- **fastp** (optional but recommended): For quality filtering
- **Network Access**: For SRA downloads from NCBI/AWS/GCP
- **Disk Space**: Substantial (SRA files + FASTQs can be 5-50GB per sample)

### System Dependencies

| Tool | Purpose | Installation |
|------|---------|--------------|
| **prefetch** | Download SRA files | `conda install -c bioconda sra-tools` or `sudo apt-get install -y sra-toolkit` |
| **fasterq-dump** | Extract FASTQs from SRA | Included in sra-tools |
| **parallel-fastq-dump** | Parallel FASTQ extraction | `conda install -c bioconda parallel-fastq-dump` |
| **fastp** | Quality filtering | `conda install -c bioconda fastp` or `sudo apt-get install -y fastp` |

**Note**: For Python packages, use `uv pip install` (primary method). System tools can use conda or system package managers.

## Output Files

### Directory Structure

**Important**: The `getfastq` step automatically creates a `getfastq/` subdirectory within the specified `out_dir`. FASTQ files are NOT placed directly in `out_dir/`, but rather in `out_dir/getfastq/{sample_id}/`.

```
out_dir/
└── getfastq/                          # ← Automatically created by amalgkit
    └── SRR12345678/
        ├── SRR12345678_1.fastq.gz     # Paired-end read 1
        ├── SRR12345678_2.fastq.gz     # Paired-end read 2
        ├── SRR12345678.fastq.gz       # Single-end (if applicable)
        ├── fastp.json                 # Quality control report
        └── fastp.html                 # HTML QC report
```

**Example**: If `out_dir: output/amalgkit/pbarbatus/fastq`, then FASTQ files will be in:
- `output/amalgkit/pbarbatus/fastq/getfastq/SRR12345678/SRR12345678_1.fastq.gz`

**Key Point**: The `getfastq/` subdirectory is created automatically by amalgkit. When configuring the `integrate` step, the `fastq_dir` parameter must point to `{out_dir}/getfastq/` (the getfastq subdirectory), not just `{out_dir}/`.

### FASTQ File Naming

- **Paired-end**: `{SRR_ID}_1.fastq.gz` and `{SRR_ID}_2.fastq.gz`
- **Single-end**: `{SRR_ID}.fastq.gz`

### Quality Reports

Each sample generates:
- **fastp.json**: Machine-readable QC metrics
- **fastp.html**: Human-readable QC report with plots

## Workflow Integration

### Position in Pipeline

```mermaid
flowchart LR
    A[select] --> B[getfastq]
    B --> C[integrate]
    C --> D[quant]
```

**getfastq** runs **after select**, **before quant**.

### Downstream Dependencies

| Step | Dependency | Description |
|------|------------|-------------|
| `integrate` | FASTQs | Adds FASTQ info to metadata |
| `quant` | FASTQs | Quantifies gene expression |

## Download Sources

### Source Priority

When multiple sources enabled, amalgkit tries in order:
1. **AWS** (fastest, US regions)
2. **GCP** (fast, global)
3. **NCBI** (fallback, slower)

### Enabling/Disabling Sources

```bash
# Use all sources (default)
--ncbi yes --aws yes --gcp yes

# NCBI only (if cloud access blocked)
--ncbi yes --aws no --gcp no

# AWS only (fastest in US)
--ncbi no --aws yes --gcp no
```

## Performance Considerations

### Runtime

**Per Sample** (typical 5-10GB SRA):
- **Download**: 5-30 minutes (depends on network speed)
- **Extraction**: 5-15 minutes (with parallel-fastq-dump)
- **Quality filtering**: 2-10 minutes (with fastp)
- **Total**: 15-60 minutes per sample

**For 100 Samples**:
- **Serial**: 25-100 hours
- **Parallel** (10 concurrent): 3-10 hours

### Disk Space

**Per Sample**:
- **SRA file**: 2-20GB
- **Raw FASTQs**: 4-40GB (uncompressed during extraction)
- **Compressed FASTQs**: 1-10GB (final output)
- **Peak usage**: ~30-70GB per sample during processing

**Important**: Use `--remove_sra yes` to delete SRA files after extraction!

### Network Usage

- **Per sample**: 2-20GB download
- **100 samples**: 200GB - 2TB total
- **Recommendation**: Use high-speed network connection

### Download Size Limits

The `--max_bp` parameter allows you to limit the download size per sample. This is useful for:
- **Testing workflows**: Limit downloads to validate pipeline without full data
- **Disk space constraints**: Control peak disk usage during processing
- **Network bandwidth management**: Reduce download time and bandwidth usage

**Important Notes**:
- After quantification, only small abundance files remain (~few MB per sample)
- Large downloads are **temporary** - quantified samples don't need re-download
- The `max_bp` limit applies per sample, not total workflow size

**Example**:
```yaml
steps:
  getfastq:
    max_bp: 10000000000  # Limit to 10GB per sample
```

### Parallelization

```bash
# Good: Moderate parallelization
--threads 8

# Better: Higher parallelization (if CPU/network allows)
--threads 16

# Optimal: Match available cores
--threads $(nproc)
```

## Common Use Cases

### 1. Download All Selected Samples

```bash
amalgkit getfastq \
  --out_dir output/amalgkit/work/fastq \
  --metadata output/amalgkit/work/metadata/metadata.tsv \
  --threads 12 \
  --remove_sra yes
```

### 2. Download Specific SRA ID

```bash
# Direct SRA download (skip metadata step)
amalgkit getfastq \
  --out_dir output/amalgkit/work \
  --id SRR12345678 \
  --threads 8
```

### 3. Download from ID List

```bash
# Create ID list
echo "SRR12345678" > sra_ids.txt
echo "SRR12345679" >> sra_ids.txt
echo "SRR12345680" >> sra_ids.txt

# Download all
amalgkit getfastq \
  --out_dir output/amalgkit/work \
  --id_list sra_ids.txt \
  --threads 12
```

### 4. HPC Array Job Processing

```bash
# Submit array job (SLURM example)
sbatch --array=1-100 download_job.sh

# In download_job.sh:
amalgkit getfastq \
  --out_dir output/amalgkit/work/fastq \
  --metadata output/work/metadata/metadata.tsv \
  --batch ${SLURM_ARRAY_TASK_ID} \
  --threads 8
```

**Effect**: Each array job downloads one sample in parallel.

### 5. Re-download Failed Samples

```bash
# Force re-download
amalgkit getfastq \
  --out_dir output/amalgkit/work \
  --metadata failed_samples.tsv \
  --redo yes \
  --threads 8
```

## METAINFORMANT Enhanced Features

The METAINFORMANT wrapper provides additional functionality beyond stock amalgkit:

### Robust SRA Download with Retries

```python
# From src/metainformant/rna/steps/getfastq.py
# Automatic retry logic for failed downloads:
# 1. Try amalgkit getfastq (uses parallel-fastq-dump)
# 2. If fails, retry with prefetch + fasterq-dump
# 3. Up to 3 attempts per sample
```

### Cloud Acceleration (METAINFORMANT-specific)

Set `accelerate: true` in your config to automatically enable all cloud sources:

```yaml
steps:
  getfastq:
    accelerate: true  # Automatically sets aws, gcp, ncbi to 'yes'
```

This is processed by `src/metainformant/rna/steps/getfastq.py` before calling amalgkit:

```python
# Extract accelerate flag and apply to amalgkit params
accelerate_enabled = bool(effective_params.pop("accelerate", False))
if accelerate_enabled:
    effective_params.setdefault("aws", "yes")
    effective_params.setdefault("gcp", "yes")
    effective_params.setdefault("ncbi", "yes")
```

### Intelligent Defaults

METAINFORMANT sets robust defaults for production workflows:

```python
# From src/metainformant/rna/steps/getfastq.py::_inject_robust_defaults()
# PFD: Auto-enabled if parallel-fastq-dump is available, otherwise disabled
# fastp: Defaults to False (set to True in config to enable quality filtering)
# accelerate: Defaults to True (enables AWS/GCP/NCBI sources)
# pfd_print/fastp_print: Default to True for better diagnostics
# prefetch_exe/pfd_exe: Auto-detected from PATH if available
```

### Progress Monitoring

The METAINFORMANT wrapper includes:
- **Live progress tracking**: File size monitoring with per-thread download rates
- **Heartbeat messages**: Every 30 seconds during long-running downloads
- **Streaming logs**: Real-time output from amalgkit and helper tools

```python
from metainformant.rna.workflow import execute_workflow, load_workflow_config

cfg = load_workflow_config("config/amalgkit_amellifera.yaml")
execute_workflow(cfg)

# Output shows:
# [20:57:10] starting step 'getfastq' -> amalgkit getfastq ...
# [20:57:40] still running step 'getfastq' (pid=12345, heartbeat #1)...
# [20:58:10] still running step 'getfastq' (pid=12345, heartbeat #2)...
# Real-time download progress from DownloadProgressMonitor
```

### Targeted Missing Sample Recovery

```python
# Automatically detects missing FASTQs and downloads only those
# See: src/metainformant/rna/steps/getfastq.py::_retry_missing_srrs()
```

## Automatic SRA-to-FASTQ Conversion

The METAINFORMANT workflow includes automatic conversion of SRA files to FASTQ format:

1. **During getfastq step**: If `amalgkit getfastq` downloads SRA files but doesn't extract FASTQ files (e.g., due to extraction failure), the workflow automatically attempts conversion using `convert_sra_to_fastq()`.

2. **For existing SRA files**: If SRA files exist but FASTQ files are missing, the workflow detects this and automatically triggers extraction.

3. **Tool selection**: The conversion process:
   - First tries `parallel-fastq-dump` (if available) - works better with local files
   - Falls back to `fasterq-dump` with `--size-check off` to prevent disk limit errors
   - Automatically compresses output FASTQ files using `pigz` or `gzip`

4. **Wrapper mechanism**: For `amalgkit` calls, a wrapper script injects `--size-check off` into `fasterq-dump` commands to prevent "disk-limit exceeded" errors. The wrapper is automatically created in `{fastq_dir}/temp/fasterq-dump` and added to PATH for `amalgkit` processes.

## Validation

After the `getfastq` step completes, validation automatically runs to verify that samples were successfully downloaded and extracted.

### Automatic Validation

Validation runs automatically after `getfastq` completes:
- Checks that FASTQ files exist for each sample
- Validates file sizes and counts
- Generates validation report: `work_dir/validation/getfastq_validation.json`

### Validation Results

Check validation results:

```bash
# View validation report
cat output/amalgkit/my_species/work/validation/getfastq_validation.json

# Run standalone validation
python3 scripts/rna/run_workflow.py config/amalgkit/my_species.yaml --validate --validate-stage extraction
```

### Common Validation Issues

**Issue**: Samples show `extraction: false` in validation

**Solutions**:
1. Check if FASTQ files exist:
   ```bash
   ls output/amalgkit/my_species/fastq/getfastq/SRR123456/
   ```
2. Verify `getfastq` step completed successfully
3. Check logs for download errors:
   ```bash
   cat output/amalgkit/my_species/work/logs/getfastq.log
   ```
4. Re-run `getfastq` step if needed

For detailed validation troubleshooting, see [Validation Guide](../VALIDATION.md).

## Troubleshooting

### Issue: getfastq reports success but produces 0 reads

**Symptoms**:
- `amalgkit getfastq` returns exit code 0 (success)
- Log shows "Sum of fastp input reads: 0 bp" and "Sum of fastp output reads: 0 bp"

**Root Cause**: fasterq-dump is failing due to disk space issues in `/tmp` (tmpfs may be full).

**Solution**: The workflow now automatically sets `TMPDIR` to the repository's `.tmp/fasterq-dump` directory (on external drive with sufficient space) before running getfastq step. This prevents "disk-limit exceeded" and "storage exhausted" errors.

**Verification**: Check that TMPDIR is set correctly:
```bash
# During getfastq execution, TMPDIR should point to repository .tmp directory
echo $TMPDIR
# Should show: /path/to/repo/.tmp/fasterq-dump
```
- No FASTQ files are created
- SRA file may be deleted after "successful" download

**Root Causes**:

1. **LITE SRA Files (Most Common)**: Some species have LITE SRA files on NCBI/GCP that are metadata-only and contain no sequence data. These files have "lite" in their URLs (e.g., `SRR123456.lite.1`). AWS typically has full SRA files.

2. **Empty/Corrupted SRA Files**: Some SRA files in the database are empty or corrupted.

3. **Conversion Failure**: The SRA-to-FASTQ conversion process fails silently.

**Solution**:

**For LITE SRA Files**:
1. **Check metadata URLs**: Inspect the `NCBI_Link`, `AWS_Link`, and `GCP_Link` columns in your metadata. If URLs contain "lite", those sources have metadata-only files.
2. **Use AWS only**: Disable NCBI and GCP, use AWS only:
   ```yaml
   steps:
     getfastq:
       aws: yes
       gcp: no   # Disable if has lite files
       ncbi: no  # Disable if has lite files
   ```
3. **Verify AWS links**: AWS links typically don't have "lite" and contain full sequence data.

**For Other Cases**: The METAINFORMANT workflow includes enhanced validation that:
1. **Immediately validates FASTQ files** after getfastq completes
2. **Detects 0-read cases** by checking for actual FASTQ files
3. **Marks downloads as failed** if no files are found, preventing false positives
4. **Provides clear error messages** distinguishing between conversion failures and 0-read cases
5. **Automatically attempts conversion** for SRA files that exist but lack FASTQ files

**Error Messages**:
- `⚠️ SRA exists but no FASTQ (conversion may have failed or produced 0 reads)` - SRA file downloaded but conversion failed
- `⚠️ getfastq succeeded but no files found (may have produced 0 reads)` - No files at all, likely 0-read case or LITE file

**Manual Check**:
```bash
# Check if sample directory has files
ls -la output/amalgkit/work/fastq/getfastq/SRR123456/

# If empty or only has .sra file, the download produced 0 reads

# Check metadata for lite files
grep "SRR123456" output/amalgkit/work/metadata/metadata.tsv | grep -i "lite"
```

**Workflow Behavior**: The workflow automatically:
- Detects these cases immediately after getfastq completes
- Attempts automatic SRA-to-FASTQ conversion for existing SRA files
- Marks the sample as failed if conversion fails (doesn't wait indefinitely)
- Continues with other samples
- Reports failed samples in the final summary

### Issue: "disk-limit exceeded" during SRA extraction

**Symptoms**:
- `fasterq-dump` fails with "disk-limit exceeded!" error
- SRA file exists but extraction fails
- Even though disk has sufficient space

**Root Cause**: `fasterq-dump` has an internal disk space check that can fail even when sufficient space is available, especially if the default temporary directory (`/tmp`) is small or full. Additionally, `/tmp` may be a tmpfs (RAM-based filesystem) with limited space.

**Solutions Implemented**:

1. **Automatic TMPDIR Configuration** (Primary Fix):
   - **Location**: `src/metainformant/rna/amalgkit.py::run_amalgkit()`
   - **Purpose**: Automatically sets `TMPDIR`, `TEMP`, and `TMP` environment variables to repository's `.tmp/fasterq-dump` directory before running getfastq step
   - **Implementation**:
     - Detects when `subcommand == "getfastq"`
     - Sets `TMPDIR` to `{repo_root}/.tmp/fasterq-dump` (on external drive with sufficient space)
     - Creates directory if it doesn't exist
     - Passes environment to subprocess calls
   - **Benefits**:
     - Prevents "disk-limit exceeded" errors when `/tmp` is full
     - Uses external drive space instead of limited tmpfs
     - Works automatically without user configuration
   - **Code Location**: `src/metainformant/rna/amalgkit.py::run_amalgkit()` (lines 299-318 and 392-406)

2. **Automatic vdb-config Repository Path Configuration** (For Prefetch Downloads):
   - **Location**: `src/metainformant/rna/workflow.py::execute_workflow()`
   - **Purpose**: Configures vdb-config repository root to point to amalgkit fastq output directory
   - **Implementation**:
     - Detects when getfastq step is in workflow
     - Gets fastq output directory from getfastq step parameters
     - Sets vdb-config repository root to `{fastq_dir}/getfastq` (where amalgkit expects SRA files)
     - This ensures `prefetch` downloads SRA files to the correct location
   - **Benefits**:
     - Prevents "SRA file download failed" errors where prefetch downloads to wrong location
     - Ensures SRA files are in the location amalgkit expects them
     - Works automatically without user configuration
   - **Code Location**: `src/metainformant/rna/workflow.py::execute_workflow()` (lines 422-460)
   - **Note**: vdb-config repository root (for prefetch) is different from VDB_CONFIG environment variable (for cache/temp)

2. **Wrapper Script for `amalgkit` Calls** (Legacy/Alternative):
   - **Location**: `{fastq_dir}/temp/fasterq-dump`
   - **Purpose**: Automatically injects `--size-check off` into all `fasterq-dump` calls made by `amalgkit`
   - **Implementation**:
     - Created automatically in `process_samples.py` when setting up download workers
     - Injected into PATH for `amalgkit` processes via environment variable
     - Wrapper script:
       ```bash
       #!/bin/bash
       exec /usr/bin/fasterq-dump --size-check off "$@"
       ```
   - **Code Location**: `src/metainformant/rna/steps/process_samples.py::_download_worker()` (lines 259-270)

2. **Direct Binary Detection for Manual Conversion**:
   - **Location**: `src/metainformant/rna/steps/getfastq.py::convert_sra_to_fastq()`
   - **Purpose**: Ensures direct calls to `convert_sra_to_fastq()` use the real `fasterq-dump` binary, not the wrapper
   - **Implementation**:
     - Checks system locations first (`/usr/bin/fasterq-dump`, `/usr/local/bin/fasterq-dump`)
     - Filters out wrapper scripts (checks if "temp" is in the path)
     - Explicitly passes `--size-check off` to the real binary
   - **Code Location**: `src/metainformant/rna/steps/getfastq.py::convert_sra_to_fastq()` (lines 578-614)

3. **Automatic SRA-to-FASTQ Conversion**:
   - **Location**: `src/metainformant/rna/steps/getfastq.py::run()`
   - **Purpose**: Automatically converts SRA files to FASTQ when SRA files exist but FASTQ files are missing
   - **Implementation**:
     - Detects samples with SRA files but no FASTQ files
     - Automatically calls `convert_sra_to_fastq()` for each sample
     - Logs success/failure for each conversion attempt
   - **Code Location**: `src/metainformant/rna/steps/getfastq.py::run()` (lines 274-295)

4. **LITE SRA File Handling**:
   - **Location**: Configuration files and workflow validation
   - **Purpose**: Prevents downloading metadata-only LITE SRA files
   - **Implementation**:
     - Configuration option to disable NCBI/GCP sources that have LITE files
     - Use AWS only when LITE files are detected
     - Validation detects 0-read cases and marks downloads as failed
   - **Code Location**: 
     - Configuration: `config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml`
     - Validation: `src/metainformant/rna/steps/process_samples.py::_download_worker()` (lines 315-327)

**Workflow Integration**:

1. **Download Phase** (`_download_worker`):
   - Creates wrapper script if it doesn't exist
   - Injects wrapper into PATH for `amalgkit getfastq` calls
   - `amalgkit` uses wrapper, which adds `--size-check off` automatically
   - Validates FASTQ files exist after download

2. **Automatic Conversion** (`getfastq.run`):
   - Detects samples with SRA but no FASTQ
   - Calls `convert_sra_to_fastq()` for each sample
   - Uses real binary (not wrapper) with explicit `--size-check off`

3. **Quantification Phase**:
   - Proceeds once FASTQ files are available
   - Deletes FASTQ files after quantification (if configured)

**Manual Fix** (if needed):
```bash
# Create wrapper manually
mkdir -p output/amalgkit/work/fastq/temp
cat > output/amalgkit/work/fastq/temp/fasterq-dump << 'EOF'
#!/bin/bash
exec /usr/bin/fasterq-dump --size-check off "$@"
EOF
chmod +x output/amalgkit/work/fastq/temp/fasterq-dump

# Or configure vdb-config repository root manually (for prefetch downloads)
# This should point to where amalgkit expects SRA files: {fastq_dir}/getfastq
vdb-config -s /repository/user/main/public/root=/path/to/amalgkit/fastq/getfastq

# Note: The workflow automatically configures vdb-config repository root to the correct location.
# Manual configuration is rarely needed.
```

**Note**: The workflow handles this automatically - manual intervention is rarely needed.

**Testing**:

```bash
# Test wrapper
output/amalgkit/pogonomyrmex_barbatus/fastq/temp/fasterq-dump --help
# Should show help with --size-check off already applied

# Test direct conversion
python3 -c "
from pathlib import Path
from metainformant.rna.steps.getfastq import convert_sra_to_fastq
success, msg, files = convert_sra_to_fastq(
    'SRR123456',
    Path('path/to/SRR123456.sra'),
    Path('output/dir'),
    threads=4
)
print(f'Success: {success}, Files: {len(files)}')
"
```

**Expected Behavior**:
1. **Wrapper for amalgkit**: All `amalgkit getfastq` calls use wrapper → no "disk-limit exceeded" errors
2. **Direct conversion**: Uses real binary with explicit flag → no duplicate flag errors
3. **Automatic recovery**: SRA files without FASTQ are automatically converted
4. **LITE file detection**: 0-read cases are detected and marked as failed

### Issue: "prefetch: command not found"

**Solutions**:
```bash
# Install SRA toolkit (system tool)
conda install -c bioconda sra-tools
# Or using system package manager:
sudo apt-get install -y sra-toolkit

# Or specify full path
--prefetch_exe /path/to/prefetch

# Note: For Python packages, use uv pip install (primary method)
```

### Issue: Download very slow

**Diagnosis**:
```bash
# Check network speed
curl -o /dev/null http://speedtest.tele2.net/100MB.zip

# Check which source is being used (check logs)
```

**Solutions**:
1. Try different download source:
   ```bash
   # Disable slow sources
   --ncbi no --aws yes --gcp yes
   ```

2. Increase bandwidth:
   ```bash
   # If on institutional network with throttling
   # Contact IT to allowlist NCBI/AWS/GCP domains
   ```

3. Download during off-peak hours

### Issue: Disk space exhausted

```
OSError: [Errno 28] No space left on device
```

**Solutions**:
1. Enable SRA removal:
   ```bash
   --remove_sra yes  # Delete SRA after extraction
   ```

2. Use external storage:
   ```bash
   --out_dir /mnt/large_disk/amalgkit/work
   ```

3. Process in smaller batches:
   ```bash
   # Download 10 samples at a time
   head -11 metadata.tsv > batch1.tsv  # Include header + 10 samples
   tail -n +2 metadata.tsv | head -10 | cat <(head -1 metadata.tsv) - > batch1.tsv
   amalgkit getfastq --metadata batch1.tsv ...
   ```

4. Check and clean temporary files:
   ```bash
   # Remove old downloads
   rm -rf output/amalgkit/work/getfastq/*/SRR*.sra
   ```

### Issue: SRA extraction fails

```
parallel-fastq-dump: error reading SRA file
```

**Solutions**:
1. Retry with fasterq-dump:
   ```bash
   --pfd no  # Use fasterq-dump instead
   ```

2. Check SRA file integrity:
   ```bash
   vdb-validate output/work/getfastq/SRR*/SRR*.sra
   ```

3. Re-download with `--redo yes`:
   ```bash
   amalgkit getfastq --redo yes --id SRR12345678
   ```

### Issue: fastp quality filtering too aggressive

**Diagnosis**:
```bash
# Check fastp.json for filtering stats
cat output/work/getfastq/SRR12345678/fastp.json | grep "total_reads"
```

**Solutions**:
1. Lower minimum read length:
   ```bash
   --min_read_length 20  # Instead of 25
   ```

2. Disable fastp:
   ```bash
   --fastp no
   ```

3. Customize fastp options:
   ```bash
   --fastp_option "-q 15 -u 50 -j /dev/null -h /dev/null"
   ```

## Best Practices

### 1. Always Remove SRA Files

```bash
# Good: Clean up as you go
--remove_sra yes

# Bad: Wastes enormous disk space
--remove_sra no
```

### 2. Use Quality Filtering

```bash
# Good: Clean, high-quality reads
--fastp yes --min_read_length 25

# Less optimal: Raw reads with adapters/low quality
--fastp no
```

### 3. Enable All Download Sources

```bash
# Good: Fast, reliable downloads
--ncbi yes --aws yes --gcp yes

# Less optimal: Slow, limited to one source
--ncbi yes --aws no --gcp no
```

### 4. Monitor Disk Space

```bash
# Before starting
df -h output/

# During download (in another terminal)
watch -n 60 'df -h output/ && du -sh output/amalgkit/*/getfastq/'

# Set up automatic cleanup
--remove_sra yes --remove_tmp yes
```

### 5. Test with Small Subset First

```bash
# Test workflow with 1-2 samples
head -3 metadata.tsv > test.tsv
amalgkit getfastq --metadata test.tsv --threads 8

# If successful, run full dataset
amalgkit getfastq --metadata metadata.tsv --threads 12
```

## Real-World Examples

### Example 1: Apis mellifera (83 Brain Samples)

```bash
amalgkit getfastq \
  --out_dir output/amalgkit/amellifera/fastq \
  --metadata output/amalgkit/amellifera/work/metadata/metadata.tsv \
  --threads 16 \
  --pfd yes \
  --fastp yes \
  --remove_sra yes \
  --min_read_length 25
```

**Result**: Downloaded 83 samples, ~350GB total, 12 hours runtime

### Example 2: Pogonomyrmex barbatus (120 Samples)

```bash
amalgkit getfastq \
  --out_dir output/amalgkit/pogonomyrmex_barbatus/fastq \
  --metadata output/amalgkit/pogonomyrmex_barbatus/work/metadata/metadata.tsv \
  --threads 12 \
  --aws yes \
  --ncbi yes \
  --gcp yes \
  --remove_sra yes
```

**Result**: Multi-source downloads with automatic retry logic, 8 hours total

## Integration with METAINFORMANT Workflow

### Automatic Execution

```python
from metainformant.rna.workflow import execute_workflow, load_workflow_config

cfg = load_workflow_config("config/amalgkit_amellifera.yaml")
execute_workflow(cfg)  # getfastq runs automatically after select
```

### Enhanced Error Recovery

The METAINFORMANT workflow adds automatic retry logic:

```python
# From src/metainformant/rna/steps/getfastq.py
# Features:
# - 3 retries per sample
# - Fallback to prefetch+fasterq-dump if parallel-fastq-dump fails
# - Targeted re-download of missing samples
# - Comprehensive error logging
```

### Configuration

```yaml
steps:
  getfastq:
    threads: 16
    pfd: yes
    fastp: yes
    remove_sra: yes
    accelerate: yes  # METAINFORMANT-specific: enables retry logic
```

## References

- **SRA Toolkit Documentation**: https://github.com/ncbi/sra-tools
- **fastp**: https://github.com/OpenGene/fastp
- **NCBI Cloud Data**: https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud/
- **METAINFORMANT Workflow**: `docs/rna/workflow.md`

## See Also

- **Previous Step**: [`select.md`](select.md) - Selecting samples for download
- **Next Step**: [`integrate.md`](integrate.md) - Integrating FASTQ info into metadata
- **Next Step**: [`quant.md`](quant.md) - Quantifying gene expression
- **Workflow Overview**: [`../amalgkit.md`](../amalgkit.md)
- **Testing**: `tests/test_rna_amalgkit_steps.py::test_getfastq_basic_execution`

---

**Last Updated**: October 29, 2025  
**AMALGKIT Version**: 0.12.19  
**Status**: ✅ Production-ready, comprehensively tested, enhanced with METAINFORMANT retry logic

