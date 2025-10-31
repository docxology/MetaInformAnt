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
  --out_dir output/amalgkit/work \
  --metadata output/amalgkit/work/metadata/pivot_qualified.tsv \
  --threads 8 \
  --pfd yes \
  --fastp yes
```

### Python API

```python
from metainformant.rna import amalgkit

result = amalgkit.getfastq(
    out_dir="output/amalgkit/work",
    metadata="output/amalgkit/work/metadata/pivot_qualified.tsv",
    threads=8,
    pfd=True,
    fastp=True
)
```

### Configuration File

```yaml
steps:
  getfastq:
    out_dir: output/amalgkit/amellifera/work
    metadata: output/amalgkit/amellifera/work/metadata/pivot_qualified.tsv
    threads: 8
    pfd: yes                # Use parallel-fastq-dump
    fastp: yes              # Quality filtering
    remove_sra: yes         # Delete SRA files after extraction
    accelerate: true        # Enable cloud acceleration (METAINFORMANT-specific)
    ncbi: yes
    aws: yes
    gcp: yes
    max_size: "50GB"        # Handle large samples
    min_size: "1MB"         # Skip tiny files
    min_read_length: 25
```

## Parameters

### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--out_dir` | PATH | `./` | Directory for intermediate and output files. |
| `--metadata` | PATH | `inferred` | Path to metadata.tsv. Default: `out_dir/metadata/metadata.tsv` |
| `--threads` | INT | `1` | Number of threads for parallel processing. |
| `--redo` | yes/no | `no` | Re-download even if FASTQs exist. |
| `--batch` | INT | `None` | Process only one SRA record (1-based index). For HPC array jobs. |
| `--entrez_email` | email | `` | Email for NCBI Entrez API. |
| `--id` | STR | `None` | Single BioProject/BioSample/SRR ID to download directly. |
| `--id_list` | PATH | `None` | File containing list of SRA IDs (one per line). |
| `--layout` | STR | `auto` | Library layout: `single`, `paired`, or `auto` (prefers paired). |
| `--max_bp` | INT | `999999999999999` | Target sequence size (bp) to download. |
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

## Input Requirements

### Prerequisites

- **Metadata Table**: Selected samples from `amalgkit select` (typically `pivot_qualified.tsv` or `pivot_selected.tsv`)
- **SRA Toolkit**: `prefetch`, `fasterq-dump`, or `parallel-fastq-dump`
- **fastp** (optional but recommended): For quality filtering
- **Network Access**: For SRA downloads from NCBI/AWS/GCP
- **Disk Space**: Substantial (SRA files + FASTQs can be 5-50GB per sample)

### System Dependencies

| Tool | Purpose | Installation |
|------|---------|--------------|
| **prefetch** | Download SRA files | `conda install -c bioconda sra-tools` |
| **fasterq-dump** | Extract FASTQs from SRA | Included in sra-tools |
| **parallel-fastq-dump** | Parallel FASTQ extraction | `conda install -c bioconda parallel-fastq-dump` |
| **fastp** | Quality filtering | `conda install -c bioconda fastp` |

## Output Files

### Directory Structure

```
out_dir/getfastq/
└── SRR12345678/
    ├── SRR12345678_1.fastq.gz         # Paired-end read 1
    ├── SRR12345678_2.fastq.gz         # Paired-end read 2
    ├── SRR12345678.fastq.gz           # Single-end (if applicable)
    ├── fastp.json                     # Quality control report
    └── fastp.html                     # HTML QC report
```

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
  --out_dir output/amalgkit/work \
  --metadata output/amalgkit/work/metadata/pivot_selected.tsv \
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
  --out_dir output/amalgkit/work \
  --metadata output/work/metadata/pivot_selected.tsv \
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
effective_params.setdefault("pfd", True)           # Use parallel-fastq-dump when available
effective_params.setdefault("max_size", "50GB")    # Handle large samples
effective_params.setdefault("min_size", "1MB")     # Skip tiny/corrupted files
effective_params.setdefault("fastp", "yes")        # Quality filtering enabled
effective_params.setdefault("remove_sra", "yes")   # Delete SRA after extraction
```

### Progress Monitoring

```python
from metainformant.rna.workflow import execute_workflow

cfg = load_workflow_config("config/amalgkit_amellifera.yaml")
execute_workflow(cfg, stream_logs=True)

# Output shows real-time progress:
# [INFO] Downloading SRR12345678... (1/100)
# [INFO] Download speed: 15MB/s
# [INFO] Time remaining: 2h 15m
```

### Targeted Missing Sample Recovery

```python
# Automatically detects missing FASTQs and downloads only those
# See: src/metainformant/rna/steps/getfastq.py::_retry_missing_srrs()
```

## Troubleshooting

### Issue: "prefetch: command not found"

**Solutions**:
```bash
# Install SRA toolkit
conda install -c bioconda sra-tools

# Or specify full path
--prefetch_exe /path/to/prefetch
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
   head -10 pivot_selected.tsv > batch1.tsv
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
head -3 pivot_selected.tsv > test.tsv
amalgkit getfastq --metadata test.tsv --threads 8

# If successful, run full dataset
amalgkit getfastq --metadata pivot_selected.tsv --threads 12
```

## Real-World Examples

### Example 1: Apis mellifera (83 Brain Samples)

```bash
amalgkit getfastq \
  --out_dir output/amalgkit/amellifera/work \
  --metadata output/amalgkit/amellifera/work/metadata/pivot_qualified.tsv \
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
  --out_dir output/amalgkit/pbarbatus/work \
  --metadata output/amalgkit/pbarbatus/work/metadata/pivot_qualified.tsv \
  --threads 12 \
  --aws yes \
  --ncbi no \
  --gcp no \
  --remove_sra yes
```

**Result**: AWS-only downloads, faster than NCBI, 8 hours total

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

