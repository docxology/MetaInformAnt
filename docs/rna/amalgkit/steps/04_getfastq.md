# amalgkit getfastq: FASTQ File Download

## Purpose

Downloads FASTQ files for RNA-seq quantification. In METAINFORMANT, the
implemented `ENADownloader` uses direct ENA HTTP retrieval as the primary
method; `StreamingPipelineOrchestrator` can fall back to SRA Toolkit when
configured and available.

## Download Method

### Primary: ENA direct retrieval (`ENADownloader`)

```python
from pathlib import Path

from metainformant.rna.retrieval.ena_downloader import ENADownloader

downloader = ENADownloader(timeout=1800, retries=3)
success, message, files = downloader.download_run(
    "<SRR_ID>", Path("output/amalgkit/<species>/work/getfastq/<SRR_ID>")
)
```

**ENA advantages:**
- Pre-compressed `.fastq.gz` files (~5 GB/sample vs ~25 GB for SRA cache)
- Existing valid gzip files are reused; corrupt or partial files are removed
- Supports 16+ concurrent workers without disk pressure
- Gzip integrity verification is built in; MD5 is an available utility, not an
  automatic comparison in `download_run()`

### Fallback: SRA Toolkit

Used automatically when an accession is absent from ENA (< 5% of cases):

```bash
# Requires sra-tools on PATH
conda install -c bioconda sra-tools
```

## Amalgkit 0.16.59 parameters

`amalgkit getfastq` is still invoked as part of the pipeline for integration with metadata and file management. Key parameters:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--out_dir` | PATH | `./` | Output directory. FASTQs land in `{out_dir}/getfastq/{SRR}/`. |
| `--metadata` | PATH | inferred | Row-per-sample metadata TSV (must have `run` column). |
| `--threads` | INT | `1` | Threads for parallel processing. |
| `--redo` | yes/no | `no` | Skip already-downloaded samples (recommended). |
| `--batch` | INT | None | Process one SRA record by 1-based index (for HPC array jobs). |
| `--max_bp` | INT | `999999999999999` | Limit target sequence extraction (downsamples or skips massive sequences). Use 10000000000 for a 10 GB cap. Default: Unlimited. |
| `--fastp` | yes/no | `yes` | Quality filtering and adapter trimming with fastp. |
| `--fastp_exe` / `--fastp_option` | PATH/text | current | Select the fastp executable and pass additional fastp arguments. |
| `--remove_sra` | yes/no | `yes` | Delete SRA files after FASTQ extraction. |
| `--remove_tmp` | yes/no | `yes` | Remove temporary extraction files after a successful sample. |
| `--fasterq_size_check` | yes/no | current | Validate available extraction space before `fasterq-dump`. |
| `--fasterq_disk_limit` / `--fasterq_disk_limit_tmp` | INT | optional | Bound extraction and temporary-disk use when required by the host. |
| `--aws` | yes/no | `yes` | Use AWS Open Data for SRA downloads (SRA fallback path). |
| `--ncbi` | yes/no | `yes` | Use NCBI directly (SRA fallback path). |
| `--gcp` | yes/no | `yes` | Use GCP (SRA fallback path). |
| `--ena` / `--ddbj` | yes/no | current | Enable the corresponding SRA source. |
| `--dump_print` | yes/no | `no` | Print extraction commands without running them. |

### Resumable rRNA filtering

Amalgkit 0.16.59 adds resumable rRNA filtering controls. They are optional
and are disabled in the Hymenoptera configurations unless a run record
documents the reference database, filter order, resource limits, and output
validation:

```yaml
steps:
  getfastq:
    rrna_filter: no
    # rrna_filter_sensitivity: 1.0
    # rrna_filter_max_seqs: 20
    # rrna_filter_chunk_spots: 5000000
    # rrna_filter_memory_limit: 32G
    # rrna_filter_jobs: 1
    # filter_order: fastp,rrna,contam
```

Do not describe a cohort as rRNA-filtered unless the run manifest records
these settings and the filtered FASTQ integrity checks pass. The resumable
implementation is useful for interrupted large-sample runs, but it does not
change the biological inclusion criteria.

### METAINFORMANT-Specific Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `accelerate` | bool | `true` | Enables all cloud sources for SRA fallback. |

## Configuration

```yaml
steps:
  getfastq:
    out_dir: output/amalgkit/<species>/fastq
    threads: 24
    # Project-wide concurrency is supplied at runtime:
    # --workers 4 --max-in-flight 8
    fastp: yes
    max_bp: 10000000000         # 10GB extraction cap / downsampling limit
    remove_sra: yes
```

## Output Structure

```
out_dir/
 getfastq/ # Created automatically by amalgkit
 SRR12345678/
 SRR12345678_1.fastq.gz # Paired-end read 1
 SRR12345678_2.fastq.gz # Paired-end read 2 (if paired)
 SRR12345678.fastq.gz # Single-end (if applicable)
```

**Note:** FASTQs are **temporary**. They are deleted immediately after successful kallisto quantification. The `abundance.tsv` file is the canonical output.

## Workflow Integration

```mermaid
flowchart LR
    A[select] --> B[getfastq]
    B --> C[integrate]
    C --> D[quant]
    D --> E[FASTQ deleted]
```

- `getfastq` runs **after select**, **before integrate/quant**
- Downstream steps consume FASTQs from `{out_dir}/getfastq/{SRR}/`

## Sample Size Limits & Downsampling

To prevent pipeline exhaustion on massive monolithic datasets, you can apply an upper limit using `max_bp`. For sequence datasets exceeding this threshold, `amalgkit` natively utilizes `fasterq-dump` and metadata filters to intelligently downsample or skip them, effectively capping the download bandwidth and final FASTQ extraction.

**Standard 10 GB Cap Framework**:
A maximum base pair configuration of 10 billion bases mathematically translates to an absolute download and extraction maximum of roughly **10 GB** per sample (yielding ~150M reads, which remains incredibly high-depth for RNA-seq quantification). 

```yaml
steps:
  getfastq:
    max_bp: 10000000000   # 10B bases limit (~10 GB max extraction cap)
```

**Customizing Limits for Your Resources**:
If you lack storage capacity and bandwidth, or conversely, if you have immense enterprise resources, you may scale this parameter up or down exactly as you wish. 

- **Aggressive Execution (Low Storage)**: `max_bp: 50000000` (50 Mbp max limit)
- **Standard Baseline**: `max_bp: 4000000000` (4 Gbp max limit)
- **High-Resource / Deep-Seq**: `max_bp: 999999999999999` (Unlimited Extraction)

Samples categorically skipped entirely (if bypasses trigger) will be logged as:

```
 SzSkip — SRR12345678 (12.2B bases, limit: 10.0B)
```

## Performance

**Per Sample** (typical 5–10 GB compressed):
- **Download**: 5–20 min (depends on network)
- **Quality filter (fastp)**: 2–5 min
- **Total**: ~10–30 min per sample

**For 100 Samples with 16 Workers**: ~1.5–3 hours total

**Disk per sample**: ~5 GB during download, 0 after quantification + cleanup.

## Troubleshooting

### Download Slow or Failing

```bash
# Check active acquisition workers
pgrep -af 'curl .*fastq|fasterq-dump' || true

# Check ENA API availability
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR14740514&result=read_run&fields=fastq_ftp" | head
```

### Sample Has No ENA Entry

The downloader automatically falls back to SRA Toolkit if ENA lacks the accession. Ensure sra-tools is installed:

```bash
fasterq-dump --version
prefetch --version
```

### Disk Space Exhausted

```bash
# Check partial/stuck downloads
find output/amalgkit -name "*.fastq.gz" -size +1G | head

# Remove partial downloads
find output/amalgkit -name "*.fastq.gz.part" -delete
find output/amalgkit -name "*.sra" -delete   # remove SRA cache files
```

### Samples Showing 0 Reads (LITE Files)

If a sample downloaded via SRA Toolkit shows 0 reads, it likely has a LITE version:
- NCBI and GCP sometimes serve metadata-only "LITE" SRA files
- ENA never has LITE files — ENA always has full sequence data
- Solution: use the implemented `ENADownloader` path through the workflow, or
  configure the SRA fallback explicitly and verify the resulting gzip files.

## See Also

- [retrieval/ena_downloader.md](../../retrieval/ena_downloader.md) — ENA download details
- [05_integrate.md](05_integrate.md) — integrate step
- [06_quant.md](06_quant.md) — quantification step
