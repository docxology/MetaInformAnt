# amalgkit getfastq: FASTQ File Download

## Purpose

Downloads FASTQ files for RNA-seq quantification. In METAINFORMANT, this step uses **ENA direct wget** as the primary download method, with SRA Toolkit as a fallback.

## Download Method

### Primary: ENA Direct wget (`download_ena.py`)

```bash
# How METAINFORMANT downloads FASTQs for each sample
python3 scripts/rna/download_ena.py --srr <SRR_ID> --out-dir output/amalgkit/<species>/fastq
```

**ENA advantages:**
- Pre-compressed `.fastq.gz` files (~5 GB/sample vs ~25 GB for SRA cache)
- 100% reliable — no LITE file issue (NCBI/GCP LITE files have zero sequence data)
- Resumable with `wget --continue`
- Supports 16+ concurrent workers without disk pressure
- MD5 + gzip integrity verification built-in

### Fallback: SRA Toolkit

Used automatically when an accession is absent from ENA (< 5% of cases):

```bash
# Requires sra-tools on PATH
conda install -c bioconda sra-tools
```

## Amalgkit getfastq Parameters

`amalgkit getfastq` is still invoked as part of the pipeline for integration with metadata and file management. Key parameters:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--out_dir` | PATH | `./` | Output directory. FASTQs land in `{out_dir}/getfastq/{SRR}/`. |
| `--metadata` | PATH | inferred | Row-per-sample metadata TSV (must have `run` column). |
| `--threads` | INT | `1` | Threads for parallel processing. |
| `--redo` | yes/no | `no` | Skip already-downloaded samples (recommended). |
| `--batch` | INT | None | Process one SRA record by 1-based index (for HPC array jobs). |
| `--max_bp` | INT | `4000000000` | Skip samples exceeding this base count. Default: 4B. |
| `--fastp` | yes/no | `yes` | Quality filtering and adapter trimming with fastp. |
| `--remove_sra` | yes/no | `yes` | Delete SRA files after FASTQ extraction. |
| `--aws` | yes/no | `yes` | Use AWS Open Data for SRA downloads (SRA fallback path). |
| `--ncbi` | yes/no | `yes` | Use NCBI directly (SRA fallback path). |
| `--gcp` | yes/no | `yes` | Use GCP (SRA fallback path). |

### METAINFORMANT-Specific Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `num_download_workers` | INT | `16` | Parallel ENA download workers (not an amalgkit param; processed by METAINFORMANT). |
| `accelerate` | bool | `true` | Enables all cloud sources for SRA fallback. |

## Configuration

```yaml
steps:
  getfastq:
    out_dir: output/amalgkit/<species>/fastq
    threads: 24
    num_download_workers: 16    # Parallel ENA workers
    fastp: yes
    max_bp: 4000000000          # Skip samples >4B bases
    remove_sra: yes
```

## Output Structure

```
out_dir/
└── getfastq/                   # Created automatically by amalgkit
    └── SRR12345678/
        ├── SRR12345678_1.fastq.gz   # Paired-end read 1
        ├── SRR12345678_2.fastq.gz   # Paired-end read 2 (if paired)
        └── SRR12345678.fastq.gz     # Single-end (if applicable)
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

## Size Filtering

Samples exceeding `max_bp` (default: 4,000,000,000 = 4B bases) are auto-skipped:

```
✗ SzSkip — SRR12345678 (7.2B bases, limit: 4.0B)
```

Skipped samples are logged with `SzSkip` in the run log. To process large samples, increase `max_bp`:

```yaml
steps:
  getfastq:
    max_bp: 10000000000   # 10B bases
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
# Check active wget workers
ps aux | grep wget | grep -v grep

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
- Solution: use `download_ena.py` directly, or set `aws: yes, ncbi: no, gcp: no` in config

## See Also

- [retrieval/ena_downloader.md](../../retrieval/ena_downloader.md) — ENA download details
- [05_integrate.md](05_integrate.md) — integrate step
- [06_quant.md](06_quant.md) — quantification step
