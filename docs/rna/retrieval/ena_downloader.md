# ENA Downloader

The `download_ena.py` script is the **primary download mechanism** for FASTQ files in METAINFORMANT. It downloads pre-compressed `.fastq.gz` files directly from the European Nucleotide Archive (ENA) using `wget`.

## Why ENA Instead of SRA Toolkit

| Factor | ENA Direct (wget) | SRA Toolkit (fasterq-dump) |
|--------|-------------------|---------------------------|
| File format | Pre-compressed `.fastq.gz` | SRA binary → extract → compress |
| Disk per sample | ~5 GB | ~25 GB (SRA cache + extracted) |
| Reliability | 100% (no LITE file issue) | ~60–80% (LITE SRA files produce 0 reads) |
| Resume support | `wget --continue` | Partial support |
| Concurrent workers | 16+ safe | Limited by SRA cache space |

## Usage

```bash
# Download a single SRR
python3 scripts/rna/download_ena.py --srr SRR14740514 --out-dir output/amalgkit/pbarbatus/fastq

# Batch download from a list
python3 scripts/rna/download_ena.py --srr-list srr_ids.txt --out-dir output/amalgkit/pbarbatus/fastq --workers 16
```

## Python API

```python
from metainformant.rna.retrieval import ena_downloader
from pathlib import Path

# Download a single file
success = ena_downloader.download_file(
    url="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/SRR123456/SRR123456_1.fastq.gz",
    dest_path=Path("output/SRR123456_1.fastq.gz"),
    expected_md5="d41d8cd98f00b204e9800998ecf8427e"
)

# Batch download with fallback
success_count, fail_count = ena_downloader.download_sra_samples(
    sra_ids=["SRR123456", "SRR123457"],
    base_out_dir=Path("output/amalgkit"),
    sort_by_size=True,
    use_fallback=True   # falls back to SRA Toolkit if ENA lacks the accession
)
```

## Integrity Verification

Every downloaded file is verified before being passed to kallisto:

1. **MD5 checksum** — compared against ENA-provided hash
2. **Gzip integrity** — `gzip -t` to catch truncated transfers

If either check fails:
1. Invalid file is deleted immediately
2. Download is retried (with exponential backoff)
3. If all retries fail, sample is marked failed and pipeline continues with the next sample

## Resumable Downloads

`wget --continue` is used by default. If a download is interrupted:
- Restarts from where it left off (no re-downloading from zero)
- Safe to re-run the script against partially downloaded files

## Fallback to SRA Toolkit

When an SRR accession has no ENA entry, the downloader falls back to `fasterq-dump`:

```bash
# Fallback is automatic when use_fallback=True (default)
# Requires SRA Toolkit on PATH:
fasterq-dump --version
```

The fallback is rare (< 5% of accessions in practice). ENA covers virtually all public SRA data.

## Output Structure

```
out_dir/
└── getfastq/
    └── <SRR>/
        ├── <SRR>_1.fastq.gz   # forward reads (paired-end)
        ├── <SRR>_2.fastq.gz   # reverse reads (paired-end)
        └── <SRR>.fastq.gz     # single-end (if applicable)
```

## See Also

- [04_getfastq.md](../amalgkit/steps/04_getfastq.md) — getfastq step documentation
- [amalgkit.md](../amalgkit/amalgkit.md) — Download architecture overview
