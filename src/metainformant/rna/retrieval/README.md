# Retrieval

Download RNA-seq FASTQ files from ENA (European Nucleotide Archive) with automatic FTP/HTTP fallback, MD5 verification, and size-sorted batch processing.

## Contents

| File | Purpose |
|------|---------|
| `ena_downloader.py` | ENA sample lookup, FASTQ download with fallback, and batch orchestration |

## Key Functions and Classes

| Symbol | Description |
|--------|-------------|
| `SampleInfo` | Dataclass holding SRA ID, file URLs, sizes, and MD5 checksums |
| `get_ena_sample_info()` | Query the ENA API for sample metadata and download links |
| `get_ena_links()` | Retrieve FTP/HTTP download URLs for an SRA accession |
| `download_file()` | Download a single file with timeout and MD5 verification |
| `download_with_fallback()` | Try ENA FTP first, fall back to HTTP or SRA toolkit |
| `sort_samples_by_size()` | Order samples by file size for efficient batch scheduling |
| `download_sra_samples()` | Batch-download multiple SRA samples with parallel workers |
| `calculate_md5()` | Compute MD5 checksum of a local file for integrity checks |
| `clean_stagnant_file()` | Remove partially downloaded files that stopped progressing |

## Usage

```python
from metainformant.rna.retrieval.ena_downloader import (
    get_ena_sample_info,
    download_sra_samples,
    sort_samples_by_size,
)

info = get_ena_sample_info("SRR11196055")
sorted_samples = sort_samples_by_size(sra_ids, max_size_bytes=5_000_000_000)
download_sra_samples(sra_ids=["SRR11196055", "SRR11196056"], output_dir="output/fastq/")
```
