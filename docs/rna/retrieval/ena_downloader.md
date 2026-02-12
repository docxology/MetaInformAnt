# ENA Downloader Documentation

The `ena_downloader` module provides robust functionality for downloading FASTQ files from the European Nucleotide Archive (ENA) and other sources. It is designed to handle large-scale data retrieval with built-in validation and error handling.

## Key Features

- **Multi-Source Fallback**: Attempts downloads from ENA FTP, ENA HTTP, and falls back to NCBI SRA (via `fasterq-dump`) if needed.
- **Integrity Verification**:
  - **MD5 Checksum**: Verifies downloaded files against ENA-provided MD5 hashes.
  - **Gzip Integrity**: Strictly verifies `.gz` files using `gzip -t` to ensure no corruption occurred during transfer.
- **Resumable Downloads**: Uses `curl -C -` to resume interrupted downloads.
- **Smart Retries**: Implements exponential backoff for network transient errors.

## Usage

### Basic Download

```python
from metainformant.rna.retrieval import ena_downloader
from pathlib import Path

# Download a single file
success = ena_downloader.download_file(
    url="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/SRR123456/SRR123456_1.fastq.gz",
    dest_path=Path("output/SRR123456_1.fastq.gz"),
    expected_md5="d41d8cd98f00b204e9800998ecf8427e"
)
```

### Batch Processing

```python
from metainformant.rna.retrieval import ena_downloader

# Download multiple samples with fallback
sra_ids = ["SRR123456", "SRR123457"]
success_count, fail_count = ena_downloader.download_sra_samples(
    sra_ids=sra_ids,
    base_out_dir=Path("output/amalgkit"),
    sort_by_size=True,
    use_fallback=True
)
```

## Integrity Logic

The module enforces strict integrity checks. If a file fails MD5 or Gzip verification:

1. The invalid file is immediately deleted (`clean_stagnant_file`).
2. The download is retried (up to `retries` limit).
3. If all retries fail, the function returns `False`.

This ensures that no corrupted data enters the downstream pipeline.
