# Retrieval

Download RNA-seq FASTQ files from ENA (European Nucleotide Archive) with automatic URL discovery and retry support.

## Contents

| File | Purpose |
|------|---------|
| `ena_downloader.py` | `ENADownloader` class for ENA FASTQ discovery and download |

## Key Classes

### ENADownloader

| Method | Description |
|--------|-------------|
| `__init__(timeout, retries)` | Configure download timeouts and retry count |
| `get_fastq_urls(sample_id)` | Query ENA Portal API for FASTQ download URLs |
| `download_run(sample_id, output_dir)` | Download all FASTQ files for an SRA run via curl |

## Usage

```python
from pathlib import Path
from metainformant.rna.retrieval.ena_downloader import ENADownloader

downloader = ENADownloader(timeout=1800, retries=3, file_workers=2)

# Discover FASTQ URLs
urls = downloader.get_fastq_urls("SRR11196055")

# Download run
success, message, files = downloader.download_run("SRR11196055", Path("output/fastq/"))
```

`file_workers` accepts 1 (default) or 2. Two overlaps distinct mate transfers.
Every mate passes gzip integrity checks; completed mates and resumable partials
are retained after failure. Reuse preserves validated file bytes and timestamps.
