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
from metainformant.rna.retrieval.ena_downloader import ENADownloader

downloader = ENADownloader(timeout=1800, retries=3)

# Discover FASTQ URLs
urls = downloader.get_fastq_urls("SRR11196055")

# Download run
success, message, files = downloader.download_run("SRR11196055", Path("output/fastq/"))
```
