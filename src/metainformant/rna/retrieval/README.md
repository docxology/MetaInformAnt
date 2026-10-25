# Retrieval

Download RNA-seq FASTQ files from ENA (European Nucleotide Archive) with automatic URL discovery and retry support.

## Contents

| File | Purpose |
|------|---------|
| `ena_downloader.py` | `ENADownloader` class for ENA FASTQ discovery and download |
| `ena_staging.py` | Bounded SSD acquisition and serialized, verified publication |
| `ena_transport.py` | Shared FIFO file/metadata request budgets across downloader instances |

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

`AMALGKIT_ENA_MAX_CONNECTIONS` and `AMALGKIT_ENA_API_CONNECTIONS` bound native
file requests (including HEAD probes) and metadata requests per origin in this
process. Zero preserves unlimited legacy behavior; explicit limits accept 1–256.
Queue waits consume the absolute deadline and retry sleeps release the lease.
`last_metadata_error` distinguishes failed discovery from an empty ENA record;
transfer errors include the native curl exit code and HTTP status when available.

`run_timeout_seconds` sets a shared ENA discovery/transfer budget across mates
and retries; zero disables the additional first-pass budget. The streaming
worker reports `deferred_download` when this budget expires, keeps the sample
pending and retains partials. A later retry phase must set the budget to zero
to use the full transfer timeout and recovery path.

Set `AMALGKIT_LOCAL_ENA_SCRATCH_DIR` to opt into SSD acquisition with an explicit
`AMALGKIT_DATA_ROOT`. Staging shares the quantification/extraction scratch
reservation ledger, requires known compressed sizes, checks available ENA MD5
values and validates layout before publication. Persisted SSD inputs require
an ENA checksum for reuse. Publication uses one slot, fsync, SHA-256 readback
and an atomic no-overwrite link on the destination filesystem. Canonical
validation witnesses are recreated after promotion. Existing bulk partials,
unknown sizes and unnumbered paired inputs retain the bulk acquisition path.
Failures retain SSD sources; canonical conflicts are refused.

`AMALGKIT_PIPELINE_VALIDATE_FASTQ_RECORDS=yes` applies full record, sequence/
quality, and mate identity/order/count checks to bulk reuse/resume as well.
The bounded chunk parser retains at most one input block and incomplete record
per stream. Raw validation witnesses use exclusively owned temporary files;
unsafe, nonregular and oversized witness objects are refused. Set
`AMALGKIT_PRESERVE_RESUMABLE_PARTS=yes` to disable legacy age-based partial
eviction, including when the extra first-pass deadline is disabled for retry.
## Metadata retries inside bounded acquisition

A nonzero `run_timeout_seconds` retains the caller's configured `api_retries`.
Both metadata requests and retry backoff consume the same absolute deadline.
Direct `get_fastq_urls` calls derive a local deadline from that run allowance
when needed; an earlier explicit deadline always wins. Explicit zero retries
still makes only one request. Backoff holds no metadata connection lease.
The local launcher selects two metadata retries while retaining its two-request
metadata cap. Transient server errors can therefore recover within the first
pass without weakening response parsing or FASTQ integrity/record validation.

