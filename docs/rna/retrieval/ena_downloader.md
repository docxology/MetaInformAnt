# ENA Downloader

MetaInformAnt provides the implemented ENA downloader as the
`ENADownloader` Python class in
`src/metainformant/rna/retrieval/ena_downloader.py`. There is no standalone
`scripts/rna/download_ena.py` entry point in this checkout.

## Python API

```python
from pathlib import Path

from metainformant.rna.retrieval.ena_downloader import ENADownloader

downloader = ENADownloader(timeout=1800, retries=3)
success, message, files = downloader.download_run(
    "SRR1234567",
    Path("output/amalgkit/pbarbatus/work/getfastq/SRR1234567"),
)
if not success:
    raise RuntimeError(message)
```

For a batch of accessions, use the implemented `download_sra_samples()`
helper. It writes the Amalgkit `getfastq/{run_id}/` layout and returns
`(success_count, failure_count)`. The higher-level
`StreamingPipelineOrchestrator` additionally provides the configured
`fasterq-dump` fallback when direct ENA retrieval fails.

## What is verified

The downloader queries the ENA Portal API, retrieves pre-compressed FASTQ
URLs, retries `curl` transfers, skips existing non-empty valid gzip files,
retains resumable `.part` files, and preserves rejected payloads with an
`.invalid` suffix. `verify_gzip_integrity()` is the transfer integrity gate
used by the downloader. `calculate_md5()` is available as a
utility, but an MD5 comparison is not automatically performed by
`download_run()` unless a caller supplies that policy separately.

## Output layout

```text
<work_dir>/getfastq/
└── <run_id>/
    ├── <run_id>_1.fastq.gz  # paired-end forward reads
    ├── <run_id>_2.fastq.gz  # paired-end reverse reads, when present
    └── <run_id>.fastq.gz     # single-end reads, when present
```

FASTQ files are temporary inputs. Successful quantification produces the
canonical `abundance.tsv` evidence and may remove FASTQ files according to the
active configuration.

## See also

- [04_getfastq.md](../amalgkit/steps/04_getfastq.md) — workflow integration
- [amalgkit.md](../amalgkit/amalgkit.md) — wrapper and monitoring overview
- [Streaming Orchestration](../ORCHESTRATION.md) — execution and fallback logic
