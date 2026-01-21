# Amalgkit Monitoring (Heartbeats + Progress Bars)

METAINFORMANT adds **heartbeat JSON files** and **optional progress bars** for long-running `amalgkit` steps executed via `metainformant.rna.amalgkit`.

This is designed for workflows that can run hours/days (notably `getfastq`, `quant`).

## What you get

- **Heartbeat JSON** written periodically so you can see if a step is alive and how far it’s progressed.
- **Progress bars** (tqdm) when an estimate of total work is available.

## Where heartbeats are written

For any amalgkit step with `out_dir=<DIR>`, heartbeats are written to:

- `<DIR>/.downloads/amalgkit-<step>.heartbeat.json`

Examples:
- `output/amalgkit/pbarbatus_test5/fastq/.downloads/amalgkit-getfastq.heartbeat.json`
- `output/amalgkit/pbarbatus_test5/quant/.downloads/amalgkit-quant.heartbeat.json`

## Per-step progress strategies

- **`getfastq`**: directory-size growth of `<out_dir>/getfastq/` (and total bytes estimated from metadata `size` column when available).
- **`metadata`**: file-count progress for key outputs under `<out_dir>/metadata/` (`metadata.tsv`, `metadata_original.tsv`, `pivot_selected.tsv`, `pivot_qualified.tsv`).
- **`quant`**: sample-count progress by counting `abundance.tsv` under `<out_dir>/quant/*/abundance.tsv`.
- **`merge`, `cstmm`, `curate`, `csca`**: file-count best-effort for key outputs; falls back to directory growth when needed.

## Heartbeat schema

Heartbeat files include fields like:

- `status`: `starting`, `running`, `completed`, `failed`
- `bytes_downloaded`: total bytes observed under the watched directory
- `progress`: a structured object when an estimate is available:
  - `type`: `directory_size` | `file_count` | `sample_count`
  - `current`, `total`, `percent`

## Controlling monitoring

All amalgkit wrappers accept:

- `monitor`: bool (default True for long-running steps)
- `show_progress`: bool (enable/disable progress bar)
- `heartbeat_interval`: seconds between heartbeat writes

Example:

```python
from metainformant.rna.amalgkit import getfastq

getfastq(
    {
        "out_dir": "output/amalgkit/demo/fastq",
        "metadata": "output/amalgkit/demo/work/metadata/metadata_selected.tsv",
        "threads": 8,
        "aws": "yes",
    },
    monitor=True,
    show_progress=True,
    heartbeat_interval=5,
)
```



