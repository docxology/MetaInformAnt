# Core Downloads (Progress + Heartbeat)

METAINFORMANT centralizes downloading logic in `metainformant.core.download` to make long-running downloads observable and robust.

## Features

- **Heartbeat files**: machine-readable JSON updated periodically during downloads.
- **Progress bars**: optional (uses `tqdm` if installed via `metainformant.core.progress`).
- **Retry + backoff**: retries with exponential backoff.
- **HTTP resume**: best-effort resume using HTTP `Range` requests (when supported).

## Heartbeat location

For a destination file `output/foo/bar.fastq.gz`, the heartbeat file is:

- `output/foo/.downloads/bar.fastq.gz.heartbeat.json`

## Python usage

```python
from metainformant.core.download import download_with_progress

result = download_with_progress(
    "https://example.com/file.fastq.gz",
    "output/rna/example/file.fastq.gz",
    heartbeat_interval=5,
    show_progress=True,
    max_retries=3,
    resume=True,
)

if not result.success:
    raise RuntimeError(result.error)
```

## Compatibility

`metainformant.core.io.download_file()` now delegates to `download_with_progress()` (with progress disabled by default) while preserving the same boolean return semantics.



