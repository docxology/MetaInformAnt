# Specification: retrieval

## 🎯 Scope
RNA data retrieval modules.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 2 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `ena_downloader.py`

## Transfer bounds

`ena_transport.ena_request_slot` shares FIFO file and metadata limits across
downloader instances by HTTP origin. Native transfer and size-probe requests
share the file budget; request waits count against the caller's deadline.
Configured limits are validated before acquisition. Release occurs on success,
failure and exceptions, before any retry sleep. Limits are process-local and
do not coordinate independent producers; exclusive campaign locks prevent those.

`ENADownloader(file_workers=1|2)` bounds concurrent file transfers per run.
Distinct output filenames are required before starting transfers. All submitted
transfers settle before fallback; gzip integrity, retry deadlines, partial resume
and retained invalid evidence apply independently to each file.

## SSD acquisition

`download_staged_ena` is an internal protocol adapter used by the streaming
engine. Its reservation covers every advertised compressed input plus bounded
invalid witnesses. A destination/URL/size/checksum/layout identity and a file
lock isolate resumable workspaces. Existing bulk inputs are preserved. Staged
files are published through `publish_verified_file`, which refuses conflicts,
source mutation and nonregular inputs, copies to a destination-side temporary,
fsyncs, verifies SHA-256, then atomically links without replacement. Sources
are removed only after every mate is published. The canonical stat-based
validation witness is regenerated on the destination volume.

`AMALGKIT_PIPELINE_DOWNLOAD_FIRST_PASS_SECONDS` bounds the ENA first pass.
The terminal `deferred_download` result remains pending and bypasses immediate
NCBI recovery. Full retry phases disable this extra budget. Local input reuse
and SRA recovery retain their separately configured timeouts.
## Bounded metadata retries

`run_timeout_seconds` must not override configured `api_retries` to zero.
Metadata requests and retry sleeps share the earlier of an explicit absolute
deadline and the local run deadline. Expired explicit deadlines issue no HTTP
request; explicit zero retries remains authoritative. Backoff cannot extend
either deadline or retain a metadata request lease. Discovery failures retain
their diagnostic reason and never constitute accepted FASTQ data.

