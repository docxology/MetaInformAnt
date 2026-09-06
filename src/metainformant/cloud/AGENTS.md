# Agent Directives: cloud

**Context**: Cloud deployment infrastructure for METAINFORMANT pipelines.

## Capabilities

GCP Compute Engine VM lifecycle management for running large-scale amalgkit RNA-seq pipelines: VM creation/teardown, SSH-based pipeline monitoring, GCS result sync, and result download (`scripts/cloud/cloud_startup.sh` and `scripts/cloud/download_results.sh` are invoked by the deployer; Docker builds and genome prep happen inside those scripts/on the VM).

## Subpackages

| File | Key Classes / Functions |
|------|------------------------|
| `cloud_config.py` | `CloudConfig` — dataclass for GCP project, zone, machine type, disk, Docker image |
| `gcp_deployer.py` | `GCPDeployer` — VM creation, monitoring, teardown via `gcloud` CLI |

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Follow REAL IMPLEMENTATION policy — all tests must use real implementations
- Use `uv` for dependency management
- Shell out to `gcloud` CLI; no Google Cloud Python SDK dependency
