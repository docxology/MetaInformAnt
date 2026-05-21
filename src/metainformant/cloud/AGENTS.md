# Agent Directives: cloud

**Context**: Cloud deployment infrastructure for METAINFORMANT pipelines.

## Capabilities

GCP Compute Engine VM lifecycle management for running large-scale amalgkit RNA-seq pipelines: VM creation/teardown, Docker container deployment, genome preparation orchestration, and result exfiltration.

## Subpackages

| File | Key Classes / Functions |
|------|------------------------|
| `cloud_config.py` | `CloudConfig` — dataclass for GCP project, zone, machine type, disk, Docker image |
| `gcp_deployer.py` | `GCPDeployer` — VM creation, monitoring, teardown via `gcloud` CLI |

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Follow NO MOCKING policy — all tests must use real implementations
- Use `uv` for dependency management
- Shell out to `gcloud` CLI; no Google Cloud Python SDK dependency
