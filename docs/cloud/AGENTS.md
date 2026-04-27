# Agent Directives: cloud

## Role

Cloud deployment and infrastructure automation. Provides API for provisioning Google Cloud Platform (GCP) VMs, installing Docker containers, uploading/downloading data, running METAINFORMANT pipelines, and collecting results.

**Primary consumers:** `scripts/cloud/deploy_gcp.py`, `scripts/cloud/prep_genomes.py`.

## Module Scope

| Submodule | Purpose |
|-----------|---------|
| `cloud_config.py` | `CloudConfig` dataclass — GCP project/zone/machine-type configuration, validation, cost estimation |
| `gcp_deployer.py` | `GCPDeployer` class — VM lifecycle (create/delete/wait), file transfer (SCP), Docker build/run/stream-logs |

## Key Source Files

- **CLI entry**: `scripts/cloud/deploy_gcp.py` — argparsing, command dispatch, CLI UX
- **Genome prep**: `scripts/cloud/prep_genomes.py` — NCBI download, indexing, GCS caching
- **Startup script**: `scripts/cloud/cloud_startup.sh` — VM init: Docker install, firewall rules, environment
- **Startup script**: `scripts/cloud/cloud_shutdown.sh` — cleanup, results compression
- **Library**: `src/metainformant/cloud/cloud_config.py` — 32-field `CloudConfig` dataclass
- **Library**: `src/metainformant/cloud/gcp_deployer.py` — ~1,000 lines GCP operations

## Related Documentation

- **Module guide**: [index.md](index.md) — Architecture, cost optimization, deployment walkthrough
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **API spec**: [SPEC.md](SPEC.md) — Type signatures, JSON schemas, error codes
- **Troubleshooting**: [TROUBLESHOOTING.md](TROUBLESHOOTING.md) — Common errors + fixes
- **Economics**: [ECONOMICS.md](ECONOMICS.md) — Cost breakdown, preemptible VM strategy
- **Architecture**: [../architecture.md](../architecture.md) — System-wide component diagram
- **RNA module**: [../rna/index.md](../rna/index.md) — Primary cloud consumer (amalgkit)

## Rules & Constraints

- GCP-native only (no AWS/Azure in this module yet)
- Thin wrapper over `gcloud` CLI — no client-library lock-in
- All VM operations idempotent where possible
- Preemptible VMs used by default for 70-80% cost savings
- Genome indices cached in `gs://metainformant-genomes/` (shared across projects)
