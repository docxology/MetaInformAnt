# Agent Directives: cloud

## Role

Cloud deployment of the amalgkit RNA-seq pipeline on Google Cloud Platform (GCP). Provides `CloudConfig` and `GCPDeployer` for creating a Compute Engine VM with a startup script, checking its status and logs over SSH, downloading results and deleting the VM.

**Primary consumer:** `scripts/cloud/deploy_gcp.py`.

## Module Scope

| Submodule | Purpose |
|-----------|---------|
| `cloud_config.py` | `CloudConfig` dataclass — GCP project/zone/machine-type and pipeline settings, validation, instance metadata, startup script path |
| `gcp_deployer.py` | `GCPDeployer` class — VM lifecycle (create/stop/start/delete), status and log tailing over SSH, `wait_for_ssh`, result download, GCS sync |

## Key Source Files

- **CLI entry**: `scripts/cloud/deploy_gcp.py` — argparse subcommands (`deploy`, `status`, `logs`, `startup-log`, `download`, `stop`, `start`, `destroy`)
- **Genome prep**: `scripts/cloud/prep_genomes.py` — downloads reference transcriptomes with `wget` and builds kallisto indexes; standalone, it does not import `metainformant.cloud`
- **Startup script**: `scripts/cloud/cloud_startup.sh` — VM init: Docker install, NVMe RAID mount, repository clone, image build, container start, optional GCS sync cron job
- **Download script**: `scripts/cloud/download_results.sh` — copies results on the VM and downloads them with `gcloud compute scp`
- **Native setup**: `scripts/cloud/vm_setup.sh` — installs the tools without Docker and starts `scripts/rna/run_all_species.py`
- **Library**: `src/metainformant/cloud/cloud_config.py` — 19-field `CloudConfig` dataclass
- **Library**: `src/metainformant/cloud/gcp_deployer.py` — `GCPDeployer`, a wrapper over the `gcloud` CLI

## Related Documentation

- **Module guide**: [index.md](index.md) — Overview, defaults, links to the guides
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O); the cloud module itself does not import them
- **API spec**: [SPEC.md](SPEC.md) — Type signatures, failure behavior, instance metadata keys, CLI and script reference
- **Troubleshooting**: [TROUBLESHOOTING.md](TROUBLESHOOTING.md) — Issues met during the first deployment and their fixes
- **Economics**: [ECONOMICS.md](ECONOMICS.md) — Spot versus standard provisioning, disk sizing, cost-related pipeline settings
- **Architecture**: [../architecture.md](../architecture.md) — System-wide component diagram
- **RNA module**: [../rna/index.md](../rna/index.md) — The pipeline the VM runs (amalgkit)

## Rules & Constraints

- GCP-native only (no AWS/Azure in this module yet)
- Thin wrapper over the `gcloud` CLI — no client-library lock-in
- `create_vm` raises on an invalid config or a failing `gcloud` call; the other lifecycle methods return `False` (or a `NOT_FOUND` status) when `gcloud` exits non-zero
- Spot provisioning is the default (`spot=True`); `--no-spot` selects on-demand pricing
- The module has no genome cache; kallisto indexes are built on the VM by the streaming orchestrator or by `prep_genomes.py`
- `download_results.sh` takes the project from `GCP_PROJECT` and falls back to the hard-coded project `cryptoptera` when it is unset; export `GCP_PROJECT`, `GCP_ZONE` and `GCP_VM` before downloading
