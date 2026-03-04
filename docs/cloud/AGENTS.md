# Agent Directives: cloud

## Role

Cloud deployment module documentation agent for METAINFORMANT.

## Scope

- `src/metainformant/cloud/` — Cloud configuration and GCP deployer
- `scripts/cloud/` — Deployment CLI, startup scripts, genome prep
- `docs/cloud/` — User-facing deployment and troubleshooting guides

## Key Components

| File | Purpose |
|------|---------|
| `cloud_config.py` | `CloudConfig` dataclass — VM type, disk, workers, spot pricing |
| `gcp_deployer.py` | VM lifecycle via `subprocess` → `gcloud` CLI (no Python SDK) |
| `deploy_gcp.py` | User-facing CLI: deploy / status / logs / download / destroy |
| `prep_genomes.py` | Download genomes + build Kallisto indices on VM |
| `cloud_startup.sh` | VM boot script — Docker install, repo clone, pipeline auto-start |

## Standards

- **Thin orchestrator pattern**: local script shells out to `gcloud`; all logic in `src/`
- **No Python GCP SDK**: uses `subprocess.run(["gcloud", ...])` exclusively
- **Real implementations only** — NO_MOCKING_POLICY applies
- **Package management**: `uv` for all Python operations
