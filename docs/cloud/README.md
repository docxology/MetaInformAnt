# Cloud Deployment Documentation

This directory covers deploying MetaInformAnt pipelines to Google Cloud Platform.

## Contents

- **[DEPLOYMENT.md](DEPLOYMENT.md)** — Quick start guide, architecture, cost estimates
- **[TROUBLESHOOTING.md](TROUBLESHOOTING.md)** — Common issues and solutions
- **[ARCHITECTURE.md](ARCHITECTURE.md)** — Technical details of the cloud module

## Key Files

| File | Purpose |
|---|---|
| `src/metainformant/cloud/cloud_config.py` | VM configuration |
| `src/metainformant/cloud/gcp_deployer.py` | VM lifecycle management |
| `scripts/cloud/deploy_gcp.py` | User-facing CLI |
| `scripts/cloud/prep_genomes.py` | Genome download + index building |
| `Dockerfile` | Container image |
