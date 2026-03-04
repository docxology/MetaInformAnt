# Cloud Deployment Module

GCP deployment for running METAINFORMANT pipelines at scale.

## Overview

The cloud module provides a **thin orchestrator** for deploying the amalgkit RNA-seq pipeline to Google Cloud Platform. A local CLI (`deploy_gcp.py`) manages the full VM lifecycle via `gcloud` — no Python GCP SDK required.

## Guides

| Guide | Description |
|-------|-------------|
| [DEPLOYMENT.md](DEPLOYMENT.md) | Step-by-step deployment walkthrough |
| [ARCHITECTURE.md](ARCHITECTURE.md) | Technical design of the cloud module |
| [TROUBLESHOOTING.md](TROUBLESHOOTING.md) | Common issues and solutions |

## Quick Start

```bash
# Install gcloud CLI
bash scripts/cloud/install_gcloud.sh
gcloud auth login
gcloud config set project YOUR_PROJECT_ID

# Deploy pipeline to GCP
python scripts/cloud/deploy_gcp.py deploy --project YOUR_PROJECT_ID

# Monitor
python scripts/cloud/deploy_gcp.py status --project YOUR_PROJECT_ID

# Download results
python scripts/cloud/deploy_gcp.py download --project YOUR_PROJECT_ID

# Tear down
python scripts/cloud/deploy_gcp.py destroy --project YOUR_PROJECT_ID
```

## Defaults

| Parameter | Value |
|-----------|-------|
| Machine type | `n2-highcpu-96` (96 vCPUs) |
| Disk | 500 GB SSD |
| Pricing | Spot (~$0.80/hr) |
| Workers | 80 |
| Max sample size | 20 GB |

**Cost estimate**: ~$7–27 for full pipeline (4–8 hours).

## Source

| File | Purpose |
|------|---------|
| `src/metainformant/cloud/cloud_config.py` | VM configuration dataclass |
| `src/metainformant/cloud/gcp_deployer.py` | VM lifecycle management |
| `scripts/cloud/deploy_gcp.py` | User-facing CLI |
| `scripts/cloud/prep_genomes.py` | Genome download + index building |

## Related Documentation

- [LINUX_TRANSFER.md](../LINUX_TRANSFER.md) — Resuming pipeline on Linux
- [RNA workflows](../rna/index.md) — RNA-seq pipeline documentation
- [Architecture](../architecture.md) — System architecture
