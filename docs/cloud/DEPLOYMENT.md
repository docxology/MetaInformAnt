# GCP Cloud Deployment Guide

Deploy the MetaInformAnt amalgkit RNA-seq pipeline to Google Cloud Platform for massively parallel processing.

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│  Local Machine                                              │
│  ┌──────────────────┐    ┌────────────────────────────┐     │
│  │ deploy_gcp.py    │───▸│ gcloud CLI                 │     │
│  │ (orchestrator)   │    │ (SSH, SCP, VM lifecycle)    │     │
│  └──────────────────┘    └────────────────────────────┘     │
│                                   │                         │
└───────────────────────────────────│─────────────────────────┘
                                    │ SSH/SCP
                                    ▼
┌─────────────────────────────────────────────────────────────┐
│  GCP VM (n2-highcpu-32, spot)                               │
│  ┌──────────────────┐    ┌────────────────────────────┐     │
│  │ prep_genomes.py  │───▸│ kallisto index (per species)│     │
│  └──────────────────┘    └────────────────────────────┘     │
│  ┌──────────────────┐    ┌────────────────────────────┐     │
│  │ run_all_species  │───▸│ streaming_orchestrator.py   │     │
│  │ (pipeline entry) │    │ (download + quant + merge)  │     │
│  └──────────────────┘    └────────────────────────────┘     │
│                                   │                         │
│                                   ▼                         │
│                          output/amalgkit/                    │
│                          (quant results, DB)                 │
└─────────────────────────────────────────────────────────────┘
```

## Quick Start

### 1. Install gcloud CLI

```bash
bash scripts/cloud/install_gcloud.sh
# Or manually:
# curl -fsSLO https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-linux-x86_64.tar.gz
# tar xzf google-cloud-cli-linux-x86_64.tar.gz && ./google-cloud-sdk/install.sh --quiet
```

### 2. Authenticate

```bash
gcloud auth login
gcloud config set project YOUR_PROJECT_ID
```

### 3. Deploy

```bash
# Deploy with defaults (n2-highcpu-32, spot, 24 workers, 32 threads)
python scripts/cloud/deploy_gcp.py deploy --project YOUR_PROJECT_ID

# Or customize:
python scripts/cloud/deploy_gcp.py deploy \
    --project YOUR_PROJECT_ID \
    --machine-type n2-highcpu-96 \
    --workers 80 --threads 96 \
    --max-gb 20.0 \
    --gcs-bucket my-results-bucket
```

### 4. Monitor

```bash
python scripts/cloud/deploy_gcp.py status  --project YOUR_PROJECT_ID
python scripts/cloud/deploy_gcp.py logs    --project YOUR_PROJECT_ID
```

### 5. Download Results

```bash
python scripts/cloud/deploy_gcp.py download --project YOUR_PROJECT_ID
```

### 6. Cleanup

```bash
python scripts/cloud/deploy_gcp.py destroy --project YOUR_PROJECT_ID
```

## Cost Estimates

| Machine Type | vCPUs | RAM | Spot/hr | On-demand/hr | Est. Total (8hr) |
|---|---|---|---|---|---|
| n2-highcpu-32 | 32 | 32 GB | ~$0.30 | ~$1.10 | $2-9 |
| n2-highcpu-64 | 64 | 64 GB | ~$0.60 | ~$2.20 | $5-18 |
| n2-highcpu-96 | 96 | 96 GB | ~$0.80 | ~$3.40 | $7-27 |

> **Note:** Default vCPU quota is 32 per project. Request a quota increase via the [GCP Console](https://console.cloud.google.com/iam-admin/quotas) for larger machines.

## File Structure

```
src/metainformant/cloud/
├── __init__.py           # Package exports
├── cloud_config.py       # VM config dataclass (machine, spot, workers...)
└── gcp_deployer.py       # VM lifecycle (create, SSH, status, download, delete)

scripts/cloud/
├── deploy_gcp.py         # CLI orchestrator (deploy/status/logs/download/destroy)
├── install_gcloud.sh     # One-click gcloud CLI installer
├── cloud_startup.sh      # VM boot script (Docker + pipeline auto-start)
├── prep_genomes.py       # Download transcriptomes + build kallisto indices
└── vm_setup.sh           # Native tool install (no Docker)

Dockerfile                # Full pipeline container image
.dockerignore             # Keeps image small
```

## Manual VM Setup (Without Docker)

If Docker builds are slow or unnecessary, install tools natively:

```bash
# SSH into VM
gcloud compute ssh metainformant-pipeline --zone us-central1-a

# Install tools
sudo apt-get install -y python3-pip r-base
sudo pip3 install --break-system-packages -e . amalgkit
wget -qO /usr/local/bin/fastp https://github.com/OpenGene/fastp/releases/download/v0.24.0/fastp && chmod +x /usr/local/bin/fastp

# Build genome indices (must run before pipeline)
python3 scripts/cloud/prep_genomes.py --threads 8

# Start pipeline
nohup python3 scripts/rna/run_all_species.py --max-gb 20.0 --workers 24 --threads 32 > output/amalgkit/pipeline.log 2>&1 &
```
