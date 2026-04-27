# Agent Directives: scripts/cloud

## Role
Thin orchestrator scripts for GCP cloud deployment and data exfiltration.

## Contents
- `deploy_gcp.py` — VM provisioning and pipeline launch
- `prep_genomes.py` — Cloud-side genome preparation
- `cloud_startup.sh` — VM bootstrap script
- `vm_setup.sh` — Environment setup on the VM
- `download_results.sh` — Exfiltrate quantification results from cloud
- `install_gcloud.sh` — One-time gcloud CLI setup
- `categorize_errors.py` — Parse and classify pipeline errors
- `print_table.py` — Format status tables for reporting

## Rules
- Scripts are thin wrappers that call `metainformant.cloud` library code
- Follow NO MOCKING policy
- Use `uv` for Python dependencies
