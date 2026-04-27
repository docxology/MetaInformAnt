# Cloud Deployment Scripts

Thin orchestrator scripts for GCP cloud deployment and data exfiltration.

## Scripts

| Script | Purpose |
|--------|---------|
| `deploy_gcp.py` | Provision GCP VM and launch pipeline |
| `prep_genomes.py` | Prepare genomes on cloud VM |
| `cloud_startup.sh` | VM bootstrap (Docker, deps) |
| `vm_setup.sh` | Environment configuration |
| `download_results.sh` | Exfiltrate results to local repo |
| `install_gcloud.sh` | One-time gcloud CLI installation |
| `categorize_errors.py` | Classify pipeline error logs |
| `print_table.py` | Format status tables |

## Usage

```bash
# Deploy a GCP VM for RNA-seq quantification
python3 scripts/cloud/deploy_gcp.py --config config/cloud/gcp_config.yaml

# Download results after completion
bash scripts/cloud/download_results.sh <instance-name>
```

## Related

- [Cloud Module](../../src/metainformant/cloud/README.md) — Library code
- [Cloud Docs](../../docs/cloud/index.md) — User guide
