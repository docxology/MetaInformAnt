# Cloud Deployment Quick Reference

Deploy METAINFORMANT pipelines to Google Cloud Platform (GCP) with automated VM lifecycle, Docker orchestration, and cost-optimized preemptible instances.

## When to Use

Use `deploy_cloud` for large-scale compute-intensive workflows (28+ species RNA-seq, 10k+ sample GWAS, genome assemblies) where local resources are insufficient—not for small test runs (<1 hour) which are faster locally.

## Table of Contents

- [Prerequisites Checklist](#prerequisites-checklist)
- [Cost Estimates](#cost-estimates)
- [Deploy Command Reference](#deploy-command-reference)
- [Flags](#flags)
- [One-Liner: Full Pipeline](#one-liner-full-pipeline)
- [Cost Optimization Cheatsheet](#cost-optimization-cheatsheet)
- [Expected Output](#expected-output)
- [Common Pitfalls](#common-pitfalls)

---

- [ ] GCP project created + billing enabled
- [ ] `gcloud` installed and authenticated (`gcloud auth login`)
- [ ] Project set as default (`gcloud config set project PROJECT_ID`)
- [ ] Required APIs enabled:
 - Compute Engine API
 - Cloud Storage API
 - Cloud Logging API
- [ ] VM quota ≥ 96 vCPUs in target region
- [ ] Service account with `roles/compute.admin` + `roles/storage.admin`

## Cost Estimates

| Pipeline | Instance | Est. Cost | Est. Time |
|----------|----------|-----------|-----------|
| RNA-seq (28 species) | n2-highcpu-96 (spot) | $7–27 | 4–8 hours |
| Genome assembly | n2-highmem-32 (on-demand) | $47 | 12 hours |
| GWAS (10k samples) | n2-highcpu-64 (spot) | $3–10 | 1–2 hours |

**Stop VM promptly** — forgetful runs = $100+/day accidentally.

## Deploy Command Reference

```bash
# Show help
python scripts/cloud/deploy_gcp.py --help

# Deploy new run
python scripts/cloud/deploy_gcp.py deploy \
 --project YOUR_PROJECT_ID \
 --config config/cloud/run_config.yaml \
 --species-list config/amellifera.txt

# Check status (all VMs)
python scripts/cloud/deploy_gcp.py status \
 --project YOUR_PROJECT_ID

# Check specific VM
python scripts/cloud/deploy_gcp.py status \
 --project YOUR_PROJECT_ID \
 --instance rna-amellifera-001

# Stream logs
python scripts/cloud/deploy_gcp.py logs \
 --project YOUR_PROJECT_ID \
 --follow

# Download completed results
python scripts/cloud/deploy_gcp.py download \
 --project YOUR_PROJECT_ID \
 --output-dir local_results/

# Tear down (IMPORTANT: stop billing!)
python scripts/cloud/deploy_gcp.py destroy \
 --project YOUR_PROJECT_ID \
 --instance rna-amellifera-001 # Specific VM
# OR
python scripts/cloud/deploy_gcp.py destroy-all \
 --project YOUR_PROJECT_ID # All project VMs (careful!)
```

## Flags

```bash
--project PROJECT # GCP project ID (required)
--zone ZONE # Compute zone (default: us-central1-a)
--machine-type TYPE # VM shape (default: n2-highcpu-96)
--preemptible BOOL # Use spot VMs (default: True)
--docker-image IMAGE # Container image (default: metainformant/amplicon:latest)
--disk-size GB # Boot disk (default: 500)
--metadata KEY=VALUE # Custom GCE metadata (repeatable)
--no-auto-delete # Prevent auto-shutdown on failure
```

## One-Liner: Full Pipeline

```bash
# Deploy → wait → download → destroy (all-in-one)
python scripts/cloud/deploy_gcp.py full-pipeline \
 --project mycproj \
 --species-list config/28_hymenoptera.txt \
 --output data/cloud_results/ \
 --auto-destroy
```

## Cost Optimization Cheatsheet

```
[OK] Use --preemptible (default) 70% savings
[OK] Batch 80+ samples per VM Better VM utilization
[OK] Set --max-runtime-hours 6 Preemptible VM limit = 24h
[OK] Shut down immediately Don't pay for idle VMs
[OK] Use regional SSD (not persistent) for temp storage (auto-deletes)
 Don't use on-demand for batch 3–4× more expensive
 Don't leave VMs running overnight Accidental $100+ bills
```

## Advanced Examples

### Custom Docker image with proprietary tools
```bash
# Build Docker image locally, push to GCR
docker build -t gcr.io/my-project/metainformant-custom:latest -f Dockerfile.custom .
docker push gcr.io/my-project/metainformant-custom:latest

# Deploy with custom image
python scripts/cloud/deploy_gcp.py deploy \
  --project my-project \
  --config config/custom.yaml \
  --docker-image gcr.io/my-project/metainformant-custom:latest \
  --machine-type n2-highmem-64
```

### Spot VM interruption handling with checkpointing
```python-snippet
# In your pipeline script, register checkpoint handlers
from metainformant.cloud.interrupt import register_preempt_handler

def checkpoint_callback(state):
    """Save intermediate results when VM is about to die"""
    save_checkpoint(state)
    upload_to_gcs(f"gs://my-bucket/checkpoints/{state.job_id}.pkl")

register_preempt_handler(checkpoint_callback)

# Later, resume from checkpoint
python scripts/rna/run_amalgkit_single.py --config config.yaml --resume-from gs://...
```

### Multi-region redundant deployment
```bash
# Deploy same pipeline in two regions simultaneously
python scripts/cloud/deploy_gcp.py deploy \
  --project my-project \
  --config config.yaml \
  --zone us-central1-a &

python scripts/cloud/deploy_gcp.py deploy \
  --project my-project \
  --config config.yaml \
  --zone europe-west4-a &
```
This creates redundant workers; combine results post-completion.

## Expected Output

### Deployment status (JSON)
```json
{
  "instance": "rna-amellifera-001",
  "status": "RUNNING",
  "zone": "us-central1-a",
  "machine_type": "n2-highcpu-96",
  "preemptible": true,
  "start_time": "2026-04-26T08:15:23Z",
  "estimated_cost_per_hour": 0.183,
  "ssh_command": "gcloud compute ssh rna-amellifera-001 --zone us-central1-a",
  "log_stream": "https://console.cloud.google.com/logs/..."
}
```

### Real-time log streaming
```
[2026-04-26 08:16:01] VM startup complete
[2026-04-26 08:16:02] Pulling Docker image metainformant/amplicon:latest...
[2026-04-26 08:16:45] Image pulled (2.1 GB)
[2026-04-26 08:16:46] Starting container...
[2026-04-26 08:16:48] Container started (PID 1234)
[2026-04-26 08:16:50] Loading reference: data/refs/amellifera/GCA_003254395.2.fna
[2026-04-26 08:17:12] Reference indexed (STAR + kallisto)
[2026-04-26 08:17:15] Starting RNA-seq pipeline for 8 samples...
[2026-04-26 08:18:01] Aligning sample_001: 87.3% mapped (45.2M reads)
...
[2026-04-26 12:34:22] Pipeline complete
[2026-04-26 12:34:25] Uploading results to gs://output/amellifera/
[2026-04-26 12:41:33] Upload complete (47.2 GB)
[2026-04-26 12:41:35] Shutting down VM in 5 minutes...
[2026-04-26 12:46:01] VM terminated
```

### Cost report after run
```
Pipeline: rna-amellifera (8 samples)
VM type: n2-highcpu-96 (preemptible)
Runtime: 4h 32m
Cost: $0.18/hr × 4.53 hr = $0.82
Storage (regional SSD, 500 GB, 5h): $0.17
Network egress (to GCS same-region): $0.00
Total: $0.99
```

## Common Pitfalls

| Problem | Likely Cause | Fix |
|---------|-------------|-----|
| `QUOTA_EXCEEDED: vCPUs` | Project quota too low for requested machine type | Request quota increase via GCP Console; reduce `--machine-type` to n2-highcpu-32; or spread across multiple smaller VMs |
| `PERMISSION_DENIED: storage.objects.create` | Service account missing Storage Admin role | Grant `roles/storage.admin` to service account; verify `GOOGLE_APPLICATION_CREDENTIALS` env var |
| `PREEMPTED` after 1 hour (before spot limit) | GCP preemptible capacity shortage in chosen zone | Retry in different zone (`--zone us-west1-b`); switch to on-demand (`--preemptible false`) for critical runs |
| `SSH connection refused` after deploy | Firewall rule missing or VM not fully booted | Wait 60-90s after `RUNNING` status; check `gcloud compute firewall-rules list`; ensure port 22 open |
| Billing surprise > $100 | VM forgotten running (destroy failed) or `--no-auto-delete` set | Set `--auto-delete` flag; create Cloud Scheduler to kill all VMs nightly; enable budget alerts in GCP billing |
| `Docker build fails: pip timeout` | Slow network during Docker build in VM | Pre-build and push to GCR locally; use `--docker-image` flag; or use `--no-docker` for local-only mode |

---

**Related:** [Cloud module docs](../cloud/index.md) | [Economics guide](../cloud/ECONOMICS.md) | [GCP Console](https://console.cloud.google.com) | [Docker best practices](../cloud/DEPLOYMENT.md)
