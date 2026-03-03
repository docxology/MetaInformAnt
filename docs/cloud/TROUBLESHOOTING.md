# GCP Cloud Deployment — Troubleshooting & Learnings

Key issues encountered and resolved during the first deployment.

## 1. CPU Quota Limits

**Problem:** `Quota 'CPUS_ALL_REGIONS' exceeded. Limit: 32.0 globally.`

**Solution:** Default GCP projects have a 32 vCPU limit. Either:
- Use `n2-highcpu-32` (fits the default quota)
- Request quota increase via [GCP Console → Quotas](https://console.cloud.google.com/iam-admin/quotas)

## 2. Genome Indices Must Be Pre-Built

**Problem:** All 22 species failed with `No genome index found`.

**Root cause:** The `streaming_orchestrator.py` assumes pre-built kallisto indices in `output/amalgkit/shared/genome/*/index/`. The orchestrator does NOT download genomes.

**Solution:** Run `scripts/cloud/prep_genomes.py` before the pipeline:
```bash
cd /opt/MetaInformAnt && python3 scripts/cloud/prep_genomes.py --threads 8
```
This reads `genome.ftp_url` and `genome.files.transcriptome_fasta` from each species YAML config, downloads from NCBI FTP, and builds kallisto indices.

## 3. amalgkit Not on PyPI

**Problem:** `pip install amalgkit` fails — no matching distribution.

**Solution:** Install from GitHub:
```bash
pip3 install git+https://github.com/kfuku52/amalgkit.git
```

## 4. kallisto index Has No `-t` Flag

**Problem:** `kallisto index -t 8` fails with `invalid option -- 't'`.

**Solution:** `kallisto index` doesn't support threading. Remove `-t` flag.

## 5. SSH Output Buffering

**Problem:** `gcloud compute ssh --command "long_running_cmd"` shows no output for minutes.

**Solution:** For long commands, use `nohup` and check logs separately:
```bash
# Start in background
gcloud compute ssh VM --command "nohup python3 script.py > /tmp/log 2>&1 &"

# Check progress later
gcloud compute ssh VM --command "tail -20 /tmp/log"
```

## 6. File Permissions on VM

**Problem:** `/opt/MetaInformAnt/` owned by root (from startup script git clone).

**Solution:** SCP files to `~/`, then `sudo cp` to `/opt/`.

## 7. Disk Space Management

**Problem:** Local disk fills to 100% from accumulated FASTQ downloads.

**Solution:**
- The pipeline auto-deletes FASTQs after successful quantification
- Monitor with `df -h`
- GCP VMs have dedicated 500 GB SSD — much more headroom

## 8. Docker vs Native Install Trade-offs

| Approach | Pros | Cons |
|---|---|---|
| **Docker** | Reproducible, isolated | Slow build (~15 min), extra complexity |
| **Native** | Faster setup, simpler | Less reproducible, version drift risk |

**Recommendation:** Use native install for development/iteration. Docker for production reproducibility.
