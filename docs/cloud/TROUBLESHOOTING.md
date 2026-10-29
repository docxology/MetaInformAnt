# GCP Cloud Deployment — Troubleshooting & Learnings

Key issues encountered and resolved during the first deployment.

> These are historical GCP lessons. The original Hymenoptera cloud
> infrastructure was decommissioned on 2026-04-15; do not treat these pages as
> evidence that a current VM or cloud dataset exists.

## 1. CPU Quota Limits

**Problem:** `Quota 'CPUS_ALL_REGIONS' exceeded. Limit: 32.0 globally.`

**Solution:** The project's vCPU limit was 32 on that deployment. Either:
- Use `n2-highcpu-32` (fits a 32 vCPU quota)
- Request quota increase via [GCP Console → Quotas](https://console.cloud.google.com/iam-admin/quotas)

`CloudConfig` defaults to `n2-highcpu-96` (96 vCPUs) and `deploy_gcp.py deploy` to `n2-standard-32` (32 vCPUs); pass `--machine-type` to match the quota.

## 2. Genome Indices

**Problem:** A species failed because its kallisto index was missing.

**Root cause:** Each species needs a non-empty kallisto index (`*.idx`) in one of the directories the orchestrator checks: `steps.quant.index_dir`, `<genome.dest_dir>/index`, `genome.index_dir`, `<work_dir>/index` or `<data root>/shared/genome/<species>/index`. When none has one, it logs `No genome index found` and builds the index through `prepare_reference_genome`, which downloads `genome.files.transcriptome_fasta` from `genome.ftp_url` with `curl` and runs `kallisto index -i`. That step fails if the download fails its gzip check or `kallisto` is not on `PATH`.

**Solution:** Make sure `kallisto` is installed, or build every index ahead of the run with `scripts/cloud/prep_genomes.py`:
```bash
cd /opt/MetaInformAnt && python3 scripts/cloud/prep_genomes.py
```
This reads `genome.ftp_url` and `genome.files.transcriptome_fasta` from each species YAML in `--config-dir` (default `config/amalgkit`), downloads them from NCBI with `wget`, and writes `<species>.idx` to `<genome.dest_dir>/index/` (default `output/amalgkit/shared/genome/<species>/index/` when the YAML sets no `dest_dir`).

## 3. amalgkit Not on PyPI

**Problem:** a generic PyPI installation does not provide the project-pinned
Amalgkit release.

**Solution:** install the pinned project extra:
```bash
uv sync --extra rna
```
The `rna` extra pins `amalgkit` to the `v0.16.60` tag, and the Docker image installs `amalgkit=0.16.60` from bioconda. `scripts/cloud/vm_setup.sh` installs amalgkit from `git+https://github.com/kfuku52/amalgkit` without a tag.

## 4. kallisto index Is Run Without Threads

**Problem:** `kallisto index -t 8` failed with `invalid option -- 't'` on the deployed kallisto build.

**Solution:** Both index builders (`prepare_reference_genome` and `scripts/cloud/prep_genomes.py`) run `kallisto index -i <index> <fasta>` with no thread option. `prep_genomes.py` accepts `--threads` but does not use it.

## 5. SSH Output Buffering

**Problem:** `gcloud compute ssh --command "long_running_cmd"` shows no output for minutes.

**Solution:** For long commands, use `nohup` and check logs separately (`VM` is a placeholder for your instance name; `GCPDeployer` runs its status and log commands over SSH with a 30 second timeout, so it is meant for short checks):
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
- By default (`AMALGKIT_RECLAIM_RAW_AFTER_QUANT=yes`) the pipeline removes a sample's FASTQ/SRA files after a verified quantification
- Monitor with `df -h`
- A GCP VM's boot disk is `pd-standard` with `disk_size_gb` (500 GB in `CloudConfig`, 1000 GB in `deploy_gcp.py deploy`; `--disk-gb` changes it)

## 8. Docker vs Native Install Trade-offs

| Approach | Pros | Cons |
|---|---|---|
| **Docker** | Reproducible, isolated | Slow image build, extra complexity |
| **Native** | Faster setup, simpler | Less reproducible, version drift risk |

**Recommendation:** Use native install for development/iteration. Docker for production reproducibility.
