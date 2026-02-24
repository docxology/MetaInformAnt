# Amalgkit FAQ

Frequently asked questions and solutions for the amalgkit RNA-seq workflow.

## 📊 Performance & Optimization

### Q: How does the pipeline process samples?

The pipeline uses **per-sample concurrency** within chunks. Each chunk of 6 samples is processed simultaneously using a `ThreadPoolExecutor`. For each sample, the steps are:

1. `getfastq` — Download and extract FASTQ files
2. `quant` — Quantify with kallisto
3. `cleanup` — Delete FASTQ files immediately

This means as soon as one sample finishes downloading, it starts quantifying while other samples continue downloading.

### Q: How many concurrent samples should I use?

**Recommended: 6 samples per chunk** (default in `run_all_species.sh --chunk-size 6`).

| Chunk Size | Threads/Sample | Disk Usage | Stability |
|------------|---------------|------------|-----------|
| 4 | 4 | ~40 GB temp | Very stable |
| **6** | **2** | **~60 GB temp** | **Stable (default)** |
| 8 | 2 | ~80 GB temp | May hit disk limits |
| 16 | 1 | ~160 GB temp | Not recommended |

### Q: How much disk space do I need?

**Minimum**: 80 GB free for 6 concurrent samples

| Component | Size |
|-----------|------|
| Per sample temp (FASTQ) | 5-15 GB |
| Per sample output (quant) | ~2 MB |
| Kallisto index | ~100-500 MB |
| Reference genome | ~200-500 MB |

### Q: What's a good processing rate?

| Rate | Assessment |
|------|------------|
| <4/hr | Slow - check network/disk |
| 4-8/hr | Normal |
| >10/hr | Excellent |

---

## ⚠️ Common Errors

### Q: "disk-limit exceeded" during extraction

**Cause**: `fasterq-dump` ran out of temp space.

**Solutions**:

1. Reduce chunk size (6 → 4) in `run_all_species.sh`
2. Free disk space
3. These samples are logged for retry on next run

### Q: "prefetch failed: path not found"

**Cause**: Race condition with multiple samples accessing same directory.

**Solution**: Reduce chunk size in `run_all_species.sh` or ensure sequential directory creation.

### Q: "prefetch failed: Lock file exists"

**Cause**: Previous download interrupted, lock file remains.

**Solution**: Delete stale lock files:

```bash
find output/amalgkit/*/fastq -name "*.lock" -delete
```

### Q: "Timeout" errors

**Cause**: Sample download/extraction exceeded time limit.

**Solutions**:

1. Network issues - retry later
2. Very large sample (>5 GB) - may need manual download

### Q: Samples stuck for hours

**Cause**: Process crashed without cleanup.

**Solution**:

```bash
# Stop the pipeline
pkill -f run_all_species

# Clean orphaned files
rm -rf output/amalgkit/*/fastq/getfastq/SRR*

# Restart
nohup bash scripts/rna/run_all_species.sh > output/amalgkit/run_all_species_incremental.log 2>&1 &
```

---

## 🔄 Recovery Procedures

### Q: How do I resume after a crash?

The pipeline automatically resumes thanks to `redo: no` in the configs:

```bash
nohup bash scripts/rna/run_all_species.sh > output/amalgkit/run_all_species_incremental.log 2>&1 &
```

It scans `quant/` for completed samples and skips them.

### Q: How do I retry failed samples?

1. Address the issue (usually disk space or network)
2. Restart the pipeline — it will retry any sample not yet in `quant/`

### Q: How do I recover from disk full?

```bash
# 1. Stop pipeline
pkill -f run_all_species

# 2. Clean all temp files
rm -rf output/amalgkit/*/fastq/getfastq/SRR*

# 3. Check disk
df -h /home

# 4. Restart with smaller chunk if needed (edit run_all_species.sh --chunk-size 4)
nohup bash scripts/rna/run_all_species.sh > output/amalgkit/run_all_species_incremental.log 2>&1 &
```

---

## 🐝 Species-Specific Notes

### Apis mellifera (Honeybee)

- **Total samples**: ~7,342
- **Estimated time**: Very long at current throughput
- **Large samples**: Some metagenomic samples (PRJNA1364028) are 5-7 GB each

### Pogonomyrmex barbatus (Harvester Ant)

- **Total samples**: 110
- **Notes**: Some samples may fail due to SRA issues, not workflow bugs

---

## 🛠️ Scripts Reference

| Script | Purpose |
|--------|---------|
| `scripts/rna/run_all_species.sh` | Main pipeline runner (sequential species, concurrent samples) |
| `scripts/rna/run_workflow.py` | Per-species workflow orchestrator |
| `scripts/package/generate_custom_summary.py` | Progress monitoring with ETAs |

### Running the Pipeline

```bash
# Run full pipeline (background)
nohup bash scripts/rna/run_all_species.sh > output/amalgkit/run_all_species_incremental.log 2>&1 &

# Check progress
.venv/bin/python scripts/package/generate_custom_summary.py

# Check which processes are active
ps -fC amalgkit | grep SRR
```

---

## 💾 Disk Cleanup

### Safe to delete

- `output/amalgkit/*/fastq/getfastq/SRR*` - temp FASTQ files (always safe after quant)

### NOT safe to delete

- `output/amalgkit/*/work/quant/` - final quantification results
- `output/amalgkit/*/work/metadata/` - sample metadata
- `output/amalgkit/*/work/index/` - kallisto index

---

## 🧬 Tissue Normalization

### Q: How do I normalize tissue metadata?

Use the tissue normalization script:

```bash
.venv/bin/python scripts/rna/normalize_tissue_metadata.py \
  --input output/amalgkit/apis_mellifera_all/work/metadata/metadata.tsv \
  --mapping config/amalgkit/tissue_mapping.yaml \
  --patches config/amalgkit/tissue_patches.yaml \
  --output output/amalgkit/apis_mellifera_all/work/metadata/metadata_normalized.tsv
```

### Q: How do I add tissue mappings?

Edit `config/amalgkit/tissue_mapping.yaml`:

```yaml
brain:
  - brain
  - Brain
  - whole brain
  - mushroom body
```

### Q: How do I patch missing tissue values?

Edit `config/amalgkit/tissue_patches.yaml`:

```yaml
samples:
  SRR12345678: brain  # From manual research
bioprojects:
  PRJEB100586: brain  # All samples in project
```

---

## 🔗 Related Resources

- [Amalgkit GitHub](https://github.com/kfuku52/amalgkit)
- [config/amalgkit/README.md](./README.md) - Configuration guide
- [config/amalgkit/AGENTS.md](./AGENTS.md) - Agent directives
