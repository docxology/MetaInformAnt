# Amalgkit FAQ

Frequently asked questions and solutions for the amalgkit RNA-seq workflow.

## 📊 Performance & Optimization

### Q: How many parallel workers should I use?

**Recommended: 4 workers** for most systems.

| Workers | Disk Usage | Stability | Rate |
|---------|------------|-----------|------|
| 2 | ~20 GB temp | Very stable | ~4/hr |
| **4** | **~40 GB temp** | **Stable** | **~8/hr** |
| 8 | ~80 GB temp | Race conditions possible | ~12/hr* |

*8 workers can cause `prefetch` lock file conflicts and disk-limit errors.

### Q: How much disk space do I need?

**Minimum**: 50 GB free for 4 workers

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

1. Reduce workers (4 → 2)
2. Free disk space
3. These samples are logged to `failed_samples.json` for retry

### Q: "prefetch failed: path not found"

**Cause**: Race condition with multiple workers accessing same directory.

**Solution**: Reduce to 4 workers or ensure sequential directory creation.

### Q: "prefetch failed: Lock file exists"

**Cause**: Previous download interrupted, lock file remains.

**Solution**: Delete stale lock files:

```bash
find output/amalgkit/*/fastq -name "*.lock" -delete
```

### Q: "Timeout" errors

**Cause**: Sample download/extraction exceeded 1-hour limit.

**Solutions**:

1. Network issues - retry later
2. Very large sample (>5 GB) - may need manual download

### Q: Samples stuck in "In Progress" for hours

**Cause**: Process crashed without cleanup.

**Solution**:

```bash
# Stop processor
pkill -f process_apis_mellifera_parallel

# Clean orphaned files
rm -rf output/amalgkit/*/fastq/getfastq/SRR*

# Restart
nohup uv run python scripts/rna/process_apis_mellifera_parallel.py &
```

---

## 🔄 Recovery Procedures

### Q: How do I resume after a crash?

The parallel processor automatically resumes:

```bash
uv run python scripts/rna/process_apis_mellifera_parallel.py
```

It scans `quant/` for completed samples and skips them.

### Q: How do I retry failed samples?

Failed samples are logged to `work/failed_samples.json`. To retry:

1. Check the file for failure reasons
2. Address the issue (usually disk space)
3. Restart the processor - it will retry any sample not in `quant/`

### Q: How do I recover from disk full?

```bash
# 1. Stop processor
pkill -f process_apis_mellifera_parallel

# 2. Clean all temp files
rm -rf output/amalgkit/*/fastq/getfastq/SRR*

# 3. Check disk
df -h /Users

# 4. Restart with fewer workers if needed
```

---

## 🐝 Species-Specific Notes

### Apis mellifera (Honeybee)

- **Total samples**: ~7,242
- **Estimated time**: 40-50 days at 4 workers
- **Large samples**: Some metagenomic samples (PRJNA1364028) are 5-7 GB each
- **Filtered**: 15 metagenomic samples excluded from processing

### Pogonomyrmex barbatus (Harvester Ant)

- **Total samples**: 110
- **Completed**: ✅ 95/110 quantified
- **Notes**: Some samples failed due to SRA issues, not workflow bugs

---

## 🛠️ Scripts Reference

| Script | Purpose |
|--------|---------|
| `process_apis_mellifera_parallel.py` | Main 4-worker parallel processor |
| `monitor_apis_mellifera.py` | Real-time progress monitoring |
| `recover_missing_parallel.py` | Retry failed downloads via ENA |
| `run_workflow.py` | Full amalgkit workflow orchestration |

### Monitoring Commands

```bash
# One-time check
uv run python scripts/rna/monitor_apis_mellifera.py

# Live monitoring (refreshes every 30s)
uv run python scripts/rna/monitor_apis_mellifera.py --watch

# Check logs
tail -50 output/amalgkit/apis_mellifera_all/work/parallel_processing.log
```

---

## 💾 Disk Cleanup

### Safe to delete

- `output/amalgkit/*/fastq/getfastq/SRR*` - temp FASTQ files (always safe after quant)
- `~/.ollama/models` - if not using local LLMs
- `~/Downloads/*.dmg` - old installers

### NOT safe to delete

- `output/amalgkit/*/work/quant/` - final quantification results
- `output/amalgkit/*/work/metadata/` - sample metadata
- `output/amalgkit/*/work/index/` - kallisto index

---

## 🔗 Related Resources

- [Amalgkit GitHub](https://github.com/kfuku52/amalgkit)
- [config/amalgkit/README.md](./README.md) - Configuration guide
- [config/amalgkit/AGENTS.md](./AGENTS.md) - Agent directives
