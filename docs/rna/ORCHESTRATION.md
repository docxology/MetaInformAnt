# RNA Workflow Orchestration

Overview of scripts for running RNA-seq workflows across species.

## Quick Links

- **[API Reference](API.md#orchestration-functions)** — Orchestration function documentation
- **[Workflow Guide](workflow.md)** — Workflow planning and execution
- **[Configuration Guide](CONFIGURATION.md)** — Configuration management
- **[Amalgkit Monitoring](amalgkit/monitoring.md)** — Monitoring active pipelines

---

## Main Orchestrator: `run_workflow.py` ⭐

**Script**: `scripts/rna/run_workflow.py`

Executes the full 11-step amalgkit pipeline for a single species.

| Feature | `run_workflow.py` |
|---------|------------------|
| Download method | ENA direct wget (100% reliable) |
| Processing mode | Per-sample: download → quantify → delete FASTQ |
| Parallel downloads | Configurable via `num_download_workers` in config |
| Auto-cleanup | Yes — FASTQs deleted after quantification |
| Status monitoring | `--status`, `--detailed` flags |

```bash
# Full end-to-end workflow
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pbarbatus.yaml

# Specific steps only
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pbarbatus.yaml \
    --steps getfastq quant merge

# Check status
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pbarbatus.yaml --status

# Cleanup downloaded but unquantified samples
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pbarbatus.yaml \
    --cleanup-unquantified
```

**Config for parallel downloads:**

```yaml
steps:
  getfastq:
    num_download_workers: 16
    threads: 24
    max_bp: 4000000000   # skip samples >4B bases
```

---

## Multi-Species: `run_all_species.py`

**Script**: `scripts/rna/run_all_species.py`

Runs `run_workflow.py` for all configured species sequentially. This is the main production orchestrator.

```bash
# Background, all 23 species, logging to file
nohup python3 scripts/rna/run_all_species.py \
  > output/amalgkit/run_all_species_incremental.log 2>&1 &

# Monitor
tail -f output/amalgkit/run_all_species_incremental.log
```

---

## Multi-Species Parallel: `orchestrate_species.py`

**Script**: `scripts/rna/orchestrate_species.py`

Alternative orchestrator that runs species through a configurable multi-species YAML.

```bash
python3 scripts/rna/orchestrate_species.py \
    --config config/amalgkit/cross_species.yaml
```

**Features:**
- Sequential/robust: one species at a time
- Resumable: skips completed steps
- Centralized logging: `output/amalgkit/logs/`

---

## Running Multiple Species in Parallel

For maximum throughput, run `run_workflow.py` in parallel for each species:

```bash
# In separate terminals or with nohup:
nohup python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_species1.yaml \
  > logs/species1.log 2>&1 &
nohup python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_species2.yaml \
  > logs/species2.log 2>&1 &
```

---

## Per-Sample Processing Workflow

The orchestrator implements an immediate per-sample processing workflow:

1. **Download** — FASTQ files downloaded via ENA wget for each sample
2. **Quantify** — Sample quantified immediately against the kallisto index
3. **Delete** — FASTQ files deleted after successful quantification

Only one sample's FASTQ files exist at any time → maximum disk efficiency.

---

## Manual Per-Sample Processing

For recovering or testing individual samples:

```python
from metainformant.rna.engine.workflow_steps import quantify_sample
from metainformant.rna.engine.sra_extraction import delete_sample_fastqs
from metainformant.rna.workflow import load_workflow_config
from metainformant.core.io import read_delimited
from pathlib import Path

cfg = load_workflow_config("config/amalgkit/amalgkit_pbarbatus.yaml")
rows = list(read_delimited(cfg.work_dir / "metadata" / "metadata.tsv", delimiter="\t"))
sample_rows = [r for r in rows if r.get("run") == "SRR14740514"]

success, message, abundance_path = quantify_sample(
    sample_id="SRR14740514",
    metadata_rows=sample_rows,
    quant_params={"out_dir": str(cfg.work_dir), "threads": 12},
    log_dir=cfg.log_dir,
)

if success and abundance_path and abundance_path.exists():
    fastq_dir = Path(cfg.per_step.get("getfastq", {}).get("out_dir", cfg.work_dir / "fastq"))
    delete_sample_fastqs("SRR14740514", fastq_dir)
```

---

## Batch Cleanup

Quantify all downloaded-but-unquantified samples and delete their FASTQs:

```python
from metainformant.rna.orchestration import cleanup_unquantified_samples
from pathlib import Path

quantified, failed = cleanup_unquantified_samples(
    Path("config/amalgkit/amalgkit_pbarbatus.yaml")
)
```

---

## Performance

| Metric | Value |
|--------|-------|
| Download speed | Fast (ENA direct wget) |
| Success rate | 100% (ENA-based) |
| Disk usage | Low (per-sample immediate cleanup) |
| Parallelism | `num_download_workers` downloads + 1 quant at a time |

---

## Troubleshooting

**Downloads failing:**
- Check `ps aux | grep wget | grep -v grep` for active workers
- Check `num_download_workers` is set in config file
- Verify network connectivity to ENA: `curl -s https://www.ebi.ac.uk/ena/portal/api/`

**Disk space issues:**
- Use `--cleanup-unquantified` to quantify and clean downloaded samples
- Use `--cleanup-partial` to remove partial downloads
- Reduce `num_download_workers` if disk is being exhausted

**Virtual environment issues:**
- Scripts auto-detect `.venv`
- Manual activation: `source .venv/bin/activate`

---

## Related Documentation

- **[CONFIGURATION.md](CONFIGURATION.md)** — Configuration management
- **[workflow.md](workflow.md)** — Workflow planning and execution
- **[amalgkit/monitoring.md](amalgkit/monitoring.md)** — Pipeline monitoring
- **[GETTING_STARTED.md](GETTING_STARTED.md)** — Setup and installation
