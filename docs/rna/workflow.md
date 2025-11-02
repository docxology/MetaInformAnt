# RNA Workflow

High-level planning and execution live in `metainformant.rna.workflow`.

## Quick Start

**For multi-species production workflows**, see **[Multi-Species Quick Start Guide](MULTI_SPECIES_QUICK_START.md)** for:
- Starting parallel workflows for multiple species
- Monitoring progress in real-time
- Troubleshooting common issues
- Performance optimization

## Plan

```python
from pathlib import Path
from metainformant.rna.workflow import AmalgkitWorkflowConfig, plan_workflow

cfg = AmalgkitWorkflowConfig(
    work_dir=Path("output/amalgkit/run1"),
    threads=12,  # Default: 12 threads for downloads and quantification
    species_list=["Apis_mellifera"]
) 
for name, params in plan_workflow(cfg):
    print(name, params)
```

## Execute

```python
from metainformant.rna.workflow import execute_workflow
codes = execute_workflow(cfg)
print(codes)
```

**Command-line execution** (production ENA workflow):
```bash
# Production workflow - recommended for reliability
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --batch-size 12 \
  --threads 12
```

**Legacy SRA-based workflow:**
```bash
# Alternative using SRA Toolkit (auto-activates venv)
python3 scripts/rna/run_multi_species.py
```

The ENA-based workflow provides:
- Direct ENA downloads with 100% reliability (vs 0% SRA Toolkit)
- Batched processing (12 samples at a time)
- Automatic resume with `wget --continue`
- Disk space management (~18 GB peak per batch)

The legacy SRA-based workflow provides:
- Automatic virtual environment activation
- Batched processing (10 samples at a time)
- SRA download optimization
- Disk space management (~20-50 GB peak usage)

### From config file

```python
from metainformant.rna.workflow import load_workflow_config, execute_workflow
cfg = load_workflow_config("config/amalgkit_pbarbatus.yaml")
codes = execute_workflow(cfg, check=True)
```

## Per-step params

```python
from metainformant.rna.configs import SpeciesProfile, AmalgkitRunLayout, build_step_params
from metainformant.rna.workflow import plan_workflow_with_params  # Note: Not exported, use plan_workflow instead

species = SpeciesProfile(name="Apis mellifera", taxon_id=7460, tissues=["brain"]) 
layout = AmalgkitRunLayout(base_dir=cfg.work_dir)
params_map = build_step_params(species, layout)
steps = plan_workflow_with_params(cfg, params_map)
```

## Artifacts

- Logs under `work_dir/logs/` per step
- JSON Lines manifest at `work_dir/amalgkit.manifest.jsonl`
- Summary JSON at `work_dir/amalgkit.report.json` and Markdown at `work_dir/amalgkit.report.md`

### RNA: Workflow

Plan and execute an `amalgkit`-based workflow using `AmalgkitWorkflowConfig`.

Functions: `plan_workflow`, `plan_workflow_with_params`, `execute_workflow`.

```mermaid
sequenceDiagram
  participant U as User
  participant W as plan_workflow
  participant X as execute_workflow
  U->>W: config
  W-->>U: [(step, params) ...]
  U->>X: config, check=True/False
  X->>X: run_amalgkit per step
  X-->>U: [return codes]
```

## Example

```python
from pathlib import Path
from metainformant.rna import AmalgkitWorkflowConfig
from metainformant.rna import workflow as wf

cfg = AmalgkitWorkflowConfig(
    work_dir=Path("./work"),
    threads=12,  # Default: 12 parallel threads
    species_list=["Apis_mellifera"]
)
steps = wf.plan_workflow(cfg)
codes = wf.execute_workflow(cfg, check=False)
```

## Monitoring Workflows

See **[Multi-Species Quick Start Guide](MULTI_SPECIES_QUICK_START.md#monitoring-progress)** for comprehensive monitoring instructions.

**Quick monitoring:**
```bash
# Real-time comprehensive monitor
python3 scripts/rna/monitor_comprehensive.py

# Continuous watch mode
watch -n 60 python3 scripts/rna/monitor_comprehensive.py

# Check running processes
ps aux | grep workflow_ena | grep -v grep
```

## How this supports meta-analysis

The workflow produces standardized, well-logged artifacts that map to common meta-analysis stages:

- Discovery/selection (`metadata`, `integrate`, `config`, `select`): defines the population of samples across studies
- Acquisition/quantification (`getfastq`, `quant`): yields per-sample abundance estimates using consistent parameters
- Aggregation/normalization (`merge`, `cstmm`): constructs study-agnostic matrices for downstream analysis
- Curation/QC (`curate`, `csca`, `sanity`): ensures input quality and comparability

Downstream, you can:

- Run DE per study and combine with p-value or effect-size methods (e.g., metaRNASeq)
- Apply batch-effect correction on the merged matrix using tools like ComBat or limma
- Perform GO/pathway enrichment using `metainformant.ontology` utilities or external R/Python libraries

All outputs default under `output/` in keeping with repository policy; override via `work_dir` in the config.

## Common Issues and Solutions

### Metadata Format Selection

The workflow intelligently selects metadata files for downstream steps. Steps like `getfastq`, `integrate`, `quant`, and `merge` require **row-per-sample format** with a `run` column containing SRA IDs.

**Issue**: `amalgkit select` creates pivot tables that lack run IDs, causing "No SRA entry found" errors.

**Solution**: The workflow automatically falls back to `metadata.filtered.tissue.tsv` or `metadata.tsv` when pivot tables are detected:

```python
# Automatic detection and fallback in workflow.py
if rows and 'run' not in rows[0]:  # Pivot format detected
    # Use metadata.filtered.tissue.tsv instead
```

### Filter Configuration

**Issue**: Overly restrictive filters can exclude 99%+ of samples.

**Solution**: Start with minimal filtering and adjust based on data quality:

```yaml
# Recommended starting point
filters:
  require_tissue: true  # Essential filter only
  
# Add stricter filters only if needed
filters:
  require_tissue: true
  min_spots: 1000000  # Add only after reviewing initial results
```

### Performance Optimization

**GetFASTQ Acceleration**: Enable parallel downloads and cloud sources for significant speedup:

```yaml
steps:
  getfastq:
    threads: 10  # Default: 10 parallel downloads
    aws: yes     # Use AWS Open Data Program
    gcp: yes     # Use Google Cloud
    ncbi: yes    # Use NCBI directly
```

**Automatic Optimizations** (handled by `run_multi_species.py`):
- SRA wrapper script: Disables conservative size checks
- Temp directory: Uses project location instead of `/tmp`
- Batched processing: 10 samples at a time (~20-50 GB peak)
- Auto-activation: Virtual environment detection and activation

### Disk Space Management

**Issue**: Large samples (100+ GB) fail with "disk-limit exceeded"

**Solution**: The script automatically handles this:
1. Creates wrapper: `output/sra_temp/fasterq-dump` with `--size-check off`
2. Sets temp directory to project location (not `/tmp`)
3. Batches processing: Downloads → Quantifies → Deletes FASTQs
4. Peak usage: ~20-50 GB (10 samples at a time)

No manual configuration needed - all handled automatically by the workflow script.

### Virtual Environment Issues

**Issue**: "amalgkit: command not found" or import errors

**Solution**: The script automatically activates virtual environment:
```bash
# Just run the script - no manual activation needed
python3 scripts/rna/run_multi_species.py
```

The script detects if it's running outside the venv and re-executes itself within it using `os.execve()`.

### SRA Download Failures

**Issue**: fasterq-dump fails with size check errors despite available disk space

**Solution**: Automatic wrapper script injection:
- Script creates: `output/sra_temp/fasterq-dump`
- Adds flag: `--size-check off`
- Symlinks tools: fastp, kallisto, seqkit
- Updates PATH automatically

All handled transparently by the workflow.
steps:
  getfastq:
    threads: 6           # Parallel processing
    pfd: yes             # Use parallel-fastq-dump
    accelerate: true     # Enable cloud mirrors (AWS/GCP)
    max_size: "50GB"     # Handle large samples
```

**Disk Space Management**: Delete FASTQ files after quantification to prevent disk exhaustion:

```yaml
steps:
  quant:
    threads: 6
    keep_fastq: no       # Delete FASTQs after processing
    redo: no             # Skip already quantified samples
```

**Expected Performance** (with optimizations):
- Download speed: 6x faster with parallel processing
- Cloud sources: 3-5x faster than NCBI only
- Disk usage: 90% reduction with FASTQ cleanup
- Total runtime: ~1-2 days for 300+ samples (vs 6-13 days without optimization)
