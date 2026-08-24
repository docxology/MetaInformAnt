# Amalgkit workflow guide

Quick-start guide for running RNA-seq analysis workflows with METAINFORMANT's amalgkit integration.

## Prerequisites

Before running workflows:

1. **Python environment**: Python 3.11+ with `uv` package manager
2. **Amalgkit CLI**: install the pinned project environment and verify `amalgkit --help`
3. **Reference assets**: transcriptome FASTA and kallisto index (see [genome_setup_guide.md](genome_setup_guide.md))
4. **External data root**: choose a mounted root with sufficient free space

## Quick Start

### 1. Verify Environment

```bash
# Check amalgkit availability
python3 -c "from metainformant.rna import check_cli_available; print(check_cli_available())"

# Check the installed Amalgkit version
amalgkit --version
```

### 2. Configure Workflow

Create or modify a YAML configuration file in `config/amalgkit/`:

```yaml
# config/amalgkit/amalgkit_<species>.yaml
species: Pogonomyrmex_barbatus
taxonomy_id: 219557
work_dir: output/amalgkit/pbarbatus/
threads: 16

genome:
  assembly_accession: GCF_000187915.1
  # ... genome configuration
```

See existing configurations for examples.

### 3. Run the current producer and downstream checkpoints

**Single-species dry run and execution**:

```bash
uv run python scripts/rna/run_all_species.py \
  --config-dir config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" \
  --species pogonomyrmex_barbatus --dry-run
uv run python scripts/rna/run_all_species.py \
  --config-dir config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" \
  --species pogonomyrmex_barbatus
bash projects/hymenoptera_amalgkit/scripts/run_all_finalization.sh \
  --data-root "$AMALGKIT_DATA_ROOT"
```

### 4. Monitor and Validate

```bash
# Check the current progress database and downstream evidence
uv run python scripts/rna/check_pipeline_status.py \
  --data-root "$AMALGKIT_DATA_ROOT" --verbose

# Check output files
ls output/amalgkit/pbarbatus/work/quant/ | wc -l
tail -20 "$AMALGKIT_DATA_ROOT"/results/full_campaign_*.log
```

## Workflow Steps

The current contract contains nine required per-species analysis stages plus
the `dataset` workspace-preparation command:

| Step | Purpose | Output |
|------|---------|--------|
| **metadata** | Retrieve SRA/ENA metadata | `work/metadata/` |
| **select** | Select samples with the pinned rules | `work/metadata/` |
| **getfastq** | Acquire FASTQ files with resume/reuse | `work/getfastq/` |
| **integrate** | Attach local read paths | `work/metadata/` |
| **quant** | Quantify with kallisto | `quant/` |
| **merge** | Merge expression data | `merged/` |
| **cstmm** | Cross-species normalization | `cstmm/` |
| **wsfilter** | Within-species filtering | `work/wsfilter/` |
| **finalize** | Final analysis tables | `work/finalize/` |
| **sanity** | Validate results | `work/sanity/` |

See [steps/README.md](steps/README.md) for detailed step documentation.

## Python API

```python-snippet
from metainformant.rna import workflow, amalgkit

# Check CLI availability
ok, help_text = amalgkit.check_cli_available()
if not ok:
    print(f"Amalgkit not available: {help_text}")

# Load and run workflow
config = workflow.load_workflow_config("config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml")
steps = workflow.plan_workflow(config)
results = workflow.execute_workflow(config)
```

## Common Issues

### "amalgkit not found"

Install amalgkit:
```bash
uv sync --extra rna  # installs the exact Amalgkit 0.16.60 project contract
```

### "No genome index"

Set up genome first:
```bash
python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_<species>.yaml
```

### "Disk space full"

- Use immediate processing mode (default): FASTQs deleted after quantification
- Set `TMPDIR` to repository temp: `export TMPDIR="$(pwd)/.tmp/bash"`
- See [DISK_SPACE_MANAGEMENT.md](../../DISK_SPACE_MANAGEMENT.md)

### "Output table missing"

Run the step-specific prerequisite validator, inspect the JSONL manifest, and
confirm that the configured input directory contains non-empty TSV tables.

## Related Documentation

- **[README.md](README.md)** - Complete amalgkit integration overview
- **[amalgkit.md](amalgkit.md)** - Detailed pipeline documentation
- **[FUNCTIONS.md](FUNCTIONS.md)** - Quick function lookup
- **[steps/README.md](steps/README.md)** — Current step guides
- **[genome_setup_guide.md](genome_setup_guide.md)** - Genome preparation

## Production Tips

1. **ENA HTTPS retrieval is the default**: valid local files are reused and
   incomplete files remain resumable `.part` files
2. **Raw reads are reclaimed after verified quantification**: disable with
   `AMALGKIT_RECLAIM_RAW_AFTER_QUANT=no`
3. **Monitor with `check_pipeline_status.py`**: current-contract counts from
   the selected data root
4. **Monitor the current log**: `tail -f "$AMALGKIT_DATA_ROOT"/results/full_campaign_*.log`
5. **Resume safely**: contract-verified quantification skips already completed
   samples, including compatible runtime drift; incomplete or unverifiable
   outputs are quarantined rather than silently re-downloaded

---

*For complete documentation, see [README.md](README.md) and the step-specific guides in `steps/`.*
