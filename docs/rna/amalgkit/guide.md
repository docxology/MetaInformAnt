# Amalgkit Workflow Guide

Quick-start guide for running RNA-seq analysis workflows with METAINFORMANT's amalgkit integration.

## Prerequisites

Before running workflows:

1. **Python Environment**: Python 3.11+ with `uv` package manager
2. **Amalgkit CLI**: `uv pip install git+https://github.com/kfuku52/amalgkit`
3. **R Environment**: R with required packages (see [R_INSTALLATION.md](R_INSTALLATION.md))
4. **Genome Files**: Transcriptome FASTA and kallisto index (see [genome_setup_guide.md](genome_setup_guide.md))

## Quick Start

### 1. Verify Environment

```bash
# Check amalgkit availability
python3 -c "from metainformant.rna import check_cli_available; print(check_cli_available())"

# Check R availability
R --version
```

### 2. Configure Workflow

Create or modify a YAML configuration file in `config/amalgkit/`:

```yaml
# config/amalgkit/amalgkit_<species>.yaml
species: Pogonomyrmex_barbatus
taxonomy_id: 219557
work_dir: output/amalgkit/pbarbatus/
threads: 12

genome:
  assembly_accession: GCF_000187915.1
  # ... genome configuration
```

See existing configurations for examples.

### 3. Run Workflow

**Full end-to-end workflow** (recommended):

```bash
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml
```

**Check status**:

```bash
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status
```

**Specific steps only**:

```bash
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps getfastq quant merge
```

**Inspect step order and commands**:

```bash
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --plan
```

### 4. Validate Results

```bash
# Verify workflow completion
bash scripts/rna/amalgkit/verify_workflow.sh pbarbatus

# Check output files
ls output/amalgkit/pbarbatus/
```

## Workflow Steps

The amalgkit pipeline includes 11 steps:

| Step | Purpose | Output |
|------|---------|--------|
| **metadata** | Download SRA metadata | `work/metadata/` |
| **integrate** | Integrate external data | `work/metadata/` |
| **config** | Generate configuration | `work/config_base/` |
| **select** | Select samples | `work/pivot_qualified.tsv` |
| **getfastq** | Download FASTQ files | `fastq/` |
| **quant** | Quantify with kallisto | `quant/` |
| **merge** | Merge expression data | `merged/` |
| **cstmm** | Cross-species normalization | `cstmm/` |
| **curate** | Quality control | `curate/` |
| **csca** | Cross-species analysis | `csca/` |
| **sanity** | Validate results | `work/sanity/` |

See [steps/README.md](steps/README.md) for detailed step documentation.

## Python API

```python
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
uv pip install git+https://github.com/kfuku52/amalgkit
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

### "R package missing"

See [R_INSTALLATION.md](R_INSTALLATION.md) and [r_packages.md](r_packages.md).

## Related Documentation

- **[README.md](README.md)** - Complete amalgkit integration overview
- **[amalgkit.md](amalgkit.md)** - Detailed pipeline documentation
- **[FUNCTIONS.md](FUNCTIONS.md)** - Quick function lookup
- **[steps/README.md](steps/README.md)** - All 11 step guides
- **[genome_setup_guide.md](genome_setup_guide.md)** - Genome preparation

## Production Tips

1. **Use direct ENA downloads**: 100% reliability vs SRA Toolkit
2. **Enable immediate processing**: Minimizes disk usage
3. **Monitor with --status**: Check progress without interrupting
4. **Verify with verify_workflow.sh**: Comprehensive validation

---

*For complete documentation, see [README.md](README.md) and the step-specific guides in `steps/`.*





