# Amalgkit Integration

METAINFORMANT wraps [amalgkit](https://github.com/kfuku52/amalgkit) — a command-line toolkit for large-scale RNA-seq meta-analysis — to orchestrate all 11 processing steps across 23 ant/bee species.

## What Amalgkit Does

Amalgkit handles the full pipeline from public database → expression matrix:

1. Downloads SRA metadata from NCBI
2. Filters samples by quality criteria
3. Downloads FASTQ files (via ENA direct wget in this project)
4. Quantifies with kallisto pseudoalignment
5. Merges per-sample abundances into expression matrices
6. Cross-species normalization (CSTMM) and correlation analysis (CSCA)

## How METAINFORMANT Uses It

```
scripts/rna/run_all_species.py
   └─ iterates species configs → run_workflow.py per species
         └─ per sample: download_ena.py (wget) → amalgkit quant → delete FASTQ
```

The orchestrator (`run_all_species.py`) drives sequential per-species execution. Each species config lives in `config/amalgkit/amalgkit_<species>.yaml`.

## Directory Structure

```
docs/rna/amalgkit/
├── README.md               ← This file
├── guide.md                ← Quick-start guide
├── amalgkit.md             ← Wrapper details and configuration
├── commands.md             ← Genome setup script reference
├── cross_species_pipeline.md  ← CSTMM + CSCA cross-species steps
├── monitoring.md           ← How to check pipeline progress
├── genome_preparation.md   ← Genome download and kallisto index build
├── genome_setup_guide.md   ← Step-by-step genome setup
├── FUNCTIONS.md            ← Python function index
├── PATH_RESOLUTION.md      ← Path resolution reference
├── R_INSTALLATION.md       ← R environment setup
├── r_packages.md           ← R package management
├── testing_coverage.md     ← Test coverage details
├── AGENTS.md               ← AI contribution notes
├── PAI.md                  ← Persistent AI context
└── steps/                  ← Per-step documentation (01–11)
```

## 11-Step Pipeline

| # | Step | Purpose | Key Output |
|---|------|---------|-----------|
| 1 | `metadata` | Fetch SRA sample metadata from NCBI | `work/metadata/metadata.tsv` |
| 2 | `config` | Generate amalgkit config files | `work/config_base/` |
| 3 | `select` | Filter samples by quality/tissue criteria | `work/metadata/pivot_qualified.tsv` |
| 4 | `getfastq` | Download FASTQ from ENA (wget) | `fastq/getfastq/<SRR>/` |
| 5 | `integrate` | Integrate local FASTQ paths into metadata | `work/metadata/metadata_integrated.tsv` |
| 6 | `quant` | Quantify with kallisto (per sample) | `work/quant/<SRR>/abundance.tsv` |
| 7 | `merge` | Merge sample abundances into matrix | `merged/merged_abundance.tsv` |
| 8 | `cstmm` | Cross-species TMM normalization | `cstmm/` |
| 9 | `curate` | Quality control, outlier removal | `curate/` |
| 10 | `csca` | Cross-species correlation analysis | `csca/` |
| 11 | `sanity` | Validate outputs | `work/sanity/` |

See [steps/README.md](steps/README.md) for detailed per-step documentation.

## Installation

```bash
# Install amalgkit into the project venv
uv pip install git+https://github.com/kfuku52/amalgkit

# Verify
amalgkit --help
```

## Quick Reference

```bash
# Run a single species end-to-end
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pbarbatus.yaml

# Run all 23 species sequentially (background)
nohup python3 scripts/rna/run_all_species.py \
  > output/amalgkit/run_all_species_incremental.log 2>&1 &

# Check status
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pbarbatus.yaml --status

# Monitor live progress
tail -f output/amalgkit/run_all_species_incremental.log
ps aux | grep wget | grep -v grep | wc -l  # active download workers
python3 scripts/rna/report_completed.py     # quant counts per species
```

## Related Documentation

- **[guide.md](guide.md)** — Quick-start guide with prerequisites
- **[steps/README.md](steps/README.md)** — All 11 steps with status table
- **[../ORCHESTRATION.md](../ORCHESTRATION.md)** — Orchestrator script documentation
- **[../CONFIGURATION.md](../CONFIGURATION.md)** — YAML config reference
- **[genome_setup_guide.md](genome_setup_guide.md)** — Reference genome preparation
