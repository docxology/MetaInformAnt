# Cross-Species RNA-seq Pipeline

This document describes the workflow for comparative analysis of RNA-seq data across multiple species using `amalgkit`.

## Prerequisites

- **Amalgkit** (v0.12.20+)
- **Python 3.10+**
- **Quantified data** for each target species (completed `amalgkit quant` and `amalgkit merge` steps)
- **Ortholog mapping** (OrthoDB or similar)

## Directory Structure

Structure the inputs as follows:

```
output/amalgkit/cross_species/
├── cstmm/              # Normalization output
├── csca/               # Correlation analysis output
├── orthologs/
│   └── orthogroups.tsv # Ortholog table (Orthogroup, Sp1, Sp2...)
└── work/
    └── merge_inputs/   # Input abundance matrices
        ├── Species_One/
        │   └── Species_One_est_counts.tsv
        └── Species_Two/
            └── Species_Two_est_counts.tsv
```

## Configuration

Create a configuration file (e.g., `amalgkit_cross_species.yaml`):

```yaml
work_dir: output/amalgkit/cross_species/work
species_list:
  - Species_One
  - Species_Two

steps:
  cstmm:
    out_dir: output/amalgkit/cross_species/cstmm
    orthogroup_table: output/amalgkit/cross_species/orthologs/orthogroups.tsv
    dir_count: output/amalgkit/cross_species/work/merge_inputs
    
  csca:
    out_dir: output/amalgkit/cross_species/csca
    orthogroup_table: output/amalgkit/cross_species/orthologs/orthogroups.tsv
    sample_group: tissue
```

## Generating Ortholog Table

Use the provided script `scripts/rna/create_ortholog_table.py` to generate the table from OrthoDB data:

```bash
# Download required files from OrthoDB
# - odb12v2_genes.tab.gz
# - odb12v2_OG2genes.tab.gz

python scripts/rna/create_ortholog_table.py \
  --species-ids 7460 144034 \
  --genes odb12v2_genes.tab.gz \
  --og2genes odb12v2_OG2genes.tab.gz \
  --output output/amalgkit/cross_species/orthologs/orthogroups.tsv
```

## Execution

Run the pipeline steps:

```bash
# Step 1: Cross-species TMM Normalization
amalgkit cstmm --config config/amalgkit/amalgkit_cross_species.yaml

# Step 2: Cross-species Correlation Analysis
amalgkit csca --config config/amalgkit/amalgkit_cross_species.yaml
```

## Outputs

- `cstmm/`: Contains normalized expression values (TMM-FPKM)
- `csca/`: Contains correlation heatmaps and PCA plots

## Download Mirror Architecture & 16-Worker Strategy

The pipeline utilizes a **Hybrid ENA / AWS Download Architecture** to maximize concurrency without exhausting local disk space.

### The Storage Constraint
Standard `fasterq-dump` (AWS/NCBI) requires ~25GB of temporary uncompressed `.sra` cache space per sample. Running 16 concurrent workers on AWS would instantly exhaust a 500GB NVMe drive. However, **EBI/ENA** provides pre-compressed `.fastq.gz` files natively via their API (only ~5GB per sample).

### The Hybrid Solution
The orchestrator (`run_all_species_parallel.py`) intercepts the download phase:
1. **Primary (ENA)**: A custom script (`scripts/rna/download_ena.py`) queries the ENA API and securely downloads the pre-compressed fastq files.
2. **Fallback (AWS)**: If a sample is missing from the EBI database, the orchestrator gracefully falls back to `amalgkit getfastq` using the AWS mirror.

By prioritizing ENA, the pipeline safely runs **16 concurrent workers** (1 CPU thread each), fully saturating the internet connection (~36+ MB/s aggregate) while keeping disk usage strictly under ~80GB. 

**Execution Command:**
Use `--max-workers 16` to launch this maximum-throughput mode:
```bash
.venv/bin/python scripts/rna/run_all_species_parallel.py --max-workers 16
```
