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
