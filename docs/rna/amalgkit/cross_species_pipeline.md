# Cross-Species RNA-seq Pipeline

This document describes the workflow for comparative analysis of RNA-seq data across multiple species using `amalgkit`'s cross-species steps (`cstmm` and `csca`).

## Prerequisites

- **Amalgkit** (v0.12.20+)
- **R** with required packages (see [R_INSTALLATION.md](R_INSTALLATION.md))
- **Quantified data** for each target species (completed `amalgkit quant` and `amalgkit merge` steps)
- **Ortholog mapping** (OrthoDB or similar)

## Directory Structure

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

## Generating the Ortholog Table

Use `scripts/rna/create_ortholog_table.py` to generate the table from OrthoDB data:

```bash
# Download required files from OrthoDB (https://www.orthodb.org/)
# - odb12v2_genes.tab.gz
# - odb12v2_OG2genes.tab.gz

python3 scripts/rna/create_ortholog_table.py \
  --species-ids 7460 144034 \
  --genes odb12v2_genes.tab.gz \
  --og2genes odb12v2_OG2genes.tab.gz \
  --output output/amalgkit/cross_species/orthologs/orthogroups.tsv
```

## Execution

Run after single-species pipelines (steps 1–7) are complete for all target species:

```bash
# Step 1: Cross-species TMM Normalization
amalgkit cstmm --config config/amalgkit/amalgkit_cross_species.yaml

# Step 2: Cross-species Correlation Analysis
amalgkit csca --config config/amalgkit/amalgkit_cross_species.yaml
```

Or via the workflow runner:

```bash
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_cross_species.yaml \
    --steps cstmm csca
```

## Outputs

- `cstmm/` — Normalized expression values (TMM-FPKM)
- `csca/` — Correlation heatmaps and PCA plots

## Download Architecture

The single-species pipelines feeding into cross-species steps use **ENA direct wget downloads**:

- `download_ena.py` queries the ENA API and downloads pre-compressed `.fastq.gz` files
- ~5 GB disk per sample vs ~25 GB for SRA Toolkit cache
- Supports 16+ concurrent workers safely

This enables the pipeline to run 16 parallel download workers (`num_download_workers: 16` in config) while keeping disk usage well below 100 GB total — compared to the ~400 GB that 16 concurrent SRA cache downloads would require.

## Multi-Species Orchestration

To process all 23 species before running cross-species steps:

```bash
# Sequential (robust): species one at a time
nohup python3 scripts/rna/run_all_species.py \
  > output/amalgkit/run_all_species_incremental.log 2>&1 &

# Or per-species in parallel via orchestrate_species.py
python3 scripts/rna/orchestrate_species.py

# Then run cross-species steps once all species complete
amalgkit cstmm --config config/amalgkit/amalgkit_cross_species.yaml
amalgkit csca  --config config/amalgkit/amalgkit_cross_species.yaml
```

## See Also

- [steps/08_cstmm.md](steps/08_cstmm.md) — cstmm step details
- [steps/10_csca.md](steps/10_csca.md) — csca step details
- [ORCHESTRATION.md](../../ORCHESTRATION.md) — Orchestrator selection
