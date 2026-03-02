# Amalgkit Pipeline Steps

Documentation for all 11 amalgkit processing steps used in METAINFORMANT.

## Step Overview

| # | Step | Purpose | Input | Output |
|---|------|---------|-------|--------|
| 1 | [metadata](01_metadata.md) | Fetch SRA sample metadata from NCBI | NCBI search string | `work/metadata/metadata.tsv` |
| 2 | [config](02_config.md) | Generate amalgkit config files | metadata.tsv | `work/config_base/` |
| 3 | [select](03_select.md) | Filter samples by quality/tissue | config_base/ | `work/metadata/pivot_qualified.tsv` |
| 4 | [getfastq](04_getfastq.md) | Download FASTQ from ENA via wget | pivot_qualified.tsv | `fastq/getfastq/<SRR>/` |
| 5 | [integrate](05_integrate.md) | Add local FASTQ paths to metadata | getfastq output + metadata | `work/metadata/metadata_integrated.tsv` |
| 6 | [quant](06_quant.md) | Quantify expression with kallisto | FASTQ files + kallisto index | `work/quant/<SRR>/abundance.tsv` |
| 7 | [merge](07_merge.md) | Merge sample abundances into matrix | quant outputs | `merged/merged_abundance.tsv` |
| 8 | [cstmm](08_cstmm.md) | Cross-species TMM normalization | merged matrices + orthologs | `cstmm/` |
| 9 | [curate](09_curate.md) | Quality control, outlier removal | cstmm output | `curate/` |
| 10 | [csca](10_csca.md) | Cross-species correlation analysis | curate output | `csca/` |
| 11 | [sanity](11_sanity.md) | Validate all outputs | all step outputs | `work/sanity/` |

## Key Notes

### Download Method (Step 4)
Downloading uses **ENA direct wget** (`scripts/rna/download_ena.py`), **not** the SRA Toolkit. This provides:
- 100% reliability (pre-compressed `.fastq.gz` from EBI servers)
- Resumable downloads (`wget --continue`)
- Smaller disk footprint (~5 GB/sample vs ~25 GB for SRA cache)

### Quantification (Step 6)
- Runs `kallisto quant` against a pre-built species index
- FASTQs are **immediately deleted** after successful quantification
- A non-empty `abundance.tsv` is the canonical proof of work
- Size-based filtering: samples >4B bases are auto-skipped (configurable with `max_bp`)

### Cross-Species Steps (8–10)
Steps 8–10 require completed data from **all target species**. Run per-species pipelines (steps 1–7) first, then run the cross-species steps together with an ortholog table.

## Execution

```bash
# End-to-end for one species
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_<species>.yaml

# Specific steps only
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_<species>.yaml \
    --steps getfastq quant merge

# All 23 species sequentially
nohup python3 scripts/rna/run_all_species.py \
  > output/amalgkit/run_all_species_incremental.log 2>&1 &
```

## Output Structure

```
output/amalgkit/<species>/
├── work/
│   ├── metadata/          # Steps 1–3 outputs
│   ├── config_base/       # Step 2 output
│   ├── index/             # Kallisto index (pre-built)
│   ├── fasta/             # Transcriptome FASTA
│   └── quant/<SRR>/       # Step 6 outputs (abundance.tsv per sample)
├── fastq/getfastq/<SRR>/  # Step 4 (temporary — deleted after quant)
├── merged/                # Step 7 output
├── cstmm/                 # Step 8 output
├── curate/                # Step 9 output
├── csca/                  # Step 10 output
└── logs/                  # Per-step logs
```

## See Also

- [Amalgkit README](../README.md) — Integration overview
- [Guide](../guide.md) — Quick-start
- [ORCHESTRATION.md](../../ORCHESTRATION.md) — Orchestrator scripts
- [genome_setup_guide.md](../genome_setup_guide.md) — Reference genome preparation
