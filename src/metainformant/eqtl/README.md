# eQTL Module

eQTL and transcriptome-variant analysis: synthetic/real input-data construction, HISAT2+bcftools variant-calling wrappers, bcftools stats parsing, and full pipeline orchestration.

## Overview

The `eqtl` module holds the reusable methods behind the eQTL integration
pipeline (GWAS variants x Amalgkit RNA-seq expression). Scripts under
`scripts/eqtl/` are thin orchestrators that import from here.

## Table of Contents

- [Architecture](#architecture)
- [Submodules](#submodules)
- [Quick Start](#quick-start)
- [Related](#related)

## Architecture

```mermaid
graph TD
    subgraph "eQTL Module"
        SYN[workflow/synthetic.py] --> |matrices| SCAN[gwas.finemapping.eqtl scan]
        VC[workflow/variant_calling.py] --> |HISAT2/samtools/bcftools| VCF[per-sample VCFs]
        VS[variant_stats.py] --> |bcftools stats| SUM[variant summaries]
        PL[pipeline.py] --> VC
        PL --> VS
    end

    SYN --> SCAN
    VCF --> VS
```

## Submodules

| Module | Purpose |
|--------|---------|
| [`workflow/synthetic.py`](workflow/synthetic.py) | Synthetic expression/genotype generators, kallisto target-ID position parsing (synthetic TSS positions), real quantification loader |
| [`workflow/variant_calling.py`](workflow/variant_calling.py) | Tool checks, FASTQ acquisition, HISAT2 index/alignment, bcftools call/filter/merge |
| [`variant_stats.py`](variant_stats.py) | bcftools stats parsing, per-sample and population summaries, allele-frequency extraction |
| [`pipeline.py`](pipeline.py) | Parameter resolution and end-to-end run orchestration |

## Quick Start

```python
from metainformant.eqtl import create_synthetic_data
from metainformant.gwas.finemapping.eqtl import cis_eqtl_scan

expr, geno, gene_pos, var_pos = create_synthetic_data(n_genes=20, n_variants=100, n_samples=30)
results = cis_eqtl_scan(expr, geno, gene_pos, var_pos, cis_window=1_000_000)
```

```bash
# Full pipeline (requires hisat2, samtools, bcftools, curl and AMALGKIT_DATA_ROOT)
uv run python scripts/eqtl/rna_snp_pipeline.py --species amellifera --n-samples 3

# Demo with synthetic data
uv run python scripts/eqtl/run_eqtl_demo.py
```

## Real and Synthetic Inputs

`load_real_expression_data()` reads real expression. `parse_gene_positions()` reads only the chromosome accession from each kallisto target ID and invents the TSS positions (evenly spaced, in the order of the IDs); real positions need a GFF/GTF annotation. `create_synthetic_genotypes()` draws random dosages and names its variants `var_<full gene id>_<i>`, built from the whole gene ID so that they are unique whenever the gene IDs are (a repeated gene ID raises `ValueError`). A scan on these inputs pairs real expression with synthetic genotypes and positions, so its associations are not biological findings: `scripts/eqtl/run_eqtl_real.py` labels its gene positions, scan, summaries and plots accordingly (`gene_positions_synthetic.tsv`, `*_synthetic_genotypes.*`).

`cis_eqtl_scan()` raises `ValueError` on a repeated ID in the gene positions, the variant positions or the genotype index. A fit that cannot be made has NaN `beta`, `se` and `pvalue`, `status == "error"` and an `error` message (also in `conditional_eqtl()` and `eqtl_effect_sizes()`). `eqtl_summary_stats()` reports `n_tests` (tests with a p-value) and `n_tests_failed`, and `write_eqtl_results()` writes missing values as `NA`.

## Tests

Zero-mocks tests: `tests/eqtl/` (synthetic frames, real bcftools parsing on
fixture text, subprocess-free where tools are unavailable).

## Related

- [eQTL Scripts](../../../scripts/eqtl/) - thin orchestrators
- [GWAS Fine-Mapping](../gwas/finemapping/) - cis/trans scans, effect sizes
- [RNA Core](../rna/) - sample_utils, ENA downloader
