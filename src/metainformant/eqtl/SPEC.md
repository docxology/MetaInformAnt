# Specification: eqtl

## Scope

eQTL and transcriptome-variant analysis for METAINFORMANT. Synthetic/real
eQTL input construction, transcriptome SNP calling (HISAT2 + bcftools),
variant statistics, and pipeline orchestration.

## Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code

## Data Structures

- **Submodules**: workflow.synthetic, workflow.variant_calling, variant_stats, pipeline
- **Key Concepts**: cis-eQTL inputs (expression genes x samples, genotype
  dosages variants x samples, gene/variant positions), HISAT2 index and
  BAM alignment, bcftools mpileup/call/filter/merge, Ti/Tv and record
  counts from bcftools stats

## API Definition

### Exports — `workflow/synthetic.py`

- `create_synthetic_data` — Synthetic expression/genotype matrices with known cis-eQTL effects
- `parse_gene_positions` — Approximate positions from kallisto target IDs
- `create_synthetic_genotypes` — Synthetic dosages matched to gene positions
- `load_real_expression_data` — Real quantification loader (TPM filter, top-gene cap)

### Exports — `workflow/variant_calling.py`

- `check_tools` — Verify hisat2/samtools/bcftools/curl availability
- `find_completed_samples` — Completed-sample discovery under amalgkit output
- `find_reference_genome` — Reference FASTA discovery
- `decompress_if_needed` — gunzip wrapper (idempotent)
- `build_hisat2_index` — Index build with marker-file resume
- `download_fastq` — Local-reuse-first FASTQ acquisition (ENA fallback)
- `align_reads` — HISAT2 alignment + samtools sort/index
- `call_variants` — bcftools mpileup+call
- `filter_variants` — bcftools quality filter
- `merge_vcfs` — Population VCF merge

### Exports — `variant_stats.py`

- `parse_bcftools_stats` — Parse SN/TS records, SNPs, indels, Ti/Tv
- `compute_sample_stats` — Per-sample VCF stats to JSON
- `compute_popgen_summary` — Population summary to JSON
- `compute_allele_frequencies` — Per-site AF extraction via bcftools query

### Exports — `pipeline.py`

- `resolve_run_parameters` — CLI/config parameter resolution
- `run_pipeline` — End-to-end transcriptome SNP-calling run
