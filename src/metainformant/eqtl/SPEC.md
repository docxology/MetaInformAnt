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

- `create_synthetic_data(n_genes=100, n_variants=500, n_samples=50, seed=42)` — Synthetic expression/genotype matrices with known cis-eQTL effects (about 30% of genes get an effect on their first cis variant). Gene IDs are `LOC<n>`, variant IDs `rs<n>`; `n_variants` must be at least `5 * n_genes`
- `parse_gene_positions(gene_ids)` — Gene table `[gene_id, chrom, tss_position]` from kallisto target IDs. Only the chromosome accession (`NC_`, `NW_` or `NT_`; `unknown` when absent) is read from the ID. The TSS positions are synthetic: `1_000_000 + 10_000 * index`, in the order of the IDs, because real coordinates need a GFF/GTF annotation
- `create_synthetic_genotypes(sample_ids, gene_positions, variants_per_gene=1, seed=42)` — Synthetic dosages matched to gene positions (MAF about 0.25). Variant IDs are `var_<full gene id>_<i>`, built from the whole gene ID so they are unique whenever the gene IDs are; a repeated gene ID in `gene_positions` raises `ValueError`
- `load_real_expression_data(quant_dir, max_samples=100, min_tpm=1.0, max_genes=1000)` — Real quantification loader (TPM filter, top-gene cap). The genotypes of an analysis built on it are synthetic unless matched genotype data are supplied

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
