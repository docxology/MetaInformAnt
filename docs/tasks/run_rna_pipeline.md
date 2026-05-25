# RNA-seq Pipeline Quick Reference

Run amalgkit RNA-seq analysis on one or many species. Integrates ENA retrieval, alignment (STAR/hisat2), quantification (Kallisto), and differential expression (DESeq2).

## When to Use

Choose `run_rna_pipeline` for bulk RNA-seq analyses across multiple species with consistent cross-species normalization—not for single-cell RNA-seq (use `singlecell/` module) or long-read transcriptomics (use `longread/`).

## Table of Contents

- [Single-Species Quick Run](#single-species-quick-run)
- [Multi-Species Parallel](#multi-species-parallel)
- [Key Output Files](#key-output-files)
- [Common Flags](#common-flags)
- [Monitoring](#monitoring)
- [Expected Output](#expected-output)
- [Common Pitfalls](#common-pitfalls)

---

```bash
# 1. Install amalgkit (if not via metainformant)
uv pip install git+https://github.com/kfuku52/amalgkit

# 2. Prepare reference genome
python3 scripts/cloud/prep_genomes.py \
  --accession GCA_003254395.2 \
  --output data/refs/amellifera/

# 3. Place FASTQ files in species-specific directory
mkdir -p data/amellifera/fastq/
# cp sample_R1.fastq.gz sample_R2.fastq.gz data/amellifera/fastq/

# 4. Write simple config
cat > config/amellifera.yaml <<EOF
species: Apis mellifera
accession: GCA_003254395.2
input_dir: data/amellifera/fastq/
output_dir: output/amellifera/
threads: 16
reference: data/refs/amellifera/GCA_003254395.2.fna
EOF

# 5. Run amalgkit
python3 scripts/rna/run_amalgkit_single.py --config config/amellifera.yaml

# 6. Monitor
tail -f output/amellifera/logs/pipeline.log
```

## Multi-Species Parallel (28 Hymenoptera)

```bash
# Edit species list file
cat > config/hymenoptera_species.txt <<EOF
GCA_003254395.2  Apis mellifera
GCA_002443255.1  Apis cerana
GCA_018340665.1  Apis dorsata
# ... add 25 more species ...
EOF

# Run orchestrated pipeline (80 workers on GCP)
python3 scripts/rna/orchestrate_species.py \
  --species-list config/hymenoptera_species.txt \
  --output output/hymenoptera_all/ \
  --cloud  # Optional: auto-deploy GCP for large runs
```

## Key Output Files

| File | Description |
|------|-------------|
| `counts.tsv` | Gene expression matrix (genes × samples) |
| `deseq2_results.tsv` | Differential expression statistics |
| `pca.tsv` + `pca.png` | Principal component analysis |
| `gene_lengths.tsv` | Transcript length annotation |
| `qc_report.html` | MultiQC aggregate quality report |
| `alignment_stats.tsv` | Read mapping rates |

## Common Flags

```bash
--config FILE      # YAML config (required)
--threads N        # CPU threads (default: detect from system)
--dry-run          # Validate config without executing
--resume           # Continue from checkpoint if interrupted
--no-docker        # Run locally without Docker container
--cloud            # Deploy to GCP (see deploy_gcp.py)
```

## Monitoring

```bash
# Check status of single run
python3 scripts/rna/status.py --output-dir output/amellifera/

# Stream logs (GCP deployment)
python3 scripts/cloud/deploy_gcp.py logs --project PROJECT_ID

# Dashboard (requires hermes agent)
hermes logs --session rna-amellifera-20250426
```

## Advanced Examples

### Custom reference with alternative spliced isoforms
```bash
# Build custom genome index with GTF annotations
python3 scripts/rna/build_reference.py \
  --fasta custom_genome.fna \
  --gtf annotation.gtf \
  --output data/refs/custom/ \
  --threads 16

# Run pipeline with custom reference
python3 scripts/rna/run_amalgkit_single.py \
  --config config/custom.yaml \
  --reference data/refs/custom/
```

### Differential expression with multiple contrasts
```python-snippet
from metainformant.rna import deseq2

# Load counts matrix
counts = deseq2.load_counts("output/amellifera/counts.tsv")

# Define multi-factor design
design = {
    'condition': ['control', 'treatment', 'control', 'treatment'],
    'batch': ['batch1', 'batch1', 'batch2', 'batch2']
}

# Run DESeq2 with batch correction
results = deseq2.run(
    counts=counts,
    design=design,
    contrast=('condition', 'treatment', 'control')
)
results.to_csv("deseq2_batch_corrected.tsv", sep='\t')
```

### Cross-species gene ortholog mapping
```python-snippet
from metainformant.rna import orthologs

# Map orthologs across 28 Hymenoptera species
ortho_map = orthologs.map_across_species(
    species_list="config/hymenoptera_species.txt",
    gene_family=" honeybee_genes",
    output="ortholog_matrix.tsv"
)
print(f"Mapped {len(ortho_map)} gene families across {ortho_map.n_species} species")
```
Expected output:
```
Mapped 15634 gene families across 28 species
```

## Expected Output

### Single-species pipeline log (tail)
```
[2026-04-26 14:02:15] Starting RNA-seq pipeline for Apis mellifera
[2026-04-26 14:02:16] Reference: data/refs/amellifera/GCA_003254395.2.fna
[2026-04-26 14:02:17] Found 8 FASTQ pairs in data/amellifera/fastq/
[2026-04-26 14:02:20] Indexing genome (STAR)... done in 4m23s
[2026-04-26 14:07:01] Aligning sample_001... 87.3% mapped (45.2M reads)
[2026-04-26 14:12:33] Quantifying with Kallisto... done
[2026-04-26 14:15:44] Building gene counts matrix: 18654 genes × 8 samples
[2026-04-26 14:18:02] Running DESeq2... normalized counts written
[2026-04-26 14:20:14] PCA plot: output/amellifera/pca.png
[2026-04-26 14:21:03] Pipeline complete: output/amellifera/
```

### MultiQC summary report excerpt
```
# MultiQC Report

## General Stats

| Sample | % Aligned | % MQ0 | % Dups | Mean Cov |
|--------|-----------|-------|--------|----------|
| S1_R1  | 87.3%     | 1.2%  | 8.7%   | 42.1x    |
| S2_R1  | 85.9%     | 1.5%  | 9.2%   | 39.8x    |
| S3_R1  | 90.1%     | 0.9%  | 7.5%   | 48.3x    |

## Alignment Rates
- Average mapping rate: 87.8% (range: 82.1% - 92.4%)
- rRNA contamination: < 0.3% (acceptable)
```

### Cross-species normalized counts preview
```
> head output/hymenoptera_all/cross_species_counts.tsv
gene_id           Amellifera  Adorsata  Acerana  ... (26 more)
Amel_004123       1245.3      1189.7    1102.4   ...
Amel_007890       890.1       892.4     745.2    ...
...
```

## Common Pitfalls

| Problem | Likely Cause | Fix |
|---------|-------------|-----|
| "Zero reads aligned for all samples" | Reference genome accession mismatch (wrong species) | Verify NCBI accession in config matches sample taxonomy; re-run `prep_genomes.py` with correct accession |
| `OutOfMemoryError: Killed` during alignment | STAR/hisat2 needs >30GB RAM for large genomes | Reduce `--threads`, switch to hisat2 (lower memory), or use `--cloud` with high-memory VM |
| "Kallisto index not found" | Reference index not built or wrong path | Run `kallisto index -i index.kallisto -g transcripts.fa`; check `reference` path in config |
| DESeq2 fails: "size factors are zero" | All genes zero in some samples (failed alignment) | Check alignment logs; exclude samples with < 20% mapping rate |
| Slow (>24h) on 28 species | Running locally on laptop; not using cloud | Add `--cloud` flag; use preemptible VMs; parallelize across species (each species runs independently) |
| "Missing fastq files" error | FASTQ filenames don't match pattern `*_R1.fastq.gz` | Rename files or set `fastq_pattern` in config; see `scripts/rna/README.md` |

---

**Related:** [Full RNA docs](../rna/index.md) | [Amalgkit manual](../rna/amalgkit/) | [Cloud deployment](../cloud/DEPLOYMENT.md) | [Quality control](../quality/index.md)
