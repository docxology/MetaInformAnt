# Data Conversion Quick Reference

Convert between common bioinformatics formats with the format helpers that ship with the package.

## When to Use

Use `data_conversion` for format translation during pipeline interoperability (e.g., FASTQ → FASTA, merging per-sample VCFs) — not for initial file creation (that's `io.read_*`/`write_*`) or schema transformations (use `pandas`). There is no single `convert` module: each conversion uses the helper that lives next to its format, and the package does not provide VCF → BED, FASTQ → BAM, BAM → BigWig, BAM → CRAM or GFF3 → BED12 converters.

## Table of Contents

- [Format support](#format-support)
- [Format support](#format-support)
- [Advanced Examples](#advanced-examples)
- [Common Pitfalls](#common-pitfalls)

---

```python-snippet
from metainformant.dna.io.fastq import convert_fastq_to_fasta

# FASTQ → FASTA
convert_fastq_to_fasta("sample_R1.fastq.gz", "sample_R1.fasta")
```

## Format support

| Format | Read | Write |
|--------|------|-------|
| FASTA | `dna.sequence.core.read_fasta()` | `dna.sequence.core.write_fasta()` |
| FASTQ | `dna.io.fastq.read_fastq()` | `dna.io.fastq.write_fastq()` |
| VCF | `dna.variation.variants.parse_vcf()` | `dna.variation.variants.write_vcf()` |
| BAM/CRAM | `longread.io.bam.read_long_read_bam()` (long reads) | — |
| BED/bedGraph/bigWig | `epigenome.analysis.tracks.load_bed_track()` / `load_bedgraph_track()` / `load_bigwig_track()` | `epigenome.analysis.tracks.save_bed_track()` / `save_bedgraph_track()` |
| GFF3 | `gwas.data.genome.parse_gff3_genes()` (gene features on chromosomes it can map to integers 1-16) | — |
| Parquet | `core.io.read_parquet()` | `core.io.write_parquet()` |
| HDF5 | `h5py.File()` | `h5py.File()` |

## Advanced Examples

### Batch conversion with progress bar
```python-snippet
from metainformant.dna.io.fastq import convert_fastq_to_fasta
from tqdm import tqdm
import glob

# Convert all FASTQ files in directory to FASTA
fastq_files = glob.glob("data/fastq/*.fastq.gz")

for fastq in tqdm(fastq_files, desc="Converting FASTQ→FASTA"):
    fasta = fastq.replace(".fastq.gz", ".fasta")
    convert_fastq_to_fasta(fastq, fasta)
```
Expected output (tqdm bar):
```
Converting FASTQ→FASTA: 100%|██████████| 28/28 [00:42<00:00,  1.50s/it]
```

### Merge multiple VCFs (cohort joint-calling)
```python-snippet
from pathlib import Path
from metainformant.gwas.data.vcf_utils import count_variants, merge_vcfs

# Merge bgzipped + indexed per-sample VCFs into a cohort VCF (runs `bcftools merge`)
merged = merge_vcfs(
    [Path("samples/sample1.vcf.gz"), Path("samples/sample2.vcf.gz")],
    "cohort.vcf.gz",
    threads=4,
)
print(f"Merged {count_variants(merged):,} variants into {merged}")
```
Expected output:
```
Merged 12,456,789 variants into cohort.vcf.gz
```

### FASTQ quality filtering (conversion + filter)
```python-snippet
from metainformant.quality.io.fastq import filter_reads

# Read, filter low-quality reads, write new FASTQ
stats = filter_reads(
    "raw_R1.fastq.gz",
    "filtered_R1.fastq.gz",
    min_quality=20,          # mean Phred score per read
    min_length=50,
    max_n_bases=0,           # drop reads that contain N
)
print(f"Kept {stats['passed_reads']:,} / {stats['total_reads']:,} reads ({stats['pass_rate']:.1f}%)")
print(f"Removed {stats['filtered_reads']:,} reads (quality < 20, length < 50 or contains N)")
```
Expected output:
```
Kept 48,231,045 / 62,451,283 reads (77.2%)
Removed 14,220,238 reads (quality < 20, length < 50 or contains N)
```

## Common Pitfalls

| Problem | Likely Cause | Fix |
|---------|-------------|-----|
| `bcftools` (or `bgzip`/`tabix`) not found when merging VCFs | `merge_vcfs`, `count_variants` and `bgzip_and_index` shell out to these tools | Install them and make sure they are on `PATH` |
| `merge_vcfs` fails on plain `.vcf` inputs | `bcftools merge` needs bgzipped, indexed inputs | Run `bgzip_and_index("sample.vcf")` from `metainformant.gwas.data.vcf_utils` on each input first |
| `Chromosome name mismatch` in track files | UCSC ("chr1") vs NCBI ("1") naming | Normalize names: `normalize_chromosome_name()` in `metainformant.gwas.data.genome` returns the integer chromosome number (1-16) or `None` if unrecognized |

---

**Related:** [Full I/O docs](../core/io.md) | [Variant tools](../dna/variants.md) | [FASTQ quality](../quality/fastq.md)
