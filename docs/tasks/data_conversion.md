# Data Conversion Quick Reference

Convert between common bioinformatics formats with a unified API.

## When to Use

Use `data_conversion` for format translation during pipeline interoperability (e.g., VCF → BED for bedtools, BAM → BigWig for IGV) — not for initial file creation (that's `io.read_*`/`write_*`) or schema transformations (use `pandas`).

## Table of Contents

- [Fast interconversion](#fast-interconversion)
- [Format support](#format-support)
- [Advanced Examples](#advanced-examples)
- [Expected Output](#expected-output)
- [Common Pitfalls](#common-pitfalls)

---

```python-snippet
from metainformant.core.io import convert

# VCF → BED
convert.vcf_to_bed("variants.vcf", "variants.bed")

# FASTQ → BAM
convert.fastq_to_bam(
    "sample_R1.fastq.gz",
    "sample_R2.fastq.gz",
    reference="genome.fasta",
    output="sample.bam"
)

# BAM → BigWig (coverage)
convert.bam_to_bigwig("sample.bam", "coverage.bw", scale=True)
```

## Format support

| Format | Read | Write |
|--------|------|-------|
| FASTA/FASTQ | `sequences.read_fasta()` | `sequences.write_fasta()` |
| VCF | `variants.read_vcf()` | `variants.write_vcf()` |
| BAM/CRAM | `io.read_bam()` | `io.write_bam()` |
| BED/GFF/GTF | `io.read_bam()` | `io.write_bed()` |
| Parquet | `pd.read_parquet()` | `df.to_parquet()` |
| HDF5 | `h5py.File()` | `h5py.File()` |

## Advanced Examples

### Batch conversion with progress bar
```python-snippet
from metainformant.core.io import convert
from tqdm import tqdm
import glob

# Convert all VCFs in directory to BED
vcf_files = glob.glob("data/vcf/*.vcf.gz")

for vcf in tqdm(vcf_files, desc="Converting VCF→BED"):
    bed = vcf.replace(".vcf.gz", ".bed")
    convert.vcf_to_bed(vcf, bed)
```
Expected output (tqdm bar):
```
Converting VCF→BED: 100%|██████████| 28/28 [00:42<00:00,  1.50s/it]
```

### BAM quality metrics extraction (before/after conversion)
```python-snippet
from metainformant.core.io import bam_metrics

# Compute per-sample stats during conversion
stats = bam_metrics.summarize("sample.bam")
print(f"Total reads: {stats.total_reads:,}")
print(f"Mapped: {stats.mapped_pct:.1f}%")
print(f"Duplicate rate: {stats.duplicate_pct:.1f}%")
print(f"Mean insert size: {stats.mean_insert_size:.0f} bp")

# Convert to CRAM (lossless compression) if reference provided
convert.bam_to_cram(
    "sample.bam",
    "sample.cram",
    reference="genome.fasta"
)
```
Expected output:
```
Total reads: 62,451,283
Mapped: 87.3%
Duplicate rate: 8.2%
Mean insert size: 342 bp
CRAM compression: 8.2 GB → 1.9 GB (4.3x)
```

### Merge multiple VCFs (cohort joint-calling)
```python-snippet
from metainformant.variants import merge

# Merge per-sample VCFs into cohort VCF
merged = merge.vcfs(
    inputs=["samples/sample1.vcf.gz", "samples/sample2.vcf.gz", "..."],
    output="cohort.vcf.gz",
    rules=[
        merge.Rule.PASS_ONLY,          # Only keep PASS variants
        merge.Rule.MIN_AC(2),          # Allele count ≥2
        merge.Rule.SAMPLE_ORDER_ALPHANUMERIC
    ]
)
print(f"Merged {merged.n_variants} variants from {merged.n_samples} samples")
```
Expected output:
```
Merged 12,456,789 variants from 8432 samples
Sites present in: [min=1, median=412, max=8432]
```

### GFF3 → BED12 (exon features)
```python-snippet
from metainformant.core.io import gff

# Extract all exon features as BED12
exons = gff.to_bed12(
    "annotation.gff3",
    feature_type="exon",
    attribute_map={"gene_id": "ID", "transcript_id": "Parent"}
)
exons.to_csv("exons.bed", sep='\t', index=False, header=False)
print(f"Written {len(exons)} exon records")
```
Expected output snippet:
```
Written 247,391 exon records
> head exons.bed
chr1    11873   12227   ENST00000456328.2_0 0   +   11873   12227   color=0,0,255
chr1    12613   12721   ENST00000456328.2_1 0   +   12613   12721   color=0,0,255
...
```

### FASTQ quality filtering (conversion + filter)
```python-snippet
from metainformant.core.io import fastq

# Read, filter low-quality reads, write new FASTQ
count = fastq.filter_quality(
    input_fastq="raw_R1.fastq.gz",
    output_fastq="filtered_R1.fastq.gz",
    min_quality=20,          # Phred score
    min_length=50,
    truncate_n=True          # Trim Ns from ends
)
print(f"Kept {count.kept:,} / {count.total:,} reads ({count.kept/count.total:.1%})")
```
Expected output:
```
Kept 48,231,045 / 62,451,283 reads (77.2%)
Removed 14,220,238 reads (quality < 20 or length < 50)
```

## Expected Output

### Conversion log (verbose)
```
[2026-04-26 10:00:01] Starting VCF → BED conversion
[2026-04-26 10:00:01] Input: raw.vcf.gz (12.3 GB, 4,123,456 variants)
[2026-04-26 10:00:02] Parsing header: 9 INFO fields, 96 samples
[2026-04-26 10:00:05] Converting 4,123,456 records...
[2026-04-26 10:01:23] Written: raw.bed (8.7 GB, 4,123,456 lines)
[2026-04-26 10:01:23] Compression: 12.3 GB → 8.7 GB (1.4x with bgzip)
[2026-04-26 10:01:24] Indexing BED with tabix...
[2026-04-26 10:01:25] Done: raw.bed.gz + raw.bed.gz.tbi
```

### BAM → BigWig coverage graph generation
```
[bigWig] Normalizing to 1× coverage (RPKM → CPM)
[bigWig] Windowsize: 10 bp, smoothing: none
[bigWig] Chromosome chr1: 248,956,422 bp → 24,895,642 bins
[bigWig] Writing coverage... [100.0%] 248.9 MB
[bigWig] Total: 19 chromosomes, 2.1 GB bigWig
[2026-04-26 10:05:01] Indexing with wigToBigWig...
[2026-04-26 10:05:02] Index written: coverage.bw
```

### Format capability matrix (show with `convert --help-formats`)
```
Supported conversions (59 total):

FROM → TO       | Lossless? | Requires ref | Notes
----------------|-----------|--------------|------------------------
VCF → BED       | Yes       | No           | Includes all INFO fields
VCF → BCF       | Yes       | No           | Binary VCF, smaller
BAM → CRAM      | Yes       | Yes          | Lossless with ref
BAM → BigWig    | No*       | No           | *coverage summary only
FASTA → FASTAQ  | N/A       | No           | Rewrite to canonical
GFF3 → BED12    | Partial   | No           | Keeps block features
```

### File size reduction benchmarks
```
Input format  | Size (GB) | Output format | Size (GB) | Ratio
--------------|-----------|---------------|-----------|------
VCF (text)    | 12.30     | BCF (binary)  |  3.87     | 3.18x
BAM           | 47.20     | CRAM (ref)    | 11.80     | 4.00x
WIG           |  8.40     | BigWig        |  0.62     | 13.5x
FASTQ (gzip)  | 15.60     | FASTQ (zstd)  | 12.10     | 1.29x
```

## Common Pitfalls

| Problem | Likely Cause | Fix |
|---------|-------------|-----|
| `Variant ID mismatch` after VCF → BED | BED 0-based; VCF 1-based coordinate shift | `convert.vcf_to_bed(..., one_based=True)` or subtract 1 from POS manually |
| `BGZF block overrun` on large VCF | Output file >4GB with bgzip (block limit) | Use plain gzip (`compression="gzip"`), or split output; bgzf supports indexing up to 512Gb with recent htslib |
| `CRAM reference required` error | CRAM uses reference-based encoding | Provide `reference=` arg; or use `--no-ref` to skip reference (larger file) |
| BAM to BigWig produces all zeros | Scaling factor wrong or coverage normalized incorrectly | Set `scale=False` for raw counts; check `--normalize cpm` vs `rpkm`; verify genome sizes file exists |
| `GFF3 attribute parse error` | Custom attribute format not standard | Provide `attribute_sep=';'` and `pair_sep=' '`; or preprocess with `gffread` |
| `MemoryError` during VCF compression | VCF parsed fully in memory before write | Use streaming: `convert.vcf_to_bcf(..., streaming=True)` or split VCF chromosomes first |
| Conversion silently drops contigs | Genome sizes file missing for bedGraph → BigWig | Provide `chrom_sizes="genome.chrom.sizes"` (from `samtools faidx`) |
| `Chromosome name mismatch` in track files | UCSC ("chr1") vs NCBI ("1") naming | Normalize names: `convert.normalize_chrom_names(..., style="ucsc")` |

---

**Related:** [Full I/O docs](../core/io.md) | [Variant tools](../dna/variants.md) | [FASTQ quality](../quality/fastq.md)
