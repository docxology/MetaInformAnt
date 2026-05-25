# DNA Analysis Quick Reference

Comprehensive DNA sequence operations for bioinformatics workflows.

## When to Use

Use `analyze_dna` tasks when working with genomic sequences, variant discovery, or phylogenetic analysis—not for RNA expression quantification or GWAS association testing.

## Table of Contents

- [Basic Operations](#basic-operations)
- [Sequence Statistics](#sequence-statistics)
- [Alignment](#alignment)
- [Phylogenetics](#phylogenetics)
- [Variant Calling](#variant-calling)
- [Expected Output](#expected-output)
- [Common Pitfalls](#common-pitfalls)

---

## Basic Operations

```python
from metainformant.dna import sequences

# Read FASTA
seqs = sequences.read_fasta("data/genome.fasta")
for name, seq in seqs.items():
    print(f"{name}: {len(seq)} bp, GC={sequences.gc_content(seq):.1f}%")

# Reverse complement
rc = sequences.reverse_complement(seqs["chr1"])

# Translate DNA → protein (6-frame)
for frame in range(3):
    protein = sequences.translate(seqs["chr1"][frame:])
```

## Sequence Statistics

```python-snippet
from metainformant.dna import composition

gc = composition.gc_content(seq)
at_ratio = composition.at_ratio(seq)
kmer_freq = composition.kmer_frequencies(seq, k=3)

# Dinucleotide bias
bias = composition.dinucleotide_bias(seq)
```

## Alignment

```python
from metainformant.dna import alignment

# Pairwise global alignment (Needleman-Wunsch)
aligned1, aligned2 = alignment.global_align(seq1, seq2, match=2, mismatch=-1, gap=-2)

# Local alignment (Smith-Waterman)
aligned1, aligned2, score = alignment.local_align(seq1, seq2, match=3, mismatch=-2, gap=-2)

# Multiple sequence alignment (ClustalW-style)
msa = alignment.multiple_align([seq1, seq2, seq3], method="clustalw")

# Compute distance matrix
dist_matrix = alignment.pairwise_distances(msa, metric="identity")
```

## Phylogenetics

```python
from metainformant.dna import phylogeny

# Build distance tree (UPGMA)
tree = phylogeny.upgma(dist_matrix)
tree.plot("tree.pdf")

# Build maximum likelihood tree (IQ-TREE wrapper)
tree = phylogeny.maximum_likelihood(msa, model="GTR+G")

# Bootstrap support
bootstraps = phylogeny.bootstrap(msa, n=1000)
```

## Variant Calling

```python-snippet
from metainformant.dna import variants

# Call SNPs from BAM
variants.call_snps(
    bam_path="sample.bam",
    reference="genome.fasta",
    output="variants.vcf"
)

# Filter VCF
variants.filter_vcf(
    input_vcf="raw.vcf",
    output_vcf="filtered.vcf",
    min_quality=30,
    min_depth=10
)

# Annotate with VEP
variants.annotate_vep("filtered.vcf", "annotated.vcf")
```

## Advanced Examples

### K-mer spectrum and repeat detection
```python-snippet
from metainformant.dna import composition

# Compute 5-mer frequencies
kmer_counts = composition.kmer_frequencies(seq, k=5)
print(f"Top 10 5-mers: {sorted(kmer_counts.items(), key=lambda x: -x[1])[:10]}")

# Detect simple repeats (microsatellites)
repeats = composition.find_repeats(seq, min_length=10)
for repeat in repeats[:5]:
    print(f"Repeat: {repeat.unit} × {repeat.count} at {repeat.position}")
```
Expected output:
```
Top 10 5-mers: [('AAAAA', 1243), ('TTTTT', 1187), ('GCGCG', 892), ...]
Repeat: CA × 12 at 14567
Repeat: AT × 8 at 89012
```

### Variant effect prediction with SnpEff
```python-snippet
from metainformant.dna import effects

# Predict functional impacts
impacts = effects.predict_effects(
    vcf="filtered.vcf",
    genome_build="GRCh38",
    annotation_db=" RefSeq"
)
impacts.to_csv("variant_effects.tsv", sep='\t')
print(f"Annotated {len(impacts)} variants")
```
Expected output:
```
Annotated 24793 variants
```

### Population genetics statistics (Tajima's D, Fst)
```python
from metainformant.dna import population

# Load multiple population VCFs
pop1 = population.load_vcf("pop1.vcf")
pop2 = population.load_vcf("pop2.vcf")

# Compute Fst between populations
fst = population.weir_cockerham_fst(pop1, pop2)
print(f"Mean Fst across genome: {fst.mean():.3f}")

# Tajima's D per window
tajima = population.tajima_d(seq, window_size=1000, step=500)
```
Expected output:
```
Mean Fst across genome: 0.142
```

## Expected Output

### Basic operations console log
```
> python analyze_dna.py
chr1: 248956422 bp, GC=41.2%
chr2: 242193529 bp, GC=40.1%
chrX: 156040895 bp, GC=39.3%
Reverse complement length check: 248956422 bp
Translation frame 0: 82985473 amino acids
```

### Alignment score report
```
Global alignment score: 142.5 (matches=78, mismatches=12, gaps=5)
Local alignment score: 89.0 (high-scoring segment: positions 234-567)
MSA dimensions: 8 sequences × 1245 positions
Pairwise identity matrix (min=0.82, max=0.99)
```

### Variant calling summary
```
[2026-04-26 10:15:23] Splitting BAM by chromosome...
[2026-04-26 10:17:45] Calling chr1: 23451 variants
[2026-04-26 10:22:11] Calling chr2: 19832 variants
[2026-04-26 10:26:03] Filtering: 41205 → 38912 (5.5% removed)
[2026-04-26 10:26:15] Annotating with VEP...
[2026-04-26 10:31:44] Complete: filtered.vcf → annotated.vcf (38912 variants)
```

## Common Pitfalls

| Problem | Likely Cause | Fix |
|---------|-------------|-----|
| `MemoryError` on large FASTA | Genome loaded fully into RAM | Use streaming: `sequences.read_fasta(path, streaming=True)` or process chromosome-by-chromosome |
| Alignment extremely slow | Quadratic O(n²) needleman-wunsch on megabase sequences | Use `local_align` for roughly [sub]strings; for full-genome alignment use `minimap2` wrapper instead |
| `VariantCallingError: No reference index` | Missing `.fai` or `.dict` for reference | Run `samtools faidx genome.fasta` and `picard CreateSequenceDictionary` |
| Low VCF variant count (near zero) | overly strict filters (MAF=0.5, depth=100) | Relax `--maf 0.01`, `--min-depth 5`, check BAM mapping rate first |
| `FileNotFoundError` for NCBI downloads | No internet or firewall blocks ENA | Set `export METAINFORMANT_OFFLINE=1` and use local genome mirror |

---

**Related:** [Sequence operations](../dna/sequences.md) | [Population genetics](../dna/population.md) | [Full DNA docs](../dna/index.md)
