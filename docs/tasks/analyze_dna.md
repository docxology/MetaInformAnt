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
from metainformant.dna.sequence.core import gc_content, read_fasta, reverse_complement
from metainformant.dna.translation import six_frame_translation, translate_dna

# Read FASTA ({sequence_id: sequence})
seqs = read_fasta("data/genome.fasta")
for name, seq in seqs.items():
    print(f"{name}: {len(seq)} bp, GC={gc_content(seq):.1%}")

# Reverse complement
rc = reverse_complement(seqs["chr1"])

# Translate DNA → protein (the three forward frames)
for frame in range(3):
    protein = translate_dna(seqs["chr1"][frame:])

# Or all six frames at once: {"+1": ..., "+2": ..., "+3": ..., "-1": ..., "-2": ..., "-3": ...}
frames = six_frame_translation(seqs["chr1"])
```

## Sequence Statistics

```python-snippet
from metainformant.dna.sequence import composition, kmer

gc = composition.gc_content(seq)
base_freq = composition.nucleotide_frequencies(seq)  # {"A": ..., "T": ..., "G": ..., "C": ...}
at_content = base_freq["A"] + base_freq["T"]
kmer_freq = kmer.kmer_frequencies(seq, k=3)  # normalised frequencies that sum to 1

# Dinucleotide frequencies
dinuc = composition.dinucleotide_frequencies(seq)
```

## Alignment

```python
import numpy as np

from metainformant.dna.alignment import distances, msa, pairwise

# Pairwise global alignment (Needleman-Wunsch); returns an AlignmentResult
global_result = pairwise.global_align(seq1, seq2, match=2, mismatch=-1, gap=-2)
aligned1, aligned2 = global_result.seq1_aligned, global_result.seq2_aligned
stats = pairwise.alignment_statistics(global_result)  # identity, matches, mismatches, gaps, ...

# Local alignment (Smith-Waterman)
local_result = pairwise.local_align(seq1, seq2, match=3, mismatch=-2, gap=-2)

# Multiple sequence alignment: {name: sequence} in, {name: aligned sequence} out.
# method="auto" uses the first of MUSCLE, MAFFT, Clustal Omega or ClustalW found on PATH,
# and a built-in progressive aligner otherwise.
aligned = msa.align_msa({"s1": seq1, "s2": seq2, "s3": seq3}, method="auto")

# Compute distance matrix (pandas DataFrame; methods: jukes_cantor, kimura, p_distance, hamming)
dist_matrix = distances.distance_matrix(aligned, method="p_distance")
identity = distances.sequence_identity_matrix(aligned)  # 1 - p_distance

print(f"Global alignment score: {global_result.score} "
      f"(matches={stats['matches']}, mismatches={stats['mismatches']}, gaps={stats['gaps']})")
print(f"Local alignment score: {local_result.score} (start positions: {local_result.start_positions})")
print(f"MSA dimensions: {len(aligned)} sequences × {len(next(iter(aligned.values())))} positions")
off_diagonal = identity.to_numpy()[np.triu_indices(len(identity), k=1)]
print(f"Pairwise identity (min={off_diagonal.min():.2f}, max={off_diagonal.max():.2f})")
```

## Phylogenetics

```python
from metainformant.dna import phylogeny

# Trees are built directly from {name: aligned sequence} and returned as a dictionary
id_to_seq = aligned

# Build distance tree (UPGMA) or neighbor-joining tree
tree = phylogeny.upgma_tree(id_to_seq)
nj_tree = phylogeny.neighbor_joining_tree(id_to_seq)

# Newick string (to_newick) and back (from_newick)
newick = phylogeny.to_newick(tree)
with open("tree.nwk", "w") as handle:
    handle.write(newick + "\n")

# Bootstrap support: returns a new tree with a "bootstrap" value (0-100) on each internal node.
# The input tree (nj_tree) is not modified.
supported = phylogeny.bootstrap_support(nj_tree, id_to_seq, n_replicates=100, method="nj")
support = {node: kids["bootstrap"] for node, kids in supported.items() if isinstance(kids, dict) and "bootstrap" in kids}
```

Maximum-likelihood tree inference (for example an IQ-TREE wrapper) is not part of the package; build those trees with an external tool and read the Newick file back with `phylogeny.from_newick`. To draw a tree, parse the Newick string with Bio.Phylo and pass it to `metainformant.visualization.genomics.trees.plot_phylo_tree` (see [visualize_results.md](visualize_results.md)).

## Variant Calling

```python-snippet
from metainformant.dna.variation import calling, variants

# Call SNPs from pileup data: one dict per position with per-base counts
# (the package does not read BAM files itself; see the note below)
pileup_data = [
    {"chrom": "chr1", "pos": 100, "ref": "C", "depth": 40, "bases": {"C": 20, "A": 20}},
    {"chrom": "chr1", "pos": 200, "ref": "G", "depth": 12, "bases": {"G": 6, "T": 6}},
]
calls = calling.call_variants_pileup(pileup_data, min_depth=10, min_qual=20.0, min_alt_freq=0.2)

# Filter the calls (depth, quality, strand bias)
passed = calling.filter_variants(calls, min_depth=10, min_qual=30.0)

# Parse, filter and write a VCF (a path ending in .gz, such as raw.vcf.gz, is read and written as gzip)
vcf_data = variants.parse_vcf("raw.vcf")
vcf_data = variants.filter_variants_by_quality(vcf_data, min_qual=30)
vcf_data = variants.filter_variants_by_maf(vcf_data, min_maf=0.01)
variants.write_vcf(vcf_data, "filtered.vcf")

print(f"Called {len(calls)} variants; {len(passed)} passed filtering")
```

For BAM input, `metainformant.eqtl.workflow.variant_calling.call_variants(bam_path, ref_fasta, output_vcf)` runs `bcftools mpileup | bcftools call` (it needs `samtools` and `bcftools` on `PATH` and an indexed BAM). The package has no VEP wrapper; for variant consequences use `predict_variant_effect` below or `metainformant.dna.annotation.functional.annotate_variants`.

## Advanced Examples

### K-mer spectrum and repeat detection
```python-snippet
from metainformant.dna.sequence import kmer

# Count 5-mers
kmer_counts = kmer.count_kmers(seq, k=5)
print(f"Top 10 5-mers: {sorted(kmer_counts.items(), key=lambda x: -x[1])[:10]}")

# Detect simple repeats (microsatellites): (unit, start, end, repeat_count) tuples, end exclusive.
# Homopolymer runs (unit length 1) are reported too; keep units of two or more bases here.
repeats = [r for r in kmer.find_microsatellites(seq, min_repeats=4, max_unit_length=6) if len(r[0]) > 1]
for unit, start, end, count in repeats[:5]:
    print(f"Repeat: {unit} × {count} at {start}")
```
Example output for a 620 bp sequence with a `CA` repeat at 400 and an `AT` repeat at 620 (values depend on your sequence):
```
Top 10 5-mers: [('ACACA', 10), ('CACAC', 8), ('ATATA', 6), ('TATAT', 6), ('CACAT', 5), ...]
Repeat: CA × 10 at 400
Repeat: AT × 8 at 620
```

### Variant effect prediction
The package does not wrap SnpEff; `predict_variant_effect` classifies a variant against a coding sequence you supply.
```python-snippet
from metainformant.dna.variation.variants import predict_variant_effect

# Predict the protein-level consequence of a single variant
cds = "ATGATCGAA"
effect = predict_variant_effect(ref="A", alt="G", position=0, coding_sequence=cds)
print(f"{effect['effect_type']}: {effect['original_aa']} -> {effect['mutated_aa']}")
```
Expected output:
```
missense: M -> V
```

### Population genetics statistics (Tajima's D, Fst)
```python
from metainformant.dna import population
from metainformant.dna.population import analysis

from metainformant.dna.sequence.core import read_fasta

# Aligned sequences (equal-length strings) for each population; Tajima's D needs at least 4 sequences
pop1_seqs = list(read_fasta("data/pop1_aligned.fasta").values())
pop2_seqs = list(read_fasta("data/pop2_aligned.fasta").values())

# Hudson's Fst between populations
fst = population.hudson_fst(pop1_seqs, pop2_seqs)
print(f"Hudson's Fst between populations: {fst:.3f}")

# Tajima's D for one set of aligned sequences, and with its two-sided p-value
tajima = population.tajimas_d(pop1_seqs)
d, p_value = analysis.calculate_tajima_d(pop1_seqs)

# There is no windowed variant: slice the alignment into windows and call tajimas_d per window
window_size, step = 1000, 500
tajima_windows = [
    (start, population.tajimas_d([s[start:start + window_size] for s in pop1_seqs]))
    for start in range(0, len(pop1_seqs[0]) - window_size + 1, step)
]
```
Example output for two simulated populations of 12 sequences (3 kb):
```
Hudson's Fst between populations: 0.002
```

## Expected Output

### Basic operations console log
Example for a FASTA file with two records (lengths and GC depend on your data):
```
> python analyze_dna.py
chr1: 900 bp, GC=51.7%
chr2: 600 bp, GC=49.5%
```

### Alignment score report
Example for three 120 bp sequences, the second and third differing from the first at up to 8 and 15 positions:
```
Global alignment score: 216.0 (matches=112, mismatches=8, gaps=0)
Local alignment score: 320.0 (start positions: (0, 0))
MSA dimensions: 3 sequences × 120 positions
Pairwise identity (min=0.82, max=0.93)
```

### Variant calling summary
```
Called 2 variants; 2 passed filtering
```

## Common Pitfalls

| Problem | Likely Cause | Fix |
|---------|-------------|-----|
| `MemoryError` on large FASTA | `read_fasta` loads the whole file into a dict | There is no streaming reader; split the FASTA per chromosome first, or iterate over the file yourself |
| Alignment extremely slow | `global_align` and `local_align` fill an n×m dynamic-programming matrix in pure Python | Use them for short sequences; for genome-scale alignment use an external aligner such as minimap2 outside the package |
| Low VCF variant count (near zero) | Overly strict filters | Relax `min_qual` and `min_maf` in `filter_variants_by_quality` / `filter_variants_by_maf` and `min_depth` in `calling.filter_variants`; check the BAM mapping rate first |
| `UnicodeDecodeError` from `variants.parse_vcf` | The file is gzip-compressed but its name does not end in a lowercase `.gz` (for example `calls.vcf.bgz`); `parse_vcf` chooses gzip from the suffix | Rename the file so that it ends in `.gz`, or decompress it first |
| `NCBIContactError: NCBI contact is required` | NCBI E-utilities clients need a contact email | Pass `email=...` or `export NCBI_EMAIL=you@example.org` |

---

**Related:** [Sequence operations](../dna/sequences.md) | [Population genetics](../dna/population.md) | [Full DNA docs](../dna/index.md)
