# Shotgun Metagenomics

Whole-genome shotgun sequencing analysis including de Bruijn graph assembly, contig scaffolding, assembly statistics, genome binning, taxonomic profiling, and community composition.

## Key Concepts

**De Bruijn graph assembly** extracts k-mers from reads to build a directed graph where (k-1)-mer nodes are connected by k-mer edges. Non-branching paths are collapsed into contigs. K-mers seen fewer than `min_kmer_coverage` times are discarded as likely errors; there is no tip removal or bubble popping. Multiple k-mer sizes improve assembly by capturing different repeat structures.

**Scaffolding** uses paired-end read linkage to group contigs, joining them with estimated gap sizes (represented as N's) based on insert size. Contigs in a scaffold are ordered by length; orientation is not inferred.

**Metagenomic binning** groups contigs into putative genomes (MAGs) using tetranucleotide frequency (TNF) composition profiles and/or log-transformed coverage (single sample or per-sample vectors). Bin quality is estimated by `assess_bin_quality`, a simplified heuristic stand-in for CheckM-style single-copy marker gene analysis: it does not search for marker genes (see below), so treat its completeness and contamination as rough placeholders.

**Taxonomic profiling** classifies reads by matching k-mers (forward strand only) against a reference index, using Lowest Common Ancestor (LCA) assignment for ambiguous k-mers (Kraken-style).

## Data Models

### Assembly

- `Contig`: `contig_id`, `sequence`, `length`, `coverage` (mean k-mer count along the contig), `gc_content` (fraction, 0--1).
- `AssemblyStats`: `total_contigs`, `total_length`, `n50`, `l50`, `n90`, `l90`, `gc_content` (fraction, 0--1), `mean_coverage`, `largest_contig`, `smallest_contig`, `mean_length`, `median_length`.
- `Scaffold`: `scaffold_id`, `contigs`, `gaps`, `sequence`, `total_length`.

### Binning

- `GenomeBin`: `bin_id`, `contig_ids`, `total_length`, `num_contigs`, `gc_content`, `mean_coverage`, `completeness`, `contamination`, `quality_score`, `marker_counts`.
- `BinningResult`: `bins`, `unbinned_contigs`, `total_contigs`, `binned_contigs`, `high_quality_bins` (never filled in by `bin_contigs` or `assess_bin_quality`; it stays 0).

### Profiling

- `TaxonProfile`: `taxon_id`, `name`, `rank`, `lineage`, `read_count`, `relative_abundance`, `kmer_count` (not filled in; stays 0).
- `CommunityProfile`: `taxa`, `total_reads`, `classified_reads`, `unclassified_reads`, `classification_rate`, `rank_summaries` (rank -> {taxon name: summed relative abundance}).
- `KmerIndex`: `kmer_to_taxon`, `taxon_lineages`, `k`, `total_kmers`, `num_taxa`.

## Function Reference

The `metainformant.metagenomics` package exposes only submodules, so import the functions from `metainformant.metagenomics.shotgun.assembly`, `metainformant.metagenomics.shotgun.binning` and `metainformant.metagenomics.shotgun.profiling`.

### Assembly

#### `assemble_contigs(reads, k_range=None, min_contig_length=200, min_kmer_coverage=2) -> List[Contig]`

Assemble reads using de Bruijn graphs at one or more k-mer sizes. Reads can be provided as `dict[str, str]` or `list[str]`. If `k_range` is None it depends on the mean read length: `[21]` below 100 bp, `[21, 33]` below 200 bp, and `[21, 33, 55]` otherwise. Every k must be at least 11 (`ValueError`); even values only log a warning. K-mers are counted on the forward strand as given, not canonicalized. Contigs from different k values are deduplicated by dropping any contig that is a substring of a longer one. The result is sorted by length, descending, with IDs renumbered `contig_000000`, `contig_000001`, and so on. An empty list is returned when no contig reaches `min_contig_length` (for example with very short or highly repetitive reads), and `calculate_assembly_stats` then raises `ValueError`.

#### `scaffold_contigs(contigs, paired_reads=None, insert_size=500, insert_std=100, min_links=3) -> List[Scaffold]`

Scaffold contigs using paired-end linkage. `paired_reads` maps pair IDs to `(forward, reverse)` sequence tuples; without it every contig becomes its own scaffold. Reads are anchored to contigs by exact 31-mers on the contig sequence as given (a mate is not reverse-complemented, so supply mates in the contigs' orientation), inter-contig links are counted, and contigs with at least `min_links` links are joined using union-find with greedy edge selection. Gaps are runs of `N` of length `max(10, insert_size - min(contig length, insert_size))`. `insert_std` is accepted but not used.

#### `calculate_assembly_stats(contigs) -> AssemblyStats`

Compute N50, L50, N90, L90, GC content, coverage, and length distribution. Raises `ValueError` for an empty contig list.

### Binning

#### `bin_contigs(contigs, coverage=None, method="composition", n_bins=None, min_contig_length=1000, seed=42) -> BinningResult`

Bin contigs (`dict` of contig ID -> sequence) into putative genomes; contigs shorter than `min_contig_length` stay unbinned. Methods: `composition` (TNF only), `coverage` (log coverage only), or `combined` (TNF plus log coverage). `coverage` maps contig IDs to one value or to a list of per-sample values; without it `combined` reduces to composition-only, and `coverage` has no features to cluster on, so the contigs typically end up in one bin. Clustering uses pure-Python k-means with k-means++ initialization (`seed`). If `n_bins` is not given it is the smaller of one bin per 3 Mbp of binned sequence and half the number of binnable contigs, raised to at least 2 and capped at the number of binnable contigs (so 1 or 2 binnable contigs give 1 or 2 bins). The returned bins carry length, GC and mean coverage only: `completeness`, `contamination` and `quality_score` stay 0.0 until `assess_bin_quality` (or `refine_bins`) is run.

#### `calculate_tetranucleotide_freq(sequence, normalize=True) -> List[float]`

Compute 256-dimensional TNF profile (alphabetical order, AAAA to TTTT) counting both strands. Raises `ValueError` for sequences shorter than 4 bases.

#### `refine_bins(bins, contigs, completeness_threshold=0.5, contamination_threshold=0.1) -> List[GenomeBin]`

Runs `assess_bin_quality` and drops bins below `completeness_threshold`. Bins whose contamination exceeds `contamination_threshold` are handled by size: a one-contig bin is dropped, a two-contig bin is kept unchanged, and a larger bin is split in two by TNF k-means, but the sub-bins are not re-assessed (their completeness is 0.0), so they are dropped unless `completeness_threshold` is 0. Retained bins are renamed `refined_bin_0000`, `refined_bin_0001`, and so on. Merging of small bins is not implemented.

#### `assess_bin_quality(bins, contigs=None) -> List[GenomeBin]`

Heuristic stand-in for CheckM-style quality assessment. It updates the bins in place and returns them, filling `completeness` (fraction of 36 bacterial single-copy marker names counted as present), `contamination` (fraction counted more than once), `marker_counts`, and `quality_score` = completeness - 5 * contamination. No marker gene is actually searched for. With `contigs`, a marker counts as present when a sequence window of the bin has enough 4-mer diversity (above 0.3 of the possible 4-mers) and 25--75% GC, which almost any non-degenerate sequence satisfies (a 3 kb bin scores completeness 1.0); markers are counted twice only for bins whose recorded `total_length` exceeds 7 Mbp. Without `contigs`, completeness is `int(36 * min(total_length / 3.5 Mbp, 1)) / 36` and contamination is 0. Use CheckM or BUSCO when quality values will be reported.

### Profiling

#### `build_kmer_index(reference_sequences, taxonomy=None, k=31) -> KmerIndex`

Build a k-mer classification index from reference genomes (forward-strand k-mers). `taxonomy` maps reference IDs to lists of `(rank, name)` tuples; without it each reference is its own taxon. Shared k-mers are resolved to the LCA of the lineages that contain them. Raises `ValueError` for empty references or `k` below 11.

#### `profile_community(reads, database=None, reference_sequences=None, reference_taxonomy=None, k=31, min_kmer_hits=2, confidence_threshold=0.5) -> CommunityProfile`

Profile the taxonomic composition of a metagenomic sample (`reads` as `dict[str, str]` or `list[str]`) by classifying reads against a k-mer index. Pass a prebuilt `database`, or `reference_sequences` (with optional `reference_taxonomy` and `k`) to build one; `k` is ignored when `database` is given. A read with fewer than `min_kmer_hits` matching k-mers is unclassified. Otherwise it is assigned to the taxon holding most of its k-mer hits if that taxon holds at least `confidence_threshold` of them, or else to the LCA of all taxa it hit. Reports per-taxon read counts and relative abundances (fractions of classified reads) and rank-level summaries.

#### `calculate_relative_abundance(profile, rank=None, min_abundance=0.0) -> Dict[str, float]`

Extract relative abundances from a profile, optionally aggregated to a specific taxonomic rank (taxa without that rank are counted as "unclassified"). Taxa below `min_abundance` are grouped as "Other".

## Usage Examples

```python
import random

from metainformant.metagenomics.shotgun.assembly import (
    assemble_contigs, calculate_assembly_stats, scaffold_contigs,
)
from metainformant.metagenomics.shotgun.binning import assess_bin_quality, bin_contigs
from metainformant.metagenomics.shotgun.profiling import build_kmer_index, profile_community

# Simulate error-free 100 bp reads (step 25) from two random 3 kb genomes
rng = random.Random(7)
genomes = {
    "genome1": "".join(rng.choices("ACGT", weights=[0.3, 0.2, 0.2, 0.3], k=3000)),
    "genome2": "".join(rng.choices("ACGT", weights=[0.2, 0.3, 0.3, 0.2], k=3000)),
}
reads = [g[i:i + 100] for g in genomes.values() for i in range(0, len(g) - 99, 25)]

# Assembly
contigs = assemble_contigs(reads, k_range=[21], min_contig_length=200)
stats = calculate_assembly_stats(contigs)
print(f"N50: {stats.n50}, contigs: {stats.total_contigs}")
scaffolds = scaffold_contigs(contigs)  # no paired reads: one scaffold per contig

# Binning
contig_dict = {c.contig_id: c.sequence for c in contigs}
coverage = {c.contig_id: c.coverage for c in contigs}
bins = bin_contigs(contig_dict, coverage, method="combined", n_bins=2)
assessed = assess_bin_quality(bins.bins, contig_dict)  # heuristic estimates, see above

# Taxonomic profiling
tax = {"genome1": [("domain", "Bacteria"), ("phylum", "Firmicutes")],
       "genome2": [("domain", "Bacteria"), ("phylum", "Proteobacteria")]}
index = build_kmer_index(genomes, tax, k=21)
profile = profile_community(reads, database=index)
print(f"Classified: {profile.classification_rate:.1%}")
```

## Configuration

This module reads no environment variables (`assembly.py`, `binning.py` and `profiling.py` define `_ENV_PREFIX = "META_"`, but no function uses it).

## Related Modules

- `metainformant.metagenomics.amplicon` -- amplicon-based analysis
- `metainformant.metagenomics.diversity` -- community diversity metrics
- `metainformant.metagenomics.functional` -- functional annotation of genes
- `metainformant.metagenomics.comparative` -- differential abundance testing
