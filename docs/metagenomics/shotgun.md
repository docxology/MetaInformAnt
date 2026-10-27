# Shotgun Metagenomics

Whole-genome shotgun sequencing analysis including de Bruijn graph assembly, contig scaffolding, assembly statistics, genome binning, taxonomic profiling, and community composition.

## Key Concepts

**De Bruijn graph assembly** extracts k-mers from reads to build a directed graph where (k-1)-mer nodes are connected by k-mer edges. Non-branching paths are collapsed into contigs. K-mers seen fewer than `min_kmer_coverage` times are discarded as likely errors; there is no tip removal or bubble popping. Multiple k-mer sizes improve assembly by capturing different repeat structures.

**Scaffolding** uses paired-end read linkage to group contigs, joining them with estimated gap sizes (represented as N's) based on insert size. Reads are anchored on both strands of every contig. Contigs in a scaffold are ordered by length and joined in the orientation they have; relative orientation is not inferred.

**Metagenomic binning** groups contigs into putative genomes (MAGs) using tetranucleotide frequency (TNF) composition profiles and/or log-transformed coverage (single sample or per-sample vectors). Bin quality is a separate step that needs single-copy marker gene annotations supplied by the caller (the package does not search the contigs for markers): `assess_bin_quality` computes CheckM-style completeness, contamination and quality score from them. Without annotations these values are `None`; they are never estimated from sequence composition or bin size.

**Taxonomic profiling** classifies reads by matching k-mers (forward strand only) against a reference index, using Lowest Common Ancestor (LCA) assignment for ambiguous k-mers (Kraken-style). K-mers shared by taxa with no common ancestor carry no taxonomic information, and a read whose matched taxa have no common ancestor counts as unclassified.

## Data Models

### Assembly

- `Contig`: `contig_id`, `sequence`, `length`, `coverage` (mean k-mer count along the contig), `gc_content` (fraction, 0--1).
- `AssemblyStats`: `total_contigs`, `total_length`, `n50`, `l50`, `n90`, `l90`, `gc_content` (fraction, 0--1), `mean_coverage`, `largest_contig`, `smallest_contig`, `mean_length`, `median_length`.
- `Scaffold`: `scaffold_id`, `contigs`, `gaps`, `sequence`, `total_length`.

### Binning

- `GenomeBin`: `bin_id`, `contig_ids`, `total_length`, `num_contigs`, `gc_content`, `mean_coverage`, `completeness`, `contamination`, `quality_score` (these three are `None` until the bin is assessed from marker gene annotations), `marker_counts` (copy number of each expected marker found in the bin; empty until assessed).
- `BinningResult`: `bins`, `unbinned_contigs`, `total_contigs`, `binned_contigs`, `high_quality_bins` (number of bins with completeness of at least 0.9 and contamination of at most 0.05; `None` when `bin_contigs` was not given `marker_genes`, because the count is then unknown; 0 when there are no bins).
- `BACTERIAL_SINGLE_COPY_MARKERS`: the default expected marker set, a tuple of 36 bacterial single-copy gene symbols (`rpoB`, `rpoC`, `rpsB`, ..., `ychF`).

### Profiling

- `TaxonProfile`: `taxon_id`, `name`, `rank`, `lineage`, `read_count`, `relative_abundance`, `kmer_count` (not filled in; stays 0).
- `CommunityProfile`: `taxa`, `total_reads`, `classified_reads`, `unclassified_reads`, `classification_rate`, `rank_summaries` (rank -> {taxon name: summed relative abundance}).
- `KmerIndex`: `kmer_to_taxon`, `taxon_lineages`, `k`, `total_kmers`, `num_taxa`.

## Function Reference

The `metainformant.metagenomics` package exposes only submodules, so import the functions from `metainformant.metagenomics.shotgun.assembly`, `metainformant.metagenomics.shotgun.binning` and `metainformant.metagenomics.shotgun.profiling`.

### Assembly

#### `assemble_contigs(reads, k_range=None, min_contig_length=200, min_kmer_coverage=2) -> List[Contig]`

Assemble reads using de Bruijn graphs at one or more k-mer sizes. Reads can be provided as `dict[str, str]` or `list[str]`. If `k_range` is None it depends on the mean read length: `[21]` below 100 bp, `[21, 33]` below 200 bp, and `[21, 33, 55]` otherwise. Every k must be at least 11 (`ValueError`); even values only log a warning. K-mers are counted on the forward strand as given, not canonicalized. Contigs from different k values are deduplicated by dropping any contig that is a substring of a longer one. The result is sorted by length, descending, with IDs renumbered `contig_000000`, `contig_000001`, and so on. The assembly is deterministic: graph start nodes are visited in sorted order, so the same reads and parameters give the same contigs with the same IDs in every process. An empty list is returned when no contig reaches `min_contig_length` (for example with very short or highly repetitive reads), and `calculate_assembly_stats` then raises `ValueError`.

#### `scaffold_contigs(contigs, paired_reads=None, insert_size=500, insert_std=100, min_links=3) -> List[Scaffold]`

Scaffold contigs using paired-end linkage. `paired_reads` maps pair IDs to `(forward, reverse)` sequence tuples; without it every contig becomes its own scaffold. Reads are anchored to contigs by exact 31-mers. Both strands of every contig are indexed, so a mate sequenced from the opposite strand (typically the reverse mate) still maps to its contig and mates can be supplied as sequenced. Inter-contig links are counted, and contigs with at least `min_links` links are joined using union-find with greedy edge selection. Contigs are listed longest first and joined in their given orientation; relative orientation is not inferred. Gaps are runs of `N` of length `max(10, insert_size - min(length of the preceding contig, insert_size))`. `insert_std` is accepted but not used.

#### `calculate_assembly_stats(contigs) -> AssemblyStats`

Compute N50, L50, N90, L90, GC content, coverage, and length distribution. Raises `ValueError` for an empty contig list.

### Binning

#### `bin_contigs(contigs, coverage=None, method="composition", n_bins=None, min_contig_length=1000, seed=42, marker_genes=None, marker_set=None) -> BinningResult`

Bin contigs (`dict` of contig ID -> sequence) into putative genomes; contigs shorter than `min_contig_length` stay unbinned. Methods: `composition` (TNF only), `coverage` (log coverage only), or `combined` (TNF plus log coverage). `coverage` maps contig IDs to one value or to a list of per-sample values; without it `combined` reduces to composition-only, and `coverage` has no features to cluster on, so the contigs typically end up in one bin. Clustering uses pure-Python k-means with k-means++ initialization (`seed`). If `n_bins` is not given it is the smaller of one bin per 3 Mbp of binned sequence and half the number of binnable contigs, raised to at least 2 and capped at the number of binnable contigs (so 1 or 2 binnable contigs give 1 or 2 bins).

The returned bins carry length, GC and mean coverage; `completeness`, `contamination` and `quality_score` are `None`, and so is `BinningResult.high_quality_bins`, unless `marker_genes` is given. With `marker_genes` (contig ID -> marker gene names found on that contig) every bin is assessed by `assess_bin_quality` (using `marker_set` as the expected markers) and `high_quality_bins` counts the bins that pass `count_high_quality_bins`.

#### `calculate_tetranucleotide_freq(sequence, normalize=True) -> List[float]`

Compute 256-dimensional TNF profile (alphabetical order, AAAA to TTTT) counting both strands. Raises `ValueError` for sequences shorter than 4 bases.

#### `refine_bins(bins, contigs, completeness_threshold=0.5, contamination_threshold=0.1, marker_genes=None, marker_set=None) -> List[GenomeBin]`

Refines bins on their marker gene completeness and contamination. `marker_genes` is required whenever `bins` is not empty (`ValueError` otherwise, since without annotations there is nothing to refine on; an empty `bins` returns `[]`). It runs `assess_bin_quality` and drops bins below `completeness_threshold`. Bins whose contamination exceeds `contamination_threshold` are handled by size: a one-contig bin is dropped, a two-contig bin is kept unchanged, and a larger bin is split in two by TNF k-means. The sub-bins are assessed with the same `marker_genes` and `marker_set` and kept when their completeness reaches `completeness_threshold` (their contamination is not checked again). Retained bins are renamed `refined_bin_0000`, `refined_bin_0001`, and so on. Merging of small bins is not implemented.

#### `assess_bin_quality(bins, contigs=None, marker_genes=None, marker_set=None) -> List[GenomeBin]`

Computes CheckM-style quality from single-copy marker gene annotations. It updates the bins in place and returns them. For a flat set of `M` expected markers:

- `completeness` = (number of expected markers present at least once in the bin) / `M`.
- `contamination` = (extra copies, the sum of `copies - 1` over the markers present) / `M`. It is 0.0 or more and exceeds 1.0 for a bin that carries several copies of most markers.
- `quality_score` = completeness - 5 * contamination (the Parks et al. convention); it is negative for a heavily contaminated bin.
- `marker_counts` holds the copy number of every expected marker found in the bin.

`marker_genes` maps a contig ID to the marker names found on that contig; a name listed twice counts as two copies, and a contig that is not in the mapping carries no markers (a warning is logged when `marker_genes` is not empty but none of its contig IDs occurs in the bins, which usually means the contig names differ). The annotations must come from a real marker search of the contigs, for example HMMER with TIGRFAM or Pfam marker profiles; this function does not search sequences, and `contigs` is accepted but unused. The expected markers `M` are `marker_set` if given and otherwise `BACTERIAL_SINGLE_COPY_MARKERS`; pass the set that matches the names in your annotations (for example a lineage-specific CheckM or BUSCO set), because annotation names outside the expected set are ignored. An empty `marker_set` raises `ValueError`.

Without `marker_genes`, `completeness`, `contamination` and `quality_score` are set to `None` for every bin and `marker_counts` is emptied; the values are never estimated from genome size or sequence composition.

#### `count_high_quality_bins(bins, min_completeness=0.9, max_contamination=0.05) -> int | None`

Counts the assessed bins with completeness of at least `min_completeness` and contamination of at most `max_contamination`. Returns `None` when any bin has no assessed quality (the count is then unknown) and 0 for an empty list.

### Profiling

#### `build_kmer_index(reference_sequences, taxonomy=None, k=31) -> KmerIndex`

Build a k-mer classification index from reference genomes (forward-strand k-mers). `taxonomy` maps reference IDs to lists of `(rank, name)` tuples; without it each reference is its own taxon. Shared k-mers are resolved to the LCA of the lineages that contain them. Raises `ValueError` for empty references or `k` below 11.

#### `profile_community(reads, database=None, reference_sequences=None, reference_taxonomy=None, k=31, min_kmer_hits=2, confidence_threshold=0.5) -> CommunityProfile`

Profile the taxonomic composition of a metagenomic sample (`reads` as `dict[str, str]` or `list[str]`) by classifying reads against a k-mer index. Pass a prebuilt `database`, or `reference_sequences` (with optional `reference_taxonomy` and `k`) to build one; `k` is ignored when `database` is given. K-mers shared by reference taxa with no common ancestor are indexed under an empty lineage and ignored when reads are classified. A read with fewer than `min_kmer_hits` informative k-mer matches is unclassified. Otherwise it is assigned to the taxon holding most of its k-mer hits if that taxon holds at least `confidence_threshold` of them, or else to the LCA of all taxa it hit; when those taxa have no common ancestor the read stays unclassified and is not counted as classified. Reports per-taxon read counts and relative abundances (fractions of classified reads) and rank-level summaries.

#### `calculate_relative_abundance(profile, rank=None, min_abundance=0.0) -> Dict[str, float]`

Extract relative abundances from a profile, optionally aggregated to a specific taxonomic rank (taxa without that rank are counted as "unclassified"). Taxa below `min_abundance` are grouped as "Other".

## Usage Examples

```python
import random

from metainformant.metagenomics.shotgun.assembly import (
    assemble_contigs, calculate_assembly_stats, scaffold_contigs,
)
from metainformant.metagenomics.shotgun.binning import (
    BACTERIAL_SINGLE_COPY_MARKERS, assess_bin_quality, bin_contigs,
)
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
print(bins.high_quality_bins)  # None: no marker gene annotations were given

# Bin quality needs marker gene annotations (contig ID -> marker names found on that contig),
# normally parsed from a marker search such as HMMER. This hand-written stand-in says the first
# contig carries every expected marker once and the second carries half of them, so the values
# below describe that stand-in, not the simulated sequences.
markers = list(BACTERIAL_SINGLE_COPY_MARKERS)
marker_genes = {contigs[0].contig_id: markers, contigs[1].contig_id: markers[:18]}
assessed = assess_bin_quality(bins.bins, marker_genes=marker_genes)
print([(b.bin_id, b.completeness, b.contamination) for b in assessed])
# [('bin_0000', 1.0, 0.0), ('bin_0001', 0.5, 0.0)]

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
