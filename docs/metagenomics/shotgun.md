# Shotgun Metagenomics

Whole-genome shotgun sequencing analysis including de Bruijn graph assembly, contig scaffolding, assembly statistics, genome binning, taxonomic profiling, and community composition.

## Key Concepts

**De Bruijn graph assembly** extracts k-mers from reads to build a directed graph where (k-1)-mer nodes are connected by k-mer edges. Non-branching paths are collapsed into contigs. Multiple k-mer sizes improve assembly by capturing different repeat structures.

**Scaffolding** uses paired-end read linkage to order and orient contigs, joining them with estimated gap sizes (represented as N's) based on insert size.

**Metagenomic binning** groups contigs into putative genomes (MAGs) using tetranucleotide frequency (TNF) composition profiles and/or differential coverage across samples. Bins are assessed for quality using single-copy marker gene analysis (CheckM-style).

**Taxonomic profiling** classifies reads by matching k-mers against a reference index, using Lowest Common Ancestor (LCA) assignment for ambiguous k-mers (Kraken-style).

## Data Models

### Assembly

- `Contig`: `contig_id`, `sequence`, `length`, `coverage`, `gc_content`.
- `AssemblyStats`: `total_contigs`, `total_length`, `n50`, `l50`, `n90`, `l90`, `gc_content`, `mean_coverage`, `largest_contig`, `smallest_contig`.
- `Scaffold`: `scaffold_id`, `contigs`, `gaps`, `sequence`, `total_length`.

### Binning

- `GenomeBin`: `bin_id`, `contig_ids`, `total_length`, `num_contigs`, `gc_content`, `mean_coverage`, `completeness`, `contamination`, `quality_score`, `marker_counts`.
- `BinningResult`: `bins`, `unbinned_contigs`, `total_contigs`, `binned_contigs`, `high_quality_bins`.

### Profiling

- `TaxonProfile`: `taxon_id`, `name`, `rank`, `lineage`, `read_count`, `relative_abundance`.
- `CommunityProfile`: `taxa`, `total_reads`, `classified_reads`, `unclassified_reads`, `classification_rate`, `rank_summaries`.
- `KmerIndex`: `kmer_to_taxon`, `taxon_lineages`, `k`, `total_kmers`, `num_taxa`.

## Function Reference

### Assembly

#### `assemble_contigs(reads, k_range=None, min_contig_length=200, min_kmer_coverage=2) -> List[Contig]`

Assemble reads using de Bruijn graphs at multiple k-mer sizes. Reads can be provided as `dict[str, str]` or `list[str]`. Default k-mer range is [21, 33, 55], adjusted by read length. Contigs from different k values are deduplicated.

#### `scaffold_contigs(contigs, paired_reads=None, insert_size=500, min_links=3) -> List[Scaffold]`

Scaffold contigs using paired-end linkage. Maps reads to contigs via k-mer anchoring, counts inter-contig links, and joins using union-find with greedy edge selection.

#### `calculate_assembly_stats(contigs) -> AssemblyStats`

Compute N50, L50, N90, L90, GC content, coverage, and length distribution.

### Binning

#### `bin_contigs(contigs, coverage=None, method="composition", n_bins=None, min_contig_length=1000) -> BinningResult`

Bin contigs into putative genomes. Methods: `composition` (TNF only), `coverage` (differential coverage only), or `combined` (recommended). Clustering uses k-means with k-means++ initialization. Number of bins estimated from total assembly size if not specified.

#### `calculate_tetranucleotide_freq(sequence, normalize=True) -> List[float]`

Compute 256-dimensional TNF profile from both strands.

#### `refine_bins(bins, contigs, completeness_threshold=0.5, contamination_threshold=0.1) -> List[GenomeBin]`

Refine bins by assessing quality with marker genes, removing low-quality bins, and splitting contaminated bins by TNF divergence.

#### `assess_bin_quality(bins, contigs=None) -> List[GenomeBin]`

CheckM-style quality assessment: estimates completeness (fraction of 36 bacterial single-copy markers found) and contamination (fraction present in multiple copies). Quality score = completeness - 5 * contamination.

### Profiling

#### `build_kmer_index(reference_sequences, taxonomy=None, k=31) -> KmerIndex`

Build a k-mer classification index from reference genomes. Shared k-mers are resolved to their LCA.

#### `profile_community(reads, database=None, reference_sequences=None, k=31, min_kmer_hits=2, confidence_threshold=0.5) -> CommunityProfile`

Profile the taxonomic composition of a metagenomic sample by classifying reads against a k-mer index. Reports per-taxon abundances and rank-level summaries.

#### `calculate_relative_abundance(profile, rank=None, min_abundance=0.0) -> Dict[str, float]`

Extract relative abundances from a profile, optionally aggregated to a specific taxonomic rank. Low-abundance taxa can be grouped as "Other".

## Usage Examples

```python
from metainformant.metagenomics import (
    assemble_contigs, calculate_assembly_stats, scaffold_contigs,
    bin_contigs, assess_bin_quality,
    build_kmer_index, profile_community,
)

# Assembly
reads = ["ATCGATCGATCGATCG" * 5, "CGATCGATCGATCGTT" * 5]
contigs = assemble_contigs(reads, k_range=[21])
stats = calculate_assembly_stats(contigs)
print(f"N50: {stats.n50}, contigs: {stats.total_contigs}")

# Binning
contig_dict = {c.contig_id: c.sequence for c in contigs}
coverage = {c.contig_id: c.coverage for c in contigs}
bins = bin_contigs(contig_dict, coverage, method="combined", n_bins=3)
assessed = assess_bin_quality(bins.bins, contig_dict)

# Taxonomic profiling
refs = {"genome1": "ATCG" * 100}
tax = {"genome1": [("domain", "Bacteria"), ("phylum", "Firmicutes")]}
index = build_kmer_index(refs, tax, k=21)
profile = profile_community(reads, database=index)
print(f"Classified: {profile.classification_rate:.1%}")
```

## Configuration

Environment variable prefix: `META_`

## Related Modules

- `metainformant.metagenomics.amplicon` -- amplicon-based analysis
- `metainformant.metagenomics.diversity` -- community diversity metrics
- `metainformant.metagenomics.functional` -- functional annotation of genes
- `metainformant.metagenomics.comparative` -- differential abundance testing
