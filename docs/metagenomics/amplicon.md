# Amplicon Metagenomics

16S/18S/ITS amplicon sequencing analysis including OTU clustering, ASV denoising, chimera detection, taxonomic classification, and paired-end read merging.

## Key Concepts

**OTU clustering** groups amplicon sequences at a similarity threshold (typically 97%) using a greedy centroid-based algorithm (VSEARCH/UCLUST-style). The most abundant or longest sequence becomes the centroid.

**ASV denoising** (DADA2-style) resolves amplicon data to single-nucleotide resolution by modeling sequencing errors. An error model estimates the probability that a less abundant sequence arose from a more abundant true sequence through sequencing errors.

**Chimera detection** identifies artificial hybrid sequences formed during PCR by the UCHIME algorithm, which tests whether a query sequence can be reconstructed as a combination of two more abundant parent sequences.

**Taxonomic classification** assigns taxonomy using naive Bayes k-mer classifiers (RDP-style) or alignment-based methods with bootstrap confidence estimation.

## Data Models

### `OTU` / `ClusteringResult`

OTU holds `centroid_id`, `centroid_sequence`, `member_ids`, `member_sequences`, and `size`. ClusteringResult wraps a list of OTUs with `threshold`, `total_sequences`, `num_otus`, and `otu_table`.

### `ASV` / `DenoisingResult`

ASV holds `sequence`, `abundance`, `sample_counts`, and `quality`. DenoisingResult wraps ASVs with `error_model`, `total_reads`, `num_asvs`, and `reads_removed`.

### `TaxonomyAssignment` / `TaxonomyNode`

TaxonomyAssignment holds `sequence_id`, `lineage` (list of (rank, name) tuples), `confidence` per rank, `method`, `best_hit_id`, and `best_hit_identity`. TaxonomyNode forms a hierarchical tree with `name`, `rank`, `count`, `children`, and `parent`.

## Function Reference

### `cluster_otus(sequences, threshold=0.97, abundance=None, sort_by="length", prefilter=True) -> ClusteringResult`

Greedy centroid-based OTU clustering with k-mer pre-filter for speed. Sequences sorted by abundance or length; each assigned to the best-matching centroid or starts a new OTU.

### `calculate_identity(seq1, seq2) -> float`

Pairwise sequence identity using Needleman-Wunsch global alignment with affine gap penalties. Returns fraction of identical positions (0.0--1.0).

### `filter_chimeras(sequences, reference_db=None, abundance=None, min_score=0.28) -> Dict[str, bool]`

UCHIME-style de novo chimera detection. Flags sequences that can be explained as a combination of two more abundant parent sequences. Returns dict mapping sequence IDs to True (chimeric) or False.

### `estimate_error_rates(quality_scores, sequences=None) -> ErrorModel`

Learn a sequencing error model from Phred quality scores. Computes position-specific transition rates (A->C, A->G, etc.) and quality-to-error mappings.

### `denoise_sequences(sequences, error_rates=None, quality_scores=None, omega_a=1e-40, min_abundance=1) -> DenoisingResult`

DADA2-style ASV denoising: dereplicates sequences, estimates error model, and merges likely error variants into true ASVs using abundance-weighted Poisson testing.

### `merge_paired_reads(forward, reverse, min_overlap=20, max_mismatch_ratio=0.2) -> Dict[str, str]`

Merge paired-end reads by finding the optimal overlap alignment. In conflict positions, selects the base with higher quality score.

### `classify_taxonomy(sequences, reference_db, reference_taxonomy=None, method="naive_bayes", confidence_threshold=0.8, k=8, bootstrap_n=100) -> List[TaxonomyAssignment]`

Classify sequences against a reference database. Methods: `naive_bayes` (RDP-style k-mer Bayesian classifier) or `blast` (k-mer overlap with consensus of top hits). Bootstrap resampling estimates per-rank confidence.

### `build_taxonomy_tree(classifications) -> TaxonomyNode`

Construct a hierarchical taxonomy tree from classification results for visualization and summary.

### `calculate_confidence(assignments, min_confidence=0.0) -> Dict[str, Dict]`

Summarize per-rank confidence statistics: mean, median, and fraction above threshold.

## Usage Examples

```python
from metainformant.metagenomics import (
    cluster_otus, filter_chimeras, denoise_sequences,
    merge_paired_reads, classify_taxonomy, build_taxonomy_tree,
)

# OTU clustering at 97% identity
sequences = {"s1": "ATCGATCGATCG", "s2": "ATCGATCGATCG", "s3": "TTTTAAAACCCC"}
result = cluster_otus(sequences, threshold=0.97)
print(f"{result.num_otus} OTUs from {result.total_sequences} sequences")

# Chimera detection
chimeras = filter_chimeras(sequences)
clean = {sid: seq for sid, seq in sequences.items() if not chimeras.get(sid)}

# ASV denoising
reads = {"r1": "ATCGATCG", "r2": "ATCGATCG", "r3": "ATCGATTG"}
asv_result = denoise_sequences(reads)
print(f"{asv_result.num_asvs} ASVs, {asv_result.reads_removed} merged")

# Paired-end merging
merged = merge_paired_reads(forward_reads, reverse_reads, min_overlap=20)

# Taxonomy classification
ref_db = {"r1": "ATCGATCGATCG", "r2": "GGGGCCCCTTTT"}
ref_tax = {"r1": [("domain", "Bacteria"), ("phylum", "Firmicutes")],
           "r2": [("domain", "Bacteria"), ("phylum", "Proteobacteria")]}
classifications = classify_taxonomy(sequences, ref_db, ref_tax)
tree = build_taxonomy_tree(classifications)
```

## Configuration

Environment variable prefix: `META_`

- `META_OTU_THRESHOLD` -- Override default OTU clustering threshold.

## Related Modules

- `metainformant.metagenomics.diversity` -- alpha/beta diversity metrics
- `metainformant.metagenomics.functional` -- functional annotation
- `metainformant.metagenomics.comparative` -- differential abundance
