# Amplicon Metagenomics

16S/18S/ITS amplicon sequencing analysis including OTU clustering, ASV denoising, chimera detection, taxonomic classification, and paired-end read merging.

## Key Concepts

**OTU clustering** groups amplicon sequences at a similarity threshold (typically 97%) using a greedy centroid-based algorithm (VSEARCH/UCLUST-style). The most abundant or longest sequence becomes the centroid.

**ASV denoising** (DADA2-style) aims for single-nucleotide resolution by modeling sequencing errors. An error model estimates the probability that a less abundant sequence arose from a more abundant true sequence through sequencing errors.

**Chimera detection** identifies artificial hybrid sequences formed during PCR with a simplified UCHIME-style test, which checks whether a query sequence can be reconstructed as a combination of two more abundant parent sequences.

**Taxonomic classification** assigns taxonomy using a naive Bayes k-mer classifier (RDP-style) with bootstrap confidence estimation, or a k-mer similarity consensus of the top reference hits (`blast`, which does not run BLAST).

## Data Models

### `OTU` / `ClusteringResult`

OTU holds `centroid_id`, `centroid_sequence`, `member_ids` (the centroid is its own first member), `member_sequences`, and `size`. ClusteringResult wraps a list of OTUs (`otus`) with `threshold`, `total_sequences`, `num_otus`, and `otu_table` (centroid ID -> number of member sequences).

### `ASV` / `DenoisingResult`

ASV holds `sequence`, `abundance`, `sample_counts`, and `quality`; `denoise_sequences` fills only `sequence` and `abundance` (`sample_counts` stays empty and `quality` stays 0.0). DenoisingResult wraps ASVs (`asvs`) with `error_model`, `total_reads`, `num_asvs`, and `reads_removed`.

### `TaxonomyAssignment` / `TaxonomyNode`

TaxonomyAssignment holds `sequence_id`, `lineage` (list of (rank, name) tuples), `confidence` (dict of rank -> score), `method`, `best_hit_id`, and `best_hit_identity`. TaxonomyNode forms a hierarchical tree with `name`, `rank`, `count`, `children` (dict of name -> node), and `parent`.

## Function Reference

The `metainformant.metagenomics` package exposes only its subpackages, so import the functions from the modules that define them:

| Module | Contents |
|--------|----------|
| `metainformant.metagenomics.amplicon.otu_clustering` | `cluster_otus`, `calculate_identity`, `filter_chimeras`, `OTU`, `ClusteringResult` |
| `metainformant.metagenomics.amplicon.asv_denoising` | `denoise_sequences`, `estimate_error_rates`, `merge_paired_reads`, `ASV`, `DenoisingResult`, `ErrorModel` |
| `metainformant.metagenomics.amplicon.taxonomy` | `classify_taxonomy`, `build_taxonomy_tree`, `calculate_confidence`, `TaxonomyAssignment`, `TaxonomyNode`, `TAXONOMIC_RANKS` |

### `cluster_otus(sequences, threshold=0.97, abundance=None, sort_by="length", prefilter=True) -> ClusteringResult`

Greedy centroid-based OTU clustering. Sequences are sorted by length (default) or, with `sort_by="abundance"` and an `abundance` dict, by abundance (ties broken by length); `abundance` is ignored otherwise. Each sequence joins the best-matching centroid whose identity is at least `threshold`, or starts a new OTU. With `prefilter=True` the alignment is skipped when the shared 8-mer fraction is below `threshold - 0.10`. This is a speed heuristic, so `prefilter=False` is slower but can merge sequences the pre-filter keeps apart. Raises `ValueError` for empty input or a `threshold` outside (0, 1]. The identity step is pure Python (roughly 0.02--0.1 s per pair of 250 bp sequences, depending on the machine), so clustering suits hundreds to a few thousand sequences, not millions; `denoise_sequences` needs no alignments and is far cheaper on the same reads. `META_OTU_THRESHOLD` can override `threshold` (see Configuration).

### `calculate_identity(seq1, seq2) -> float`

Pairwise sequence identity using Needleman-Wunsch global alignment with affine gap penalties (match +2, mismatch -1, gap open -5, gap extend -1; full dynamic-programming matrices, no banding). Returns identical positions divided by alignment length including gap columns (0.0--1.0). Raises `ValueError` if either input is an empty string; `-` characters are then removed (a sequence made only of `-` is not rejected: it scores 0.0, or 1.0 against another all-`-` sequence).

### `filter_chimeras(sequences, reference_db=None, abundance=None, min_divergence=1.0, min_score=0.28) -> Dict[str, bool]`

Simplified UCHIME-style chimera detection (not the UCHIME program). Each query is tested against its two most similar parents using a chimeric model built at 10% breakpoints, and is flagged when the best score is at least `min_score`. Without `reference_db`, sequences are processed by descending `abundance` (descending length if `abundance` is None), the first two are never flagged, and parents are the non-chimeric sequences accepted so far; with `reference_db`, the reference sequences are the parents. Sequences shorter than 10 nt are never flagged. `min_divergence` is accepted but not used. Returns dict mapping sequence IDs to True (chimeric) or False.

### `estimate_error_rates(quality_scores, sequences=None) -> ErrorModel`

Learn a sequencing error model. `quality_scores` is a list with one list of Phred scores per read. If `sequences` is given with one entry per read, the 4x4 base transition matrix is estimated from each read's differences to the per-position consensus base; otherwise it is derived from the mean Phred score. The `ErrorModel` holds `transition_rates` (rows and columns A, C, G, T; one matrix for all positions), `quality_error_rates` (Phred 0--41 -> 10^(-Q/10)), and `positions_modeled` (longest read length). Raises `ValueError` for empty `quality_scores`.

### `denoise_sequences(sequences, error_rates=None, quality_scores=None, abundance=None, omega_a=1e-40, min_abundance=1, band_size=16) -> DenoisingResult`

Simplified DADA2-style ASV denoising. Dereplicates identical sequences (or takes per-sequence counts from `abundance`), drops those below `min_abundance`, and visits the rest by descending abundance. A sequence is merged into the first earlier ASV whose abundance-weighted Poisson test gives p > `omega_a`, provided the lengths differ by at most `band_size`; otherwise it becomes a new ASV. For sequences of up to 200 reads the upper-tail probability is computed as `1 - CDF`, which cannot resolve values below about 1e-16, so a rare variant can be merged even when its true tail probability is far below `omega_a`. `error_rates` is an `ErrorModel`; if None and `quality_scores` (dict of read ID -> Phred list) is given, `estimate_error_rates` is called with the sequences as well, so the transition matrix comes from per-position consensus differences and the Phred values themselves are not used; without `quality_scores` a default Phred-30 model is used. `reads_removed` counts the reads merged into other ASVs. Raises `ValueError` for empty input.

### `merge_paired_reads(forward, reverse, min_overlap=20, max_mismatch_ratio=0.2, quality_forward=None, quality_reverse=None) -> Dict[str, str]`

Merge paired-end reads by overlap. The reverse read is reverse-complemented internally, so pass it as sequenced. The largest overlap of at least `min_overlap` bases with a mismatch ratio of at most `max_mismatch_ratio` is used. At mismatching overlap positions the base with the higher Phred score is kept when `quality_forward` and `quality_reverse` (dicts of read ID -> list of Phred ints) are both given; otherwise the forward base is kept. Only pairs that merge appear in the result; raises `ValueError` if the two dicts share no read ID.

### `classify_taxonomy(sequences, reference_db, reference_taxonomy=None, method="naive_bayes", confidence_threshold=0.8, k=8, bootstrap_n=100) -> List[TaxonomyAssignment]`

Classify sequences against a reference database (`reference_db`: reference ID -> sequence). `reference_taxonomy` maps reference IDs to lists of (rank, name) tuples; if None, each lineage is parsed from the reference ID split on `;` (for example `Bacteria;Firmicutes;Bacilli`, assigned to domain, phylum, class, ... in order). Ranks are `domain, phylum, class, order, family, genus, species`; a rank missing from a reference lineage counts as `unclassified`. Methods:

- `naive_bayes`: RDP-style k-mer Bayesian classifier using `k`; per-rank confidence is the fraction of `bootstrap_n` k-mer bootstrap resamples that reproduce the best taxon at that rank. It always picks the best-scoring taxon, so a sequence unrelated to every reference can still get confidence 1.0 (against the two 32 bp toy references in `GETTING_STARTED.md`, nearly all unrelated random sequences received the first reference's lineage at confidence 1.0, because they match neither reference and ties go to the first reference listed; with unequal reference lengths the shorter one tends to win instead). The resampling is unseeded (seed Python's `random` for repeatable confidences).
- `blast`: a k-mer approach, not BLAST. It scores every reference by 11-mer Jaccard similarity and takes a score-weighted vote among the top 5 hits; confidence is the vote share, and `best_hit_id` / `best_hit_identity` are filled (`best_hit_identity` is the top hit's 11-mer Jaccard score, not a sequence identity). A sequence that shares no 11-mer with any reference is left `unclassified`. `k` and `bootstrap_n` are ignored.

Ranks with confidence below `confidence_threshold` are reported as `unclassified`. Raises `ValueError` for empty input or an unknown method.

### `build_taxonomy_tree(classifications) -> TaxonomyNode`

Construct a hierarchical taxonomy tree from classification results for visualization and summary. The root is a `TaxonomyNode` named `root` (rank `root`); a lineage stops at its first `unclassified` rank, which becomes an `unclassified_<rank>` child.

### `calculate_confidence(assignments, min_confidence=0.0) -> Dict[str, Dict]`

Summarize per-rank confidence. Returns `{rank: {"mean", "median", "fraction_confident", "n"}}`, where `fraction_confident` is the fraction of assignments with confidence of at least `min_confidence`.

## Usage Examples

```python
from metainformant.metagenomics.amplicon.asv_denoising import denoise_sequences, merge_paired_reads
from metainformant.metagenomics.amplicon.otu_clustering import cluster_otus, filter_chimeras
from metainformant.metagenomics.amplicon.taxonomy import build_taxonomy_tree, classify_taxonomy

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

# Paired-end merging: a 62 bp amplicon read as a 40 bp forward read and a
# 42 bp reverse read (20 bp overlap); the reverse read is given as sequenced
forward_reads = {"r1": "ATCGGCTAAGCTTGCAGGATCCTAGGCTAACGTTAGCCGA"}
reverse_reads = {"r1": "TAGCTAGCCGATCGATCCGTAATCGGCTAACGTTAGCCTAGG"}
merged = merge_paired_reads(forward_reads, reverse_reads, min_overlap=20)
print(merged["r1"])  # the 62 bp amplicon

# Taxonomy classification
ref_db = {"r1": "ATCGATCGATCG", "r2": "GGGGCCCCTTTT"}
ref_tax = {"r1": [("domain", "Bacteria"), ("phylum", "Firmicutes")],
           "r2": [("domain", "Bacteria"), ("phylum", "Proteobacteria")]}
classifications = classify_taxonomy(sequences, ref_db, ref_tax)
tree = build_taxonomy_tree(classifications)
```

## Configuration

The only environment variable this package reads is `META_OTU_THRESHOLD`, in `cluster_otus`. When it is set to a value that parses as a float, it replaces the `threshold` argument on every call, even when the caller passes one explicitly, and `ClusteringResult.threshold` reports the value actually used. The environment value is not range-checked (the (0, 1] check applies only to the argument); an unparsable value logs a warning and the argument is used. No other amplicon function reads environment variables (the `_ENV_PREFIX = "META_"` constants in the other modules are unused).

## Related Modules

- `metainformant.metagenomics.diversity` -- alpha/beta diversity metrics
- `metainformant.metagenomics.functional` -- functional annotation
- `metainformant.metagenomics.comparative` -- differential abundance
