# Amplicon

Amplicon-based metagenomics analysis for 16S/ITS marker gene sequencing, providing OTU clustering, ASV denoising, and taxonomic classification.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports otu_clustering, asv_denoising, taxonomy submodules |
| `otu_clustering.py` | Greedy centroid-based OTU clustering and chimera detection |
| `asv_denoising.py` | DADA2-style ASV denoising with error model estimation |
| `taxonomy.py` | Naive Bayes and k-mer similarity ("blast") taxonomic classification |

## Key Functions

| Function | Description |
|----------|-------------|
| `otu_clustering.cluster_otus()` | Cluster sequences into OTUs at a given identity threshold |
| `otu_clustering.filter_chimeras()` | Detect and remove chimeric sequences (UCHIME-style) |
| `otu_clustering.calculate_identity()` | Compute pairwise sequence identity |
| `asv_denoising.denoise_sequences()` | Denoise amplicon sequences into ASVs |
| `asv_denoising.estimate_error_rates()` | Estimate a base-substitution error model from Phred quality scores |
| `asv_denoising.merge_paired_reads()` | Merge overlapping paired-end reads |
| `taxonomy.classify_taxonomy()` | Classify sequences to taxonomic lineages |
| `taxonomy.build_taxonomy_tree()` | Build hierarchical taxonomy tree from assignments |
| `taxonomy.calculate_confidence()` | Summarize per-rank confidence of classifications (mean, median, fraction confident) |

## Usage

```python
from metainformant.metagenomics.amplicon import otu_clustering, asv_denoising, taxonomy

otus = otu_clustering.cluster_otus(sequences, threshold=0.97)
asv_result = asv_denoising.denoise_sequences(sequences, quality_scores=qualities)
asv_sequences = {f"asv_{i}": asv.sequence for i, asv in enumerate(asv_result.asvs)}
assignments = taxonomy.classify_taxonomy(asv_sequences, reference_db, reference_taxonomy)
```

## Behavior notes

- `cluster_otus(threshold=None)` uses `META_OTU_THRESHOLD` when it is set and 0.97 otherwise; an explicit `threshold` wins. The k-mer pre-filter skips only alignments that provably cannot reach the threshold, so it never changes the result.
- `filter_chimeras(min_divergence=...)` is a margin in percentage points of identity (UCHIME `mindiv`); a negative or NaN value raises `ValueError`.
- `denoise_sequences` merges a variant into an ASV when the Poisson upper tail `poisson.sf(abundance - 1, lambda)` exceeds `omega_a`; Phred `quality_scores` set the error probability of each mismatch and the `quality` of each ASV.
- `classify_taxonomy(method="naive_bayes")` abstains (`unclassified`) at the first rank whose bootstrap confidence is below `confidence_threshold` and at every rank under it, and for queries with no known k-mer; its bootstrap is seeded (`seed=42`).
- SciPy is optional: without it the Poisson and chi-square tails fall back to exact pure-Python functions in `metagenomics/_special.py`.
