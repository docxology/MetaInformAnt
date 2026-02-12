# Amplicon

Amplicon-based metagenomics analysis for 16S/ITS marker gene sequencing, providing OTU clustering, ASV denoising, and taxonomic classification.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports otu_clustering, asv_denoising, taxonomy submodules |
| `otu_clustering.py` | Greedy centroid-based OTU clustering and chimera detection |
| `asv_denoising.py` | DADA2-style ASV denoising with error model estimation |
| `taxonomy.py` | Naive Bayes and alignment-based taxonomic classification |

## Key Functions

| Function | Description |
|----------|-------------|
| `otu_clustering.cluster_otus()` | Cluster sequences into OTUs at a given identity threshold |
| `otu_clustering.filter_chimeras()` | Detect and remove chimeric sequences (UCHIME-style) |
| `otu_clustering.calculate_identity()` | Compute pairwise sequence identity |
| `asv_denoising.denoise_sequences()` | Denoise amplicon sequences into ASVs |
| `asv_denoising.estimate_error_rates()` | Estimate position-specific error model from quality scores |
| `asv_denoising.merge_paired_reads()` | Merge overlapping paired-end reads |
| `taxonomy.classify_taxonomy()` | Classify sequences to taxonomic lineages |
| `taxonomy.build_taxonomy_tree()` | Build hierarchical taxonomy tree from assignments |
| `taxonomy.calculate_confidence()` | Bootstrap confidence estimation for classifications |

## Usage

```python
from metainformant.metagenomics.amplicon import otu_clustering, asv_denoising, taxonomy

otus = otu_clustering.cluster_otus(sequences, threshold=0.97)
asvs = asv_denoising.denoise_sequences(sequences, qualities)
assignments = taxonomy.classify_taxonomy(asvs, reference_db)
```
