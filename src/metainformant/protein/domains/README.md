# Protein Domains

Protein domain detection, scanning, and family classification. Identifies structural domains (signal peptides, transmembrane helices, coiled-coils, zinc fingers, leucine zippers) and classifies proteins into families based on domain composition.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Re-exports `classification` and `detection` modules |
| `detection.py` | Domain scanning, signal peptide, transmembrane, and motif detection |
| `classification.py` | Family classification, domain similarity, and enrichment analysis |

## Key Functions

| Function | Description |
|----------|-------------|
| `detection.scan_domains()` | Scan sequence against domain profiles by pattern; a hit has a heuristic `score` and `e_value` None (a pattern hit has no score distribution) |
| `detection.detect_signal_peptide()` | Predict N-terminal signal peptides |
| `detection.predict_transmembrane()` | Predict transmembrane helices via hydrophobicity |
| `detection.find_coiled_coils()` | Identify coiled-coil regions using heptad scoring |
| `detection.identify_zinc_fingers()` | Detect C2H2 zinc finger motifs |
| `detection.detect_leucine_zipper()` | Find leucine zipper dimerization motifs |
| `detection.compute_domain_architecture()` | Summarize full domain architecture |
| `classification.classify_protein_family()` | Classify protein into known family by domains |
| `classification.compute_domain_similarity()` | Compute pairwise domain architecture similarity |
| `classification.cluster_by_domain()` | Cluster proteins by shared domain content |
| `classification.domain_enrichment()` | One-sided Fisher test of domain over-representation; proteins are counted once per domain |

## Usage

```python
from metainformant.protein.domains import detection, classification

domains = detection.scan_domains(sequence, profiles)
family = classification.classify_protein_family(domains, sequence)
```
