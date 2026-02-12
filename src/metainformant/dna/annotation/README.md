# DNA Annotation

Gene prediction and functional annotation tools for DNA sequences, including ORF finding, regulatory element detection, variant classification, and domain identification.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `gene_prediction` and `functional` |
| `gene_prediction.py` | Re-exports from `gene_finding` and `gene_annotation` submodules |
| `gene_finding.py` | ORF prediction, genetic code tables, IUPAC handling |
| `gene_annotation.py` | Coding region annotation, CpG islands, splice sites, codon usage |
| `functional.py` | Variant classification, impact prediction, conservation scoring, domain identification |

## Key Functions

| Function | Description |
|----------|-------------|
| `predict_orfs()` | Predict open reading frames in a DNA sequence |
| `annotate_coding_regions()` | Classify sequence regions as coding or non-coding |
| `annotate_cpg_islands()` | Identify CpG islands in a DNA sequence |
| `find_regulatory_elements()` | Detect known regulatory motifs |
| `find_splice_sites()` | Predict donor/acceptor splice sites |
| `compute_codon_usage()` | Calculate codon usage statistics |
| `mask_repeats()` | Mask repetitive elements in a sequence |
| `annotate_variants()` | Classify variants (synonymous, nonsynonymous, nonsense, frameshift) |
| `predict_variant_impact()` | Score variant impact using substitution matrices |
| `compute_conservation_score()` | Compute positional conservation from alignments |
| `identify_functional_domains()` | Identify protein domains via motif scanning |

## Usage

```python
from metainformant.dna.annotation import gene_prediction, functional

orfs = gene_prediction.predict_orfs(sequence)
regions = gene_prediction.annotate_coding_regions(sequence)
variants = functional.annotate_variants(ref_seq, alt_seq, position)
impact = functional.predict_variant_impact(ref_aa, alt_aa)
```
