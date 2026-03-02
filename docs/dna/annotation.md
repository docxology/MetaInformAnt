### DNA: Feature Annotation

The `metainformant.dna.annotation` package provides tools for predicting genes, annotating sequence features, and predicting the functional impact of variants.

```mermaid
flowchart TD
  seq[DNA Sequence] --> find[Gene Finding/Prediction]
  seq --> annot[Feature Annotation]
  find --> orf[ORFs]
  annot --> cpg[CpG Islands]
  annot --> mask[Repeat Masking]
  annot --> reg[Regulatory Elements]
  annot --> splice[Splice Sites]
  
  var[Variants] --> func[Functional Annotation]
  func --> class[Effect Classification]
  func --> impact[Impact Prediction]
```

---

## Gene Prediction (`dna.annotation.gene_prediction`)

`gene_prediction.py` integrates tools from `gene_finding.py` and `gene_annotation.py`.

### `predict_orfs`

Find all open reading frames across all 6 reading frames.

```python
from metainformant.dna.annotation import predict_orfs

orfs = predict_orfs("ATGCGTAAATGATAG", min_length=10)
# Returns list of dicts with frame, start, end, length, sequence, and translated protein
```

### `annotate_coding_regions`

Classify regions of a sequence as coding or non-coding using hexamer scoring.

```python
from metainformant.dna.annotation import annotate_coding_regions

result = annotate_coding_regions(sequence)
# Returns {"regions": [...], "coding_fraction": 0.45, ...}
```

---

## Sequence Feature Annotation

### `find_regulatory_elements`

Search for regulatory elements (promoters, enhancers) using IUPAC-aware consensus motifs.

```python
from metainformant.dna.annotation import find_regulatory_elements

elements = find_regulatory_elements(sequence, elements=["TATA_box", "CAAT_box"])
# Returns list of matched elements with their positions and scores
```

### `annotate_cpg_islands`

Find CpG islands using sliding window criteria (GC content, observed/expected CpG ratio).

```python
from metainformant.dna.annotation import annotate_cpg_islands

islands = annotate_cpg_islands(sequence, min_length=200, min_gc=0.5, min_obs_exp=0.6)
```

### `compute_codon_usage`

Calculate codon frequencies, RSCU, and the Codon Adaptation Index (CAI).

```python
from metainformant.dna.annotation import compute_codon_usage

usage = compute_codon_usage(coding_sequence)
# Returns {"codon_frequencies": {...}, "rscu": {...}, "cai": 0.85, ...}
```

### `mask_repeats`

Identify and soft-mask tandem repeats and simple sequence repeats.

```python
from metainformant.dna.annotation import mask_repeats

masked = mask_repeats(sequence, min_length=10)
# Returns sequence with lowercase repeats: "ACGTacgtacgtacgtACGT"
```

### `find_splice_sites`

Identify potential GT-AG splice donor and acceptor sites.

```python
from metainformant.dna.annotation import find_splice_sites

sites = find_splice_sites(sequence)
```

---

## Functional Variant Annotation (`dna.annotation.functional`)

### `annotate_variants`

Classify variants against a reference coding sequence (e.g., synonymous, nonsynonymous, nonsense, frameshift).

```python
from metainformant.dna.annotation.functional import annotate_variants

annotated = annotate_variants(variants_list, reference_sequence)
```

### `predict_variant_impact`

Score variant impact using BLOSUM62 matrices, Grantham physicochemical distances, and conservation.

```python
from metainformant.dna.annotation.functional import predict_variant_impact

impact = predict_variant_impact(annotated_variant, context={"conservation": 0.8})
# Returns {"impact_score": 0.85, "impact_category": "probably_damaging", ...}
```

### `compute_conservation_score`

Compute Shannon entropy-based conservation at an alignment position.

```python
from metainformant.dna.annotation.functional import compute_conservation_score

score = compute_conservation_score(alignment_list, position=10)
# Returns 0.0 (highly variable) to 1.0 (perfectly conserved)
```
