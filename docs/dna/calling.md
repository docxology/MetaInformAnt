### DNA: Variant Calling

The `metainformant.dna.variation.calling` module provides tools for calling variants from pileup data, genotyping using Bayesian models, filtering, and merging call sets.

```mermaid
flowchart TD
  pileup[Pileup Data] --> call[call_variants_pileup]
  call --> gt[Bayesian Genotyping]
  gt --> filter[filter_variants]
  filter --> merge[merge_variant_calls\nConsensus]
  merge --> annot[Context Annotation]
  merge --> stats[Variant Stats]
```

---

## Variant Calling Pipeline

### `call_variants_pileup`

Call SNP variants from pileup data using a simple Bayesian genotyping model. Filters by depth, base quality, and alternate allele frequency.

```python
from metainformant.dna.variation.calling import call_variants_pileup

pileup = [
    {"chrom": "chr1", "pos": 100, "ref": "A", "depth": 50, "bases": {"A": 40, "G": 10}}
]

variants = call_variants_pileup(pileup, min_depth=10, min_qual=20.0, min_alt_freq=0.1)
# Returns list of dicts with 'chrom', 'pos', 'ref', 'alt', 'genotype', etc.
```

### `genotype_variants`

Compute maximum likelihood genotype from allele counts using a binomial/multinomial likelihood model given an expected sequencing error rate.

```python
from metainformant.dna.variation.calling import genotype_variants

result = genotype_variants(allele_counts={"ref": 40, "alt": 10}, ploidy=2, error_rate=0.01)
# Returns {"genotype": "0/1", "genotype_quality": 99.0, ...}
```

### `filter_variants`

Apply multi-criteria filters to variant calls (e.g., depth, quality, strand bias). Adds a `filter` field ("PASS", or reasons like "LowDepth;LowQual").

```python
from metainformant.dna.variation.calling import filter_variants

passed = filter_variants(variants, min_depth=20, min_qual=30.0)
```

### `merge_variant_calls`

Merge variant calls from multiple callers (or replicates) requiring a minimum number of callers to agree.

```python
from metainformant.dna.variation.calling import merge_variant_calls

caller1 = [...]
caller2 = [...]

consensus = merge_variant_calls([caller1, caller2], min_callers=2)
```

---

## Variant Analysis

### `compute_variant_stats`

Compute summary statistics for a set of variant calls.

```python
from metainformant.dna.variation.calling import compute_variant_stats

stats = compute_variant_stats(consensus)
# Returns {"total_variants": 1000, "ti_tv_ratio": 2.1, "het_hom_ratio": 1.5, ...}
```

### `annotate_variant_context`

Add sequence context and trinucleotide context to variant calls. Converts mutations to the standard 96-channel SBS framework (pyrimidine context) used for mutational signature analysis.

```python
from metainformant.dna.variation.calling import annotate_variant_context

annotated = annotate_variant_context(variants, reference_sequence, window=5)
# Adds fields: "upstream_context", "downstream_context", "trinucleotide_mutation", "sbs_channel"
```
