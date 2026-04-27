# Star Allele Matching Engine

PharmVar standardizes haplotype labels (e.g. `*1`, `*4`, `*2A`) — this module
implements the matching algorithm.

## Matching algorithm (detailed walk‑through)

Given observed variant set V, gene G with alleles A₁…Aₙ. Each allele Aᵢ has
canonical defining variants Dᵢ (dbSNP rsIDs) from PharmVar.

```
Step 1 — Filter: keep alleles where Dᵢ ⊆ V
Step 2 — Score each kept allele:  cov = |V ∩ Dᵢ| / |Dᵢ|         (coverage)
Step 3 — Sort descending by (cov, activity, lexicographic name)
Step 4 — Return top‑k alleles (k = 1 fast; 2–3 balanced; all minimal accurate)
```

Complexity O(n × m) — n alleles, m avg variants per allele. For CYP2D6 (n≈140)
worst case ~140 × 9 = 1 260 comparisons → <1 ms thanks to Cython hash tables.

## Duplication / CNV handling (CYP2D6 specific)

CNV marker rsIDs: `rs28371685`, `rs529216` (and a few others). If any present in V
→ duplication inferred. Matching then permits multi‑copy star alleles like `*1x2`
or `*1x3`. `determine_diplotype()` incorporates copy number into activity
calculation (× copy count).

If no evidence, copy number defaults to 2 even if one allele is a deletion (`*5`,
activity 0.0). This conservative assumption avoids over‑calling duplication.

## Activity score reference

Per‑gene activity values (excerpt):

| Gene | *1 | *2 | *4 | *10 | *41 | *5 (del) |
|------|----|----|----|-----|-----|----------|
| CYP2D6 | 1.0 | 1.0 | 0.0 | 0.25 | 0.5 | 0.0 |
| CYP2C19 | 1.0 | 1.0 | 0.0 | 0.5 | — | — |
| TPMT | 1.0 | 1.0 | 0.0 | — | — | — |

Complete tables in `src/metainformant/pharmacogenomics/alleles/activity.py`.

## Custom gene / allele addition

```python
import json
from metainformant.pharmacogenomics.alleles.star_allele import load_allele_definitions

cyp2e1 = {
  "gene": "CYP2E1",
  "star_alleles": [
    {"name":"*1","variants":["rs2031920"],"function":"normal","activity":1.0},
    {"name":"*5B","variants":["rs3813867"],"function":"decreased","activity":0.5},
  ],
}
with open('/tmp/cyp2e1.json','w') as fh: json.dump(cyp2e1, fh, indent=2)
load_allele_definitions('CYP2E1', filepath='/tmp/cyp2e1.json')
```

Now `call_star_alleles(variants, gene='CYP2E1')` uses your table.

## Debugging mismatches

If expected alleles not returned:

```python
from metainformant.pharmacogenomics.alleles.star_allele import _BUILTIN_ALLELE_DEFINITIONS
# Which variants do you have vs which alleles expect?
alleles = _BUILTIN_ALLELE_DEFINITIONS['CYP2D6']
for a in alleles:
    if a.defining_variants.issubset(your_variants):
        print('MATCH:', a.name, a.defining_variants)
```

Common pitfall: using old rsIDs (e.g. `rs3918` for `*2` was replaced by newer
rs number). Always reference PharmVar's current variant list per allele.

## Performance internals

Matching uses Cython’s `HashTable` for `set` operations on variant IDs (strings).
Activity lookup is a plain Python dict. The whole function runs without GIL
release, but still <10 ms for worst‑case genes. Future versions will move the
inner loop to Cython for additional speed.

## Known edge cases

- **Allele with zero defining variants** — treated as wildcard (matches always).
  Used for conceptual full‑gene deletions when no marker exists (rare; currently
  only `*5` in CYP2D6).
- **Shared variants across alleles** — e.g. `*1` and `*2` share some markers. Both
  pass filter; ranking by higher coverage decides, then activity.
- **Tie on coverage & activity** — lexical order ( `*2` before `*4` ) provides
  deterministic output; graph‑based tie resolution planned for v1.0.
