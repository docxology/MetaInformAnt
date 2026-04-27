# ACMG Variant Classification: Pharmacogenomics

Applied to pharmacogenes to classify pathogenicity of rare variants lacking population
frequency or functional data. Criteria adapted from ACMG/AMP 2015 guidelines with
gene‑specific tweaks.

## Criteria table (PGx‑adapted)

| Code | Full name | Strength | When it applies |
|------|-----------|----------|-----------------|
| PVS1 | Null variant | Very strong | Nonsense, frameshift, canonical splice, whole‑gene deletion |
| PS1 | Same amino‑acid change as established pathogenic | Strong | Missense identical to known pathogenic |
| PM1 | Located in mutational hot‑spot / critical domain | Moderate | Variant in well‑established functional domain |
| PM2 | Absent/rare in population databases | Moderate | gnomAD AF < gene‑specific threshold |
| PP3 | Multiple computational lines of evidence | Supporting | CADD>20, SIFT<0.05, PolyPhen>0.95 |
| BA1 | Allele frequency ≥ 5 % in population | Stand‑alone | Benign by prevalence |

PGx genes often have lower variability, so PM2 thresholds are stricter (e.g. CYP2D6
AF < 0.0005).

## Full example

```python
from metainformant.pharmacogenomics.clinical.pathogenicity import (
    classify_variant_acmg,
    ACMGClassification,
    ACMGCriteria,
)

variant = {
    'rsid': 'rs4244285',          # CYP2C19 *2 loss‑of‑function
    'gene': 'CYP2C19',
    'gnomad_af': 0.0002,
    'clinvar_sig': 'Pathogenic',
    'sift': 0.0,
    'polyphen': 1.0,
    'cadd_score': 25.0,
    'provean_score': -5.0,
}

result = classify_variant_acmg(
    variant,
    gene='CYP2C19',
    criteria_applied=[ACMGCriteria.PVS1, ACMGCriteria.PM2, ACMGCriteria.PP3],
)
print(result['classification'].value)  # PATHOGENIC or LIKELY_PATHOGENIC
print(result['evidence_summary'])      # list of applied criteria with weights
```

## Gene‑specific thresholds

| Gene | PM2 AF cutoff | PVS1 exceptions |
|------|---------------|-----------------|
| CYP2D6 | < 0.0005 | Deletion `*5` always PVS1 |
| CYP2C19 | < 0.001 | Splice variants require RNAseq confirmation |
| TPMT | < 0.002 | Missense must have functional assay support |
| NUDT15 | < 0.001 | |

These live in `src/metainformant/pharmacogenomics/clinical/pathogenicity.py` in
`_GENE_THRESHOLDS`.

## Output fields

```python
{
  'classification': ACMGClassification.PATHOGENIC,
  'activity_score': 0.0,
  'evidence_summary': [
    {'code':'PVS1','weight':1.0,'reason':'nonsense variant'},
    {'code':'PM2','weight':0.5,'reason':'gnomAD AF=2e-4'},
    {'code':'PP3','weight':0.5,'reason':'CADD=25, SIFT=0.0'}
  ],
  'recommendation': 'Consider functional assay for confirmation'
}
```

## Combining with phenotype

Once a variant is classified, its effect on activity enters the star-allele table.
If you add a novel LoF variant via custom JSON, you can use ACMG classification
as evidence to justify the activity assignment (0.0).

## References

- Richards et al. (2015) *Genet Med* 17:405 — standards and guidelines.
- CPIC allele function criteria (2024) — gene‑specific mappings.
