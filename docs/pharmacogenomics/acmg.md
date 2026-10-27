# ACMG Variant Classification: Pharmacogenomics

`clinical.pathogenicity` classifies variants with the ACMG/AMP 2015 five-tier system
(Richards et al., *Genet Med* 17:405). It is meant for pharmacogene variants that lack
population-frequency or functional data. The functions work on a plain dictionary of evidence
that you supply; nothing is fetched from ClinVar or gnomAD.

## Criteria the module evaluates

`apply_acmg_criteria(variant_data)` evaluates every criterion independently and returns
`{criterion: bool}`. The criteria most often decisive, and the `variant_data` keys they read:

| Code | Strength | Met when |
|------|----------|----------|
| PVS1 | Very strong | `consequence` is `frameshift`, `nonsense`, `stop_gained`, `splice_donor`, `splice_acceptor`, `start_lost`, `initiator_codon` or `exon_deletion` |
| PS1 | Strong | `same_aa_pathogenic` is true |
| PS2 | Strong | `de_novo` and `paternity_confirmed` are true |
| PS3 | Strong | `functional_data` has `result` of `damaging`, `loss_of_function` or `deleterious` with `quality` `strong` or `well-established` |
| PS4 | Strong | `odds_ratio` > 5 and `case_count` ≥ 5 |
| PM1 | Moderate | `in_hotspot` or `in_functional_domain` is true |
| PM2 | Moderate | the highest `allele_frequency` is below `pm2_threshold` (default 0.0001) |
| PM4 | Moderate | `consequence` is `inframe_deletion`, `inframe_insertion` or `stop_lost` |
| PP1 | Supporting | `segregation_data['lod_score']` ≥ 1.5 |
| PP3 | Supporting | at least two `computational_predictions` and more than half are deleterious (`damaging`, `deleterious`, `probably_damaging`, `disease_causing`, `pathogenic`, `high`) |
| BA1 | Stand-alone benign | the highest `allele_frequency` exceeds 0.05 |
| BS1 | Strong benign | the highest `allele_frequency` exceeds `disease_af_threshold` (default 0.01) and BA1 is not met |
| BP4 | Supporting benign | at least two predictions and more than half are benign |
| BP7 | Supporting benign | `consequence` contains `synonymous` and `splice_impact` is not true |

`allele_frequency` is a float or a `{population: frequency}` dict. The remaining criteria (PM3,
PM5, PM6, PP2, PP4, PP5, BS2 to BS4, BP1 to BP3, BP5, BP6) read boolean flags or small dicts in
the same way; see the `apply_acmg_criteria` source for the exact keys. There is no gene-specific
frequency table: the PM2 and BS1 cut-offs are the per-variant `pm2_threshold` and
`disease_af_threshold` keys.

## Full example

```python
from metainformant.pharmacogenomics.clinical.pathogenicity import (
    classify_variant_acmg,
    apply_acmg_criteria,
    ACMGClassification,
)

# Illustrative evidence, not a real gnomAD/ClinVar record
variant = {
    'rsid': 'rs0000000',
    'gene': 'CYP2C19',
    'consequence': 'frameshift',
    'allele_frequency': 0.00002,  # a float, or {population: frequency}
    'computational_predictions': {'sift': 'deleterious', 'polyphen': 'probably_damaging', 'cadd': 'damaging'},
}

criteria = apply_acmg_criteria(variant)  # {criterion: bool}; here PVS1, PM2 and PP3 are True
result = classify_variant_acmg(variant)  # evaluates the criteria itself
print(result['classification'].value)    # Pathogenic
print(result['criteria_met'])            # ['PM2', 'PP3', 'PVS1']

# Or pass pre-evaluated criteria instead of evaluating them from the variant data
result = classify_variant_acmg(variant, evidence={'PVS1': True, 'PM2': True})
print(result['classification'] is ACMGClassification.LIKELY_PATHOGENIC)  # True
```

`classify_variant_acmg(variant, evidence=None)` takes no gene or criteria-list argument: the
criteria come from the variant dict or from `evidence`.

## Output fields

```python
{
  'classification': ACMGClassification.PATHOGENIC,   # Pathogenic, Likely Pathogenic,
                                                     # Uncertain Significance, Likely Benign, Benign
  'criteria_met': ['PM2', 'PP3', 'PVS1'],            # sorted
  'criteria_details': {'PVS1': 'Null variant (...)', ...},   # description of each met criterion
  'score_summary': {'very_strong_pathogenic': 1, 'strong_pathogenic': 0, 'moderate_pathogenic': 1,
                    'supporting_pathogenic': 1, 'standalone_benign': 0, 'strong_benign': 0,
                    'supporting_benign': 0},
  'confidence': 'high',                              # 'high' (3 or more criteria), 'moderate' (2), 'low'
}
```

`aggregate_evidence(criteria)` is the combining step: it applies the rules of Richards et al.
2015, Table 5, to the counts of very strong, strong, moderate and supporting pathogenic criteria
and of stand-alone, strong and supporting benign criteria, and returns an `ACMGClassification`.
Anything the rules do not cover is `Uncertain Significance`.

## Related helpers

- `query_clinvar(variant_id)` looks a variant up (by rsID or HGVS) in a small ClinVar-style
  record set bundled in the module (6 pharmacogene variants) and returns the record or `None`.
  It makes no network request.
- `check_gnomad_frequency(variant, population=None, threshold=0.05)` compares a frequency you
  supply with the BA1 (5%) and BS1 cut-offs and reports which population has the highest one.
- `visualization.plots.plot_acmg_criteria(criteria_results, output_path)` draws which criteria
  are met.

## Combining with phenotype

The classification is not connected to star-allele calling or to activity values. If you add a
novel loss-of-function allele through a custom definition file, you may cite its ACMG
classification as the evidence for giving it an activity value of 0.0, but you assign that
value yourself.

## References

- Richards et al. (2015) *Genet Med* 17:405: standards and guidelines for the interpretation
  of sequence variants.
