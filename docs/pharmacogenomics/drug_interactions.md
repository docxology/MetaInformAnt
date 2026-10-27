# Drug–Drug & Drug–Gene Interactions

Computes severity‑weighted interaction scores from CYP inhibition/induction and
PGx phenotype, returning actionable alerts.

## Core API

```python-snippet
from metainformant.pharmacogenomics.clinical.drug_interaction import analyze_drug_gene_interactions
from metainformant.pharmacogenomics.interaction.drug_interactions import (
    polypharmacy_risk,
    cyp_inhibition_prediction,
    predict_drug_interaction,
)

# Drug–gene interactions for a CYP2C9 poor metabolizer
interactions = analyze_drug_gene_interactions(
    drugs=['warfarin','amiodarone','simvastatin','fluconazole'],
    genotype_data={'CYP2C9': {'phenotype': 'Poor Metabolizer'}},
)
for ix in interactions:
    if ix.severity.value in ('Major', 'Moderate'):
        print(f"ALERT {ix.drug}: {ix.recommendation}")

# Drug–drug interaction for a pair of medications
pair = predict_drug_interaction('warfarin', 'amiodarone')
print(pair['severity'], pair['recommendation'])
```

## Severity levels

| Severity | Definition | Example |
|----------|------------|---------|
| major | Avoid the combination or adjust/substitute | fluconazole + warfarin; PM CYP2D6 + codeine |
| moderate | Consider monitoring or dose reduction | amiodarone + warfarin |
| minor | Low clinical impact | caffeine + cimetidine |
| none | No documented interaction | — |

Drug–drug severities are lowercase strings (`'major'`, `'moderate'`, `'minor'`, `'none'`); drug–gene
severities are `InteractionSeverity` values (`'Major'`, `'Moderate'`, `'Minor'`, `'None'`).

## CYP inhibition taxonomy

| Category | AUC multiplier | Examples |
|----------|----------------|----------|
| Strong inhibitor | >5× | ketoconazole, itraconazole, ritonavir |
| Moderate inhibitor | 2–5× | amiodarone, fluoxetine, fluconazole |
| Weak inhibitor | 1.25–2× | cimetidine, trimethoprim |
| Strong inducer | <0.2× | rifampin, carbamazepine, phenytoin |
| Moderate inducer | 0.2–0.5× | efavirenz, modafinil |

Inhibitor and inducer tables for CYP1A2, 2C9, 2C19, 2D6 and 3A4 are built into
`interaction/drug_interactions.py`; entries are labelled `strong`, `moderate` or `weak`.

## CYP inhibition profile

```python
profile = cyp_inhibition_prediction('ketoconazole')
print(profile['summary'])   # 'ketoconazole is a strong inhibitor of CYP3A4'
for effect in profile['affected_enzymes']:
    print(effect['enzyme'], effect['effect_type'], effect['inhibition_type'])  # CYP3A4 inhibitor strong
```

Each call looks one drug up in the built-in inhibitor and inducer tables and also lists the
CYP enzymes it is a substrate of (`profile['substrate_of']`).

## Polypharmacy composite risk

```python
risk = polypharmacy_risk(patient_meds)
print(f"Risk score = {risk['risk_score']:.2f}")
print(f"Major interactions = {risk['severity_counts']['major']}")
```

Scoring: pairwise severities are weighted `minor=0.15, moderate=0.5, major=1.0`, 0.1 is added for each
extra drug competing for the same CYP enzyme, and the sum is divided by the number of drug pairs
(maximum 1.0). At least two medications are required.

## Report inclusion

`generate_clinical_report()` runs `analyze_drug_gene_interactions()` for the requested drugs (every
drug with a built-in CPIC guideline when `drugs=None`; none when `drugs=[]`), using the phenotypes it
classified from each gene's diplotype, and adds an `interaction_summary` (counts per severity and an
overall risk level) and prioritized `clinical_actions` for the major and moderate findings. See
[reporting.md](reporting.md).

## Extending the knowledge base

The built-in interaction database is a plain dictionary keyed by the sorted, lowercased drug
pair. To add an interaction, extend a copy and pass it as `interaction_db`:

```python-snippet
from metainformant.pharmacogenomics.interaction.drug_interactions import (
    default_interaction_database,
    polypharmacy_risk,
)

db = default_interaction_database()
db[('my_drug', 'warfarin')] = {
    'severity': 'major',
    'mechanism': 'CYP2C9 inhibition',
    'description': 'my_drug inhibits CYP2C9 and raises warfarin exposure',
    'evidence_level': 'B',
    'recommendation': 'Reduce the warfarin dose and monitor INR',
}
risk = polypharmacy_risk(['my_drug', 'warfarin'], interaction_db=db)
```

`predict_drug_interaction` and `polypharmacy_risk` accept `interaction_db`, and
`cyp_inhibition_prediction` accepts a custom `cyp_profiles` dictionary
(`{'inhibitors': {cyp: [{'drug', 'type', 'potency'}]}, 'inducers': {...}}`).

Run tests (`pytest tests/pharmacogenomics/test_drug_interactions_depth.py`) to validate.

## Pitfalls

- **Double counting** — If drug A inhibits and drug B induces the same CYP, their effects
  on that enzyme may offset each other, but individual pairwise severities still apply;
  both alerts show.
- **Phenotype labels** — `analyze_drug_gene_interactions` reads each gene's `phenotype` key and
  matches it to the CPIC guideline's phenotype names: pass a `MetabolizerPhenotype`, its full name
  (`'Poor Metabolizer'`) or the two-letter abbreviation (`'PM'`). The bare token `'poor'` matches
  nothing and the gene yields no entry (the dose-adjustment table in `metabolism` does accept it;
  see `MetabolizerPhenotype.from_string` to normalise any label first).
- **Missing CYP entry** — If a drug is not in the internal tables, it is silently
  ignored (contributes no severity). Add it with a custom `interaction_db` as above.
