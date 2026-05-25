# Drug–Drug & Drug–Gene Interactions

Computes severity‑weighted interaction scores from CYP inhibition/induction and
PGx phenotype, returning actionable alerts.

## Core API

```python-snippet
from metainformant.pharmacogenomics.interaction.drug_interactions import (
    analyze_drug_gene_interactions,
    polypharmacy_risk,
    cyp_inhibition_prediction,
)

interactions = analyze_drug_gene_interactions(
    gene='CYP2C9',
    phenotype='Poor Metabolizer',
    drugs=['warfarin','amiodarone','simvastatin','fluconazole'],
)
for ix in interactions:
    if ix['severity'] in ('contraindicated','major'):
        print(f"ALERT {ix['drug']}: {ix['recommendation']}")
```

## Severity levels

| Severity | Definition | Example |
|----------|------------|---------|
| contraindicated | Absolutely avoid combination | PM CYP2C9 + fluconazole + warfarin |
| major | Strong recommendation to adjust/substitute | IM CYP2D6 + codeine |
| moderate | Consider monitoring or dose reduction | NM + simvastatin |
| minor | Low clinical impact | Any + caffeine |

## CYP inhibition taxonomy

| Category | AUC multiplier | Examples |
|----------|----------------|----------|
| Strong inhibitor | >5× | ketoconazole, itraconazole, ritonavir |
| Moderate inhibitor | 2–5× | amiodarone, fluoxetine, fluconazole |
| Weak inhibitor | 1.25–2× | cimetidine, trimethoprim |
| Strong inducer | <0.2× | rifampin, carbamazepine, phenytoin |
| Moderate inducer | 0.2–0.5× | efavirenz, modafinil |

Each CYP enzyme (2C9, 2C19, 2D6, 3A4) has its own inhibition table loaded from
YAML in `interaction/data/`.

## Net CYP inhibition score

```python
score = cyp_inhibition_prediction(['ketoconazole','rifampin'], enzyme='CYP3A4')
print(score['category'])   # 'neutral' — strong inhibitor balances strong inducer
print(score['mechanism'])  # 'mixed_inhibition_induction'
```

The algorithm sums log₂ fold‑change from each drug (inhibitor +1, inducer −1) then
bins the net score.

## Polypharmacy composite risk

```python
risk = polypharmacy_risk(
    drugs=patient_meds,
    patient_phenotype='Normal Metabolizer',
    gene='CYP2C9',
)
print(f"Risk score = {risk['risk_score']:.2f}")
print(f"High‑severity count = {risk['high_severity_count']}")
```

Scoring: `minor=1, moderate=2, major=3, contraindicated=4`; weighted sum over
drugs/pairs.

## Report inclusion

`generate_clinical_report()` automatically runs `analyze_drug_gene_interactions()`
for each gene and adds a **DDI Warnings** section showing only `contraindicated`
and `major` findings (configurable via `ddi_severity_threshold`).

## Extending the knowledge base

Drug definitions live in `src/metainformant/pharmacogenomics/interaction/data/` as
YAML. To add a new drug:

```yaml
# my_drug.yaml
name: my_drug
inhibits: [CYP2C19, CYP3A4]
inducer: []
therapeutic_index: narrow   # 'narrow' or 'wide'
pgx_class: anticoagulant
```

Then register in `_DRUG_REGISTRY` inside `drug_interactions.py`:

```python-snippet
from metainformant.pharmacogenomics.interaction.data.my_drug import MY_DRUG
_DRUG_REGISTRY['my_drug'] = MY_DRUG
```

Run tests (`pytest tests/pharmacogenomics/interaction/`) to validate.

## Pitfalls

- **Double counting** — If drug A inhibits and drug B induces the same CYP, the net
  score may be neutral, but individual pairwise severities still apply; both alerts
  show.
- **Phenotype mismatches** — `phenotype` must match one of `{'Poor','PM',…}`; use
  `MetabolizerPhenotype.POOR.value` for safety.
- **Missing CYP entry** — If a drug is not in the internal tables, it is silently
  ignored (contributes no severity). Add it via YAML as above.
