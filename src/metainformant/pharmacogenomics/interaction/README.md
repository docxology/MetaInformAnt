# Interaction

Drug-drug interaction prediction and polypharmacy risk assessment for pharmacogenomics, providing CYP enzyme inhibition/induction profiling and interaction severity scoring.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports drug_interactions submodule |
| `drug_interactions.py` | DDI prediction, polypharmacy risk, CYP inhibition profiling |

## Key Functions

| Function | Description |
|----------|-------------|
| `drug_interactions.predict_drug_interaction()` | Predict interaction between two drugs from built-in database |
| `drug_interactions.polypharmacy_risk()` | Assess polypharmacy risk from multiple concurrent medications |
| `drug_interactions.cyp_inhibition_prediction()` | Predict CYP enzyme inhibition/induction potential |
| `drug_interactions.default_interaction_database()` | Get built-in drug interaction database |

## Usage

```python
from metainformant.pharmacogenomics.interaction import drug_interactions

interaction = drug_interactions.predict_drug_interaction("fluoxetine", "codeine")
risk = drug_interactions.polypharmacy_risk(["warfarin", "aspirin", "omeprazole"])
cyp = drug_interactions.cyp_inhibition_prediction("fluconazole")
```
