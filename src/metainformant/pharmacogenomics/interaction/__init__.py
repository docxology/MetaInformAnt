"""Drug-drug interaction analysis subpackage for pharmacogenomics.

Provides drug-drug interaction prediction, polypharmacy risk assessment,
and CYP enzyme inhibition/induction profiling.
"""

from __future__ import annotations

from .drug_interactions import (
    cyp_inhibition_prediction,
    default_interaction_database,
    polypharmacy_risk,
    predict_drug_interaction,
)

__all__ = [
    "predict_drug_interaction",
    "polypharmacy_risk",
    "cyp_inhibition_prediction",
    "default_interaction_database",
]
