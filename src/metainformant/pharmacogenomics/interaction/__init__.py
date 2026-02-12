"""Drug-drug interaction analysis subpackage for pharmacogenomics.

Provides drug-drug interaction prediction, polypharmacy risk assessment,
and CYP enzyme inhibition/induction profiling."""
from __future__ import annotations

from . import drug_interactions

__all__ = ['drug_interactions']
