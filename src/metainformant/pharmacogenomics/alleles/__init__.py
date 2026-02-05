"""Pharmacogenomics allele calling subpackage.

Provides star allele calling, diplotype determination, activity scoring,
and metabolizer phenotype prediction for pharmacogenes.
"""

from __future__ import annotations

from . import star_allele
from . import diplotype
from . import phenotype

__all__ = ["star_allele", "diplotype", "phenotype"]
