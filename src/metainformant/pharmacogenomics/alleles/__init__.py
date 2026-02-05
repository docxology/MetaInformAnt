"""Pharmacogenomics allele calling subpackage.

Provides star allele calling, diplotype determination, activity scoring,
and metabolizer phenotype prediction for pharmacogenes.
"""

from __future__ import annotations

from . import diplotype, phenotype, star_allele

__all__ = ["star_allele", "diplotype", "phenotype"]
