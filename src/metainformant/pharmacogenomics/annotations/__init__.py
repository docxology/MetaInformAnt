"""Pharmacogenomics annotation subpackage.

Provides CPIC guideline lookups, PharmGKB annotation queries,
and FDA drug label parsing for pharmacogenomic biomarker information.
"""

from __future__ import annotations

from . import cpic
from . import pharmgkb
from . import drug_labels

__all__ = ["cpic", "pharmgkb", "drug_labels"]
