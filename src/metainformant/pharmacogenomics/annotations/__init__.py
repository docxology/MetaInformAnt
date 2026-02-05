"""Pharmacogenomics annotation subpackage.

Provides CPIC guideline lookups, PharmGKB annotation queries,
and FDA drug label parsing for pharmacogenomic biomarker information.
"""

from __future__ import annotations

from . import cpic, drug_labels, pharmgkb

__all__ = ["cpic", "pharmgkb", "drug_labels"]
