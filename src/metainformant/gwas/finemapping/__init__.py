"""GWAS statistical fine-mapping submodule.

Provides credible set analysis, colocalization, SuSiE regression,
and conditional analysis for identifying causal variants from GWAS signals."""
from __future__ import annotations

from . import colocalization, credible_sets, eqtl

__all__ = ['colocalization', 'credible_sets', 'eqtl']
