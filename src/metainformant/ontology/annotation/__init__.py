"""Annotation bridge: GWAS hits -> genes -> GO annotations."""
from __future__ import annotations

from .annotate import (
    gwas_hits_to_genes,
    genes_to_go_annotations,
    rank_genes_by_pvalue,
    build_background_from_vcf_genes,
)

__all__ = [
    "gwas_hits_to_genes",
    "genes_to_go_annotations",
    "rank_genes_by_pvalue",
    "build_background_from_vcf_genes",
]
