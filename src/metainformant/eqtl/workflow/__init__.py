"""eQTL cross-domain workflow adapters.

Composes rna-seq inputs (quantification discovery, ENA retrieval) into
the eQTL input-construction and transcriptome variant-calling workflows;
cross-domain composition is sanctioned here (see
``scripts/quality/check_module_boundaries.py``).
"""

from __future__ import annotations

from . import synthetic, variant_calling  # noqa: F401
from .synthetic import (
    create_synthetic_data,
    create_synthetic_genotypes,
    load_real_expression_data,
    parse_gene_positions,
)
from .variant_calling import (
    align_reads,
    build_hisat2_index,
    call_variants,
    check_tools,
    decompress_if_needed,
    download_fastq,
    filter_variants,
    find_completed_samples,
    find_reference_genome,
    merge_vcfs,
)

__all__ = [
    "align_reads",
    "build_hisat2_index",
    "call_variants",
    "check_tools",
    "create_synthetic_data",
    "create_synthetic_genotypes",
    "decompress_if_needed",
    "download_fastq",
    "filter_variants",
    "find_completed_samples",
    "find_reference_genome",
    "load_real_expression_data",
    "merge_vcfs",
    "parse_gene_positions",
    "synthetic",
    "variant_calling",
]
