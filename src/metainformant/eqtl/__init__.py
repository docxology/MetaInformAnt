"""eQTL and transcriptome-variant analysis module for METAINFORMANT.

This module provides reusable eQTL analysis components: synthetic and real
expression/genotype data preparation, transcriptome SNP-calling pipeline
methods (HISAT2 alignment + bcftools variant calling), and bcftools stats
parsing. Scripts under ``scripts/eqtl/`` are thin orchestrators that import
from this module.

Example:
    >>> from metainformant.eqtl.synthetic import create_synthetic_data
    >>> expr, geno, gene_pos, var_pos = create_synthetic_data(
    ...     n_genes=5, n_variants=25, n_samples=12
    ... )
    >>> expr.shape
    (5, 12)
"""

from __future__ import annotations

from . import pipeline, synthetic, variant_calling, variant_stats
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
from .variant_stats import (
    compute_popgen_summary,
    compute_sample_stats,
    parse_bcftools_stats,
)

__all__ = [
    "align_reads",
    "build_hisat2_index",
    "call_variants",
    "check_tools",
    "compute_popgen_summary",
    "compute_sample_stats",
    "create_synthetic_data",
    "create_synthetic_genotypes",
    "decompress_if_needed",
    "download_fastq",
    "filter_variants",
    "find_completed_samples",
    "find_reference_genome",
    "load_real_expression_data",
    "merge_vcfs",
    "parse_bcftools_stats",
    "parse_gene_positions",
    "pipeline",
    "synthetic",
    "variant_calling",
    "variant_stats",
]
