"""
eQTL Analysis Wrapper.

This module provides high-level interfaces for running expression Quantitative
Trait Locus (eQTL) analysis, linking genotype data with RNA-seq expression levels.
"""

from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

from metainformant.core.utils.logging import get_logger
from metainformant.gwas.analysis.association import association_test_linear

logger = get_logger(__name__)

#: cis window (bp each side of the TSS) used when positions are given without an explicit ``cis_window``.
DEFAULT_CIS_WINDOW = 1_000_000

_GENE_POSITION_COLUMNS = ("gene_id", "chrom", "tss_position")
_VARIANT_POSITION_COLUMNS = ("variant_id", "chrom", "position")


def _validated_positions(table: pd.DataFrame, columns: Tuple[str, ...], name: str) -> pd.DataFrame:
    """Return ``table`` after checking its required columns and that its identifier column is unique."""
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise ValueError(f"{name} is missing required column(s) {missing}; expected columns {list(columns)}")
    duplicated = table[columns[0]][table[columns[0]].duplicated()].unique().tolist()
    if duplicated:
        raise ValueError(f"{name} lists {columns[0]} more than once: {duplicated[:5]}")
    return table


def _cis_candidates(
    transcripts: List[Any],
    variants: List[Any],
    gene_positions: pd.DataFrame,
    variant_positions: pd.DataFrame,
    cis_window: int,
) -> Dict[Any, List[Tuple[Any, Optional[int]]]]:
    """Map each transcript to the variants within ``cis_window`` bp of its TSS on the same chromosome.

    Both sides of the window are included (``|position - tss| <= cis_window``), as in
    :func:`metainformant.gwas.finemapping.eqtl.cis_eqtl_scan`. The variants of a transcript keep the order of
    the genotype columns and come with their signed distance to the TSS (variant position minus TSS). A
    transcript or variant that has no row in the position tables cannot be placed: transcripts get no
    candidates and variants are never candidates; both are counted in a warning.
    """
    genes = _validated_positions(gene_positions, _GENE_POSITION_COLUMNS, "gene_positions")
    sites = _validated_positions(variant_positions, _VARIANT_POSITION_COLUMNS, "variant_positions")

    tss_by_gene = dict(zip(genes["gene_id"], zip(genes["chrom"].astype(str), genes["tss_position"])))
    site_by_variant = dict(zip(sites["variant_id"], zip(sites["chrom"].astype(str), sites["position"])))

    placed: List[Tuple[Any, str, Any]] = [
        (variant, site_by_variant[variant][0], site_by_variant[variant][1])
        for variant in variants
        if variant in site_by_variant
    ]
    if len(placed) < len(variants):
        logger.warning(f"{len(variants) - len(placed)} variant(s) have no row in variant_positions; not tested")

    candidates: Dict[Any, List[Tuple[Any, Optional[int]]]] = {}
    unplaced_transcripts = 0
    for transcript in transcripts:
        if transcript not in tss_by_gene:
            unplaced_transcripts += 1
            candidates[transcript] = []
            continue
        chrom, tss = tss_by_gene[transcript]
        candidates[transcript] = [
            (variant, int(position - tss))
            for variant, variant_chrom, position in placed
            if variant_chrom == chrom and abs(position - tss) <= cis_window
        ]
    if unplaced_transcripts:
        logger.warning(f"{unplaced_transcripts} transcript(s) have no row in gene_positions; not tested")
    return candidates


def run_eqtl_analysis(
    genotype_matrix: pd.DataFrame,
    expression_matrix: pd.DataFrame,
    covariates: Optional[pd.DataFrame] = None,
    cis_window: Optional[int] = None,
    max_transcripts: Optional[int] = None,
    *,
    variant_positions: Optional[pd.DataFrame] = None,
    gene_positions: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """
    Run basic linear eQTL analysis for all transcripts.

    Without position tables every variant is tested against every transcript. With
    ``variant_positions`` and ``gene_positions`` the scan is restricted to cis pairs: a variant is
    tested against a transcript only when both lie on the same chromosome and the variant is within
    ``cis_window`` bp of the transcript's TSS (either side). A ``cis_window`` is meaningful only
    with positions, and is rejected without them rather than silently ignored.

    Note: A full all-vs-all eQTL scan is computationally expensive.
    This is a simplified wrapper for verification of the data pipeline.

    Args:
        genotype_matrix: (Samples x Variants) - 0/1/2 dosage.
        expression_matrix: (Samples x Transcripts) - Normalized expression.
        covariates: (Samples x Covariates) - Optional.
        cis_window: Half-width of the cis window in base pairs around each TSS. Needs
            ``variant_positions`` and ``gene_positions``; defaults to ``DEFAULT_CIS_WINDOW`` (1 Mb)
            when they are given and this is omitted.
        max_transcripts: Optional cap for smoke/demo runs. If None, scan all transcripts.
        variant_positions: Optional DataFrame with columns [variant_id, chrom, position], one row per
            genotype column (``variant_id`` matches the column labels). Give it together with
            ``gene_positions``.
        gene_positions: Optional DataFrame with columns [gene_id, chrom, tss_position], one row per
            expression column (``gene_id`` matches the column labels). Give it together with
            ``variant_positions``.

    Returns:
        pd.DataFrame: Summary statistics (variant, transcript, beta, p_value, ...). A cis scan adds a
        ``distance`` column: the variant position minus the TSS, in bp.

    Raises:
        ValueError: If only one of the position tables is given, ``cis_window`` is given without
            positions or is negative, a position table lacks a required column or repeats an
            identifier, or ``max_transcripts`` is below 1.
    """
    if (variant_positions is None) != (gene_positions is None):
        raise ValueError("variant_positions and gene_positions must be given together")
    cis_scan = variant_positions is not None and gene_positions is not None
    if cis_window is not None and not cis_scan:
        raise ValueError("cis_window needs variant_positions and gene_positions to place variants relative to genes")
    if cis_window is not None and cis_window < 0:
        raise ValueError(f"cis_window must be non-negative, got {cis_window}")

    results = []

    # Align samples
    common_samples = genotype_matrix.index.intersection(expression_matrix.index)

    if covariates is not None:
        common_samples = common_samples.intersection(covariates.index)

    if len(common_samples) == 0:
        logger.error("No common samples between genotype and expression matrices.")
        return pd.DataFrame()

    logger.info(f"Running eQTL analysis on {len(common_samples)} samples.")

    # Subset data
    G = genotype_matrix.loc[common_samples]
    E = expression_matrix.loc[common_samples]
    C = covariates.loc[common_samples] if covariates is not None else None

    # For this verification implementation, we run a naive loop
    # In production, use TensorQTL or MatrixeQTL approaches

    if max_transcripts is not None and max_transcripts < 1:
        raise ValueError("max_transcripts must be a positive integer or None")

    transcripts = E.columns if max_transcripts is None else E.columns[:max_transcripts]

    # The (variant, distance) pairs each transcript is tested against; distance is None for an all-vs-all scan
    candidates: Dict[Any, List[Tuple[Any, Optional[int]]]]
    if variant_positions is not None and gene_positions is not None:
        window = DEFAULT_CIS_WINDOW if cis_window is None else cis_window
        candidates = _cis_candidates(list(transcripts), list(G.columns), gene_positions, variant_positions, window)
    else:
        candidates = {transcript: [(variant, None) for variant in G.columns] for transcript in transcripts}

    # Covariate rows are invariant across the variant/transcript loops
    covariate_rows = C.T.values.tolist() if C is not None else None

    for transcript_id in transcripts:
        phenotype = E[transcript_id]

        for variant_id, distance in candidates[transcript_id]:
            try:
                assoc_stats = association_test_linear(
                    genotypes=G[variant_id].astype(float).tolist(),
                    phenotypes=phenotype.astype(float).tolist(),
                    covariates=covariate_rows,
                )

                assoc_stats["transcript"] = transcript_id
                assoc_stats["variant"] = variant_id
                if distance is not None:
                    assoc_stats["distance"] = distance
                results.append(assoc_stats)

            except Exception as e:
                logger.warning(f"Failed eQTL for {variant_id}/{transcript_id}: {e}")

    if not results:
        return pd.DataFrame()

    return pd.DataFrame(results)
