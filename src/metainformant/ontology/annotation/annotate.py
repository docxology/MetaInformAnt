"""Annotation bridge: GWAS hits -> gene lists -> GO annotations.

Converts GWAS association results and gene annotations (from step 8g) into
gene lists and ranked gene lists suitable for ORA and GSEA, and fetches GO
annotations from either the QuickGO REST API or a local OBO gene-annotation
file.

Typical usage::

    from metainformant.ontology.annotation.annotate import (
        gwas_hits_to_genes,
        genes_to_go_annotations,
        rank_genes_by_pvalue,
    )

    genes = gwas_hits_to_genes(association_results, top_n=20, gene_annotations=annotations)
    gene_sets = genes_to_go_annotations(genes, taxon_id=7460, source="api")
    ranked = rank_genes_by_pvalue(association_results, gene_annotations=annotations)
"""

from __future__ import annotations

import math
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Gene extraction from GWAS results
# ---------------------------------------------------------------------------


def gwas_hits_to_genes(
    association_results: list[dict[str, Any]],
    *,
    top_n: int = 20,
    gene_annotations: list[dict[str, Any]] | None = None,
    pvalue_threshold: float = 0.05,
) -> list[str]:
    """Extract gene names from GWAS top hits.

    Derives gene names from step-8g Ensembl annotations when available,
    falling back to SNP IDs for unannotated hits.

    Args:
        association_results: List of GWAS result dicts with at least
            ``p_value``, ``snp``, ``chrom``, ``pos`` keys.
        top_n: Maximum number of top hits to extract genes from.
        gene_annotations: List of annotation dicts from step 8g
            (``annotate_top_hits_ensembl``), each with ``snp``,
            ``nearest_gene``, ``nearby_genes``.
        pvalue_threshold: Soft threshold; hits above this are deprioritised
            (all top_n are still taken from sorted results).

    Returns:
        Deduplicated list of gene identifiers.

    Examples:
        >>> genes = gwas_hits_to_genes(results, top_n=10)
        >>> len(genes) <= 10
        True
    """
    # Sort by p-value ascending
    sorted_hits = sorted(association_results, key=lambda r: r.get("p_value", 1.0))

    # Build SNP -> nearest_gene map from step-8g annotations
    snp_to_genes: dict[str, list[str]] = {}
    for ann in gene_annotations or []:
        snp = ann.get("snp", "")
        nearby = ann.get("nearby_genes", [])
        nearest = ann.get("nearest_gene", "")
        genes_for_snp: list[str] = []
        if nearest:
            genes_for_snp.append(nearest)
        for g in nearby:
            name = g.get("gene_name", "")
            if name and name not in genes_for_snp:
                genes_for_snp.append(name)
        if snp and genes_for_snp:
            snp_to_genes[snp] = genes_for_snp

    seen: set[str] = set()
    genes: list[str] = []

    for hit in sorted_hits[:top_n]:
        snp = hit.get("snp", "")

        if snp in snp_to_genes:
            for gname in snp_to_genes[snp]:
                if gname not in seen:
                    seen.add(gname)
                    genes.append(gname)
        # If no gene annotation, skip this hit — no positional labels, no SNP-id placeholders.
        # Real gene annotations come from step 8g (NCBI/Ensembl API).

    logger.info(
        "gwas_hits_to_genes: %d unique genes from top %d hits " "(pvalue_threshold=%.2e)",
        len(genes),
        top_n,
        pvalue_threshold,
    )
    return genes


def rank_genes_by_pvalue(
    association_results: list[dict[str, Any]],
    *,
    gene_annotations: list[dict[str, Any]] | None = None,
    min_neg_log_p: float = 0.0,
) -> list[tuple[str, float]]:
    """Build a ranked gene list for GSEA from GWAS p-values.

    Computes ``-log10(p_value)`` as the rank metric and maps SNPs to gene
    names using step-8g annotations.  SNPs without gene annotations are
    represented by their genomic position label.

    Args:
        association_results: GWAS result dicts.
        gene_annotations: Step-8g Ensembl annotation list.
        min_neg_log_p: Minimum ``-log10(p)`` to include (filters near-null).

    Returns:
        List of ``(gene_name, -log10_p)`` tuples, sorted descending.

    Examples:
        >>> ranked = rank_genes_by_pvalue(results)
        >>> ranked[0][1] >= ranked[-1][1]
        True
    """
    snp_to_gene: dict[str, str] = {}
    for ann in gene_annotations or []:
        snp = ann.get("snp", "")
        nn = ann.get("nearest_gene", "")
        if snp and nn:
            snp_to_gene[snp] = nn

    seen: dict[str, float] = {}

    for hit in association_results:
        p = hit.get("p_value", 1.0)
        # Clamp to avoid -inf
        p_clamped = max(p, 1e-300)
        neg_log_p = -math.log10(p_clamped)

        if neg_log_p < min_neg_log_p:
            continue

        snp = hit.get("snp", "")
        gene = snp_to_gene.get(snp, "")

        # Only include SNPs that have a real gene annotation from step 8g.
        # No positional chrom:pos placeholders — those are not valid gene IDs for GSEA.
        if not gene:
            continue

        # If multiple SNPs map to same gene, take the most significant
        if gene not in seen or neg_log_p > seen[gene]:
            seen[gene] = neg_log_p

    ranked = sorted(seen.items(), key=lambda x: x[1], reverse=True)
    logger.info(
        "rank_genes_by_pvalue: %d unique genes ranked, top=%s (%.2f)",
        len(ranked),
        ranked[0][0] if ranked else "N/A",
        ranked[0][1] if ranked else 0.0,
    )
    return ranked


# ---------------------------------------------------------------------------
# Gene -> GO annotation fetching
# ---------------------------------------------------------------------------


def genes_to_go_annotations(
    gene_list: list[str],
    *,
    taxon_id: int = 7460,
    source: str = "api",
    go_aspects: list[str] | None = None,
    annotation_file: str | None = None,
    rate_limit_s: float = 0.2,
) -> dict[str, set[str]]:
    """Fetch GO annotations for a list of genes.

    Args:
        gene_list: Gene identifiers (symbols or Ensembl IDs).
        taxon_id: NCBI taxon ID (default ``7460`` — *Apis mellifera*).
        source: ``"api"`` (QuickGO REST, default) or ``"gaf"`` (local GAF file).
        go_aspects: Optional list of GO aspects to restrict
            (``"biological_process"``, ``"molecular_function"``,
            ``"cellular_component"``).
        annotation_file: Path to a GAF file (required when ``source="gaf"``).
        rate_limit_s: Pause between API calls (``source="api"`` only).

    Returns:
        ``dict[go_id, set[gene_name]]`` suitable for ORA/GSEA.

    Examples:
        >>> go_anns = genes_to_go_annotations(["SOD1", "CAT"], taxon_id=7460)
        >>> isinstance(go_anns, dict)
        True
    """
    if source == "api":
        from metainformant.ontology.core.go_api import build_gene_set_from_api

        return build_gene_set_from_api(
            gene_list,
            taxon_id=taxon_id,
            go_aspects=go_aspects,
            rate_limit_s=rate_limit_s,
        )

    elif source == "gaf":
        return _load_gaf_annotations(
            gene_list=gene_list,
            annotation_file=annotation_file,
            go_aspects=go_aspects,
        )

    else:
        raise ValueError(f"Unknown annotation source '{source}'. Use 'api' or 'gaf'.")


def _load_gaf_annotations(
    *,
    gene_list: list[str],
    annotation_file: str | None,
    go_aspects: list[str] | None,
) -> dict[str, set[str]]:
    """Load GO annotations from a Gene Annotation Format (GAF) file.

    Internal fallback for when API access is unavailable.

    Args:
        gene_list: Genes to annotate (subset filter).
        annotation_file: Path to GAF 2.x file.
        go_aspects: Optional aspect filter (P/C/F short codes accepted).

    Returns:
        ``dict[go_id, set[gene_name]]``.
    """
    if not annotation_file:
        logger.warning("GAF source requested but no annotation_file provided")
        return {}

    import gzip
    import pathlib

    gene_set = set(gene_list)
    aspect_map = {
        "biological_process": "P",
        "molecular_function": "F",
        "cellular_component": "C",
        "BP": "P",
        "MF": "F",
        "CC": "C",
    }
    allowed_aspects: set[str] | None = None
    if go_aspects:
        allowed_aspects = {aspect_map.get(a, a) for a in go_aspects}

    gene_sets: dict[str, set[str]] = {}
    path = pathlib.Path(annotation_file)

    opener = gzip.open if path.suffix == ".gz" else open

    try:
        with opener(path, "rt") as fh:
            for line in fh:
                if line.startswith("!"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 10:
                    continue
                gene_symbol = parts[2]
                go_id = parts[4]
                aspect = parts[8]  # P / C / F

                if gene_symbol not in gene_set:
                    continue
                if allowed_aspects and aspect not in allowed_aspects:
                    continue

                gene_sets.setdefault(go_id, set()).add(gene_symbol)

        logger.info(
            "GAF annotation: %d GO terms for %d genes from %s",
            len(gene_sets),
            len(gene_set),
            annotation_file,
        )
    except Exception as exc:
        logger.warning("Failed to load GAF file %s: %s", annotation_file, exc)

    return gene_sets


# ---------------------------------------------------------------------------
# Background gene set helpers
# ---------------------------------------------------------------------------


def build_background_from_vcf_genes(
    gene_annotations: list[dict[str, Any]],
) -> list[str]:
    """Collect all genes seen in a set of GWAS hit annotations.

    Used to build an unbiased background gene set for ORA.

    Args:
        gene_annotations: Ensembl step-8g annotation list.

    Returns:
        Deduplicated list of all annotated gene names.

    Examples:
        >>> bg = build_background_from_vcf_genes(annotations)
        >>> isinstance(bg, list)
        True
    """
    seen: set[str] = set()
    for ann in gene_annotations:
        nn = ann.get("nearest_gene", "")
        if nn:
            seen.add(nn)
        for g in ann.get("nearby_genes", []):
            name = g.get("gene_name", "")
            if name:
                seen.add(name)
    result = sorted(seen)
    logger.info("build_background_from_vcf_genes: %d background genes", len(result))
    return result
