"""Gene-set enrichment analysis for GWAS hits.

Implements Fisher's exact test (from scipy.stats) for enrichment of
GWAS hits in gene sets. Also provides a QuickGO REST API interface for
GO term enrichment, and a wrapper for computing pathway summaries from
annotated GWAS results.

Key analyses:
  - Overrepresentation analysis: Fisher's exact test (one-sided)
  - Multiple testing correction: Benjamini-Hochberg FDR
  - GO term lookup via QuickGO REST API (EBI)

References:
  Subramanian et al. (2005) PNAS — GSEA
  Boyle et al. (2004) Bioinformatics — GO enrichment
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

import requests
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

try:
    from scipy.stats import fisher_exact, chi2  # noqa: F401

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

QUICKGO_BASE = "https://www.ebi.ac.uk/QuickGO/services"
QUICKGO_TIMEOUT = 10


def fisher_enrichment_test(
    hit_genes: List[str],
    gene_set: List[str],
    background_genes: List[str],
) -> Dict[str, Any]:
    """One-sided Fisher's exact test for gene-set enrichment.

    Tests: are GWAS hit genes enriched in `gene_set` relative to the
    background (post-QC genome-wide tested genes)?

    2×2 table:
                     |  In set  |  Not in set  |
    GWAS hits        |    a     |      b       |
    Not hits         |    c     |      d       |

    Args:
        hit_genes: Gene names near GWAS significant variants.
        gene_set: Names of genes in the pathway/ontology being tested.
        background_genes: All genes testable (whole genome, post-QC).

    Returns:
        Dict with:
            odds_ratio, p_value (one-tailed), n_overlap, n_hit, n_set,
            n_background, fold_enrichment, genes_in_set
    """
    hit_set = set(hit_genes)
    gs = set(gene_set)
    bg = set(background_genes)

    # Constrain everything to background
    hit_bg = hit_set & bg
    gs_bg = gs & bg

    a = len(hit_bg & gs_bg)  # hits in gene set
    b = len(hit_bg - gs_bg)  # hits not in gene set
    c = len((bg - hit_bg) & gs_bg)  # background in gene set, not hit
    d = len(bg - hit_bg - gs_bg)  # background not in gene set, not hit

    if HAS_SCIPY:
        odds_ratio_arr, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
        odds_ratio = float(odds_ratio_arr)
    else:
        # Manual odds ratio + hypergeometric p via pure Python
        odds_ratio = (a * d) / (b * c) if (b * c) > 0 else float("inf")
        p_value = _hypergeometric_sf(a, len(hit_bg), len(gs_bg), len(bg))

    expected = len(hit_bg) * len(gs_bg) / max(1, len(bg))
    fold_enrichment = a / expected if expected > 0 else 0.0

    return {
        "n_overlap": a,
        "n_hit": len(hit_bg),
        "n_set": len(gs_bg),
        "n_background": len(bg),
        "odds_ratio": round(odds_ratio, 4),
        "p_value": p_value,
        "fold_enrichment": round(fold_enrichment, 4),
        "genes_in_set": sorted(hit_bg & gs_bg),
    }


def bh_correction(
    p_values: List[float],
    alpha: float = 0.05,
) -> Tuple[List[float], List[bool]]:
    """Benjamini-Hochberg FDR correction.

    Args:
        p_values: Raw p-values.
        alpha: FDR threshold.

    Returns:
        Tuple of (q_values, significant_flags).
    """
    n = len(p_values)
    if n == 0:
        return [], []
    order = sorted(range(n), key=lambda i: p_values[i])
    q_values = [1.0] * n
    max_q = 1.0
    for rank, idx in enumerate(reversed(order)):
        q = p_values[idx] * n / (n - rank)
        max_q = min(max_q, q)
        q_values[idx] = max_q
    significant = [q <= alpha for q in q_values]
    return q_values, significant


def go_term_info(go_id: str) -> Dict[str, Any]:
    """Fetch GO term name and definition from QuickGO REST API.

    Args:
        go_id: GO term ID e.g. "GO:0006915".

    Returns:
        Dict with name, definition, aspect. Empty dict on failure.
    """
    url = f"{QUICKGO_BASE}/ontology/go/terms/{go_id}"
    try:
        resp = requests.get(
            url, timeout=QUICKGO_TIMEOUT, headers={"Accept": "application/json"}
        )
        resp.raise_for_status()
        data = resp.json()
        results = data.get("results", [{}])
        if results:
            r = results[0]
            return {
                "go_id": go_id,
                "name": r.get("name", ""),
                "definition": r.get("definition", {}).get("text", ""),
                "aspect": r.get("aspect", ""),
            }
    except Exception as e:
        logger.debug(f"QuickGO API error for {go_id}: {e}")
    return {}


def pathway_enrichment_from_annotations(
    annotated_hits: List[Dict[str, Any]],
    all_annotated_genes: List[str],
    gene_sets: Optional[Dict[str, List[str]]] = None,
) -> List[Dict[str, Any]]:
    """Run Fisher enrichment tests across multiple gene sets.

    Args:
        annotated_hits: Output of annotate_top_hits_ensembl() — each entry
            has 'nearby_genes' list with gene_name.
        all_annotated_genes: Background genes (all tested genes with annotation).
        gene_sets: Dict of {pathway_name: [gene_names]}. If None, uses a minimal
            built-in set (metabolic, stress response, immune).

    Returns:
        List of enrichment results sorted by p-value, with BH q-values added.
    """
    if gene_sets is None:
        # Minimal built-in gene sets as fallback
        gene_sets = {
            "Metabolic_process_broad": [
                "ACO1",
                "ALDH2",
                "CYP6A2",
                "FASN",
                "GSTA1",
                "HADH",
                "IDH1",
                "LDHA",
                "MDH2",
                "PGK1",
                "PFKM",
                "PKM",
                "SDH",
            ],
            "Stress_response": [
                "HSP70",
                "HSP90",
                "HSPA5",
                "DNAJB1",
                "ATF4",
                "DDIT3",
                "SOD1",
                "SOD2",
                "CAT",
                "GPX1",
                "PRX2",
                "TXN",
            ],
            "Immune_signaling": [
                "RELISH",
                "DORSAL",
                "STAT92E",
                "TOLL",
                "MYD88",
                "IRAK1",
                "TRAF6",
                "NFKB1",
                "TNF",
                "IL6",
                "IL1B",
            ],
        }

    # Collect hit genes from annotations
    hit_genes = set()
    for ann in annotated_hits:
        for g in ann.get("nearby_genes", []):
            name = g.get("gene_name", "")
            if name:
                hit_genes.add(name)

    if not hit_genes or not all_annotated_genes:
        logger.info("Pathway enrichment: no hit genes or background to test")
        return []

    results = []
    for pathway_name, gene_list in gene_sets.items():
        r = fisher_enrichment_test(
            hit_genes=list(hit_genes),
            gene_set=gene_list,
            background_genes=all_annotated_genes,
        )
        r["pathway"] = pathway_name
        r["n_gene_set_total"] = len(gene_list)
        results.append(r)

    # BH correction
    p_vals = [r["p_value"] for r in results]
    q_vals, sig_flags = bh_correction(p_vals)
    for r, q, sig in zip(results, q_vals, sig_flags):
        r["q_value"] = round(q, 6)
        r["significant_fdr05"] = sig

    results.sort(key=lambda r: r["p_value"])
    n_sig = sum(1 for r in results if r["significant_fdr05"])
    logger.info(
        f"Pathway enrichment: {len(results)} sets tested, "
        f"{n_sig} significant at FDR<0.05"
    )
    return results


# ── Pure-Python hypergeometric fallback ────────────────────────────────────────


def _hypergeometric_sf(k: int, n: int, K: int, N: int) -> float:
    """P(X >= k) for hypergeometric(N, K, n) i.e. one-tailed enrichment p."""
    # P(X = x) = C(K,x)*C(N-K,n-x)/C(N,n)
    p_val = 0.0
    for x in range(k, min(n, K) + 1):
        p_val += _log_hypergeom_pmf(x, n, K, N)
    return min(1.0, math.exp(p_val))  # type: ignore


def _log_hypergeom_pmf(k: int, n: int, K: int, N: int) -> float:
    return _log_comb(K, k) + _log_comb(N - K, n - k) - _log_comb(N, n)


def _log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return _log_factorial(n) - _log_factorial(k) - _log_factorial(n - k)


def _log_factorial(n: int) -> float:
    if n <= 1:
        return 0.0
    return sum(math.log(i) for i in range(2, n + 1))
