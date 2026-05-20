"""QuickGO REST API client for live Gene Ontology annotation.

Retrieves GO term metadata and per-gene annotations directly from the
EMBL-EBI QuickGO REST API without requiring a local OBO file.  All
requests are rate-limited and session-cached.

API base: https://www.ebi.ac.uk/QuickGO/services/

Typical usage::

    from metainformant.ontology.core.go_api import (
        fetch_go_term,
        fetch_gene_go_annotations,
        build_gene_set_from_api,
    )

    term = fetch_go_term("GO:0008150")
    annotations = build_gene_set_from_api(["AMel_GB43025", "VG10"], taxon_id=7460)
"""

from __future__ import annotations

import time
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

_QUICKGO_BASE = "https://www.ebi.ac.uk/QuickGO/services"

# Session-level in-memory cache
_term_cache: dict[str, dict[str, Any]] = {}
_annotation_cache: dict[str, list[dict[str, Any]]] = {}

_DEFAULT_RATE_LIMIT_S = 0.2  # respect QuickGO rate limits


# ---------------------------------------------------------------------------
# Term metadata
# ---------------------------------------------------------------------------


def fetch_go_term(go_id: str, *, rate_limit_s: float = _DEFAULT_RATE_LIMIT_S) -> dict[str, Any]:
    """Fetch metadata for a single GO term from QuickGO REST API.

    Args:
        go_id: GO term ID in ``GO:XXXXXXX`` format.
        rate_limit_s: Minimum delay between API calls.

    Returns:
        Dictionary with keys ``id``, ``name``, ``aspect``, ``definition``,
        ``synonyms``, ``is_obsolete``.  Returns empty dict if the term
        cannot be retrieved.

    Examples:
        >>> term = fetch_go_term("GO:0008150")
        >>> term["name"]
        'biological_process'
    """
    if go_id in _term_cache:
        return _term_cache[go_id]

    try:
        import json
        import urllib.request

        url = f"{_QUICKGO_BASE}/ontology/go/terms/{go_id}"
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        time.sleep(rate_limit_s)
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode())

        results = data.get("results", [])
        if not results:
            return {}

        r = results[0]
        term = {
            "id": r.get("id", go_id),
            "name": r.get("name", ""),
            "aspect": r.get("aspect", ""),
            "definition": r.get("definition", {}).get("text", ""),
            "synonyms": [s.get("name", "") for s in r.get("synonyms", [])],
            "is_obsolete": r.get("isObsolete", False),
        }
        _term_cache[go_id] = term
        logger.debug("Fetched GO term %s: %s (%s)", go_id, term["name"], term["aspect"])
        return term
    except Exception as exc:
        logger.warning("fetch_go_term(%s) failed: %s", go_id, exc)
        return {}


def fetch_go_terms_batch(
    go_ids: list[str],
    *,
    rate_limit_s: float = _DEFAULT_RATE_LIMIT_S,
) -> dict[str, dict[str, Any]]:
    """Fetch metadata for multiple GO terms.

    Args:
        go_ids: List of GO term IDs.
        rate_limit_s: Pause between API calls.

    Returns:
        Dictionary mapping GO term ID to metadata dict.

    Examples:
        >>> terms = fetch_go_terms_batch(["GO:0008150", "GO:0003674"])
        >>> terms["GO:0003674"]["name"]
        'molecular_function'
    """
    result: dict[str, dict[str, Any]] = {}
    for gid in go_ids:
        info = fetch_go_term(gid, rate_limit_s=rate_limit_s)
        if info:
            result[gid] = info
    logger.info("Fetched metadata for %d/%d GO terms", len(result), len(go_ids))
    return result


# ---------------------------------------------------------------------------
# Per-gene annotations
# ---------------------------------------------------------------------------


def fetch_gene_go_annotations(
    gene_id: str,
    *,
    taxon_id: int = 7460,
    go_aspects: list[str] | None = None,
    rate_limit_s: float = _DEFAULT_RATE_LIMIT_S,
    max_results: int = 200,
) -> list[dict[str, Any]]:
    """Fetch GO annotations for a gene from QuickGO.

    Args:
        gene_id: Gene product ID or gene symbol (e.g. ``"AMel_GB43025"``).
        taxon_id: NCBI taxon ID.  Default ``7460`` (*Apis mellifera*).
        go_aspects: Optional filter list: ``["biological_process",
            "molecular_function", "cellular_component"]``.
        rate_limit_s: Pause between API calls.
        max_results: Maximum annotations to retrieve.

    Returns:
        List of annotation dicts with keys ``gene_id``, ``go_id``,
        ``go_name``, ``aspect``, ``evidence_code``, ``reference``.

    Examples:
        >>> anns = fetch_gene_go_annotations("SOD1", taxon_id=7460)
        >>> anns[0]["go_id"]
        'GO:...'
    """
    cache_key = f"{gene_id}:{taxon_id}"
    if cache_key in _annotation_cache:
        return _annotation_cache[cache_key]

    try:
        import json
        import urllib.parse
        import urllib.request

        params: dict[str, str] = {
            "taxonId": str(taxon_id),
            "geneProductId": gene_id,
            "limit": str(min(max_results, 200)),
            "page": "1",
        }
        if go_aspects:
            # QuickGO uses "aspect" filter
            params["aspect"] = ",".join(go_aspects)

        query = urllib.parse.urlencode(params)
        url = f"{_QUICKGO_BASE}/annotation/search?{query}"
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        time.sleep(rate_limit_s)
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode())

        annotations: list[dict[str, Any]] = []
        for ann in data.get("results", []):
            annotations.append(
                {
                    "gene_id": gene_id,
                    "go_id": ann.get("goId", ""),
                    "go_name": ann.get("goName", ""),
                    "aspect": ann.get("goAspect", ""),
                    "evidence_code": ann.get("evidenceCode", ""),
                    "reference": ann.get("reference", ""),
                }
            )

        _annotation_cache[cache_key] = annotations
        logger.debug("Gene %s (taxon %s): %d GO annotations", gene_id, taxon_id, len(annotations))
        return annotations

    except Exception as exc:
        logger.warning("fetch_gene_go_annotations(%s) failed: %s", gene_id, exc)
        return []


# ---------------------------------------------------------------------------
# Strategy A: Per-gene-product annotation lookup (requires UniProt/Ensembl IDs)
# ---------------------------------------------------------------------------


def build_gene_set_from_api(
    gene_list: list[str],
    *,
    taxon_id: int = 7460,
    go_aspects: list[str] | None = None,
    rate_limit_s: float = _DEFAULT_RATE_LIMIT_S,
    min_evidence_codes: list[str] | None = None,
) -> dict[str, set[str]]:
    """Build a GO-term → gene-set mapping via live QuickGO per-gene API.

    This function queries QuickGO for each gene in ``gene_list`` using the
    gene product ID directly.  It works when ``gene_list`` contains recognised
    QuickGO gene product IDs (UniProtKB accessions, Ensembl IDs such as
    ``ENSAMEG00000012345``, or MGI/FlyBase identifiers).

    For Apis mellifera GWAS workflows where only NCBI gene symbols are
    available, use :func:`build_taxon_go_gene_sets` instead (taxon-level
    reverse lookup).

    Args:
        gene_list: Gene product IDs recognisable by QuickGO.
        taxon_id: NCBI taxon ID (used as filter alongside ``geneProductId``).
        go_aspects: Optional aspect filter.
        rate_limit_s: Pause between API calls.
        min_evidence_codes: Optional evidence code filter (e.g. ``["EXP", "IDA"]``).

    Returns:
        ``dict[go_id, set[gene_id]]``.

    Examples:
        >>> gene_sets = build_gene_set_from_api(["UniProtKB:Q8MKI0"], taxon_id=7460)
        >>> isinstance(gene_sets, dict)
        True
    """
    allowed_ev = set(min_evidence_codes) if min_evidence_codes else None
    gene_sets: dict[str, set[str]] = {}

    for gene in gene_list:
        anns = fetch_gene_go_annotations(
            gene,
            taxon_id=taxon_id,
            go_aspects=go_aspects,
            rate_limit_s=rate_limit_s,
        )
        for ann in anns:
            if allowed_ev and ann.get("evidence_code", "") not in allowed_ev:
                continue
            go_id = ann.get("go_id", "")
            if go_id:
                gene_sets.setdefault(go_id, set()).add(gene)

    total_terms = len(gene_sets)
    total_genes = len(gene_list)
    logger.info(
        "build_gene_set_from_api: %d genes -> %d GO terms (taxon %d)",
        total_genes,
        total_terms,
        taxon_id,
    )
    return gene_sets


# ---------------------------------------------------------------------------
# Strategy B: Taxon-wide GO annotation reverse-lookup (works with any gene names)
# ---------------------------------------------------------------------------


def build_taxon_go_gene_sets(
    taxon_id: int = 7460,
    *,
    go_aspects: list[str] | None = None,
    max_pages: int = 50,
    results_per_page: int = 200,
    rate_limit_s: float = _DEFAULT_RATE_LIMIT_S,
    min_evidence_codes: list[str] | None = None,
) -> dict[str, set[str]]:
    """Build comprehensive GO-term → gene-set mapping for an entire taxon.

    Paginates through all QuickGO annotations for a given taxon to build
    the complete mapping of GO terms to annotated gene symbols.  This is
    the recommended approach when working with NCBI gene symbols rather than
    UniProt IDs, since QuickGO gene product IDs often use UniProtKB accessions.

    The function is intended for use as the background gene-set reference in
    ORA.  Filter the returned dict to terms whose gene set overlaps with your
    hit gene list for enrichment testing.

    Args:
        taxon_id: NCBI taxon ID (default ``7460`` — *Apis mellifera*).
        go_aspects: Optional list of GO aspects to include
            (``"biological_process"``, ``"molecular_function"``,
            ``"cellular_component"``).
        max_pages: Maximum number of pages to retrieve (cap on total results;
            one page = ``results_per_page`` annotations).
        results_per_page: Annotations per page (max 200 for QuickGO).
        rate_limit_s: Pause between API calls.
        min_evidence_codes: Optional evidence code filter.

    Returns:
        ``dict[go_id, set[gene_product_id]]`` — complete GO-term gene-set map
        for the taxon.  Gene names are QuickGO gene product IDs (UniProtKB
        accessions).

    Examples:
        >>> gene_sets = build_taxon_go_gene_sets(taxon_id=7460, max_pages=2)
        >>> len(gene_sets) > 0
        True
    """
    import json
    import urllib.parse
    import urllib.request

    allowed_ev = set(min_evidence_codes) if min_evidence_codes else None
    gene_sets: dict[str, set[str]] = {}
    n_annotations = 0
    n_pages = 0

    # Build aspect filter
    _ASPECT_SHORT = {
        "biological_process": "biological_process",
        "molecular_function": "molecular_function",
        "cellular_component": "cellular_component",
        "BP": "biological_process",
        "MF": "molecular_function",
        "CC": "cellular_component",
    }

    base_params: dict[str, str] = {
        "taxonId": str(taxon_id),
        "geneProductType": "protein",
        "limit": str(min(results_per_page, 200)),
    }
    if go_aspects:
        normalised = [_ASPECT_SHORT.get(a, a) for a in go_aspects]
        base_params["goAspect"] = ",".join(normalised)

    for page in range(1, max_pages + 1):
        params = dict(base_params)
        params["page"] = str(page)
        query = urllib.parse.urlencode(params)
        url = f"{_QUICKGO_BASE}/annotation/search?{query}"

        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            time.sleep(rate_limit_s)
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = json.loads(resp.read().decode())
        except Exception as exc:
            logger.warning("build_taxon_go_gene_sets page %d failed: %s", page, exc)
            break

        results = data.get("results", [])
        if not results:
            break  # exhausted

        for ann in results:
            ev = ann.get("evidenceCode", "")
            if allowed_ev and ev not in allowed_ev:
                continue
            go_id = ann.get("goId", "")
            gene_product = ann.get("geneProductId", "")
            if go_id and gene_product:
                gene_sets.setdefault(go_id, set()).add(gene_product)
                n_annotations += 1

        n_pages += 1
        page_info = data.get("pageInfo", {})
        total_pages = page_info.get("total", 1)
        if page >= total_pages:
            break

    logger.info(
        "build_taxon_go_gene_sets: taxon %d → %d GO terms from %d annotations " "(%d pages fetched)",
        taxon_id,
        len(gene_sets),
        n_annotations,
        n_pages,
    )
    return gene_sets


# ---------------------------------------------------------------------------
# Strategy C: Fetch GO annotations for a list of UniProt accessions from NCBI mapping
# ---------------------------------------------------------------------------


def map_symbols_to_uniprot(
    ncbi_gene_symbols: list[str],
    taxon_id: int = 7460,
    *,
    rate_limit_s: float = _DEFAULT_RATE_LIMIT_S,
) -> dict[str, list[str]]:
    """Map NCBI gene symbols to UniProtKB accessions using the UniProt search API.

    Args:
        ncbi_gene_symbols: List of gene symbols (e.g. ``["LOC411641", "NRT1"]``).
        taxon_id: NCBI taxon ID.
        rate_limit_s: Pause between batch API calls.

    Returns:
        ``dict[gene_symbol, [uniprot_accession, ...]]``.  Symbols that cannot
        be mapped are omitted.

    Examples:
        >>> mapping = map_symbols_to_uniprot(["LOC411641"], taxon_id=7460)
        >>> isinstance(mapping, dict)
        True
    """
    import json
    import urllib.parse
    import urllib.request

    if not ncbi_gene_symbols:
        return {}

    mapping: dict[str, list[str]] = {}

    for sym in ncbi_gene_symbols:
        query = urllib.parse.urlencode(
            {"query": f"({sym}) AND (taxonomy_id:{taxon_id})", "format": "json", "fields": "accession"}
        )
        url = f"https://rest.uniprot.org/uniprotkb/search?{query}"

        try:
            req = urllib.request.Request(url)
            time.sleep(rate_limit_s)
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = json.loads(resp.read().decode())

            accessions = [res.get("primaryAccession") for res in data.get("results", []) if res.get("primaryAccession")]
            if accessions:
                mapping[sym] = accessions
        except Exception as exc:
            logger.warning("map_symbols_to_uniprot failed for %s: %s", sym, exc)

    logger.info(
        "map_symbols_to_uniprot: %d/%d symbols mapped (taxon %d)",
        len(mapping),
        len(ncbi_gene_symbols),
        taxon_id,
    )

    return mapping


# ---------------------------------------------------------------------------
# Utility: annotate a list of GO IDs with names
# ---------------------------------------------------------------------------


def annotate_go_ids(
    go_ids: list[str],
    *,
    rate_limit_s: float = _DEFAULT_RATE_LIMIT_S,
) -> dict[str, str]:
    """Map GO IDs to human-readable names.

    Args:
        go_ids: List of GO IDs.
        rate_limit_s: Pause between API calls.

    Returns:
        ``dict[go_id, name]``.

    Examples:
        >>> names = annotate_go_ids(["GO:0008150"])
        >>> names["GO:0008150"]
        'biological_process'
    """
    names: dict[str, str] = {}
    for go_id in go_ids:
        info = fetch_go_term(go_id, rate_limit_s=rate_limit_s)
        if info:
            names[go_id] = info.get("name", go_id)
        else:
            names[go_id] = go_id
    return names
