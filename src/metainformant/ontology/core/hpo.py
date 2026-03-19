"""Human Phenotype Ontology (HPO) client and GWAS trait mapper.

Queries the HPO REST API (https://hpo.jax.org/api/) to look up phenotype
terms, search by free-text, and map GWAS phenotype labels to HPO IDs.

No local OBO file is required.

Typical usage::

    from metainformant.ontology.core.hpo import map_phenotype_to_hpo, fetch_hpo_term

    hp_ids = map_phenotype_to_hpo("honey production")
    for hp_id in hp_ids:
        term = fetch_hpo_term(hp_id)
        print(hp_id, term["name"])
"""

from __future__ import annotations

import time
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

_OLS_BASE = "https://www.ebi.ac.uk/ols4/api"
_DEFAULT_RATE_LIMIT_S = 0.2

# Session-level caches
_term_cache: dict[str, dict[str, Any]] = {}
_search_cache: dict[str, list[dict[str, Any]]] = {}


# ---------------------------------------------------------------------------
# Term lookup
# ---------------------------------------------------------------------------


def fetch_hpo_term(hp_id: str, *, rate_limit_s: float = _DEFAULT_RATE_LIMIT_S) -> dict[str, Any]:
    """Fetch metadata for an HPO term using EBI OLS API.

    Args:
        hp_id: HPO term ID in ``HP:XXXXXXX`` format.
        rate_limit_s: Pause before API call.

    Returns:
        Dict with ``id``, ``name``, ``definition``, ``synonyms``,
        ``category``.  Returns empty dict on failure.

    Examples:
        >>> term = fetch_hpo_term("HP:0000118")
        >>> term["name"]
        'Phenotypic abnormality'
    """
    if hp_id in _term_cache:
        return _term_cache[hp_id]

    try:
        import urllib.request
        import json

        # OLS uses short form (e.g. HP_0000118)
        short_form = hp_id.replace(":", "_")
        url = f"{_OLS_BASE}/ontologies/hp/terms?short_form={short_form}"
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        time.sleep(rate_limit_s)
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode())
        
        embedded = data.get("_embedded", {})
        terms = embedded.get("terms", [])
        if not terms:
            return {}
        
        term_data = terms[0]
        
        definition = ""
        if term_data.get("description"):
            definition = term_data["description"][0]
            
        synonyms = term_data.get("synonyms") or []

        term = {
            "id": hp_id,
            "name": term_data.get("label", ""),
            "definition": definition,
            "synonyms": synonyms,
            "category": "", # OLS doesn't directly map 'category' in the same way
        }
        _term_cache[hp_id] = term
        logger.debug("Fetched HPO term %s: %s", hp_id, term["name"])
        return term
    except Exception as exc:
        logger.warning("fetch_hpo_term(%s) failed: %s", hp_id, exc)
        return {}


# ---------------------------------------------------------------------------
# Free-text search
# ---------------------------------------------------------------------------


def search_hpo_terms(
    query_str: str,
    *,
    max_results: int = 10,
    rate_limit_s: float = _DEFAULT_RATE_LIMIT_S,
) -> list[dict[str, Any]]:
    """Search HPO terms by free text using EBI OLS API.

    Args:
        query_str: Search query (e.g. ``"honey production"``).
        max_results: Maximum results to return.
        rate_limit_s: Pause before API call.

    Returns:
        List of result dicts with ``id``, ``name``, ``definition``,
        ``score`` (relevance).

    Examples:
        >>> results = search_hpo_terms("cardiac arrhythmia")
        >>> results[0]["name"]
        'arrhythmia'
    """
    cache_key = f"{query_str}:{max_results}"
    if cache_key in _search_cache:
        return _search_cache[cache_key]

    try:
        import urllib.request
        import urllib.parse
        import json

        params = {
            "q": query_str, 
            "ontology": "hp", 
            "rows": str(max_results)
        }
        url = f"{_OLS_BASE}/search?{urllib.parse.urlencode(params)}"
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        time.sleep(rate_limit_s)
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode())

        results = []
        for hit in data.get("response", {}).get("docs", []):
            results.append({
                "id": hit.get("obo_id", ""),
                "name": hit.get("label", ""),
                "definition": hit.get("description", [""])[0] if hit.get("description") else "",
                "score": hit.get("score", 0.0),
            })

        _search_cache[cache_key] = results
        logger.info("HPO search '%s': %d results", query_str, len(results))
        return results

    except Exception as exc:
        logger.warning("search_hpo_terms('%s') failed: %s", query_str, exc)
        return []


# ---------------------------------------------------------------------------
# Phenotype-to-HPO mapper
# ---------------------------------------------------------------------------


def map_phenotype_to_hpo(
    phenotype_label: str,
    *,
    max_results: int = 5,
    rate_limit_s: float = _DEFAULT_RATE_LIMIT_S,
) -> list[str]:
    """Map a GWAS phenotype label to HPO term IDs.

    Also attempts synonym expansion using known agricultural/quantitative
    trait mappings before falling back to live API search.

    Args:
        phenotype_label: Human-readable phenotype (e.g. ``"honey yield"``).
        max_results: Maximum HPO terms to return.
        rate_limit_s: Pause before API calls.

    Returns:
        List of HP IDs (may be empty if no match found).

    Examples:
        >>> ids = map_phenotype_to_hpo("body mass index")
        >>> ids[0].startswith("HP:")
        True
    """
    # Local synonym map for common quantitative/agricultural traits
    _SYNONYM_MAP: dict[str, list[str]] = {
        "honey yield": ["HP:0000118"],       # root — no direct HPO equivalent; use for context
        "honey production": ["HP:0000118"],
        "body weight": ["HP:0004324"],        # increased body weight
        "body mass index": ["HP:0002860"],    # BMI-related
        "blood pressure": ["HP:0000822"],     # hypertension
        "height": ["HP:0000256"],             # macrocephaly / stature
        "stature": ["HP:0004322"],
        "glucose": ["HP:0011013"],            # blood glucose
        "cholesterol": ["HP:0003077"],        # hypercholesterolaemia
    }

    label_lower = phenotype_label.lower().strip()

    for key, hp_ids in _SYNONYM_MAP.items():
        if key in label_lower:
            logger.info(
                "map_phenotype_to_hpo: '%s' matched synonym -> %s",
                phenotype_label,
                hp_ids,
            )
            return hp_ids

    # Fall back to live API search
    results = search_hpo_terms(phenotype_label, max_results=max_results, rate_limit_s=rate_limit_s)
    hp_ids_found = [r["id"] for r in results if r.get("id", "").startswith("HP:")]

    if hp_ids_found:
        logger.info(
            "map_phenotype_to_hpo: '%s' -> %d HP IDs via API search",
            phenotype_label,
            len(hp_ids_found),
        )
    else:
        logger.warning(
            "map_phenotype_to_hpo: no HPO terms found for '%s'", phenotype_label
        )

    return hp_ids_found


# ---------------------------------------------------------------------------
# Batch term description
# ---------------------------------------------------------------------------


def describe_hpo_terms(
    hp_ids: list[str],
    *,
    rate_limit_s: float = _DEFAULT_RATE_LIMIT_S,
) -> list[dict[str, Any]]:
    """Fetch descriptions for a list of HPO IDs.

    Args:
        hp_ids: List of HP IDs.
        rate_limit_s: Pause between API calls.

    Returns:
        List of term metadata dicts (skips failed lookups).

    Examples:
        >>> terms = describe_hpo_terms(["HP:0000001", "HP:0004324"])
        >>> terms[0]["name"]
        'All'
    """
    terms = []
    for hp_id in hp_ids:
        info = fetch_hpo_term(hp_id, rate_limit_s=rate_limit_s)
        if info:
            terms.append(info)
    logger.info("describe_hpo_terms: %d/%d fetched", len(terms), len(hp_ids))
    return terms
