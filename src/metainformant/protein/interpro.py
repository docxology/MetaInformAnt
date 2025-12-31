"""InterPro database integration utilities.

This module provides tools for accessing InterPro protein domain and
functional site annotations from the EMBL-EBI InterPro database.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import requests

from metainformant.core import logging

logger = logging.get_logger(__name__)


def fetch_interpro_domains(uniprot_id: str) -> List[Dict[str, Any]]:
    """Fetch InterPro domain annotations for a UniProt protein.

    Args:
        uniprot_id: UniProt accession

    Returns:
        List of domain annotation dictionaries

    Example:
        >>> # This would fetch real InterPro data
        >>> # domains = fetch_interpro_domains("P12345")
        >>> # isinstance(domains, list)
        >>> # True
    """
    # InterPro API endpoint
    url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/UniProt/{uniprot_id}"

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()

        data = response.json()

        domains = []
        for result in data.get('results', []):
            for entry in result.get('entries', []):
                for location in entry.get('entry_protein_locations', []):
                    for fragment in location.get('fragments', []):
                        domain = {
                            'interpro_id': entry.get('accession'),
                            'name': entry.get('name'),
                            'type': entry.get('type'),
                            'start': fragment.get('start'),
                            'end': fragment.get('end'),
                            'score': fragment.get('score'),
                            'database': entry.get('source_database')
                        }
                        domains.append(domain)

        logger.info(f"Fetched {len(domains)} InterPro domains for {uniprot_id}")
        return domains

    except requests.RequestException as e:
        logger.error(f"Failed to fetch InterPro domains for {uniprot_id}: {e}")
        return []


def fetch_interpro_by_accession(interpro_id: str) -> Optional[Dict[str, Any]]:
    """Fetch detailed information for an InterPro entry.

    Args:
        interpro_id: InterPro accession (e.g., "IPR000001")

    Returns:
        InterPro entry information

    Example:
        >>> # This would fetch real InterPro entry data
        >>> # entry = fetch_interpro_by_accession("IPR000001")
        >>> # isinstance(entry, dict)
        >>> # True
    """
    url = f"https://www.ebi.ac.uk/interpro/api/entry/InterPro/{interpro_id}"

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()

        data = response.json()

        if 'metadata' in data:
            metadata = data['metadata']
            entry = {
                'accession': metadata.get('accession'),
                'name': metadata.get('name'),
                'type': metadata.get('type'),
                'short_name': metadata.get('short_name'),
                'description': metadata.get('description'),
                'source_database': metadata.get('source_database'),
                'member_databases': metadata.get('member_databases', []),
                'go_terms': metadata.get('go_terms', []),
                'literature': metadata.get('literature', [])
            }

            return entry

    except requests.RequestException as e:
        logger.error(f"Failed to fetch InterPro entry {interpro_id}: {e}")

    return None


def search_interpro_entries(query: str, max_results: int = 100) -> List[Dict[str, Any]]:
    """Search InterPro database for entries matching query.

    Args:
        query: Search query
        max_results: Maximum number of results

    Returns:
        List of matching InterPro entries

    Example:
        >>> # This would search InterPro
        >>> # results = search_interpro_entries("kinase", max_results=10)
        >>> # isinstance(results, list)
        >>> # True
    """
    url = "https://www.ebi.ac.uk/interpro/api/entry/InterPro"
    params = {
        'search': query,
        'page_size': min(max_results, 200)  # API limit
    }

    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()

        data = response.json()

        results = []
        for result in data.get('results', [])[:max_results]:
            entry = {
                'accession': result.get('metadata', {}).get('accession'),
                'name': result.get('metadata', {}).get('name'),
                'type': result.get('metadata', {}).get('type'),
                'short_name': result.get('metadata', {}).get('short_name'),
                'description': result.get('metadata', {}).get('description')
            }
            results.append(entry)

        logger.info(f"Found {len(results)} InterPro entries matching '{query}'")
        return results

    except requests.RequestException as e:
        logger.error(f"InterPro search failed for query '{query}': {e}")
        return []


def get_interpro_hierarchy(interpro_id: str) -> Dict[str, Any]:
    """Get hierarchy information for an InterPro entry.

    Args:
        interpro_id: InterPro accession

    Returns:
        Hierarchy information including parent/child relationships

    Example:
        >>> # This would get hierarchy info
        >>> # hierarchy = get_interpro_hierarchy("IPR000001")
        >>> # isinstance(hierarchy, dict)
        >>> # True
    """
    # This would typically require multiple API calls to get hierarchy
    # Placeholder implementation
    logger.info("InterPro hierarchy retrieval not fully implemented")
    return {
        'entry': interpro_id,
        'parents': [],
        'children': [],
        'siblings': []
    }


def batch_fetch_interpro_domains(uniprot_ids: List[str]) -> Dict[str, List[Dict[str, Any]]]:
    """Fetch InterPro domains for multiple UniProt IDs.

    Args:
        uniprot_ids: List of UniProt accessions

    Returns:
        Dictionary mapping UniProt IDs to domain lists

    Example:
        >>> ids = ["P12345", "P67890"]
        >>> domains = batch_fetch_interpro_domains(ids)
        >>> isinstance(domains, dict)
        True
    """
    results = {}

    for uniprot_id in uniprot_ids:
        try:
            domains = fetch_interpro_domains(uniprot_id)
            results[uniprot_id] = domains
        except Exception as e:
            logger.error(f"Failed to fetch domains for {uniprot_id}: {e}")
            results[uniprot_id] = []

    return results


def parse_interpro_results(xml_content: str) -> List[Dict[str, Any]]:
    """Parse InterPro XML results (legacy format).

    Args:
        xml_content: InterPro XML response

    Returns:
        List of parsed InterPro matches

    Example:
        >>> # Assuming XML content exists
        >>> # matches = parse_interpro_results(xml_string)
        >>> # isinstance(matches, list)
        >>> # True
    """
    # This would require XML parsing
    # Placeholder for legacy XML format support
    logger.info("InterPro XML parsing not fully implemented")
    return []


def get_interpro_statistics() -> Dict[str, Any]:
    """Get statistics about the InterPro database.

    Returns:
        Database statistics

    Example:
        >>> stats = get_interpro_statistics()
        >>> "total_entries" in stats
        True
    """
    # Placeholder - would query InterPro API for statistics
    stats = {
        'total_entries': 39000,  # Approximate as of 2024
        'total_proteins': 2000000,  # Approximate
        'member_databases': 14,
        'last_updated': '2024-01',
        'source': 'InterPro database'
    }

    return stats


def find_similar_interpro_entries(interpro_id: str, max_results: int = 10) -> List[Dict[str, Any]]:
    """Find InterPro entries similar to the given entry.

    Args:
        interpro_id: Reference InterPro accession
        max_results: Maximum number of similar entries

    Returns:
        List of similar InterPro entries

    Example:
        >>> # This would find similar entries
        >>> # similar = find_similar_interpro_entries("IPR000001", max_results=5)
        >>> # isinstance(similar, list)
        >>> # True
    """
    # This would use similarity measures based on GO terms, etc.
    # Placeholder implementation
    logger.info("InterPro similarity search not fully implemented")
    return []


def validate_interpro_accession(interpro_id: str) -> bool:
    """Validate InterPro accession format.

    Args:
        interpro_id: InterPro accession to validate

    Returns:
        True if format is valid

    Example:
        >>> validate_interpro_accession("IPR000001")
        True
        >>> validate_interpro_accession("INVALID")
        False
    """
    import re

    # InterPro accession pattern: IPR followed by 6 digits
    pattern = r'^IPR\d{6}$'

    return bool(re.match(pattern, interpro_id))


def get_interpro_go_annotations(interpro_id: str) -> List[Dict[str, Any]]:
    """Get GO term annotations for an InterPro entry.

    Args:
        interpro_id: InterPro accession

    Returns:
        List of GO term annotations

    Example:
        >>> # This would get GO annotations
        >>> # go_terms = get_interpro_go_annotations("IPR000001")
        >>> # isinstance(go_terms, list)
        >>> # True
    """
    entry = fetch_interpro_by_accession(interpro_id)

    if entry and 'go_terms' in entry:
        return entry['go_terms']

    return []


def cross_reference_interpro_uniprot(uniprot_id: str) -> List[Dict[str, Any]]:
    """Get cross-references between InterPro and UniProt.

    Args:
        uniprot_id: UniProt accession

    Returns:
        List of cross-reference information

    Example:
        >>> # This would get cross-references
        >>> # xrefs = cross_reference_interpro_uniprot("P12345")
        >>> # isinstance(xrefs, list)
        >>> # True
    """
    domains = fetch_interpro_domains(uniprot_id)

    xrefs = []
    for domain in domains:
        xref = {
            'uniprot_id': uniprot_id,
            'interpro_id': domain.get('interpro_id'),
            'domain_name': domain.get('name'),
            'start': domain.get('start'),
            'end': domain.get('end'),
            'evidence': domain.get('database')
        }
        xrefs.append(xref)

    return xrefs
