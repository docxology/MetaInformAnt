"""UniProt database integration utilities.

This module provides tools for accessing and parsing UniProt protein data,
including sequence retrieval, annotations, and functional information.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import requests

from metainformant.core import logging

logger = logging.get_logger(__name__)


def fetch_uniprot_record(uniprot_id: str) -> Dict[str, Any]:
    """Fetch complete UniProt record for a protein.

    Args:
        uniprot_id: UniProt accession or entry name

    Returns:
        Dictionary containing UniProt record data

    Raises:
        requests.RequestException: If API request fails

    Example:
        >>> # This would fetch real UniProt data
        >>> # record = fetch_uniprot_record("P12345")
        >>> # "sequence" in record
        >>> # True
    """
    base_url = "https://www.uniprot.org/uniprotkb"
    url = f"{base_url}/{uniprot_id}.json"

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()

        data = response.json()

        # Extract relevant information
        record = {
            'accession': data.get('primaryAccession'),
            'entry_name': data.get('uniProtkbId'),
            'protein_name': data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value'),
            'organism': data.get('organism', {}).get('scientificName'),
            'taxon_id': data.get('organism', {}).get('taxonId'),
            'sequence': data.get('sequence', {}).get('value'),
            'length': data.get('sequence', {}).get('length'),
            'gene_name': data.get('genes', [{}])[0].get('geneName', {}).get('value') if data.get('genes') else None,
            'function': data.get('comments', [{}])[0].get('texts', [{}])[0].get('value') if data.get('comments') else None,
            'subcellular_location': _extract_subcellular_location(data),
            'domains': _extract_domains(data),
            'ptms': _extract_ptms(data)
        }

        logger.info(f"Fetched UniProt record for {uniprot_id}")
        return record

    except requests.RequestException as e:
        logger.error(f"Failed to fetch UniProt record for {uniprot_id}: {e}")
        raise


def _extract_subcellular_location(data: Dict[str, Any]) -> List[str]:
    """Extract subcellular location information."""
    locations = []

    comments = data.get('comments', [])
    for comment in comments:
        if comment.get('commentType') == 'SUBCELLULAR LOCATION':
            for location in comment.get('subcellularLocations', []):
                location_name = location.get('location', {}).get('value')
                if location_name:
                    locations.append(location_name)

    return locations


def _extract_domains(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Extract domain information."""
    domains = []

    features = data.get('features', [])
    for feature in features:
        if feature.get('type') == 'Domain':
            domain = {
                'name': feature.get('description'),
                'start': feature.get('location', {}).get('start', {}).get('value'),
                'end': feature.get('location', {}).get('end', {}).get('value')
            }
            domains.append(domain)

    return domains


def _extract_ptms(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Extract post-translational modification information."""
    ptms = []

    features = data.get('features', [])
    for feature in features:
        if feature.get('type') in ['Modified residue', 'Glycosylation', 'Disulfide bond']:
            ptm = {
                'type': feature.get('type'),
                'description': feature.get('description'),
                'position': feature.get('location', {}).get('position', {}).get('value') or
                           feature.get('location', {}).get('start', {}).get('value')
            }
            ptms.append(ptm)

    return ptms


def fetch_uniprot_fasta(uniprot_id: str) -> Optional[str]:
    """Fetch protein sequence in FASTA format from UniProt.

    Args:
        uniprot_id: UniProt accession

    Returns:
        FASTA sequence string, or None if not found

    Example:
        >>> # This would fetch real FASTA data
        >>> # fasta = fetch_uniprot_fasta("P12345")
        >>> # fasta.startswith(">")
        >>> # True
    """
    url = f"https://www.uniprot.org/uniprotkb/{uniprot_id}.fasta"

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()

        return response.text

    except requests.RequestException as e:
        logger.error(f"Failed to fetch FASTA for {uniprot_id}: {e}")
        return None


def parse_uniprot_fasta_header(header: str) -> Dict[str, str]:
    """Parse UniProt FASTA header to extract metadata.

    Args:
        header: FASTA header line (without >)

    Returns:
        Dictionary with parsed header information

    Example:
        >>> header = "sp|P12345|PROT_HUMAN Protein name OS=Homo sapiens"
        >>> parsed = parse_uniprot_fasta_header(header)
        >>> parsed['accession'] == 'P12345'
        True
    """
    # UniProt FASTA header format: db|accession|entry_name description OS=organism
    parsed = {
        'database': '',
        'accession': '',
        'entry_name': '',
        'description': '',
        'organism': '',
        'gene_name': ''
    }

    if not header:
        return parsed

    # Split by spaces, but handle OS= and GN= specially
    parts = header.split()

    if len(parts) >= 1:
        # First part: db|accession|entry_name
        id_part = parts[0]
        id_components = id_part.split('|')
        if len(id_components) >= 3:
            parsed['database'] = id_components[0]
            parsed['accession'] = id_components[1]
            parsed['entry_name'] = id_components[2]

    # Extract description and organism
    description_parts = []
    i = 1
    while i < len(parts):
        part = parts[i]
        if part.startswith('OS='):
            # Organism
            organism_parts = []
            while i < len(parts) and not parts[i].startswith('GN=') and not parts[i].startswith('PE='):
                organism_parts.append(parts[i])
                i += 1
            parsed['organism'] = ' '.join(organism_parts).replace('OS=', '')
        elif part.startswith('GN='):
            # Gene name
            parsed['gene_name'] = part.replace('GN=', '')
            i += 1
        else:
            description_parts.append(part)
            i += 1

    parsed['description'] = ' '.join(description_parts)

    return parsed


def get_uniprot_annotations(uniprot_id: str) -> List[Dict[str, Any]]:
    """Get functional annotations for a UniProt entry.

    Args:
        uniprot_id: UniProt accession

    Returns:
        List of annotation dictionaries

    Example:
        >>> # This would get real annotations
        >>> # annotations = get_uniprot_annotations("P12345")
        >>> # isinstance(annotations, list)
        >>> # True
    """
    try:
        record = fetch_uniprot_record(uniprot_id)
        annotations = []

        # Extract GO terms
        if 'go_terms' in record:  # Would need to parse from API response
            for go_term in record.get('go_terms', []):
                annotations.append({
                    'type': 'GO',
                    'id': go_term.get('id'),
                    'name': go_term.get('name'),
                    'aspect': go_term.get('aspect')
                })

        # Extract keywords
        if 'keywords' in record:  # Would need to parse from API response
            for keyword in record.get('keywords', []):
                annotations.append({
                    'type': 'Keyword',
                    'name': keyword
                })

        return annotations

    except Exception as e:
        logger.error(f"Failed to get annotations for {uniprot_id}: {e}")
        return []


def search_uniprot_proteins(query: str, max_results: int = 100) -> List[Dict[str, Any]]:
    """Search UniProt database for proteins matching query.

    Args:
        query: Search query
        max_results: Maximum number of results

    Returns:
        List of matching protein records

    Example:
        >>> # This would search UniProt
        >>> # results = search_uniprot_proteins("kinase AND human", max_results=10)
        >>> # isinstance(results, list)
        >>> # True
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        'query': query,
        'format': 'json',
        'size': min(max_results, 500)  # API limit
    }

    try:
        response = requests.get(base_url, params=params, timeout=30)
        response.raise_for_status()

        data = response.json()

        results = []
        for result in data.get('results', [])[:max_results]:
            protein = {
                'accession': result.get('primaryAccession'),
                'entry_name': result.get('uniProtkbId'),
                'protein_name': result.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value'),
                'organism': result.get('organism', {}).get('scientificName'),
                'sequence_length': result.get('sequence', {}).get('length'),
                'gene_name': result.get('genes', [{}])[0].get('geneName', {}).get('value') if result.get('genes') else None
            }
            results.append(protein)

        logger.info(f"Found {len(results)} proteins matching query: {query}")
        return results

    except requests.RequestException as e:
        logger.error(f"UniProt search failed for query '{query}': {e}")
        return []


def get_uniprot_taxonomy_info(taxon_id: int) -> Optional[Dict[str, Any]]:
    """Get taxonomy information from UniProt.

    Args:
        taxon_id: NCBI taxonomy ID

    Returns:
        Taxonomy information dictionary

    Example:
        >>> # This would get taxonomy info
        >>> # taxonomy = get_uniprot_taxonomy_info(9606)
        >>> # taxonomy['scientific_name'] == 'Homo sapiens'
        >>> # True
    """
    url = f"https://www.uniprot.org/taxonomy/{taxon_id}.json"

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()

        data = response.json()

        taxonomy = {
            'taxon_id': data.get('taxonId'),
            'scientific_name': data.get('scientificName'),
            'common_name': data.get('commonName'),
            'rank': data.get('rank'),
            'lineage': data.get('lineage', []),
            'parent_id': data.get('parent', {}).get('taxonId')
        }

        return taxonomy

    except requests.RequestException as e:
        logger.error(f"Failed to get taxonomy info for {taxon_id}: {e}")
        return None


def batch_fetch_uniprot_records(uniprot_ids: List[str]) -> Dict[str, Dict[str, Any]]:
    """Fetch multiple UniProt records in batch.

    Args:
        uniprot_ids: List of UniProt accessions

    Returns:
        Dictionary mapping IDs to records

    Example:
        >>> ids = ["P12345", "P67890"]
        >>> records = batch_fetch_uniprot_records(ids)
        >>> isinstance(records, dict)
        True
    """
    results = {}

    for uniprot_id in uniprot_ids:
        try:
            record = fetch_uniprot_record(uniprot_id)
            results[uniprot_id] = record
        except Exception as e:
            logger.error(f"Failed to fetch record for {uniprot_id}: {e}")
            results[uniprot_id] = None

    return results


def validate_uniprot_accession(accession: str) -> bool:
    """Validate UniProt accession format.

    Args:
        accession: UniProt accession to validate

    Returns:
        True if format is valid

    Example:
        >>> validate_uniprot_accession("P12345")
        True
        >>> validate_uniprot_accession("INVALID")
        False
    """
    import re

    # UniProt accession patterns
    patterns = [
        r'^[A-Z]\d{5}$',      # P12345
        r'^[A-Z]\d{9}$',      # A0A1234567
        r'^[A-Z]\d{3}[A-Z]\d{2}$',  # P123A45 (rare)
    ]

    return any(re.match(pattern, accession) for pattern in patterns)
