"""Entrez database utilities for NCBI data retrieval.

This module provides utilities for accessing NCBI databases through the
Entrez programming interface, including GenBank record retrieval,
sequence downloads, and metadata queries.
"""

from __future__ import annotations

import time
from typing import Any, Dict, List, Optional

import requests

from metainformant.core import logging

logger = logging.get_logger(__name__)


class EntrezClient:
    """Client for NCBI Entrez utilities."""

    def __init__(self, api_key: Optional[str] = None, email: Optional[str] = None):
        """Initialize Entrez client.

        Args:
            api_key: NCBI API key for higher rate limits
            email: Email address (required by NCBI)
        """
        if not email:
            logger.warning("No email provided. NCBI requires email for Entrez usage.")

        self.api_key = api_key
        self.email = email
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.session = requests.Session()

        # Set headers
        headers = {'User-Agent': 'metainformant/0.2.0'}
        if email:
            headers['From'] = email
        self.session.headers.update(headers)

    def search(self, db: str, query: str, max_results: int = 100) -> Dict[str, Any]:
        """Perform Entrez search.

        Args:
            db: Database to search (nucleotide, protein, etc.)
            query: Search query
            max_results: Maximum number of results

        Returns:
            Search results dictionary
        """
        params = {
            'db': db,
            'term': query,
            'retmax': min(max_results, 10000),  # NCBI limit
            'retmode': 'json',
            'usehistory': 'y'
        }

        if self.api_key:
            params['api_key'] = self.api_key
        if self.email:
            params['email'] = self.email

        try:
            response = self.session.get(f"{self.base_url}/esearch.fcgi", params=params)
            response.raise_for_status()

            return response.json()

        except requests.RequestException as e:
            logger.error(f"Entrez search failed: {e}")
            return {}

    def fetch(self, db: str, ids: List[str], rettype: str = 'fasta',
              retmode: str = 'text') -> str:
        """Fetch records from Entrez database.

        Args:
            db: Database name
            ids: List of IDs to fetch
            rettype: Retrieval type (fasta, gb, etc.)
            retmode: Retrieval mode (text, xml, etc.)

        Returns:
            Fetched data as string
        """
        params = {
            'db': db,
            'id': ','.join(ids),
            'rettype': rettype,
            'retmode': retmode
        }

        if self.api_key:
            params['api_key'] = self.api_key
        if self.email:
            params['email'] = self.email

        # Respect NCBI rate limits (3 requests/second without API key)
        if not self.api_key:
            time.sleep(0.4)

        try:
            response = self.session.get(f"{self.base_url}/efetch.fcgi", params=params)
            response.raise_for_status()

            return response.text

        except requests.RequestException as e:
            logger.error(f"Entrez fetch failed: {e}")
            return ""

    def summary(self, db: str, ids: List[str]) -> Dict[str, Any]:
        """Get summary information for records.

        Args:
            db: Database name
            ids: List of IDs

        Returns:
            Summary data dictionary
        """
        params = {
            'db': db,
            'id': ','.join(ids),
            'retmode': 'json'
        }

        if self.api_key:
            params['api_key'] = self.api_key
        if self.email:
            params['email'] = self.email

        try:
            response = self.session.get(f"{self.base_url}/esummary.fcgi", params=params)
            response.raise_for_status()

            return response.json()

        except requests.RequestException as e:
            logger.error(f"Entrez summary failed: {e}")
            return {}

    def link(self, db: str, ids: List[str], linkname: str) -> Dict[str, Any]:
        """Find linked records in other databases.

        Args:
            db: Source database
            ids: Source IDs
            linkname: Link name (e.g., nucleotide_protein)

        Returns:
            Link data dictionary
        """
        params = {
            'db': db,
            'id': ','.join(ids),
            'linkname': linkname,
            'retmode': 'json'
        }

        if self.api_key:
            params['api_key'] = self.api_key
        if self.email:
            params['email'] = self.email

        try:
            response = self.session.get(f"{self.base_url}/elink.fcgi", params=params)
            response.raise_for_status()

            return response.json()

        except requests.RequestException as e:
            logger.error(f"Entrez link failed: {e}")
            return {}


def search_genbank(query: str, max_results: int = 100,
                  api_key: Optional[str] = None, email: Optional[str] = None) -> List[Dict[str, Any]]:
    """Search GenBank nucleotide database.

    Args:
        query: Search query
        max_results: Maximum results
        api_key: NCBI API key
        email: Email address

    Returns:
        List of GenBank records
    """
    client = EntrezClient(api_key=api_key, email=email)

    search_result = client.search('nucleotide', query, max_results)

    if 'esearchresult' not in search_result:
        return []

    id_list = search_result['esearchresult'].get('idlist', [])
    if not id_list:
        return []

    # Get summaries
    summary_result = client.summary('nucleotide', id_list[:max_results])

    records = []
    if 'result' in summary_result:
        for uid in id_list[:max_results]:
            if uid in summary_result['result']:
                record = summary_result['result'][uid]
                records.append({
                    'id': uid,
                    'accession': record.get('accessionversion', ''),
                    'title': record.get('title', ''),
                    'organism': record.get('organism', ''),
                    'length': int(record.get('slen', 0)),
                    'moltype': record.get('moltype', ''),
                    'created': record.get('createdate', ''),
                    'updated': record.get('updatedate', '')
                })

    return records


def fetch_genbank_record(accession: str, api_key: Optional[str] = None,
                        email: Optional[str] = None) -> Optional[str]:
    """Fetch complete GenBank record.

    Args:
        accession: Accession number
        api_key: NCBI API key
        email: Email address

    Returns:
        GenBank record as string
    """
    client = EntrezClient(api_key=api_key, email=email)
    return client.fetch('nucleotide', [accession], rettype='gb', retmode='text')


def fetch_fasta_sequence(accession: str, api_key: Optional[str] = None,
                        email: Optional[str] = None) -> Optional[str]:
    """Fetch sequence in FASTA format.

    Args:
        accession: Accession number
        api_key: NCBI API key
        email: Email address

    Returns:
        FASTA sequence as string
    """
    client = EntrezClient(api_key=api_key, email=email)
    fasta_data = client.fetch('nucleotide', [accession], rettype='fasta', retmode='text')

    if not fasta_data:
        return None

    # Parse FASTA to extract just the sequence
    lines = fasta_data.strip().split('\n')
    if len(lines) < 2:
        return None

    # Skip header and join sequence lines
    sequence = ''.join(lines[1:]).replace('\n', '')
    return sequence


def get_sequence_features(accession: str, api_key: Optional[str] = None,
                         email: Optional[str] = None) -> List[Dict[str, Any]]:
    """Get sequence features from GenBank record.

    Args:
        accession: Accession number
        api_key: NCBI API key
        email: Email address

    Returns:
        List of sequence features
    """
    record = fetch_genbank_record(accession, api_key, email)
    if not record:
        return []

    # Basic parsing of GenBank features (simplified)
    features = []

    # Look for FEATURES section
    lines = record.split('\n')
    in_features = False

    for line in lines:
        if line.startswith('FEATURES'):
            in_features = True
            continue
        elif in_features and line.startswith('ORIGIN'):
            break

        if in_features and line.strip():
            # Parse feature lines (simplified)
            if not line.startswith(' '):
                # Feature type line
                parts = line.strip().split()
                if len(parts) >= 2:
                    feature_type = parts[0]
                    location = parts[1]
                    features.append({
                        'type': feature_type,
                        'location': location,
                        'qualifiers': {}
                    })

    return features


def search_protein_records(query: str, max_results: int = 100,
                          api_key: Optional[str] = None,
                          email: Optional[str] = None) -> List[Dict[str, Any]]:
    """Search NCBI protein database.

    Args:
        query: Search query
        max_results: Maximum results
        api_key: NCBI API key
        email: Email address

    Returns:
        List of protein records
    """
    client = EntrezClient(api_key=api_key, email=email)

    search_result = client.search('protein', query, max_results)

    if 'esearchresult' not in search_result:
        return []

    id_list = search_result['esearchresult'].get('idlist', [])
    if not id_list:
        return []

    # Get summaries
    summary_result = client.summary('protein', id_list[:max_results])

    records = []
    if 'result' in summary_result:
        for uid in id_list[:max_results]:
            if uid in summary_result['result']:
                record = summary_result['result'][uid]
                records.append({
                    'id': uid,
                    'accession': record.get('accessionversion', ''),
                    'name': record.get('name', ''),
                    'organism': record.get('organism', ''),
                    'length': int(record.get('slen', 0)),
                    'created': record.get('createdate', ''),
                    'updated': record.get('updatedate', '')
                })

    return records


def link_nucleotide_to_protein(accession: str, api_key: Optional[str] = None,
                              email: Optional[str] = None) -> List[str]:
    """Find protein records linked to a nucleotide sequence.

    Args:
        accession: Nucleotide accession
        api_key: NCBI API key
        email: Email address

    Returns:
        List of linked protein accessions
    """
    client = EntrezClient(api_key=api_key, email=email)

    # Search for the nucleotide record first
    search_result = client.search('nucleotide', accession, 1)

    if 'esearchresult' not in search_result:
        return []

    id_list = search_result['esearchresult'].get('idlist', [])
    if not id_list:
        return []

    # Get links to protein database
    link_result = client.link('nucleotide', id_list, 'nucleotide_protein')

    protein_ids = []
    if 'linksets' in link_result:
        for linkset in link_result['linksets']:
            if 'linksetdbs' in linkset:
                for linksetdb in linkset['linksetdbs']:
                    if linksetdb.get('linkname') == 'nucleotide_protein':
                        protein_ids.extend(linksetdb.get('links', []))

    return protein_ids



