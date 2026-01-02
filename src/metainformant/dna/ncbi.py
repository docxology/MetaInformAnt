"""NCBI API integration utilities.

This module provides tools for interacting with NCBI databases including
Entrez, BLAST, and other NCBI services for genomic data retrieval.
"""

from __future__ import annotations

import time
from typing import Any, Dict, List, Optional

import requests

from metainformant.core import logging

logger = logging.get_logger(__name__)


class NCBIClient:
    """Client for NCBI API interactions."""

    def __init__(self, api_key: Optional[str] = None, email: Optional[str] = None):
        """Initialize NCBI client.

        Args:
            api_key: NCBI API key for higher rate limits
            email: Email address for NCBI requests
        """
        self.api_key = api_key
        self.email = email
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.session = requests.Session()

        # Set headers
        headers = {'User-Agent': 'metainformant/0.2.0'}
        if email:
            headers['From'] = email
        self.session.headers.update(headers)

    def search_nucleotide(self, query: str, max_results: int = 100) -> List[Dict[str, Any]]:
        """Search NCBI Nucleotide database.

        Args:
            query: Search query
            max_results: Maximum number of results

        Returns:
            List of search results
        """
        # First, search using esearch
        search_params = {
            'db': 'nucleotide',
            'term': query,
            'retmax': min(max_results, 1000),  # NCBI limit
            'retmode': 'json',
            'usehistory': 'y'
        }

        if self.api_key:
            search_params['api_key'] = self.api_key

        try:
            response = self.session.get(f"{self.base_url}/esearch.fcgi", params=search_params)
            response.raise_for_status()

            search_data = response.json()

            if 'esearchresult' not in search_data:
                return []

            id_list = search_data['esearchresult'].get('idlist', [])
            webenv = search_data['esearchresult'].get('webenv')
            query_key = search_data['esearchresult'].get('querykey')

            if not id_list:
                return []

            # Fetch summaries using esummary
            summary_params = {
                'db': 'nucleotide',
                'id': ','.join(id_list[:max_results]),
                'retmode': 'json'
            }

            if webenv and query_key:
                summary_params.update({'WebEnv': webenv, 'query_key': query_key})

            if self.api_key:
                summary_params['api_key'] = self.api_key

            summary_response = self.session.get(f"{self.base_url}/esummary.fcgi", params=summary_params)
            summary_response.raise_for_status()

            summary_data = summary_response.json()

            results = []
            if 'result' in summary_data:
                for uid in id_list[:max_results]:
                    if uid in summary_data['result']:
                        result = summary_data['result'][uid]
                        results.append({
                            'id': uid,
                            'accession': result.get('accessionversion', ''),
                            'title': result.get('title', ''),
                            'organism': result.get('organism', ''),
                            'length': result.get('slen', 0),
                            'moltype': result.get('moltype', ''),
                            'created': result.get('createdate', ''),
                            'updated': result.get('updatedate', '')
                        })

            return results

        except requests.RequestException as e:
            logger.error(f"NCBI search failed: {e}")
            return []

    def fetch_sequence(self, accession: str) -> Optional[str]:
        """Fetch DNA sequence by accession.

        Args:
            accession: Sequence accession

        Returns:
            DNA sequence string, or None if not found
        """
        try:
            # Use efetch to get sequence
            fetch_params = {
                'db': 'nucleotide',
                'id': accession,
                'rettype': 'fasta',
                'retmode': 'text'
            }

            if self.api_key:
                fetch_params['api_key'] = self.api_key

            response = self.session.get(f"{self.base_url}/efetch.fcgi", params=fetch_params)
            response.raise_for_status()

            content = response.text

            # Parse FASTA format
            lines = content.strip().split('\n')
            if len(lines) < 2:
                return None

            # Skip header line and join sequence
            sequence = ''.join(lines[1:]).replace('\n', '')

            return sequence

        except requests.RequestException as e:
            logger.error(f"Failed to fetch sequence {accession}: {e}")
            return None

    def search_protein(self, query: str, max_results: int = 100) -> List[Dict[str, Any]]:
        """Search NCBI Protein database.

        Args:
            query: Search query
            max_results: Maximum number of results

        Returns:
            List of protein search results
        """
        # Similar to nucleotide search but for protein database
        search_params = {
            'db': 'protein',
            'term': query,
            'retmax': min(max_results, 1000),
            'retmode': 'json',
            'usehistory': 'y'
        }

        if self.api_key:
            search_params['api_key'] = self.api_key

        try:
            response = self.session.get(f"{self.base_url}/esearch.fcgi", params=search_params)
            response.raise_for_status()

            search_data = response.json()

            if 'esearchresult' not in search_data:
                return []

            id_list = search_data['esearchresult'].get('idlist', [])

            if not id_list:
                return []

            # Fetch summaries
            summary_params = {
                'db': 'protein',
                'id': ','.join(id_list[:max_results]),
                'retmode': 'json'
            }

            if self.api_key:
                summary_params['api_key'] = self.api_key

            summary_response = self.session.get(f"{self.base_url}/esummary.fcgi", params=summary_params)
            summary_response.raise_for_status()

            summary_data = summary_response.json()

            results = []
            if 'result' in summary_data:
                for uid in id_list[:max_results]:
                    if uid in summary_data['result']:
                        result = summary_data['result'][uid]
                        results.append({
                            'id': uid,
                            'accession': result.get('accessionversion', ''),
                            'name': result.get('name', ''),
                            'organism': result.get('organism', ''),
                            'length': result.get('slen', 0),
                            'created': result.get('createdate', ''),
                            'updated': result.get('updatedate', '')
                        })

            return results

        except requests.RequestException as e:
            logger.error(f"NCBI protein search failed: {e}")
            return []

    def get_taxonomy_info(self, tax_id: int) -> Optional[Dict[str, Any]]:
        """Get taxonomy information for a taxonomic ID.

        Args:
            tax_id: NCBI taxonomy ID

        Returns:
            Taxonomy information dictionary
        """
        try:
            summary_params = {
                'db': 'taxonomy',
                'id': str(tax_id),
                'retmode': 'json'
            }

            if self.api_key:
                summary_params['api_key'] = self.api_key

            response = self.session.get(f"{self.base_url}/esummary.fcgi", params=summary_params)
            response.raise_for_status()

            data = response.json()

            if 'result' in data and str(tax_id) in data['result']:
                result = data['result'][str(tax_id)]
                return {
                    'tax_id': tax_id,
                    'scientific_name': result.get('scientificname', ''),
                    'common_name': result.get('commonname', ''),
                    'rank': result.get('rank', ''),
                    'division': result.get('division', ''),
                    'lineage': result.get('lineage', ''),
                }

            return None

        except requests.RequestException as e:
            logger.error(f"Failed to get taxonomy info for {tax_id}: {e}")
            return None


def search_nucleotide(query: str, max_results: int = 100, api_key: Optional[str] = None) -> List[Dict[str, Any]]:
    """Convenience function for nucleotide database search.

    Args:
        query: Search query
        max_results: Maximum results to return
        api_key: NCBI API key

    Returns:
        List of search results
    """
    client = NCBIClient(api_key=api_key)
    return client.search_nucleotide(query, max_results)


def fetch_sequence(accession: str, api_key: Optional[str] = None) -> Optional[str]:
    """Convenience function for sequence fetching.

    Args:
        accession: Sequence accession
        api_key: NCBI API key

    Returns:
        DNA sequence string
    """
    client = NCBIClient(api_key=api_key)
    return client.fetch_sequence(accession)


def search_protein(query: str, max_results: int = 100, api_key: Optional[str] = None) -> List[Dict[str, Any]]:
    """Convenience function for protein database search.

    Args:
        query: Search query
        max_results: Maximum results to return
        api_key: NCBI API key

    Returns:
        List of protein search results
    """
    client = NCBIClient(api_key=api_key)
    return client.search_protein(query, max_results)


def get_taxonomy_info(tax_id: int, api_key: Optional[str] = None) -> Optional[Dict[str, Any]]:
    """Convenience function for taxonomy information.

    Args:
        tax_id: NCBI taxonomy ID
        api_key: NCBI API key

    Returns:
        Taxonomy information
    """
    client = NCBIClient(api_key=api_key)
    return client.get_taxonomy_info(tax_id)



