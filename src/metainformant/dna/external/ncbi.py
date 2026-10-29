"""NCBI API integration utilities.

This module provides tools for interacting with NCBI databases including
Entrez, BLAST, and other NCBI services for genomic data retrieval.

A request that fails raises :class:`NCBIRequestError` (a
``requests.RequestException``); it no longer returns an empty list or ``None``,
which a caller could not tell from "no hits".  Fields NCBI leaves out of a
record are ``None``, not an invented ``""`` or ``0``.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import requests

from metainformant.core.ncbi import resolve_ncbi_contact
from metainformant.core.utils import logging

from ._http import NCBIRequestError, get_json_object, get_text, resolve_api_key
from .entrez import (
    _esearch_ids,
    _esummary_documents,
    _require_positive_int,
    parse_fasta_sequence,
    summarize_nucleotide_doc,
    summarize_protein_doc,
    summarize_taxonomy_doc,
)

logger = logging.get_logger(__name__)

_MAX_SEARCH_RESULTS = 1000  # most records one search is asked for; esummary ids travel in a single GET URL

__all__ = [
    "NCBIClient",
    "NCBIRequestError",
    "fetch_sequence",
    "get_taxonomy_info",
    "search_nucleotide",
    "search_protein",
]


class NCBIClient:
    """Client for NCBI API interactions."""

    def __init__(
        self,
        api_key: Optional[str] = None,
        email: Optional[str] = None,
        *,
        allow_anonymous: bool = False,
        timeout: float = 30.0,
    ):
        """Initialize NCBI client.

        Args:
            api_key: NCBI API key for higher rate limits. If omitted, ``NCBI_API_KEY`` is used.
            email: Explicit contact email. If omitted, ``NCBI_EMAIL`` is used.
            allow_anonymous: Explicitly permit a request without an email.
            timeout: Per-request timeout in seconds.
        """
        if timeout <= 0:
            raise ValueError("timeout must be greater than zero")

        contact = resolve_ncbi_contact(email, allow_anonymous=allow_anonymous)
        self.api_key = resolve_api_key(api_key)
        self.email = contact.email
        self.contact_mode = contact.mode
        self.timeout = timeout
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.session = requests.Session()

        # Set headers
        headers = {"User-Agent": "metainformant/0.4.0"}
        if contact.email:
            headers["From"] = contact.email
        self.session.headers.update(headers)

    def _credentials(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Add the API key and contact email to ``params`` when configured."""
        if self.api_key:
            params["api_key"] = self.api_key
        if self.email:
            params["email"] = self.email
        return params

    def search_nucleotide(self, query: str, max_results: int = 100) -> List[Dict[str, Any]]:
        """Search NCBI Nucleotide database.

        Args:
            query: Search query
            max_results: Maximum number of results, from 1 to 1000

        Returns:
            List of search results (an empty list means the search matched nothing)

        Raises:
            ValueError: If ``max_results`` is not an integer from 1 to 1000.
            NCBIRequestError: If a request fails.
        """
        # First, search using esearch
        search_params: dict[str, str | int] = {
            "db": "nucleotide",
            "term": query,
            "retmax": _require_positive_int(max_results, "max_results", maximum=_MAX_SEARCH_RESULTS),
            "retmode": "json",
            "usehistory": "y",
        }
        search_data = get_json_object(
            self.session,
            f"{self.base_url}/esearch.fcgi",
            operation="NCBI esearch",
            params=self._credentials(search_params),
            timeout=self.timeout,
        )
        id_list = _esearch_ids(search_data, "NCBI esearch")[:max_results]
        if not id_list:
            return []

        # Fetch summaries using esummary
        summary_params: Dict[str, Any] = {"db": "nucleotide", "id": ",".join(id_list), "retmode": "json"}
        esearch_block = search_data["esearchresult"]
        webenv = esearch_block.get("webenv")
        query_key = esearch_block.get("querykey")
        if webenv and query_key:
            summary_params.update({"WebEnv": webenv, "query_key": query_key})

        summary_data = get_json_object(
            self.session,
            f"{self.base_url}/esummary.fcgi",
            operation="NCBI esummary",
            params=self._credentials(summary_params),
            timeout=self.timeout,
        )
        return [summarize_nucleotide_doc(uid, doc) for uid, doc in _esummary_documents(summary_data, id_list)]

    def fetch_sequence(self, accession: str) -> Optional[str]:
        """Fetch DNA sequence by accession.

        Args:
            accession: Sequence accession

        Returns:
            DNA sequence string, or None if NCBI returned an empty record

        Raises:
            NCBIRequestError: If the request fails or the response is not a single FASTA record.
        """
        # Use efetch to get sequence
        fetch_params: Dict[str, Any] = {"db": "nucleotide", "id": accession, "rettype": "fasta", "retmode": "text"}
        content = get_text(
            self.session,
            f"{self.base_url}/efetch.fcgi",
            operation="NCBI efetch",
            params=self._credentials(fetch_params),
            timeout=self.timeout,
        )
        return parse_fasta_sequence(content)

    def search_protein(self, query: str, max_results: int = 100) -> List[Dict[str, Any]]:
        """Search NCBI Protein database.

        Args:
            query: Search query
            max_results: Maximum number of results, from 1 to 1000

        Returns:
            List of protein search results (an empty list means the search matched nothing)

        Raises:
            ValueError: If ``max_results`` is not an integer from 1 to 1000.
            NCBIRequestError: If a request fails.
        """
        # Similar to nucleotide search but for protein database
        search_params: dict[str, str | int] = {
            "db": "protein",
            "term": query,
            "retmax": _require_positive_int(max_results, "max_results", maximum=_MAX_SEARCH_RESULTS),
            "retmode": "json",
            "usehistory": "y",
        }
        search_data = get_json_object(
            self.session,
            f"{self.base_url}/esearch.fcgi",
            operation="NCBI esearch",
            params=self._credentials(search_params),
            timeout=self.timeout,
        )
        id_list = _esearch_ids(search_data, "NCBI esearch")[:max_results]
        if not id_list:
            return []

        # Fetch summaries
        summary_params: Dict[str, Any] = {"db": "protein", "id": ",".join(id_list), "retmode": "json"}
        summary_data = get_json_object(
            self.session,
            f"{self.base_url}/esummary.fcgi",
            operation="NCBI esummary",
            params=self._credentials(summary_params),
            timeout=self.timeout,
        )
        return [summarize_protein_doc(uid, doc) for uid, doc in _esummary_documents(summary_data, id_list)]

    def get_taxonomy_info(self, tax_id: int) -> Optional[Dict[str, Any]]:
        """Get taxonomy information for a taxonomic ID.

        Args:
            tax_id: NCBI taxonomy ID

        Returns:
            Taxonomy information dictionary, or None if NCBI has no such taxonomy ID

        Raises:
            NCBIRequestError: If the request fails.
        """
        summary_params: Dict[str, Any] = {"db": "taxonomy", "id": str(tax_id), "retmode": "json"}
        data = get_json_object(
            self.session,
            f"{self.base_url}/esummary.fcgi",
            operation="NCBI esummary",
            params=self._credentials(summary_params),
            timeout=self.timeout,
        )

        result = data.get("result")
        if isinstance(result, dict) and isinstance(result.get(str(tax_id)), dict):
            return summarize_taxonomy_doc(tax_id, result[str(tax_id)])

        return None


def search_nucleotide(
    query: str,
    max_results: int = 100,
    api_key: Optional[str] = None,
    email: Optional[str] = None,
    *,
    allow_anonymous: bool = False,
    timeout: float = 30.0,
) -> List[Dict[str, Any]]:
    """Convenience function for nucleotide database search.

    Args:
        query: Search query
        max_results: Maximum results to return, from 1 to 1000
        api_key: NCBI API key
        email: Explicit contact email; otherwise ``NCBI_EMAIL``.
        allow_anonymous: Explicitly permit an anonymous request.
        timeout: Per-request timeout in seconds.

    Returns:
        List of search results

    Raises:
        ValueError: If ``max_results`` is not an integer from 1 to 1000.
        NCBIRequestError: If a request fails.
    """
    client = NCBIClient(api_key=api_key, email=email, allow_anonymous=allow_anonymous, timeout=timeout)
    return client.search_nucleotide(query, max_results)


def fetch_sequence(
    accession: str,
    api_key: Optional[str] = None,
    email: Optional[str] = None,
    *,
    allow_anonymous: bool = False,
    timeout: float = 30.0,
) -> Optional[str]:
    """Convenience function for sequence fetching.

    Args:
        accession: Sequence accession
        api_key: NCBI API key
        email: Explicit contact email; otherwise ``NCBI_EMAIL``.
        allow_anonymous: Explicitly permit an anonymous request.
        timeout: Per-request timeout in seconds.

    Returns:
        DNA sequence string

    Raises:
        NCBIRequestError: If the request fails or the response is not a single FASTA record.
    """
    client = NCBIClient(api_key=api_key, email=email, allow_anonymous=allow_anonymous, timeout=timeout)
    return client.fetch_sequence(accession)


def search_protein(
    query: str,
    max_results: int = 100,
    api_key: Optional[str] = None,
    email: Optional[str] = None,
    *,
    allow_anonymous: bool = False,
    timeout: float = 30.0,
) -> List[Dict[str, Any]]:
    """Convenience function for protein database search.

    Args:
        query: Search query
        max_results: Maximum results to return, from 1 to 1000
        api_key: NCBI API key
        email: Explicit contact email; otherwise ``NCBI_EMAIL``.
        allow_anonymous: Explicitly permit an anonymous request.
        timeout: Per-request timeout in seconds.

    Returns:
        List of protein search results

    Raises:
        ValueError: If ``max_results`` is not an integer from 1 to 1000.
        NCBIRequestError: If a request fails.
    """
    client = NCBIClient(api_key=api_key, email=email, allow_anonymous=allow_anonymous, timeout=timeout)
    return client.search_protein(query, max_results)


def get_taxonomy_info(
    tax_id: int,
    api_key: Optional[str] = None,
    email: Optional[str] = None,
    *,
    allow_anonymous: bool = False,
    timeout: float = 30.0,
) -> Optional[Dict[str, Any]]:
    """Convenience function for taxonomy information.

    Args:
        tax_id: NCBI taxonomy ID
        api_key: NCBI API key
        email: Explicit contact email; otherwise ``NCBI_EMAIL``.
        allow_anonymous: Explicitly permit an anonymous request.
        timeout: Per-request timeout in seconds.

    Returns:
        Taxonomy information, or None if NCBI has no such taxonomy ID

    Raises:
        NCBIRequestError: If the request fails.
    """
    client = NCBIClient(api_key=api_key, email=email, allow_anonymous=allow_anonymous, timeout=timeout)
    return client.get_taxonomy_info(tax_id)
