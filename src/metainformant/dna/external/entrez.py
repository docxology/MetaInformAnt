"""Entrez database utilities for NCBI data retrieval.

This module provides utilities for accessing NCBI databases through the
Entrez programming interface, including GenBank record retrieval,
sequence downloads, and metadata queries.

A request that fails (transport error, HTTP error status, response that is not
the expected JSON) raises :class:`NCBIRequestError`, a
``requests.RequestException``.  The clients used to log the failure and return
``{}``, ``""`` or ``[]``, which a caller could not tell from "no hits".
"""

from __future__ import annotations

import time
from io import StringIO
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import requests

from metainformant.core.ncbi import resolve_ncbi_contact
from metainformant.core.utils import logging

from ._http import NCBIRequestError, get_json_object, get_text, resolve_api_key

logger = logging.get_logger(__name__)

__all__ = [
    "EntrezClient",
    "NCBIRequestError",
    "fetch_fasta_sequence",
    "fetch_genbank_record",
    "get_genome_from_ncbi",
    "get_sequence_features",
    "link_nucleotide_to_protein",
    "parse_fasta_sequence",
    "search_genbank",
    "search_protein_records",
    "summarize_nucleotide_doc",
    "summarize_protein_doc",
    "summarize_taxonomy_doc",
]


_MAX_SEARCH_RESULTS = 10000  # most records one esearch call is asked for


def _require_positive_int(value: object, name: str, maximum: Optional[int] = None) -> int:
    """Return ``value`` if it is an integer from 1 to ``maximum``, else raise ``ValueError``.

    A result limit above what a client retrieves used to be cut to that maximum without a word.
    """
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer, got {value!r}")
    if maximum is not None and value > maximum:
        raise ValueError(f"{name} must be at most {maximum}, the most this client retrieves in one search; got {value}")
    return value


def _optional_int(value: object, field: str) -> Optional[int]:
    """An integer field of an esummary document; absent is ``None``, never ``0``."""
    if value is None or value == "":
        return None
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    if isinstance(value, str) and value.strip().isdigit():
        return int(value)
    raise NCBIRequestError(f"NCBI returned an unexpected value for {field}: {value!r}")


def summarize_nucleotide_doc(uid: str, doc: Mapping[str, Any]) -> Dict[str, Any]:
    """One nucleotide esummary document as a record; fields NCBI omits are ``None``, not ``""`` or ``0``."""
    return {
        "id": uid,
        "accession": doc.get("accessionversion"),
        "title": doc.get("title"),
        "organism": doc.get("organism"),
        "length": _optional_int(doc.get("slen"), "slen"),
        "moltype": doc.get("moltype"),
        "created": doc.get("createdate"),
        "updated": doc.get("updatedate"),
    }


def summarize_protein_doc(uid: str, doc: Mapping[str, Any]) -> Dict[str, Any]:
    """One protein esummary document as a record; fields NCBI omits are ``None``, not ``""`` or ``0``."""
    return {
        "id": uid,
        "accession": doc.get("accessionversion"),
        "name": doc.get("name"),
        "organism": doc.get("organism"),
        "length": _optional_int(doc.get("slen"), "slen"),
        "created": doc.get("createdate"),
        "updated": doc.get("updatedate"),
    }


def summarize_taxonomy_doc(tax_id: int, doc: Mapping[str, Any]) -> Dict[str, Any]:
    """One taxonomy esummary document as a record; fields NCBI omits are ``None``, not ``""``."""
    return {
        "tax_id": tax_id,
        "scientific_name": doc.get("scientificname"),
        "common_name": doc.get("commonname"),
        "rank": doc.get("rank"),
        "division": doc.get("division"),
        "lineage": doc.get("lineage"),
    }


def parse_fasta_sequence(text: str) -> Optional[str]:
    """Return the sequence of a single-record FASTA body, or ``None`` if the record is empty.

    Raises:
        NCBIRequestError: If the body is not FASTA (an error message served with HTTP 200 would otherwise
            be joined into a "sequence") or holds more than one record.
    """
    lines = text.strip().splitlines()
    if not lines:
        return None
    if not lines[0].startswith(">"):
        raise NCBIRequestError(f"NCBI returned something other than FASTA: {lines[0][:80]!r}")
    if any(line.startswith(">") for line in lines[1:]):
        raise NCBIRequestError("NCBI returned several FASTA records where one was expected")
    sequence = "".join(line.strip() for line in lines[1:])
    return sequence or None


def _esearch_ids(result: Mapping[str, Any], operation: str = "Entrez esearch") -> List[str]:
    """The UIDs of an esearch JSON response.

    Raises:
        NCBIRequestError: If the response has no ``esearchresult`` or NCBI reports an ``ERROR`` in it.
    """
    block = result.get("esearchresult")
    if not isinstance(block, Mapping):
        raise NCBIRequestError(f"{operation} failed: response has no 'esearchresult'")
    if block.get("ERROR"):
        raise NCBIRequestError(f"{operation} failed: {block['ERROR']}")
    return [str(uid) for uid in block.get("idlist") or []]


def _esummary_documents(
    summary: Mapping[str, Any], ids: Sequence[str], operation: str = "Entrez esummary"
) -> List[Tuple[str, Mapping[str, Any]]]:
    """The ``(uid, document)`` pairs of an esummary JSON response, in the order of ``ids``.

    Raises:
        NCBIRequestError: If the response has no ``result`` object.
    """
    result = summary.get("result")
    if not isinstance(result, Mapping):
        detail = summary.get("esummaryresult")
        raise NCBIRequestError(f"{operation} failed: response has no 'result'" + (f" ({detail})" if detail else ""))
    return [(uid, result[uid]) for uid in ids if isinstance(result.get(uid), Mapping)]


class EntrezClient:
    """Client for NCBI Entrez utilities."""

    def __init__(
        self,
        api_key: Optional[str] = None,
        email: Optional[str] = None,
        *,
        allow_anonymous: bool = False,
        timeout: float = 30.0,
    ):
        """Initialize Entrez client.

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

    def search(self, db: str, query: str, max_results: int = 100) -> Dict[str, Any]:
        """Perform Entrez search.

        Args:
            db: Database to search (nucleotide, protein, etc.)
            query: Search query
            max_results: Maximum number of results, from 1 to 10000

        Returns:
            Search results dictionary

        Raises:
            ValueError: If ``max_results`` is not an integer from 1 to 10000.
            NCBIRequestError: If the request fails.
        """
        params: dict[str, str | int] = {
            "db": db,
            "term": query,
            "retmax": _require_positive_int(max_results, "max_results", maximum=_MAX_SEARCH_RESULTS),
            "retmode": "json",
            "usehistory": "y",
        }

        if self.api_key:
            params["api_key"] = self.api_key
        if self.email:
            params["email"] = self.email

        return get_json_object(
            self.session,
            f"{self.base_url}/esearch.fcgi",
            operation="Entrez esearch",
            params=params,
            timeout=self.timeout,
        )

    def fetch(self, db: str, ids: List[str], rettype: str = "fasta", retmode: str = "text") -> str:
        """Fetch records from Entrez database.

        Args:
            db: Database name
            ids: List of IDs to fetch
            rettype: Retrieval type (fasta, gb, etc.)
            retmode: Retrieval mode (text, xml, etc.)

        Returns:
            Fetched data as string

        Raises:
            NCBIRequestError: If the request fails.
        """
        params = {"db": db, "id": ",".join(ids), "rettype": rettype, "retmode": retmode}

        if self.api_key:
            params["api_key"] = self.api_key
        if self.email:
            params["email"] = self.email

        # Respect NCBI rate limits (3 requests/second without API key)
        if not self.api_key:
            time.sleep(0.4)

        return get_text(
            self.session,
            f"{self.base_url}/efetch.fcgi",
            operation="Entrez efetch",
            params=params,
            timeout=self.timeout,
        )

    def summary(self, db: str, ids: List[str]) -> Dict[str, Any]:
        """Get summary information for records.

        Args:
            db: Database name
            ids: List of IDs

        Returns:
            Summary data dictionary

        Raises:
            NCBIRequestError: If the request fails.
        """
        params: dict[str, str] = {"db": db, "id": ",".join(ids), "retmode": "json"}

        if self.api_key:
            params["api_key"] = self.api_key
        if self.email:
            params["email"] = self.email

        return get_json_object(
            self.session,
            f"{self.base_url}/esummary.fcgi",
            operation="Entrez esummary",
            params=params,
            timeout=self.timeout,
        )

    def link(self, db: str, ids: List[str], linkname: str) -> Dict[str, Any]:
        """Find linked records in other databases.

        Args:
            db: Source database
            ids: Source IDs
            linkname: Link name (e.g., nucleotide_protein)

        Returns:
            Link data dictionary

        Raises:
            NCBIRequestError: If the request fails.
        """
        params: dict[str, str] = {"db": db, "id": ",".join(ids), "linkname": linkname, "retmode": "json"}

        if self.api_key:
            params["api_key"] = self.api_key
        if self.email:
            params["email"] = self.email

        return get_json_object(
            self.session,
            f"{self.base_url}/elink.fcgi",
            operation="Entrez elink",
            params=params,
            timeout=self.timeout,
        )


def search_genbank(
    query: str, max_results: int = 100, api_key: Optional[str] = None, email: Optional[str] = None
) -> List[Dict[str, Any]]:
    """Search GenBank nucleotide database.

    Args:
        query: Search query
        max_results: Maximum results, from 1 to 10000
        api_key: NCBI API key
        email: Email address

    Returns:
        List of GenBank records (an empty list means the search matched nothing)

    Raises:
        ValueError: If ``max_results`` is not an integer from 1 to 10000.
        NCBIRequestError: If a request fails.
    """
    client = EntrezClient(api_key=api_key, email=email)

    id_list = _esearch_ids(client.search("nucleotide", query, max_results))
    if not id_list:
        return []

    # Get summaries
    summary_result = client.summary("nucleotide", id_list[:max_results])
    return [
        summarize_nucleotide_doc(uid, doc) for uid, doc in _esummary_documents(summary_result, id_list[:max_results])
    ]


def fetch_genbank_record(accession: str, api_key: Optional[str] = None, email: Optional[str] = None) -> Optional[str]:
    """Fetch complete GenBank record.

    Args:
        accession: Accession number
        api_key: NCBI API key
        email: Email address

    Returns:
        GenBank record as string

    Raises:
        NCBIRequestError: If the request fails.
    """
    client = EntrezClient(api_key=api_key, email=email)
    return client.fetch("nucleotide", [accession], rettype="gb", retmode="text")


def fetch_fasta_sequence(accession: str, api_key: Optional[str] = None, email: Optional[str] = None) -> Optional[str]:
    """Fetch sequence in FASTA format.

    Args:
        accession: Accession number
        api_key: NCBI API key
        email: Email address

    Returns:
        The sequence, or ``None`` if NCBI returned an empty record

    Raises:
        NCBIRequestError: If the request fails or the response is not a single FASTA record.
    """
    client = EntrezClient(api_key=api_key, email=email)
    return parse_fasta_sequence(client.fetch("nucleotide", [accession], rettype="fasta", retmode="text"))


def get_genome_from_ncbi(accession: str, email: Optional[str] = None, api_key: Optional[str] = None) -> Any:
    """Fetch and parse one nucleotide record from NCBI as a ``SeqRecord``.

    This compatibility-level convenience function keeps the public DNA API
    aligned with the integration workflow while delegating transport to the
    canonical :class:`EntrezClient` implementation.

    Args:
        accession: NCBI nucleotide accession, optionally with a version.
        email: Contact email sent to NCBI.
        api_key: Optional NCBI API key.

    Returns:
        A Biopython ``SeqRecord`` parsed from the FASTA response.

    Raises:
        ValueError: If ``accession`` is empty.
        RuntimeError: If NCBI returns no sequence data.
        NCBIRequestError: If the request itself fails.
        ImportError: If Biopython is not installed.
    """
    if not accession or not accession.strip():
        raise ValueError("NCBI accession must not be empty")

    fasta_data = EntrezClient(api_key=api_key, email=email).fetch(
        "nucleotide", [accession.strip()], rettype="fasta", retmode="text"
    )
    if not fasta_data.strip():
        raise RuntimeError(f"NCBI returned no FASTA sequence for {accession}")

    try:
        from Bio import SeqIO
    except ImportError as exc:
        raise ImportError("Biopython is required to parse NCBI genome records") from exc

    try:
        return SeqIO.read(StringIO(fasta_data), "fasta")
    except Exception as exc:
        raise RuntimeError(f"Failed to parse NCBI FASTA record for {accession}: {exc}") from exc


def get_sequence_features(
    accession: str, api_key: Optional[str] = None, email: Optional[str] = None
) -> List[Dict[str, Any]]:
    """Get sequence features from GenBank record.

    Args:
        accession: Accession number
        api_key: NCBI API key
        email: Email address

    Returns:
        List of sequence features

    Raises:
        NCBIRequestError: If the request fails.
    """
    record = fetch_genbank_record(accession, api_key, email)
    if not record:
        return []

    return _parse_feature_table(record)


def _parse_feature_table(text: str) -> List[Dict[str, Any]]:
    """Parse the FEATURES section of a GenBank flat-file record.

    Feature rows are indented by exactly five spaces (feature key followed
    by its location); qualifier rows are indented deeper and start with
    ``/``. Continuation lines of a wrapped qualifier value are joined onto
    that value.

    Args:
        text: Full GenBank record text

    Returns:
        List of ``{"type", "location", "qualifiers"}`` dicts, where
        ``qualifiers`` maps qualifier names to their unquoted string values.
    """
    features: List[Dict[str, Any]] = []
    current: Optional[Dict[str, Any]] = None
    in_features = False

    for line in text.split("\n"):
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if in_features and (line.startswith("ORIGIN") or line.startswith("//")):
            break
        if not in_features or not line.strip():
            continue

        stripped = line.strip()
        if line.startswith("     ") and not line.startswith("      "):
            # Feature row: feature key (columns 6+) followed by its location
            parts = stripped.split(None, 1)
            if len(parts) == 2:
                current = {"type": parts[0], "location": parts[1], "qualifiers": {}}
                features.append(current)
        elif current is not None:
            if stripped.startswith("/"):
                name, _, value = stripped[1:].partition("=")
                current["qualifiers"][name.strip()] = value.strip()
            elif current["qualifiers"]:
                # Continuation of a wrapped qualifier value
                last_name = next(reversed(current["qualifiers"]))
                current["qualifiers"][last_name] += " " + stripped

    # Unquote qualifier values once, after any wrapped continuation lines
    for feature in features:
        for name, value in feature["qualifiers"].items():
            if value.startswith('"'):
                closing = value.rfind('"')
                feature["qualifiers"][name] = value[1:closing] if closing > 0 else value[1:]

    return features


def search_protein_records(
    query: str, max_results: int = 100, api_key: Optional[str] = None, email: Optional[str] = None
) -> List[Dict[str, Any]]:
    """Search NCBI protein database.

    Args:
        query: Search query
        max_results: Maximum results, from 1 to 10000
        api_key: NCBI API key
        email: Email address

    Returns:
        List of protein records (an empty list means the search matched nothing)

    Raises:
        ValueError: If ``max_results`` is not an integer from 1 to 10000.
        NCBIRequestError: If a request fails.
    """
    client = EntrezClient(api_key=api_key, email=email)

    id_list = _esearch_ids(client.search("protein", query, max_results))
    if not id_list:
        return []

    # Get summaries
    summary_result = client.summary("protein", id_list[:max_results])
    return [summarize_protein_doc(uid, doc) for uid, doc in _esummary_documents(summary_result, id_list[:max_results])]


def link_nucleotide_to_protein(accession: str, api_key: Optional[str] = None, email: Optional[str] = None) -> List[str]:
    """Find protein records linked to a nucleotide sequence.

    Args:
        accession: Nucleotide accession
        api_key: NCBI API key
        email: Email address

    Returns:
        List of linked protein record UIDs (numeric NCBI identifiers;
        elink returns numeric UIDs, not accessions)

    Raises:
        NCBIRequestError: If a request fails.
    """
    client = EntrezClient(api_key=api_key, email=email)

    # Search for the nucleotide record first
    id_list = _esearch_ids(client.search("nucleotide", accession, 1))
    if not id_list:
        return []

    # Get links to protein database
    link_result = client.link("nucleotide", id_list, "nucleotide_protein")

    protein_ids = []
    if "linksets" in link_result:
        for linkset in link_result["linksets"]:
            if "linksetdbs" in linkset:
                for linksetdb in linkset["linksetdbs"]:
                    if linksetdb.get("linkname") == "nucleotide_protein":
                        protein_ids.extend(linksetdb.get("links", []))

    return protein_ids
