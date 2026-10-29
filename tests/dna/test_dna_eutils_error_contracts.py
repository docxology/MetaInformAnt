"""Offline tests for the E-utilities clients (``dna.external.entrez`` and ``dna.external.ncbi``).

Two defects are covered.  A failed call used to be logged and answered with ``{}``, ``""``, ``[]`` or
``None``, so an outage looked like "no hits"; and a field missing from a record was filled with ``""`` or
``0``, so a record of unknown length looked like an empty sequence.

Failures are produced for real, without a stub: a client whose base URL has no scheme is refused by
``requests`` before any I/O, and a client pointed at the discard port of the local host gets a refused
connection.  The record shapes in the pure-function tests follow NCBI's E-utilities documentation; they
were not fetched for this change, because only Datasets endpoints were in scope for live checks.
"""

from __future__ import annotations

from typing import Any, Callable, Iterator

import pytest
import requests

from metainformant.dna.external import entrez, ncbi

CLIENT_KEY = "unit-test-key-0123"
CLIENT_EMAIL = "tests@example.org"


@pytest.fixture
def entrez_client() -> Iterator[entrez.EntrezClient]:
    client = entrez.EntrezClient(api_key=CLIENT_KEY, email=CLIENT_EMAIL, timeout=2)
    client.base_url = "not-a-url"
    yield client
    client.session.close()


@pytest.fixture
def general_client() -> Iterator[ncbi.NCBIClient]:
    client = ncbi.NCBIClient(api_key=CLIENT_KEY, email=CLIENT_EMAIL, timeout=2)
    client.base_url = "not-a-url"
    yield client
    client.session.close()


# --------------------------------------------------------------------------------------
# A failed call raises
# --------------------------------------------------------------------------------------


def test_primary_client_failures_raise_instead_of_returning_empty(entrez_client: entrez.EntrezClient) -> None:
    calls: list[Callable[[], Any]] = [
        lambda: entrez_client.search("nucleotide", "phiX174"),
        lambda: entrez_client.fetch("nucleotide", ["NC_001422.1"]),
        lambda: entrez_client.summary("nucleotide", ["9626372"]),
        lambda: entrez_client.link("nucleotide", ["9626372"], "nucleotide_protein"),
    ]
    for call in calls:
        with pytest.raises(requests.RequestException):
            call()


def test_general_client_failures_raise_instead_of_returning_empty(general_client: ncbi.NCBIClient) -> None:
    calls: list[Callable[[], Any]] = [
        lambda: general_client.search_nucleotide("phiX174", max_results=3),
        lambda: general_client.fetch_sequence("NC_001422.1"),
        lambda: general_client.search_protein("phiX174", max_results=3),
        lambda: general_client.get_taxonomy_info(10847),
    ]
    for call in calls:
        with pytest.raises(requests.RequestException):
            call()


@pytest.fixture
def dead_proxy(monkeypatch: pytest.MonkeyPatch) -> None:
    """Send every outbound connection of this process to a closed local port: a real failure, no stub."""
    for name in ("http_proxy", "https_proxy", "all_proxy", "HTTP_PROXY", "HTTPS_PROXY", "ALL_PROXY"):
        monkeypatch.setenv(name, "http://127.0.0.1:9")
    for name in ("no_proxy", "NO_PROXY"):
        monkeypatch.setenv(name, "")
    monkeypatch.delenv("NCBI_API_KEY", raising=False)


def test_module_functions_raise_when_the_service_cannot_be_reached(dead_proxy: None) -> None:
    """Every convenience function built its own client and returned ``[]``, ``None`` or ``""`` on failure."""
    calls: list[Callable[[], Any]] = [
        lambda: entrez.search_genbank("phiX174", 3, api_key=CLIENT_KEY, email=CLIENT_EMAIL),
        lambda: entrez.search_protein_records("phiX174", 3, api_key=CLIENT_KEY, email=CLIENT_EMAIL),
        lambda: entrez.fetch_fasta_sequence("NC_001422.1", api_key=CLIENT_KEY, email=CLIENT_EMAIL),
        lambda: entrez.fetch_genbank_record("NC_001422.1", api_key=CLIENT_KEY, email=CLIENT_EMAIL),
        lambda: entrez.get_sequence_features("NC_001422.1", api_key=CLIENT_KEY, email=CLIENT_EMAIL),
        lambda: entrez.link_nucleotide_to_protein("NC_001422.1", api_key=CLIENT_KEY, email=CLIENT_EMAIL),
        lambda: ncbi.search_nucleotide("phiX174", 3, api_key=CLIENT_KEY, email=CLIENT_EMAIL),
        lambda: ncbi.search_protein("phiX174", 3, api_key=CLIENT_KEY, email=CLIENT_EMAIL),
        lambda: ncbi.fetch_sequence("NC_001422.1", api_key=CLIENT_KEY, email=CLIENT_EMAIL),
        lambda: ncbi.get_taxonomy_info(10847, api_key=CLIENT_KEY, email=CLIENT_EMAIL),
    ]
    for call in calls:
        with pytest.raises(requests.RequestException):
            call()


def test_refused_connection_is_reported_without_the_key_or_contact(entrez_client: entrez.EntrezClient) -> None:
    """A real transport failure: nothing listens on the discard port.  The text must not carry credentials."""
    entrez_client.base_url = "http://127.0.0.1:9"

    with pytest.raises(requests.RequestException) as caught:
        entrez_client.search("nucleotide", "phiX174")

    message = str(caught.value)
    assert CLIENT_KEY not in message and CLIENT_EMAIL not in message
    assert "Entrez esearch failed" in message
    assert caught.value.__cause__ is None and caught.value.__suppress_context__
    assert getattr(caught.value, "status_code", "missing") is None


# --------------------------------------------------------------------------------------
# Configuration the documentation promises
# --------------------------------------------------------------------------------------


def test_api_key_falls_back_to_the_environment_variable(monkeypatch: pytest.MonkeyPatch) -> None:
    """docs/dna/ncbi.md tells users to set ``NCBI_API_KEY``; the clients used to ignore it."""
    monkeypatch.setenv("NCBI_API_KEY", " env-key ")
    clients = [
        entrez.EntrezClient(email=CLIENT_EMAIL),
        ncbi.NCBIClient(email=CLIENT_EMAIL),
        entrez.EntrezClient(api_key="explicit-key", email=CLIENT_EMAIL),
    ]
    try:
        assert [client.api_key for client in clients] == ["env-key", "env-key", "explicit-key"]
    finally:
        for client in clients:
            client.session.close()

    monkeypatch.setenv("NCBI_API_KEY", "   ")
    blank = entrez.EntrezClient(email=CLIENT_EMAIL)
    monkeypatch.delenv("NCBI_API_KEY")
    unset = ncbi.NCBIClient(email=CLIENT_EMAIL)
    try:
        assert blank.api_key is None and unset.api_key is None
    finally:
        blank.session.close()
        unset.session.close()


# --------------------------------------------------------------------------------------
# Arguments
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("max_results", [0, -1, True, 2.5, "10", None])
def test_result_limits_must_be_positive_integers(
    max_results: Any, entrez_client: entrez.EntrezClient, general_client: ncbi.NCBIClient
) -> None:
    """Before this change a limit of 0 or -1 became ``retmax=0`` or a negative slice of the hits."""
    for call in (
        lambda: entrez_client.search("nucleotide", "x", max_results),
        lambda: general_client.search_nucleotide("x", max_results),
        lambda: general_client.search_protein("x", max_results),
    ):
        with pytest.raises(ValueError, match="max_results"):
            call()


def test_result_limits_above_the_client_maximum_are_refused_not_truncated(
    entrez_client: entrez.EntrezClient, general_client: ncbi.NCBIClient
) -> None:
    """Asking for 5000 records used to return at most 1000 (or 10000) without saying so."""
    with pytest.raises(ValueError, match="at most 10000"):
        entrez_client.search("nucleotide", "x", 10001)
    for call in (
        lambda: general_client.search_nucleotide("x", 1001),
        lambda: general_client.search_protein("x", 1001),
    ):
        with pytest.raises(ValueError, match="at most 1000"):
            call()


def test_result_limits_at_the_client_maximum_are_sent(
    entrez_client: entrez.EntrezClient, general_client: ncbi.NCBIClient
) -> None:
    """The boundary is accepted: the call proceeds to the (unreachable) service and fails there."""
    with pytest.raises(requests.RequestException):
        entrez_client.search("nucleotide", "x", 10000)
    with pytest.raises(requests.RequestException):
        general_client.search_nucleotide("x", 1000)


# --------------------------------------------------------------------------------------
# Records: absent stays absent
# --------------------------------------------------------------------------------------


def test_nucleotide_record_leaves_missing_fields_none() -> None:
    record = entrez.summarize_nucleotide_doc("9626372", {"accessionversion": "NC_001422.1", "slen": 5386})

    assert record == {
        "id": "9626372",
        "accession": "NC_001422.1",
        "title": None,
        "organism": None,
        "length": 5386,
        "moltype": None,
        "created": None,
        "updated": None,
    }


def test_missing_length_is_none_not_zero() -> None:
    assert entrez.summarize_nucleotide_doc("1", {})["length"] is None
    assert entrez.summarize_protein_doc("1", {"slen": ""})["length"] is None
    assert entrez.summarize_protein_doc("1", {"slen": "312"})["length"] == 312
    assert entrez.summarize_protein_doc("1", {"slen": 0})["length"] == 0  # a reported zero is kept


@pytest.mark.parametrize("length", ["n/a", 3.5, True, [], {}])
def test_unparseable_length_raises(length: Any) -> None:
    with pytest.raises(requests.RequestException, match="slen"):
        entrez.summarize_nucleotide_doc("1", {"slen": length})


def test_protein_and_taxonomy_records_leave_missing_fields_none() -> None:
    protein = entrez.summarize_protein_doc("5", {"name": "capsid", "organism": "Escherichia phage phiX174"})
    taxon = entrez.summarize_taxonomy_doc(7460, {"scientificname": "Apis mellifera"})

    assert protein["accession"] is None and protein["length"] is None and protein["name"] == "capsid"
    assert taxon == {
        "tax_id": 7460,
        "scientific_name": "Apis mellifera",
        "common_name": None,
        "rank": None,
        "division": None,
        "lineage": None,
    }


# --------------------------------------------------------------------------------------
# Responses that are not what was asked for
# --------------------------------------------------------------------------------------


def test_esearch_ids_reads_hits_and_zero_hits() -> None:
    assert entrez._esearch_ids({"esearchresult": {"count": "2", "idlist": ["11", "12"]}}) == ["11", "12"]
    assert entrez._esearch_ids({"esearchresult": {"count": "0", "idlist": []}}) == []


@pytest.mark.parametrize(
    "payload",
    [
        {},
        {"esearchresult": []},
        {"esearchresult": {"ERROR": "Invalid db name specified: nucleotid"}},
        {"error": "rate"},
    ],
)
def test_esearch_without_a_result_block_or_with_an_error_raises(payload: dict[str, Any]) -> None:
    with pytest.raises(requests.RequestException, match="Entrez esearch failed"):
        entrez._esearch_ids(payload)


def test_esummary_documents_keep_order_and_skip_unknown_uids() -> None:
    result = {"result": {"uids": ["2", "1"], "1": {"title": "one"}, "2": {"title": "two"}}}

    assert entrez._esummary_documents(result, ["1", "2", "3"]) == [("1", {"title": "one"}), ("2", {"title": "two"})]


def test_esummary_without_a_result_object_raises_with_the_service_text() -> None:
    with pytest.raises(requests.RequestException, match="Invalid uid 12345"):
        entrez._esummary_documents({"esummaryresult": ["Invalid uid 12345 at position=0"]}, ["12345"])
    with pytest.raises(requests.RequestException, match="no 'result'"):
        entrez._esummary_documents({}, ["1"])


def test_fasta_body_is_parsed_to_a_sequence() -> None:
    assert entrez.parse_fasta_sequence(">rec description\nACGT\nTTGA\n") == "ACGTTTGA"
    assert entrez.parse_fasta_sequence(">rec\r\nACGT\r\nTT\r\n") == "ACGTTT"
    assert entrez.parse_fasta_sequence("\n  >rec\nACGT\n\n") == "ACGT"


def test_empty_fasta_records_give_none() -> None:
    assert entrez.parse_fasta_sequence("") is None
    assert entrez.parse_fasta_sequence("\n\n") is None
    assert entrez.parse_fasta_sequence(">rec only a header\n") is None


@pytest.mark.parametrize(
    "body",
    [
        "Error: Failed to retrieve sequence\nPlease try again later\n",
        "<?xml version='1.0'?>\n<ERROR>Invalid uid</ERROR>\n",
        ">first\nACGT\n>second\nTTGA\n",
    ],
)
def test_bodies_that_are_not_one_fasta_record_raise(body: str) -> None:
    """An error text with two or more lines used to be joined into a 'sequence'."""
    with pytest.raises(requests.RequestException):
        entrez.parse_fasta_sequence(body)
