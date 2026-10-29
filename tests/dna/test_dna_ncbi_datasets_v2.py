"""Offline tests for the NCBI Datasets v2 port of the genome functions.

Everything here runs without an external network. Paging uses real local HTTP; parsers are exercised
on responses saved from the live service
(``tests/data/ncbi_datasets/``, provenance in its README); a few edge cases are derived from a saved report
by deleting or editing one field, and say so.  Query builders are checked against the exact URLs that were
sent to the service to capture the fixtures.  Archive tests build real zip files with ``zipfile``.
Failure behaviour uses an unreachable base URL, which ``requests`` rejects before any I/O.
"""

from __future__ import annotations

import copy
import json
import stat
import threading
import zipfile
from contextlib import contextmanager
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Iterator
from urllib.parse import parse_qs, urlsplit

import pytest
import requests

from metainformant.core.ncbi import NCBIContactError
from metainformant.dna.external import genomes, ncbi_datasets
from metainformant.dna.external._http import (
    NCBIRequestError,
    describe_http_error,
    describe_transport_error,
    redact_secrets,
)
from metainformant.dna.external.ncbi_datasets import (
    AssemblyNotFoundError,
    NCBIDatasetsClient,
    NCBIDatasetsError,
    ReportPage,
)

FIXTURES = Path(__file__).resolve().parent.parent / "data" / "ncbi_datasets"
AMEL = "GCF_003254395.2"
BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"


def load(name: str) -> Any:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


def amel_report() -> dict[str, Any]:
    report: dict[str, Any] = load(f"dataset_report_accession_{AMEL}.json")["reports"][0]
    return report


def phix_report() -> dict[str, Any]:
    report: dict[str, Any] = load("dataset_report_accession_GCF_000819615.1.json")["reports"][0]
    return report


@pytest.fixture
def unreachable_client() -> Iterator[NCBIDatasetsClient]:
    """A real client whose base URL has no scheme: ``requests`` refuses it before any network I/O."""
    client = NCBIDatasetsClient(email="tests@example.org")
    client.base_url = "not-a-url"
    yield client
    client.close()


# --------------------------------------------------------------------------------------
# Assembly reports
# --------------------------------------------------------------------------------------


def test_honey_bee_report_parses_to_metadata() -> None:
    metadata = ncbi_datasets.parse_assembly_report(amel_report())

    assert metadata == {
        "accession": AMEL,
        "organism_name": "Apis mellifera",
        "tax_id": 7460,
        "assembly_name": "Amel_HAv3.1",
        "assembly_level": "Chromosome",
        "assembly_status": "current",
        "assembly_method": "FALCON v. 0.5.0; Arcs v. 1.0.1; Links v. 1.8.5; BioNano Solve v. 3.1",
        "refseq_category": "reference genome",
        "source_database": "SOURCE_DATABASE_REFSEQ",
        "paired_accession": "GCA_003254395.2",
        "genome_size": 225234541,
        "gc_percent": 32.5,
        "chromosome_count": 16,
        "contig_count": 227,
        "contig_n50": 5382476,
        "scaffold_count": 176,
        "scaffold_n50": 13619445,
        "release_date": "2018-09-10",
        "submitter": "Uppsala University",
    }
    # NCBI serializes uint64 fields as strings; the metadata carries a number.
    assert isinstance(amel_report()["assembly_stats"]["total_sequence_length"], str)
    assert isinstance(metadata["genome_size"], int)


def test_old_contig_and_scaffold_field_names_do_not_exist_in_reports() -> None:
    """``contig_count`` / ``scaffold_count`` were read from fields the schema does not have."""
    stats = amel_report()["assembly_stats"]
    assert "contig_count" not in stats and "scaffold_count" not in stats
    assert stats["number_of_contigs"] == 227 and stats["number_of_scaffolds"] == 176


def test_sparse_report_keeps_absent_fields_none() -> None:
    """Phage phiX174's report has no submitter, assembly method or RefSeq category."""
    metadata = ncbi_datasets.parse_assembly_report(phix_report())

    assert metadata["organism_name"] == "Escherichia phage phiX174"
    assert metadata["genome_size"] == 5386
    assert metadata["gc_percent"] == 45.0
    assert metadata["submitter"] is None
    assert metadata["assembly_method"] is None
    assert metadata["refseq_category"] is None


def test_missing_organism_object_gives_none_not_unknown() -> None:
    """Derived from the honey bee report by deleting ``organism``."""
    report = amel_report()
    del report["organism"]

    metadata = ncbi_datasets.parse_assembly_report(report)

    assert metadata["organism_name"] is None
    assert metadata["tax_id"] is None


def test_report_without_accession_is_rejected() -> None:
    with pytest.raises(NCBIDatasetsError, match="without an accession"):
        ncbi_datasets.parse_assembly_report({})


def test_malformed_numeric_fields_are_rejected() -> None:
    """Derived from the honey bee report: a size that is not digits is a schema violation, not a zero."""
    report = amel_report()
    report["assembly_stats"]["total_sequence_length"] = "about 225 Mb"
    with pytest.raises(NCBIDatasetsError, match="total_sequence_length"):
        ncbi_datasets.parse_assembly_report(report)

    report = amel_report()
    report["assembly_stats"]["gc_percent"] = "32.5"
    with pytest.raises(NCBIDatasetsError, match="gc_percent"):
        ncbi_datasets.parse_assembly_report(report)

    report = amel_report()
    report["assembly_info"] = ["not", "an", "object"]
    with pytest.raises(NCBIDatasetsError, match="assembly_info"):
        ncbi_datasets.parse_assembly_report(report)


def test_honey_bee_report_is_identical_across_endpoints() -> None:
    """Capture cross-check: accession lookup, reference-only lookup and RefSeq lookup agree."""
    by_accession = load(f"dataset_report_accession_{AMEL}.json")
    by_reference = load("dataset_report_taxon_Apis_mellifera_reference_only.json")
    by_source = load("dataset_report_taxon_Apis_mellifera_refseq_page1.json")

    assert by_accession["reports"] == by_reference["reports"] == by_source["reports"]


# --------------------------------------------------------------------------------------
# Pages: empty results, errors, paging
# --------------------------------------------------------------------------------------


def test_page_with_reports_and_totals() -> None:
    page = ncbi_datasets.parse_report_page(load(f"dataset_report_accession_{AMEL}.json"))

    assert len(page.reports) == 1 and page.reports[0]["accession"] == AMEL
    assert page.total_count == 1
    assert page.next_page_token is None
    assert page.errors == () and page.warnings == ()


@pytest.mark.parametrize(
    ("fixture", "content_type", "total_count"),
    [
        ("dataset_report_accession_GCF_999999999.9_not_found.json", None, None),
        ("dataset_report_taxon_unknown_name.json", None, None),
        ("dataset_report_taxon_2886930_reference_only_empty.json", "ASSM_ACC", None),
        ("dataset_report_taxon_Apis_mellifera_refseq_page2_empty.json", None, 1),
    ],
)
def test_empty_bodies_are_empty_pages(fixture: str, content_type: str | None, total_count: int | None) -> None:
    """Unknown accession, unknown taxon, no reference genome and an exhausted page are results, not errors."""
    page = ncbi_datasets.parse_report_page(load(fixture))

    assert page == ReportPage(reports=(), total_count=total_count, content_type=content_type)


def test_error_messages_are_surfaced_not_dropped() -> None:
    page = ncbi_datasets.parse_report_page(load("dataset_report_bioproject_PRJNA000000_error_message.json"))

    assert page.reports == ()
    assert len(page.errors) == 1
    assert "NO_ASSEMBLIES_FOR_BIOPROJECTS" in page.errors[0]
    assert "No genomes found for one or more of the requested accessions" in page.errors[0]
    assert "PRJNA000000" in page.errors[0]


def test_warning_messages_are_kept_separate_from_errors() -> None:
    page = ncbi_datasets.parse_report_page(
        {"messages": [{"warning": {"reason": "replaced", "message": "see GCF_1.2"}}]}
    )

    assert page.errors == ()
    assert page.warnings == ("replaced: see GCF_1.2",)


@pytest.mark.parametrize(
    "payload",
    [
        [],
        "reports",
        {"reports": {}},
        {"reports": [1]},
        {"reports": [{"accession": "x"}], "total_count": "many"},
        {"total_count": True},
        {"next_page_token": 5},
        {"messages": {}},
        {"messages": [{"error": "text"}]},
    ],
)
def test_unexpected_page_shapes_are_rejected(payload: object) -> None:
    with pytest.raises(NCBIDatasetsError):
        ncbi_datasets.parse_report_page(payload)


def test_paging_follows_tokens_while_reports_arrive() -> None:
    first = ncbi_datasets.parse_report_page(load("dataset_report_taxon_Apis_mellifera_assm_acc_page1.json"))
    second = ncbi_datasets.parse_report_page(load("dataset_report_taxon_Apis_mellifera_assm_acc_page2.json"))

    assert first.total_count == second.total_count == 13
    assert [r["accession"] for r in first.reports] == ["GCA_003254395.2"]
    assert [r["accession"] for r in second.reports] == [AMEL]

    token = ncbi_datasets.continuation_token(first)
    assert token == first.next_page_token
    assert ncbi_datasets.continuation_token(second, {token}) == second.next_page_token != token


def test_paging_stops_on_empty_trailing_page() -> None:
    """A one-result query with ``page_size=1`` still answers with a token; its page is empty."""
    last = ncbi_datasets.parse_report_page(load("dataset_report_taxon_Apis_mellifera_refseq_page1.json"))
    trailing = ncbi_datasets.parse_report_page(load("dataset_report_taxon_Apis_mellifera_refseq_page2_empty.json"))

    assert last.total_count == 1 and len(last.reports) == 1
    assert ncbi_datasets.continuation_token(last) == last.next_page_token  # cannot know yet
    assert trailing.reports == () and trailing.total_count == 1
    assert ncbi_datasets.continuation_token(trailing) is None
    # A token on a page with nothing in it ends the walk too.
    assert ncbi_datasets.continuation_token(ReportPage(next_page_token="abc")) is None


def test_paging_stops_on_repeated_token() -> None:
    page = ncbi_datasets.parse_report_page(load("dataset_report_taxon_Apis_mellifera_assm_acc_page1.json"))

    assert page.next_page_token is not None
    assert ncbi_datasets.continuation_token(page, {page.next_page_token}) is None


def test_default_page_holds_twenty_reports_of_many() -> None:
    page = ncbi_datasets.parse_report_page(load("dataset_report_taxon_2886930_assm_acc_page1.json"))

    assert len(page.reports) == ncbi_datasets.DEFAULT_PAGE_SIZE == 20
    assert page.total_count == 150
    assert ncbi_datasets.continuation_token(page) == page.next_page_token


@contextmanager
def serve_report_pages(
    report_pages: list[list[dict[str, Any]]], content_type: str = "ASSM_ACC"
) -> Iterator[tuple[str, list[str]]]:
    """Serve real saved report rows over local HTTP with deterministic page tokens."""
    payloads: dict[str, bytes] = {}
    for index, reports in enumerate(report_pages):
        payload: dict[str, Any] = {
            "reports": reports,
            "total_count": sum(len(page) for page in report_pages),
            "content_type": content_type,
        }
        if index + 1 < len(report_pages):
            payload["next_page_token"] = f"page-{index + 1}"
        payloads["" if index == 0 else f"page-{index}"] = json.dumps(payload).encode()
    requested_tokens: list[str] = []

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:
            token = parse_qs(urlsplit(self.path).query).get("page_token", [""])[0]
            requested_tokens.append(token)
            body = payloads.get(token)
            if body is None:
                self.send_response(404)
                self.send_header("Content-Length", "0")
                self.end_headers()
                return
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, format: str, *args: object) -> None:
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, kwargs={"poll_interval": 0.01}, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{server.server_port}", requested_tokens
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)
        assert not thread.is_alive()


@pytest.mark.parametrize("content_type", ["COMPLETE", "ASSM_ACC", "PAIRED_ACC"])
@pytest.mark.parametrize("location", ["within_page", "across_pages"])
@pytest.mark.parametrize("entrypoint", ["collect_reports", "iter_pages"])
def test_http_paging_rejects_repeated_assembly_accessions(content_type: str, location: str, entrypoint: str) -> None:
    """Derive an inconsistent response by repeating one real saved assembly row."""
    row = load("dataset_report_taxon_2886930_assm_acc_page1.json")["reports"][0]
    pages = [[row, row]] if location == "within_page" else [[row], [row]]
    request = ncbi_datasets.build_taxon_report_request("Apis mellifera", returned_content=content_type)

    with serve_report_pages(pages, content_type) as (url, requested), NCBIDatasetsClient(timeout=2) as client:
        client.base_url = url
        with pytest.raises(NCBIDatasetsError, match=f"assembly accession {row['accession']}.*more than once"):
            if entrypoint == "collect_reports":
                client.collect_reports(request)
            else:
                list(client.iter_pages(request))

    assert requested == ([""] if location == "within_page" else ["", "page-1"])


@pytest.mark.parametrize("entrypoint", ["collect_reports", "iter_pages"])
@pytest.mark.parametrize("identities", ["paired", "different_versions", "sequence_rows"])
def test_http_paging_preserves_distinct_report_identities(entrypoint: str, identities: str) -> None:
    """Paired databases and accession versions are distinct; sequence rows have their own keys."""
    if identities == "sequence_rows":
        rows = load(f"sequence_reports_{AMEL}_page1.json")["reports"][:2]
        request = ncbi_datasets.build_sequence_report_request(AMEL)
    else:
        rows = load("dataset_report_taxon_2886930_assm_acc_page1.json")["reports"][:2]
        if identities == "different_versions":
            # Derive a distinct revision from the saved accession, without
            # claiming that revision exists on the public service.
            rows[1] = rows[0] | {"accession": rows[0]["accession"].split(".")[0] + ".2"}
        else:
            rows = [row | {"paired_accession": rows[1 - index]["accession"]} for index, row in enumerate(rows)]
        request = ncbi_datasets.build_taxon_report_request("Apis mellifera", returned_content="PAIRED_ACC")

    with serve_report_pages([[rows[0]], [rows[1]]]) as (url, requested), NCBIDatasetsClient(timeout=2) as client:
        client.base_url = url
        if entrypoint == "collect_reports":
            collection = client.collect_reports(request)
            assert collection.total_count == 2
            actual = collection.reports
        else:
            actual = [report for page in client.iter_pages(request) for report in page.reports]

    assert actual == rows
    assert requested == ["", "page-1"]


def test_http_paging_limit_does_not_request_an_unneeded_later_page() -> None:
    row = load("dataset_report_taxon_2886930_assm_acc_page1.json")["reports"][0]
    request = ncbi_datasets.build_taxon_report_request("Apis mellifera")
    with serve_report_pages([[row], [row]]) as (url, requested), NCBIDatasetsClient(timeout=2) as client:
        client.base_url = url
        collection = client.collect_reports(request, limit=1)

    assert collection.reports == [row] and collection.total_count == 2
    assert requested == [""]


# --------------------------------------------------------------------------------------
# Query builders
# --------------------------------------------------------------------------------------


def test_queries_reproduce_the_urls_that_captured_the_fixtures() -> None:
    token_from_page1 = load("dataset_report_taxon_Apis_mellifera_refseq_page1.json")["next_page_token"]
    amel_taxon = f"{BASE}/genome/taxon/Apis%20mellifera/dataset_report"
    refseq_page1 = ncbi_datasets.build_taxon_report_request("Apis mellifera", assembly_source="refseq", page_size=1)
    assm_acc_page1 = ncbi_datasets.build_taxon_report_request(
        "Apis mellifera", returned_content="ASSM_ACC", page_size=1
    )
    token_from_acc_page1 = load("dataset_report_taxon_Apis_mellifera_assm_acc_page1.json")["next_page_token"]

    assert ncbi_datasets.build_accession_report_request(AMEL).url() == f"{BASE}/genome/accession/{AMEL}/dataset_report"
    assert refseq_page1.url() == f"{amel_taxon}?page_size=1&filters.assembly_source=refseq"
    assert (
        refseq_page1.with_page_token(token_from_page1).url()
        == f"{amel_taxon}?page_size=1&filters.assembly_source=refseq&page_token={token_from_page1}"
    )
    assert assm_acc_page1.url() == f"{amel_taxon}?page_size=1&returned_content=ASSM_ACC"
    assert (
        assm_acc_page1.with_page_token(token_from_acc_page1).url()
        == f"{amel_taxon}?page_size=1&returned_content=ASSM_ACC&page_token={token_from_acc_page1}"
    )
    assert (
        ncbi_datasets.build_sequence_report_request(AMEL, page_size=3).url()
        == f"{BASE}/genome/accession/{AMEL}/sequence_reports?page_size=3"
    )


def test_taxon_query_carries_every_filter_it_is_given() -> None:
    request = ncbi_datasets.build_taxon_report_request(
        "Apis mellifera",
        reference_only=True,
        assembly_source="refseq",
        exclude_paired_reports=True,
        returned_content="PAIRED_ACC",
        page_size=50,
    )

    assert request.path == "/genome/taxon/Apis%20mellifera/dataset_report"
    assert dict(request.params) == {
        "page_size": "50",
        "filters.reference_only": "true",
        "filters.assembly_source": "refseq",
        "filters.exclude_paired_reports": "true",
        "returned_content": "PAIRED_ACC",
    }


def test_taxon_query_defaults_send_only_page_size() -> None:
    request = ncbi_datasets.build_taxon_report_request("Apis mellifera")

    assert dict(request.params) == {"page_size": "20"}
    assert request.with_page_token(None) is request
    assert request.with_page_token("") is request
    assert dict(request.with_page_token("tok").params) == {"page_size": "20", "page_token": "tok"}


def test_taxon_name_cannot_change_the_endpoint() -> None:
    path = urlsplit(ncbi_datasets.build_taxon_report_request("../gene/id/1").url()).path
    assert path == "/datasets/v2/genome/taxon/..%2Fgene%2Fid%2F1/dataset_report"

    url = ncbi_datasets.build_taxon_report_request("Apis?page_size=1000#x").url()
    assert urlsplit(url).path == "/datasets/v2/genome/taxon/Apis%3Fpage_size%3D1000%23x/dataset_report"
    assert parse_qs(urlsplit(url).query) == {"page_size": ["20"]}


def test_a_comma_in_a_taxon_name_would_query_several_taxa() -> None:
    """The service splits the path segment on commas, percent-encoded or not, so the builder refuses them."""
    page = ncbi_datasets.parse_report_page(load("dataset_report_taxon_two_names_joined_by_comma.json"))

    assert [report["accession"] for report in page.reports] == [AMEL, "GCF_055798465.1"]
    assert page.total_count == 2
    with pytest.raises(ValueError, match="single organism"):
        ncbi_datasets.build_taxon_report_request("Apis mellifera,Pogonomyrmex barbatus")


def test_taxon_accepts_taxonomy_id() -> None:
    assert ncbi_datasets.build_taxon_report_request(7460).path == "/genome/taxon/7460/dataset_report"
    assert ncbi_datasets.build_taxon_report_request(" 7460 ").path == "/genome/taxon/7460/dataset_report"


@pytest.mark.parametrize(
    "taxon", ["", "   ", "Apis mellifera, Pogonomyrmex barbatus", "Apis\x00mellifera", True, None, 1.5]
)
def test_taxon_input_validation(taxon: Any) -> None:
    with pytest.raises(ValueError):
        ncbi_datasets.normalize_taxon(taxon)


@pytest.mark.parametrize("page_size", [0, -1, 1001, True, "5", 2.0, None])
def test_page_size_limits(page_size: Any) -> None:
    with pytest.raises(ValueError, match="page_size"):
        ncbi_datasets.build_taxon_report_request("Apis mellifera", page_size=page_size)
    with pytest.raises(ValueError, match="page_size"):
        ncbi_datasets.build_sequence_report_request(AMEL, page_size=page_size)


def test_page_size_bounds_are_accepted() -> None:
    for page_size in (1, ncbi_datasets.MAX_PAGE_SIZE):
        request = ncbi_datasets.build_taxon_report_request("Apis mellifera", page_size=page_size)
        assert request.params["page_size"] == str(page_size)


def test_assembly_source_and_content_validation() -> None:
    for source in ("all", "refseq", "genbank"):
        ncbi_datasets.build_taxon_report_request("Apis mellifera", assembly_source=source)
    with pytest.raises(ValueError, match="assembly_source"):
        ncbi_datasets.build_taxon_report_request("Apis mellifera", assembly_source="RefSeq")
    with pytest.raises(ValueError, match="returned_content"):
        ncbi_datasets.build_taxon_report_request("Apis mellifera", returned_content="FULL")


@pytest.mark.parametrize(
    ("accession", "valid"),
    [
        (AMEL, True),
        ("GCA_000001405.29", True),
        ("GCF_003254395", False),  # the service reports versionless accessions as invalid
        ("gcf_003254395.2", False),
        ("GCF_003254395.2\n", False),
        (" GCF_003254395.2", False),
        ("GCF_00325439.2", False),
        ("NC_000001.11", False),
        ("chr1", False),
        ("../genome/accession/GCF_003254395.2", False),
        ("", False),
        (None, False),
    ],
)
def test_accession_validation(accession: Any, valid: bool) -> None:
    if valid:
        assert ncbi_datasets.normalize_assembly_accession(accession) == accession
    else:
        with pytest.raises(ValueError, match="Invalid genome accession"):
            ncbi_datasets.normalize_assembly_accession(accession)


@pytest.mark.parametrize("include", [None, [], "", ["", " "], ()])
def test_package_query_requests_hydrated_genomic_fasta_by_default(include: Any) -> None:
    request = ncbi_datasets.build_package_request(AMEL, include)

    assert request.path == f"/genome/accession/{AMEL}/download"
    assert dict(request.params) == {
        "include_annotation_type": ["GENOME_FASTA"],
        "hydrated": "FULLY_HYDRATED",
    }
    assert request.url() == (
        f"{BASE}/genome/accession/{AMEL}/download" "?include_annotation_type=GENOME_FASTA&hydrated=FULLY_HYDRATED"
    )


def test_package_query_matches_urls_echoed_by_the_service() -> None:
    """The service's own ``download_summary`` echoes the ``download`` URL for an include selection."""
    gff = load(f"download_summary_{AMEL}_genome_gff3.json")["hydrated"]["url"]
    mine = ncbi_datasets.build_package_request(AMEL, ["genome", "gff3"]).url()
    assert gff == mine

    many = load(f"download_summary_{AMEL}_rna_protein_cds_gtf_gbff_seqreport.json")["hydrated"]["url"]
    mine = ncbi_datasets.build_package_request(AMEL, "rna,protein,cds,gtf,gbff,seq-report").url()
    assert many == mine

    # The service's default summary omits file selection. The download endpoint
    # then returns only metadata, despite the CLI's default including a genome.
    default = load(f"download_summary_{AMEL}_default.json")["hydrated"]["url"]
    assert parse_qs(urlsplit(default).query) == {"hydrated": ["FULLY_HYDRATED"]}


def test_include_aliases_agree_with_service_command_lines() -> None:
    """Every ``--include`` name the service prints maps to the file type it was echoed with."""
    summary = load(f"download_summary_{AMEL}_rna_protein_cds_gtf_gbff_seqreport.json")["hydrated"]
    echoed_types = parse_qs(urlsplit(summary["url"]).query)["include_annotation_type"]
    cli_names = summary["cli_download_command_line"].split("--include ")[1].split(",")

    assert sorted(ncbi_datasets.resolve_include_types(cli_names)) == sorted(echoed_types)

    gff = load(f"download_summary_{AMEL}_genome_gff3.json")["hydrated"]
    assert ncbi_datasets.resolve_include_types(gff["cli_download_command_line"].split("--include ")[1].split(",")) == [
        "GENOME_GFF",
        "GENOME_FASTA",
    ]


def test_legacy_include_parameter_changes_nothing_on_the_service() -> None:
    """The old code sent ``include=genome,gff3``; the service answered as if nothing had been asked."""
    default = load(f"download_summary_{AMEL}_default.json")
    legacy = load(f"download_summary_{AMEL}_legacy_include_ignored.json")
    proper = load(f"download_summary_{AMEL}_genome_gff3.json")

    assert legacy == default
    assert proper != default
    assert "include" not in ncbi_datasets.build_package_request(AMEL, ["genome", "gff3"]).params


def test_include_names_accept_cli_and_service_spellings() -> None:
    assert ncbi_datasets.resolve_include_types("genome,gff3") == ["GENOME_FASTA", "GENOME_GFF"]
    assert ncbi_datasets.resolve_include_types(["GENOME_GFF", "genome", "gff3"]) == ["GENOME_GFF", "GENOME_FASTA"]
    assert ncbi_datasets.resolve_include_types(["seq_report", "SEQUENCE_REPORT", " Seq-Report "]) == ["SEQUENCE_REPORT"]
    assert ncbi_datasets.resolve_include_types(["gbff", ""]) == ["GENOME_GBFF"]
    assert ncbi_datasets.resolve_include_types(None) == []
    assert ncbi_datasets.resolve_include_types([]) == []
    assert sorted(ncbi_datasets.INCLUDE_FILE_TYPES) == sorted(
        {
            "GENOME_FASTA",
            "GENOME_GFF",
            "GENOME_GTF",
            "GENOME_GBFF",
            "RNA_FASTA",
            "PROT_FASTA",
            "CDS_FASTA",
            "SEQUENCE_REPORT",
        }
    )


@pytest.mark.parametrize("include", [["bogus"], ["none"], "genome,fastq", [1], ["DEFAULT"]])
def test_unknown_include_name_is_an_error(include: Any) -> None:
    with pytest.raises(ValueError):
        ncbi_datasets.resolve_include_types(include)


# --------------------------------------------------------------------------------------
# Reference selection
# --------------------------------------------------------------------------------------


def test_reference_selection_prefers_the_labelled_reference() -> None:
    assert ncbi_datasets.select_reference_assembly([phix_report(), amel_report()]) == amel_report()
    assert ncbi_datasets.select_reference_assembly([amel_report(), phix_report()]) == amel_report()


def test_reference_selection_of_nothing_is_none() -> None:
    assert ncbi_datasets.select_reference_assembly([]) is None
    assert ncbi_datasets.select_reference_assembly([{"source_database": "SOURCE_DATABASE_REFSEQ"}]) is None


def test_reference_selection_ranks_superseded_below_current() -> None:
    """Derived from the honey bee report by marking it superseded."""
    superseded = amel_report()
    superseded["assembly_info"]["assembly_status"] = "previous"

    assert ncbi_datasets.select_reference_assembly([superseded, phix_report()]) == phix_report()
    assert ncbi_datasets.select_reference_assembly([superseded]) == superseded


def test_reference_selection_on_abbreviated_reports_prefers_the_refseq_copy() -> None:
    page = ncbi_datasets.parse_report_page(load("dataset_report_taxon_2886930_assm_acc_page1.json"))

    chosen = ncbi_datasets.select_reference_assembly(page.reports)

    assert chosen is not None and chosen["accession"] == "GCF_000819615.1"


def test_reference_selection_breaks_ties_by_representative_and_date() -> None:
    """Derived from the honey bee report: representative < reference, and newer release wins a tie."""
    representative = amel_report()
    representative["assembly_info"]["refseq_category"] = "representative genome"
    newer = amel_report()
    newer["accession"] = "GCF_999999998.1"
    newer["assembly_info"]["release_date"] = "2030-01-01"

    assert ncbi_datasets.select_reference_assembly([representative, amel_report()]) == amel_report()
    assert ncbi_datasets.select_reference_assembly([amel_report(), newer]) == newer


# --------------------------------------------------------------------------------------
# Sequence reports
# --------------------------------------------------------------------------------------


def test_sequence_rows_are_keyed_by_accession() -> None:
    page = ncbi_datasets.parse_report_page(load(f"sequence_reports_{AMEL}_page1.json"))

    assert page.total_count == 177 and len(page.reports) == 3  # the rest lives on later pages
    assert ncbi_datasets.chromosome_lengths_from_sequence_reports(page.reports) == {
        "NC_037638.1": 27754200,
        "NC_037639.1": 16089512,
        "NC_037640.1": 13619445,
    }


def test_a_walk_that_lost_reports_is_not_returned_as_complete() -> None:
    """Page 1 of the sequence report holds 3 of the 177 sequences NCBI announces for the assembly."""
    page = ncbi_datasets.parse_report_page(load(f"sequence_reports_{AMEL}_page1.json"))
    partial = ncbi_datasets.ReportCollection(list(page.reports), page.total_count)

    with pytest.raises(NCBIDatasetsError, match="announced 177 sequences.*3 were received"):
        ncbi_datasets.require_complete(partial, "sequences")
    # With a limit, only the part asked for has to arrive.
    assert ncbi_datasets.require_complete(partial, "sequences", limit=3) == partial.reports
    with pytest.raises(NCBIDatasetsError):
        ncbi_datasets.require_complete(partial, "sequences", limit=4)
    # No announced total, nothing to compare against; and a complete walk passes.
    assert ncbi_datasets.require_complete(ncbi_datasets.ReportCollection(partial.reports, None), "sequences")
    assert ncbi_datasets.require_complete(ncbi_datasets.ReportCollection(partial.reports, 3), "sequences")
    assert ncbi_datasets.require_complete(ncbi_datasets.ReportCollection([], None), "sequences") == []


def test_sequence_rows_can_be_keyed_by_chromosome_name() -> None:
    page = ncbi_datasets.parse_report_page(load(f"sequence_reports_{AMEL}_page1.json"))

    assert ncbi_datasets.chromosome_lengths_from_sequence_reports(page.reports, names="chromosome") == {
        "LG1": 27754200,
        "LG2": 16089512,
        "LG3": 13619445,
    }


def test_chromosome_naming_never_lets_one_sequence_overwrite_another() -> None:
    """Derived from a saved row: scaffolds and repeated names fall back to the accession.

    The Entrez code this replaces keyed ``chromosome N`` titles by name, so a chromosome's patch sequences
    replaced the chromosome itself.
    """
    rows = load(f"sequence_reports_{AMEL}_page1.json")["reports"]
    scaffold = copy.deepcopy(rows[0]) | {"role": "unlocalized-scaffold", "refseq_accession": "NW_000000001.1"}
    twin = copy.deepcopy(rows[0]) | {"refseq_accession": "NC_999999999.1", "length": 5}
    unplaced = copy.deepcopy(rows[0]) | {
        "role": "unplaced-scaffold",
        "chr_name": "Un",
        "refseq_accession": "NW_000000002.1",
    }

    lengths = ncbi_datasets.chromosome_lengths_from_sequence_reports(
        [rows[0], scaffold, twin, unplaced], names="chromosome"
    )

    assert lengths == {
        "LG1": 27754200,  # the chromosome keeps its name
        "NW_000000001.1": 27754200,  # unlocalized scaffold shares the chromosome name: accession
        "NC_999999999.1": 5,  # second assembled molecule called LG1: accession
        "NW_000000002.1": 27754200,  # unplaced: accession
    }


def test_sequence_rows_without_accession_length_or_unique_accession_are_rejected() -> None:
    row = load(f"sequence_reports_{AMEL}_page1.json")["reports"][0]

    no_accession = {k: v for k, v in row.items() if k not in ("refseq_accession", "genbank_accession")}
    with pytest.raises(NCBIDatasetsError, match="without an accession"):
        ncbi_datasets.chromosome_lengths_from_sequence_reports([no_accession])

    no_length = {k: v for k, v in row.items() if k != "length"}
    with pytest.raises(NCBIDatasetsError, match="without a length"):
        ncbi_datasets.chromosome_lengths_from_sequence_reports([no_length])

    with pytest.raises(NCBIDatasetsError, match="more than once"):
        ncbi_datasets.chromosome_lengths_from_sequence_reports([row, row])

    genbank_only = {k: v for k, v in row.items() if k != "refseq_accession"}
    assert ncbi_datasets.chromosome_lengths_from_sequence_reports([genbank_only]) == {"CM009931.2": 27754200}


def test_sequence_naming_mode_is_validated() -> None:
    with pytest.raises(ValueError, match="names"):
        ncbi_datasets.chromosome_lengths_from_sequence_reports([], names="ucsc")


# --------------------------------------------------------------------------------------
# Error text
# --------------------------------------------------------------------------------------


def test_error_description_uses_the_service_message() -> None:
    """``rpcStatus`` is the documented error object: ``code``, ``message``, ``details``."""
    message = describe_http_error(
        400, "Bad Request", '{"code": 3, "message": "invalid include_annotation_type", "details": []}', "GET /x"
    )

    assert message == "GET /x failed: HTTP 400 Bad Request: invalid include_annotation_type"


def test_error_description_handles_plain_text_empty_and_long_bodies() -> None:
    assert (
        describe_http_error(404, "Not Found", "no such\n  resource ", "op")
        == "op failed: HTTP 404 Not Found: no such resource"
    )
    assert describe_http_error(503, None, "", "op") == "op failed: HTTP 503"
    long = describe_http_error(500, "Server Error", "x" * 5000, "op")
    assert len(long) < 400 and long.endswith("...")


def test_credentials_are_redacted_from_error_text() -> None:
    url = "https://eutils.example/esearch.fcgi?db=nucleotide&api_key=SECRET123&email=me%40example.org&term=x"
    assert "SECRET123" not in redact_secrets(url) and "me%40example.org" not in redact_secrets(url)
    assert "db=nucleotide" in redact_secrets(url) and "term=x" in redact_secrets(url)

    failure = requests.ConnectionError(f"HTTPSConnectionPool: Max retries exceeded with url: {url} (Caused by ...)")
    assert "SECRET123" not in describe_transport_error(failure)
    assert "SECRET123" not in describe_http_error(401, "Unauthorized", f"bad key in {url}", "op")
    assert redact_secrets("API_KEY=abc") == "API_KEY=***"


# --------------------------------------------------------------------------------------
# Package archives
# --------------------------------------------------------------------------------------

PACKAGE_MEMBERS = {
    "README.md": b"package readme",
    "md5sum.txt": b"checksums",
    "ncbi_dataset/data/dataset_catalog.json": b"{}",
    f"ncbi_dataset/data/{AMEL}/{AMEL}_Amel_HAv3.1_genomic.fna": bytes(range(256)),
    f"ncbi_dataset/data/{AMEL}/genomic.gff": b"gff payload",
}


def build_package_zip(path: Path) -> Path:
    """A real zip with the layout NCBI documents for a genome package (opaque payload bytes)."""
    report_line = json.dumps(amel_report()).encode()
    members = {"ncbi_dataset/data/assembly_data_report.jsonl": report_line + b"\n", **PACKAGE_MEMBERS}
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, payload in members.items():
            archive.writestr(name, payload)
    return path


def test_package_zip_extracts_with_the_documented_layout(tmp_path: Path) -> None:
    archive = build_package_zip(tmp_path / f"{AMEL}.zip")
    destination = tmp_path / AMEL

    written = ncbi_datasets.extract_package_zip(archive, destination)

    expected = {"ncbi_dataset/data/assembly_data_report.jsonl", *PACKAGE_MEMBERS}
    assert {p.relative_to(destination.resolve()).as_posix() for p in written} == expected
    assert (destination / "README.md").read_bytes() == b"package readme"
    genomic = destination / "ncbi_dataset" / "data" / AMEL / f"{AMEL}_Amel_HAv3.1_genomic.fna"
    assert genomic.read_bytes() == bytes(range(256))
    report_line = (destination / "ncbi_dataset/data/assembly_data_report.jsonl").read_text().splitlines()[0]
    assert json.loads(report_line)["accession"] == AMEL


def test_validate_genome_files_sees_files_in_the_package_layout(tmp_path: Path) -> None:
    """Files sit three levels down and the FASTA ends in ``.fna``; the old check looked in neither place."""
    destination = tmp_path / AMEL
    ncbi_datasets.extract_package_zip(build_package_zip(tmp_path / "p.zip"), destination)

    results = genomes.validate_genome_files(destination)

    assert results == {"valid": True, "fasta_files": 1, "gff_files": 1, "total_files": 6, "issues": []}


def test_validate_genome_files_recognises_gzip_and_flat_layouts(tmp_path: Path) -> None:
    (tmp_path / "a.fa").write_bytes(b"x")
    (tmp_path / "b.fasta").write_bytes(b"x")
    (tmp_path / "c.fna.gz").write_bytes(b"x")
    (tmp_path / "d.gff3").write_bytes(b"x")
    (tmp_path / "e.GFF.GZ").write_bytes(b"x")
    (tmp_path / "f.txt").write_bytes(b"x")

    results = genomes.validate_genome_files(tmp_path)

    assert (results["valid"], results["fasta_files"], results["gff_files"], results["total_files"]) == (True, 3, 2, 6)


def test_validate_genome_files_flags_a_missing_fasta_and_a_missing_directory(tmp_path: Path) -> None:
    (tmp_path / "assembly_data_report.jsonl").write_bytes(b"{}")

    results = genomes.validate_genome_files(tmp_path)

    assert results["valid"] is False and results["issues"] == ["No FASTA files found"]
    assert genomes.validate_genome_files(tmp_path / "absent") == {"valid": False, "error": "Directory does not exist"}
    assert genomes.validate_genome_files(tmp_path / "assembly_data_report.jsonl") == {
        "valid": False,
        "error": "Path is not a directory",
    }


def symlink_member(name: str) -> zipfile.ZipInfo:
    info = zipfile.ZipInfo(name)
    info.external_attr = (stat.S_IFLNK | 0o777) << 16
    return info


@pytest.mark.parametrize(
    "unsafe_name",
    [
        pytest.param("../evil.txt", id="parent"),
        pytest.param("ncbi_dataset/../../evil.txt", id="nested-parent"),
        pytest.param("/abs/evil.txt", id="absolute"),
        pytest.param("C:/evil.txt", id="drive-letter"),
        pytest.param("ncbi_dataset\\..\\..\\evil.txt", id="backslash-parent"),
    ],
)
def test_package_zip_refuses_unsafe_paths_before_writing_anything(tmp_path: Path, unsafe_name: str) -> None:
    archive = tmp_path / "bad.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        zf.writestr("README.md", b"fine")  # a valid member before the unsafe one must not be written either
        zf.writestr(unsafe_name, b"payload")
    destination = tmp_path / "out"

    with pytest.raises(NCBIDatasetsError, match="Unsafe path"):
        ncbi_datasets.extract_package_zip(archive, destination)

    assert not destination.exists()
    assert not (tmp_path / "evil.txt").exists()


def test_package_zip_refuses_symlink_members(tmp_path: Path) -> None:
    archive = tmp_path / "link.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        zf.writestr(symlink_member("ncbi_dataset/link"), "/etc/passwd")
    destination = tmp_path / "out"

    with pytest.raises(NCBIDatasetsError, match="Symlink"):
        ncbi_datasets.extract_package_zip(archive, destination)

    assert not destination.exists()


def test_package_zip_that_is_cut_short_or_not_a_zip_raises(tmp_path: Path) -> None:
    whole = build_package_zip(tmp_path / "whole.zip").read_bytes()
    truncated = tmp_path / "truncated.zip"
    truncated.write_bytes(whole[: len(whole) // 2])
    html = tmp_path / "error.zip"
    html.write_bytes(b"<html><body>Service unavailable</body></html>")

    for broken in (truncated, html):
        with pytest.raises(NCBIDatasetsError, match="complete zip archive"):
            ncbi_datasets.extract_package_zip(broken, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_package_zip_extraction_merges_into_an_existing_directory(tmp_path: Path) -> None:
    destination = tmp_path / AMEL
    destination.mkdir()
    (destination / "keep.txt").write_bytes(b"kept")

    ncbi_datasets.extract_package_zip(build_package_zip(tmp_path / "p.zip"), destination)

    assert (destination / "keep.txt").read_bytes() == b"kept"
    assert (destination / "md5sum.txt").exists()


# --------------------------------------------------------------------------------------
# Public functions: arguments are checked before anything is sent
# --------------------------------------------------------------------------------------


def test_listing_rejects_bad_taxa_and_sources(unreachable_client: NCBIDatasetsClient) -> None:
    with pytest.raises(ValueError):
        genomes.list_genome_assemblies("", client=unreachable_client)
    with pytest.raises(ValueError, match="assembly_source"):
        genomes.list_genome_assemblies("Apis mellifera", assembly_source="ensembl", client=unreachable_client)


def test_package_rejects_bad_arguments_before_touching_the_disk(
    tmp_path: Path, unreachable_client: NCBIDatasetsClient
) -> None:
    target = tmp_path / "out"

    with pytest.raises(ValueError, match="Invalid genome accession"):
        genomes.download_genome_package("NC_000001", target, client=unreachable_client)
    with pytest.raises(ValueError, match="Unsupported genome package file type"):
        genomes.download_genome_package(AMEL, target, include=["bogus"], client=unreachable_client)

    assert not target.exists()


def test_reference_lookup_validates_include_before_any_query(
    tmp_path: Path, unreachable_client: NCBIDatasetsClient
) -> None:
    with pytest.raises(ValueError, match="Unsupported genome package file type"):
        genomes.download_reference_genome("Apis mellifera", tmp_path, include=["bogus"], client=unreachable_client)
    with pytest.raises(ValueError):
        genomes.download_reference_genome("", tmp_path, client=unreachable_client)


def test_chromosome_lengths_keep_the_contact_policy(
    monkeypatch: pytest.MonkeyPatch, unreachable_client: NCBIDatasetsClient
) -> None:
    monkeypatch.delenv("NCBI_EMAIL", raising=False)

    with pytest.raises(NCBIContactError):
        genomes.get_chromosome_lengths(AMEL, client=unreachable_client)
    with pytest.raises(ValueError, match="names"):
        genomes.get_chromosome_lengths(AMEL, email="tests@example.org", names="ucsc", client=unreachable_client)
    with pytest.raises(ValueError, match="Invalid genome accession"):
        genomes.get_chromosome_lengths("chr1", email="tests@example.org", client=unreachable_client)


# --------------------------------------------------------------------------------------
# Public functions: a failed query raises, it no longer returns an empty result
# --------------------------------------------------------------------------------------


def test_listing_failure_raises_instead_of_returning_empty(unreachable_client: NCBIDatasetsClient) -> None:
    with pytest.raises(requests.RequestException):
        genomes.list_genome_assemblies("Apis mellifera", max_results=5, client=unreachable_client)


def test_metadata_failure_raises_instead_of_returning_an_error_dict(unreachable_client: NCBIDatasetsClient) -> None:
    with pytest.raises(requests.RequestException):
        genomes.get_genome_metadata(AMEL, client=unreachable_client)


def test_reference_lookup_failure_raises_instead_of_returning_none(
    tmp_path: Path, unreachable_client: NCBIDatasetsClient
) -> None:
    with pytest.raises(requests.RequestException):
        genomes.download_reference_genome("Apis mellifera", tmp_path, client=unreachable_client)
    assert list(tmp_path.iterdir()) == []


def test_chromosome_lengths_failure_raises_instead_of_returning_empty(unreachable_client: NCBIDatasetsClient) -> None:
    with pytest.raises(requests.RequestException):
        genomes.get_chromosome_lengths(AMEL, email="tests@example.org", client=unreachable_client)


def test_package_failure_raises_and_leaves_no_partial_zip(
    tmp_path: Path, unreachable_client: NCBIDatasetsClient
) -> None:
    with pytest.raises(NCBIDatasetsError):
        genomes.download_genome_package(AMEL, tmp_path, client=unreachable_client)

    assert list(tmp_path.glob("*.zip")) == [] and not (tmp_path / AMEL).exists()


def test_best_effort_failure_names_the_cause_and_leaves_nothing(
    tmp_path: Path, unreachable_client: NCBIDatasetsClient
) -> None:
    with pytest.raises(RuntimeError, match="All genome download methods failed.*NCBI Datasets API"):
        genomes.download_genome_package_best_effort(AMEL, tmp_path, ["genome"], client=unreachable_client)

    assert not (tmp_path / AMEL).exists()


def test_failures_subclass_the_library_exception_type(unreachable_client: NCBIDatasetsClient) -> None:
    request = ncbi_datasets.build_accession_report_request(AMEL)

    with pytest.raises(NCBIDatasetsError) as caught:
        unreachable_client.get_json(request)

    assert isinstance(caught.value, NCBIRequestError) and isinstance(caught.value, requests.RequestException)
    assert caught.value.status_code is None
    assert "GET /genome/accession/" in str(caught.value)


def test_missing_assembly_is_a_lookup_error_not_a_transport_failure() -> None:
    assert issubclass(AssemblyNotFoundError, LookupError)
    assert not issubclass(AssemblyNotFoundError, requests.RequestException)


# --------------------------------------------------------------------------------------
# Client
# --------------------------------------------------------------------------------------


def test_client_sends_key_and_contact_as_headers_only(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    monkeypatch.delenv("NCBI_EMAIL", raising=False)

    with NCBIDatasetsClient(api_key="KEY-123", email="tests@example.org") as client:
        assert client.session.headers["api-key"] == "KEY-123"
        assert client.session.headers["From"] == "tests@example.org"
        assert client.session.headers["User-Agent"].startswith("metainformant/")
        request = ncbi_datasets.build_taxon_report_request("Apis mellifera")
        assert "KEY-123" not in request.url(client.base_url)
        assert "tests@example.org" not in request.url(client.base_url)

    with NCBIDatasetsClient() as anonymous:
        assert "api-key" not in anonymous.session.headers
        assert "From" not in anonymous.session.headers  # a contact is never invented


def test_client_reads_key_and_contact_from_the_environment(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("NCBI_API_KEY", " ENV-KEY ")
    monkeypatch.setenv("NCBI_EMAIL", "env@example.org")

    with NCBIDatasetsClient() as client:
        assert client.api_key == "ENV-KEY" and client.email == "env@example.org"
        assert client.session.headers["api-key"] == "ENV-KEY"


@pytest.mark.parametrize("timeout", [0, -1, -0.5])
def test_client_rejects_non_positive_timeouts(timeout: float) -> None:
    with pytest.raises(ValueError, match="timeout"):
        NCBIDatasetsClient(timeout=timeout)


def test_collect_reports_rejects_a_zero_limit(unreachable_client: NCBIDatasetsClient) -> None:
    with pytest.raises(ValueError, match="limit"):
        unreachable_client.collect_reports(ncbi_datasets.build_accession_report_request(AMEL), limit=0)
