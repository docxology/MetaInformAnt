"""NCBI Datasets v2 REST API: request builders, response parsers and a small client.

Deciding *what to ask* and *how to read the answer* is done by pure functions,
so it can be tested on saved responses without a network: the
``build_*_request`` functions return a :class:`DatasetsRequest`, and the
``parse_*`` functions turn a decoded JSON body into plain Python values.
:class:`NCBIDatasetsClient` is the only piece that talks to the network, and it
is shared by every genome function in :mod:`.genomes`.

Endpoints used (OpenAPI specification:
https://www.ncbi.nlm.nih.gov/datasets/docs/v2/openapi3/openapi3.docs.yaml):

``GET /genome/accession/{accessions}/dataset_report``
    Assembly report for one assembly accession.
``GET /genome/taxon/{taxons}/dataset_report``
    Assembly reports for a taxon (name or NCBI Taxonomy ID); paged.
``GET /genome/accession/{accession}/sequence_reports``
    One row per sequence (chromosome, scaffold, ...) of an assembly; paged.
``GET /genome/accession/{accessions}/download``
    Genome data package (zip).

Response conventions that the parsers rely on (observed on the live service
and described by the OpenAPI schema ``v2reportsAssemblyDataReportPage``):

* An accession or taxon with nothing to report answers HTTP 200 with a body
  such as ``{}`` or ``{"content_type": "ASSM_ACC"}``: there is no ``reports``
  key and no ``total_count``.  Empty is not an error.
* ``total_sequence_length`` and the other ``uint64`` fields are JSON strings;
  the smaller counters are JSON numbers.
* A page can carry ``next_page_token`` although no further results exist; the
  follow-up page is then ``{"total_count": N}`` with no ``reports``.
* Errors are HTTP 4xx/5xx with a ``text/plain`` body, or, for some lookups,
  HTTP 200 with ``messages[].error`` entries.
"""

from __future__ import annotations

import os
import re
import shutil
import stat
import tempfile
import zipfile
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any, Collection, Dict, Iterable, Iterator, List, Mapping, NamedTuple, Optional, Sequence, Set, Tuple
from urllib.parse import quote

import requests

from metainformant import __version__

from ._http import NCBIRequestError, checked_get, describe_transport_error, get_json_object, resolve_api_key

__all__ = [
    "ASSEMBLY_SOURCES",
    "DATASETS_API_BASE_URL",
    "DEFAULT_PAGE_SIZE",
    "INCLUDE_FILE_TYPES",
    "MAX_PAGE_SIZE",
    "SEQUENCE_NAME_MODES",
    "AssemblyNotFoundError",
    "DatasetsRequest",
    "NCBIDatasetsClient",
    "NCBIDatasetsError",
    "ReportCollection",
    "ReportPage",
    "build_accession_report_request",
    "build_package_request",
    "build_sequence_report_request",
    "build_taxon_report_request",
    "chromosome_lengths_from_sequence_reports",
    "client_scope",
    "continuation_token",
    "extract_package_zip",
    "normalize_assembly_accession",
    "normalize_taxon",
    "parse_assembly_report",
    "parse_report_page",
    "require_complete",
    "resolve_include_types",
    "select_reference_assembly",
]

DATASETS_API_BASE_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2"
DEFAULT_PAGE_SIZE = 20  # server default for report endpoints
MAX_PAGE_SIZE = 1000  # documented maximum for report endpoints
ASSEMBLY_SOURCES = ("all", "refseq", "genbank")
SEQUENCE_NAME_MODES = ("accession", "chromosome")
_RETURNED_CONTENT = ("COMPLETE", "ASSM_ACC", "PAIRED_ACC")
_MAX_PAGES = 1000  # upper bound on pages followed in one walk; a safeguard, not a tuning knob
_DOWNLOAD_CHUNK_BYTES = 1024 * 1024

_ASSEMBLY_ACCESSION = re.compile(r"GC[FA]_[0-9]{9}\.[0-9]+")

# ``datasets download genome --include`` names -> ``include_annotation_type`` values.  The mapping was
# cross-checked against the ``cli_download_command_line`` that ``/download_summary`` echoes for each
# value, e.g. GENOME_FASTA + GENOME_GFF -> ``--include gff3,genome``.
_CLI_INCLUDE_NAMES: Dict[str, str] = {
    "genome": "GENOME_FASTA",
    "gff3": "GENOME_GFF",
    "gtf": "GENOME_GTF",
    "gbff": "GENOME_GBFF",
    "rna": "RNA_FASTA",
    "protein": "PROT_FASTA",
    "cds": "CDS_FASTA",
    "seq-report": "SEQUENCE_REPORT",
}
INCLUDE_FILE_TYPES = tuple(_CLI_INCLUDE_NAMES.values())
_INCLUDE_LOOKUP: Dict[str, str] = {
    **_CLI_INCLUDE_NAMES,
    **{api_name.lower().replace("_", "-"): api_name for api_name in INCLUDE_FILE_TYPES},
}

_REFERENCE_CATEGORY_RANK = {"reference genome": 2, "representative genome": 1}
_REFSEQ_SOURCE = "SOURCE_DATABASE_REFSEQ"


class NCBIDatasetsError(NCBIRequestError):
    """NCBI Datasets failed, reported an error, or returned a response that cannot be used."""


class AssemblyNotFoundError(LookupError):
    """NCBI Datasets answered, but has no report for the requested assembly."""


# --------------------------------------------------------------------------------------
# Requests
# --------------------------------------------------------------------------------------


class DatasetsRequest(NamedTuple):
    """A GET request: an already percent-encoded ``path`` plus its query ``params``."""

    path: str
    params: Mapping[str, str | List[str]]

    def with_page_token(self, token: Optional[str]) -> DatasetsRequest:
        """Return a copy that asks for the page after ``token`` (unchanged when ``token`` is empty)."""
        if not token:
            return self
        return DatasetsRequest(self.path, {**self.params, "page_token": token})

    def url(self, base_url: str = DATASETS_API_BASE_URL) -> str:
        """The full URL, including the encoded query string, exactly as ``requests`` will send it."""
        prepared = requests.Request("GET", base_url + self.path, params=dict(self.params)).prepare()
        return str(prepared.url)


def normalize_assembly_accession(accession: str) -> str:
    """Return ``accession`` if it is a versioned GenBank/RefSeq assembly accession, else raise ``ValueError``.

    The Datasets API rejects accessions without a version (``GCF_003254395`` is reported as invalid),
    so a version is required here.
    """
    if not isinstance(accession, str) or not _ASSEMBLY_ACCESSION.fullmatch(accession):
        raise ValueError(f"Invalid genome accession: {accession}")
    return accession


def normalize_taxon(taxon: str | int) -> str:
    """Return a single taxon name or NCBI Taxonomy ID, or raise ``ValueError``.

    The API treats commas in the path as list separators, so a name containing one would silently
    query several taxa.
    """
    if isinstance(taxon, bool) or not isinstance(taxon, (str, int)):
        raise ValueError(f"taxon must be an organism name or NCBI Taxonomy ID, got {type(taxon).__name__}")
    text = str(taxon).strip()
    if not text:
        raise ValueError("taxon must be a non-empty organism name or NCBI Taxonomy ID")
    if "," in text:
        raise ValueError(f"taxon must name a single organism (commas separate taxa in the API): {text!r}")
    if any(ord(char) < 32 or ord(char) == 127 for char in text):
        raise ValueError(f"taxon contains control characters: {text!r}")
    return text


def _check_page_size(page_size: int) -> int:
    if isinstance(page_size, bool) or not isinstance(page_size, int) or not 1 <= page_size <= MAX_PAGE_SIZE:
        raise ValueError(f"page_size must be an integer from 1 to {MAX_PAGE_SIZE}, got {page_size!r}")
    return page_size


def build_taxon_report_request(
    taxon: str | int,
    *,
    reference_only: bool = False,
    assembly_source: str = "all",
    exclude_paired_reports: bool = False,
    returned_content: str = "COMPLETE",
    page_size: int = DEFAULT_PAGE_SIZE,
    page_token: Optional[str] = None,
) -> DatasetsRequest:
    """Build ``GET /genome/taxon/{taxon}/dataset_report``.

    Args:
        taxon: Scientific or common name, or NCBI Taxonomy ID, at any rank.  Child taxa are included.
        reference_only: Limit to reference genome assemblies (``filters.reference_only``).
        assembly_source: ``"all"``, ``"refseq"`` (GCF_) or ``"genbank"`` (GCA_).
        exclude_paired_reports: For GenBank/RefSeq pairs return only the RefSeq copy.
        returned_content: ``"COMPLETE"`` reports, or abbreviated ``"ASSM_ACC"`` / ``"PAIRED_ACC"``.
        page_size: Reports per page, 1 to 1000 (server default 20).
        page_token: Token from the previous page's ``next_page_token``.
    """
    if assembly_source not in ASSEMBLY_SOURCES:
        raise ValueError(f"assembly_source must be one of {ASSEMBLY_SOURCES}, got {assembly_source!r}")
    if returned_content not in _RETURNED_CONTENT:
        raise ValueError(f"returned_content must be one of {_RETURNED_CONTENT}, got {returned_content!r}")
    params: Dict[str, str | List[str]] = {"page_size": str(_check_page_size(page_size))}
    if reference_only:
        params["filters.reference_only"] = "true"
    if assembly_source != "all":
        params["filters.assembly_source"] = assembly_source
    if exclude_paired_reports:
        params["filters.exclude_paired_reports"] = "true"
    if returned_content != "COMPLETE":
        params["returned_content"] = returned_content
    if page_token:
        params["page_token"] = page_token
    # Percent-encode everything (including "/") so an organism name cannot change the endpoint.
    path = f"/genome/taxon/{quote(normalize_taxon(taxon), safe='')}/dataset_report"
    return DatasetsRequest(path, params)


def build_accession_report_request(accession: str) -> DatasetsRequest:
    """Build ``GET /genome/accession/{accession}/dataset_report`` for one assembly accession."""
    return DatasetsRequest(f"/genome/accession/{normalize_assembly_accession(accession)}/dataset_report", {})


def build_sequence_report_request(accession: str, *, page_size: int = MAX_PAGE_SIZE) -> DatasetsRequest:
    """Build ``GET /genome/accession/{accession}/sequence_reports`` (one row per sequence of the assembly)."""
    path = f"/genome/accession/{normalize_assembly_accession(accession)}/sequence_reports"
    return DatasetsRequest(path, {"page_size": str(_check_page_size(page_size))})


def resolve_include_types(include: Sequence[str] | str | None) -> List[str]:
    """Translate file-type names to ``include_annotation_type`` values, in order and without duplicates.

    Accepts the ``datasets`` CLI names (``genome``, ``gff3``, ``gtf``, ``gbff``, ``rna``, ``protein``,
    ``cds``, ``seq-report``) and the API names (``GENOME_FASTA``, ``GENOME_GFF``, ...), either as a
    sequence or as one comma-separated string.  ``None`` or an empty selection returns ``[]``;
    :func:`build_package_request` supplies the genomic FASTA default explicitly.

    Raises:
        ValueError: If a name is not a supported file type.  The API itself answers HTTP 400 for an
            unknown ``include_annotation_type`` and silently ignores unknown parameter names.
    """
    if include is None:
        return []
    tokens = include.split(",") if isinstance(include, str) else list(include)
    resolved: List[str] = []
    for token in tokens:
        if not isinstance(token, str):
            raise ValueError(f"File type names must be strings, got {type(token).__name__}")
        key = token.strip().lower().replace("_", "-")
        if not key:
            continue
        api_name = _INCLUDE_LOOKUP.get(key)
        if api_name is None:
            supported = ", ".join(sorted(_CLI_INCLUDE_NAMES))
            raise ValueError(f"Unsupported genome package file type {token!r}; choose from: {supported}")
        if api_name not in resolved:
            resolved.append(api_name)
    return resolved


def build_package_request(accession: str, include: Sequence[str] | str | None = None) -> DatasetsRequest:
    """Build a fully hydrated package request, defaulting explicitly to genomic FASTA.

    The REST endpoint returns only metadata if file selection is omitted; the
    CLI's genomic FASTA default must therefore be supplied by this adapter.
    """
    params: Dict[str, str | List[str]] = {
        "include_annotation_type": resolve_include_types(include) or ["GENOME_FASTA"],
        "hydrated": "FULLY_HYDRATED",
    }
    return DatasetsRequest(f"/genome/accession/{normalize_assembly_accession(accession)}/download", params)


# --------------------------------------------------------------------------------------
# Responses
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class ReportPage:
    """One decoded page of a report endpoint (``v2reportsAssemblyDataReportPage`` / ``v2SequenceReportPage``)."""

    reports: Tuple[Mapping[str, Any], ...] = ()
    total_count: Optional[int] = None
    next_page_token: Optional[str] = None
    content_type: Optional[str] = None
    errors: Tuple[str, ...] = ()
    warnings: Tuple[str, ...] = ()


class ReportCollection(NamedTuple):
    """Reports gathered across pages, with the ``total_count`` NCBI announced (``None`` if it did not)."""

    reports: List[Mapping[str, Any]]
    total_count: Optional[int]


def require_complete(
    collected: ReportCollection, description: str, *, limit: Optional[int] = None
) -> List[Mapping[str, Any]]:
    """Return the collected reports, or raise if NCBI announced more of them than were received.

    ``total_count`` is sent on the first page (sequence reports) or on every page (assembly reports); a walk
    that ends with fewer reports than that, or than ``limit`` if that is smaller, lost some on the way and must
    not be returned as if it were the whole answer.
    """
    if collected.total_count is None:
        return collected.reports
    expected = collected.total_count if limit is None else min(limit, collected.total_count)
    if len(collected.reports) < expected:
        raise NCBIDatasetsError(
            f"NCBI Datasets announced {collected.total_count} {description} but {len(collected.reports)} were received"
        )
    return collected.reports


def _optional_int(value: object, field: str) -> Optional[int]:
    """An integer field; ``uint64`` values arrive as digit strings."""
    if value is None or value == "":
        return None
    if isinstance(value, bool):
        raise NCBIDatasetsError(f"Unexpected value for {field}: {value!r}")
    if isinstance(value, int):
        return value
    if isinstance(value, str) and value.isdigit():
        return int(value)
    raise NCBIDatasetsError(f"Unexpected value for {field}: {value!r}")


def _optional_float(value: object, field: str) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise NCBIDatasetsError(f"Unexpected value for {field}: {value!r}")
    return float(value)


def _optional_str(value: object, field: str) -> Optional[str]:
    if value is None:
        return None
    if not isinstance(value, str):
        raise NCBIDatasetsError(f"Unexpected value for {field}: {value!r}")
    return value


def _mapping(value: object, field: str) -> Mapping[str, Any]:
    """A nested object; an absent one is empty, a non-object is a schema violation."""
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise NCBIDatasetsError(f"Expected an object for {field}, got {type(value).__name__}")
    return value


def _message_text(entry: Mapping[str, Any]) -> str:
    """One ``messages[].error`` / ``messages[].warning`` object as a readable line."""
    codes = [str(value) for key, value in entry.items() if key.endswith(("_error_code", "_warning_code"))]
    details = [str(entry[key]) for key in ("reason", "message") if entry.get(key)]
    text = ": ".join(codes[:1] + details) or "unspecified message"
    identifiers = entry.get("invalid_identifiers") or []
    if entry.get("unrecognized_identifier"):
        identifiers = [*identifiers, entry["unrecognized_identifier"]]
    if identifiers:
        text += f" (identifiers: {', '.join(str(item) for item in identifiers)})"
    return text


def parse_report_page(payload: object) -> ReportPage:
    """Decode one page of a report endpoint.

    Missing ``reports`` / ``total_count`` / ``next_page_token`` keys mean "none", not an error: NCBI
    answers an unknown accession with ``{}`` and a taxon without matching assemblies with
    ``{"content_type": "ASSM_ACC"}``.  ``messages[].error`` entries are returned in ``errors`` for the
    caller to raise; they are never dropped.

    Raises:
        NCBIDatasetsError: If the payload does not have the documented shape.
    """
    if not isinstance(payload, Mapping):
        raise NCBIDatasetsError(f"Expected a JSON object from NCBI Datasets, got {type(payload).__name__}")

    raw_reports = payload.get("reports")
    if raw_reports is None:
        raw_reports = []
    if not isinstance(raw_reports, list) or not all(isinstance(item, Mapping) for item in raw_reports):
        raise NCBIDatasetsError("Expected 'reports' to be a list of objects")

    token = payload.get("next_page_token")
    if token is not None and not isinstance(token, str):
        raise NCBIDatasetsError(f"Unexpected value for next_page_token: {token!r}")

    errors: List[str] = []
    warnings: List[str] = []
    raw_messages = payload.get("messages")
    if raw_messages is not None:
        if not isinstance(raw_messages, list) or not all(isinstance(item, Mapping) for item in raw_messages):
            raise NCBIDatasetsError("Expected 'messages' to be a list of objects")
        for message in raw_messages:
            if "error" in message:
                errors.append(_message_text(_mapping(message["error"], "messages[].error")))
            if "warning" in message:
                warnings.append(_message_text(_mapping(message["warning"], "messages[].warning")))

    return ReportPage(
        reports=tuple(raw_reports),
        total_count=_optional_int(payload.get("total_count"), "total_count"),
        next_page_token=token or None,
        content_type=_optional_str(payload.get("content_type"), "content_type"),
        errors=tuple(errors),
        warnings=tuple(warnings),
    )


def continuation_token(page: ReportPage, seen_tokens: Collection[str] = ()) -> Optional[str]:
    """Return the ``page_token`` that fetches the page after ``page``, or ``None`` when paging is over.

    A page with no reports ends the walk even if it carries a token, and so does a token that has been
    used already: NCBI attaches ``next_page_token`` to a page that already holds the last result (a taxon
    with one assembly queried with ``page_size=1`` answers with a token whose page is
    ``{"total_count": 1}``), and a repeated token must not loop forever.
    """
    if not page.reports or not page.next_page_token or page.next_page_token in seen_tokens:
        return None
    return page.next_page_token


def parse_assembly_report(report: Mapping[str, Any]) -> Dict[str, Any]:
    """Flatten one ``v2reportsAssemblyDataReport`` into the metadata dictionary of this package.

    Fields NCBI does not report are ``None``; nothing is invented.  ``genome_size`` is converted from
    the ``uint64`` string NCBI sends, and the contig/scaffold counts come from ``number_of_contigs`` /
    ``number_of_scaffolds`` (there is no ``contig_count`` or ``scaffold_count`` in the schema).

    Raises:
        NCBIDatasetsError: If the report has no accession or a field has the wrong type.
    """
    accession = _optional_str(report.get("accession"), "accession")
    if not accession:
        raise NCBIDatasetsError("NCBI Datasets returned an assembly report without an accession")
    organism = _mapping(report.get("organism"), "organism")
    info = _mapping(report.get("assembly_info"), "assembly_info")
    stats = _mapping(report.get("assembly_stats"), "assembly_stats")
    return {
        "accession": accession,
        "organism_name": _optional_str(organism.get("organism_name") or organism.get("sci_name"), "organism_name"),
        "tax_id": _optional_int(organism.get("tax_id"), "tax_id"),
        "assembly_name": _optional_str(info.get("assembly_name"), "assembly_name"),
        "assembly_level": _optional_str(info.get("assembly_level"), "assembly_level"),
        "assembly_status": _optional_str(info.get("assembly_status"), "assembly_status"),
        "assembly_method": _optional_str(info.get("assembly_method"), "assembly_method"),
        "refseq_category": _optional_str(info.get("refseq_category"), "refseq_category"),
        "source_database": _optional_str(report.get("source_database"), "source_database"),
        "paired_accession": _optional_str(report.get("paired_accession"), "paired_accession"),
        "genome_size": _optional_int(stats.get("total_sequence_length"), "total_sequence_length"),
        "gc_percent": _optional_float(stats.get("gc_percent"), "gc_percent"),
        "chromosome_count": _optional_int(stats.get("total_number_of_chromosomes"), "total_number_of_chromosomes"),
        "contig_count": _optional_int(stats.get("number_of_contigs"), "number_of_contigs"),
        "contig_n50": _optional_int(stats.get("contig_n50"), "contig_n50"),
        "scaffold_count": _optional_int(stats.get("number_of_scaffolds"), "number_of_scaffolds"),
        "scaffold_n50": _optional_int(stats.get("scaffold_n50"), "scaffold_n50"),
        "release_date": _optional_str(info.get("release_date"), "release_date"),
        "submitter": _optional_str(info.get("submitter"), "submitter"),
    }


def select_reference_assembly(reports: Iterable[Mapping[str, Any]]) -> Optional[Mapping[str, Any]]:
    """Pick the best reference assembly among ``reports``, or ``None`` if there are none.

    Callers query with ``filters.reference_only``, so every report is a reference by NCBI's definition
    and the ranking only breaks ties: current assemblies before superseded ones, ``reference genome``
    before ``representative genome``, RefSeq before GenBank, then the newest release date and the
    highest accession for a deterministic result.  Reports without ``assembly_info`` (abbreviated
    ``ASSM_ACC`` reports) are ranked by source and accession alone.
    """

    def rank(report: Mapping[str, Any]) -> Tuple[int, int, int, str, str]:
        raw_info = report.get("assembly_info")
        info: Mapping[str, Any] = raw_info if isinstance(raw_info, Mapping) else {}
        status = info.get("assembly_status")
        category = str(info.get("refseq_category") or "").strip().lower()
        accession = str(report.get("accession") or "")
        is_refseq = report.get("source_database") == _REFSEQ_SOURCE or accession.startswith("GCF_")
        return (
            0 if status not in (None, "current") else 1,
            _REFERENCE_CATEGORY_RANK.get(category, 0),
            1 if is_refseq else 0,
            str(info.get("release_date") or ""),
            accession,
        )

    candidates = [report for report in reports if report.get("accession")]
    return max(candidates, key=rank) if candidates else None


def chromosome_lengths_from_sequence_reports(
    reports: Iterable[Mapping[str, Any]], *, names: str = "accession"
) -> Dict[str, int]:
    """Map sequence names to lengths in base pairs from ``sequence_reports`` rows.

    Args:
        reports: Rows of ``v2reportsSequenceInfo`` (all pages).
        names: ``"accession"`` keys every sequence by its RefSeq accession (GenBank accession when there
            is none).  ``"chromosome"`` keys assembled molecules by NCBI's chromosome name (``1``, ``X``,
            ``LG1``, ``MT``) and everything else (unlocalized or unplaced scaffolds, which share or lack a
            chromosome name) by accession; a chromosome name that is already taken also falls back to
            the accession.

    Raises:
        ValueError: If ``names`` is not one of the two modes.
        NCBIDatasetsError: If a row has no accession or length, or an accession occurs twice.
    """
    if names not in SEQUENCE_NAME_MODES:
        raise ValueError(f"names must be one of {SEQUENCE_NAME_MODES}, got {names!r}")
    lengths: Dict[str, int] = {}
    accessions_seen: Set[str] = set()
    for row in reports:
        accession = _optional_str(row.get("refseq_accession") or row.get("genbank_accession"), "sequence accession")
        if not accession:
            raise NCBIDatasetsError("NCBI Datasets returned a sequence report without an accession")
        if accession in accessions_seen:
            raise NCBIDatasetsError(f"NCBI Datasets returned sequence {accession} more than once")
        accessions_seen.add(accession)
        length = _optional_int(row.get("length"), f"length of {accession}")
        if length is None:
            raise NCBIDatasetsError(f"NCBI Datasets returned sequence {accession} without a length")
        key = accession
        chromosome = _optional_str(row.get("chr_name"), "chr_name")
        if names == "chromosome" and row.get("role") == "assembled-molecule" and chromosome and chromosome != "Un":
            if chromosome not in lengths and chromosome not in accessions_seen:
                key = chromosome
        lengths[key] = length
    return lengths


# --------------------------------------------------------------------------------------
# Package extraction
# --------------------------------------------------------------------------------------


def extract_package_zip(zip_path: str | Path, destination: str | Path) -> List[Path]:
    """Extract a genome data package into ``destination`` and return the files written.

    Every member is checked before anything is written: absolute paths, ``..`` components, drive
    letters, symlinks and paths that resolve outside ``destination`` are refused.  The package keeps
    the layout NCBI documents, e.g. ``ncbi_dataset/data/<accession>/<accession>_<name>_genomic.fna``.

    Raises:
        NCBIDatasetsError: If the file is not a zip archive or a member is unsafe.
    """
    zip_path = Path(zip_path)
    destination = Path(destination)
    try:
        archive = zipfile.ZipFile(zip_path)
    except zipfile.BadZipFile as exc:
        raise NCBIDatasetsError(
            f"NCBI Datasets did not return a complete zip archive ({zip_path.name}): {exc}"
        ) from None

    with archive:
        root = destination.resolve()
        plan: List[Tuple[zipfile.ZipInfo, Path]] = []
        for member in archive.infolist():
            name = member.filename.replace("\\", "/")
            relative = PurePosixPath(name)
            if relative.is_absolute() or ".." in relative.parts or (len(name) > 1 and name[1] == ":"):
                raise NCBIDatasetsError(f"Unsafe path in genome package: {member.filename!r}")
            if stat.S_ISLNK((member.external_attr >> 16) & 0o177777):
                raise NCBIDatasetsError(f"Symlink in genome package is not allowed: {member.filename!r}")
            target = (root / relative).resolve()
            if not target.is_relative_to(root):
                raise NCBIDatasetsError(f"Genome package member escapes the destination: {member.filename!r}")
            plan.append((member, target))

        destination.mkdir(parents=True, exist_ok=True)
        written: List[Path] = []
        for member, target in plan:
            if member.is_dir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(member) as source, target.open("wb") as sink:
                shutil.copyfileobj(source, sink, _DOWNLOAD_CHUNK_BYTES)
            written.append(target)
    return written


# --------------------------------------------------------------------------------------
# Client
# --------------------------------------------------------------------------------------


class NCBIDatasetsClient:
    """Session-based client for the NCBI Datasets v2 REST API.

    Args:
        api_key: NCBI API key, sent in the ``api-key`` header (never in the URL).  Without a key the
            service allows about 5 requests per second, with one 10.  Falls back to ``NCBI_API_KEY``.
        email: Optional contact, sent in the ``From`` header.  Falls back to ``NCBI_EMAIL``.  The
            Datasets API does not require one, and none is ever invented.
        timeout: Per-request timeout in seconds, for connecting and for each read.
    """

    def __init__(self, *, api_key: Optional[str] = None, email: Optional[str] = None, timeout: float = 60.0) -> None:
        if timeout <= 0:
            raise ValueError("timeout must be greater than zero")
        self.api_key = resolve_api_key(api_key)
        self.email = (email or os.environ.get("NCBI_EMAIL", "")).strip() or None
        self.timeout = timeout
        self.base_url = DATASETS_API_BASE_URL
        self.session = requests.Session()
        headers = {"User-Agent": f"metainformant/{__version__}"}
        if self.email:
            headers["From"] = self.email
        if self.api_key:
            headers["api-key"] = self.api_key
        self.session.headers.update(headers)

    def close(self) -> None:
        """Release the underlying connection pool."""
        self.session.close()

    def __enter__(self) -> NCBIDatasetsClient:
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.close()

    def get_json(self, request: DatasetsRequest) -> Dict[str, Any]:
        """GET a JSON endpoint and return the decoded object.

        Raises:
            NCBIDatasetsError: On a transport failure, an HTTP error status or a body that is not JSON.
        """
        return get_json_object(
            self.session,
            self.base_url + request.path,
            operation=f"NCBI Datasets GET {request.path}",
            params=request.params,
            timeout=self.timeout,
            headers={"Accept": "application/json"},
            error_class=NCBIDatasetsError,
        )

    def fetch_page(self, request: DatasetsRequest) -> ReportPage:
        """GET one page of a report endpoint; ``messages[].error`` entries raise instead of vanishing."""
        page = parse_report_page(self.get_json(request))
        if page.errors:
            raise NCBIDatasetsError(f"NCBI Datasets GET {request.path} reported: {'; '.join(page.errors)}")
        return page

    def iter_pages(self, request: DatasetsRequest) -> Iterator[ReportPage]:
        """Yield report pages, rejecting repeated assembly identities before yielding a page.

        GenBank and RefSeq paired accessions remain distinct identities. Sequence
        report rows use their own identifiers and are not checked as assemblies.
        """
        seen_tokens: Set[str] = set()
        seen_accessions: Set[str] = set()
        token: Optional[str] = None
        for _ in range(_MAX_PAGES):
            page = self.fetch_page(request.with_page_token(token))
            if request.path.endswith("/dataset_report"):
                for report in page.reports:
                    accession = report.get("accession")
                    if isinstance(accession, str) and accession:
                        if accession in seen_accessions:
                            raise NCBIDatasetsError(
                                f"NCBI Datasets GET {request.path} returned assembly accession "
                                f"{accession} more than once; reports are inconsistent"
                            )
                        seen_accessions.add(accession)
            yield page
            token = continuation_token(page, seen_tokens)
            if token is None:
                return
            seen_tokens.add(token)
        raise NCBIDatasetsError(f"NCBI Datasets GET {request.path} did not finish paging within {_MAX_PAGES} pages")

    def collect_reports(self, request: DatasetsRequest, *, limit: Optional[int] = None) -> ReportCollection:
        """Collect reports across pages, stopping as soon as ``limit`` reports are in hand."""
        if limit is not None and limit < 1:
            raise ValueError("limit must be at least 1")
        reports: List[Mapping[str, Any]] = []
        total_count: Optional[int] = None
        for index, page in enumerate(self.iter_pages(request)):
            if index == 0:
                total_count = page.total_count
            for report in page.reports:
                reports.append(report)
                if limit is not None and len(reports) >= limit:
                    return ReportCollection(reports, total_count)
        return ReportCollection(reports, total_count)

    def download_package(self, request: DatasetsRequest, destination: str | Path) -> Path:
        """Stream a package in 1 MiB chunks and atomically publish it on completion.

        Downloads use an owned temporary sibling of the destination. A failed
        request, interrupted stream or publication failure preserves any existing
        destination and removes only that owned partial file.

        Raises:
            NCBIDatasetsError: On a transport failure, an HTTP error status or an interrupted stream.
        """
        destination = Path(destination)
        operation = f"NCBI Datasets package download {request.path}"
        partial: Path | None = None
        try:
            with checked_get(
                self.session,
                self.base_url + request.path,
                operation=operation,
                params=request.params,
                timeout=self.timeout,
                stream=True,
                error_class=NCBIDatasetsError,
            ) as response:
                with tempfile.NamedTemporaryFile(
                    mode="wb",
                    prefix=f".{destination.name}.",
                    suffix=".part",
                    dir=destination.parent,
                    delete=False,
                ) as sink:
                    partial = Path(sink.name)
                    for chunk in response.iter_content(chunk_size=_DOWNLOAD_CHUNK_BYTES):
                        if chunk:
                            sink.write(chunk)
            partial.replace(destination)
            partial = None
        except NCBIDatasetsError:
            raise
        except requests.RequestException as exc:  # the connection broke while the body was streaming
            raise NCBIDatasetsError(f"{operation} failed: {describe_transport_error(exc)}") from None
        finally:
            if partial is not None:
                partial.unlink(missing_ok=True)
        return destination


@contextmanager
def client_scope(
    client: Optional[NCBIDatasetsClient] = None, *, email: Optional[str] = None
) -> Iterator[NCBIDatasetsClient]:
    """Yield ``client``, or a temporary client (sending ``email`` as contact) that is closed on exit.

    A client passed in by the caller is never closed here.
    """
    if client is not None:
        yield client
        return
    owned = NCBIDatasetsClient(email=email)
    try:
        yield owned
    finally:
        owned.close()
