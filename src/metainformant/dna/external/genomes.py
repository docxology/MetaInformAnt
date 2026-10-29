"""Genome data retrieval and processing utilities.

This module provides tools for downloading genome assemblies, validating
accessions, and retrieving genome metadata from the NCBI Datasets v2 REST API.
The request construction and response parsing live in :mod:`.ncbi_datasets`
as pure functions; every function here shares one client, and a failed request
raises instead of returning an empty result.
"""

from __future__ import annotations

import os
import re
from pathlib import Path, PurePosixPath
from typing import Any, Dict, List, Optional, Sequence
from urllib.parse import urlparse

from metainformant.core.ncbi import resolve_ncbi_contact
from metainformant.core.utils import logging

from .ncbi_datasets import (
    MAX_PAGE_SIZE,
    SEQUENCE_NAME_MODES,
    AssemblyNotFoundError,
    NCBIDatasetsClient,
    NCBIDatasetsError,
    build_accession_report_request,
    build_package_request,
    build_sequence_report_request,
    build_taxon_report_request,
    chromosome_lengths_from_sequence_reports,
    client_scope,
    extract_package_zip,
    normalize_assembly_accession,
    parse_assembly_report,
    require_complete,
    resolve_include_types,
    select_reference_assembly,
)

logger = logging.get_logger(__name__)

# How many reference candidates ``download_reference_genome`` examines.  A species has a handful at
# most; a higher rank asked for by mistake would otherwise page through thousands of reports.
_REFERENCE_CANDIDATES = 100

_FASTA_SUFFIXES = (".fna", ".fa", ".fasta")
_GFF_SUFFIXES = (".gff", ".gff3")


def download_genome_package(
    accession: str,
    output_dir: str | Path,
    include: Sequence[str] | str | None = None,
    *,
    client: Optional[NCBIDatasetsClient] = None,
) -> Path:
    """Download a genome data package from the NCBI Datasets v2 API and extract it.

    The zip is streamed to ``<output_dir>/<accession>.zip`` (it is never held in memory), extracted to
    ``<output_dir>/<accession>/`` and deleted.  The extracted files keep the layout NCBI documents:
    ``ncbi_dataset/data/assembly_data_report.jsonl`` and, per assembly,
    ``ncbi_dataset/data/<accession>/<accession>_<assembly name>_genomic.fna`` plus any requested annotation.

    Args:
        accession: Versioned assembly accession (e.g., "GCF_003254395.2" for Apis mellifera Amel_HAv3.1).
        output_dir: Output directory for downloaded files
        include: File types to download, as ``datasets`` CLI names (``genome``, ``gff3``, ``gtf``,
            ``gbff``, ``rna``, ``protein``, ``cds``, ``seq-report``) or API names (``GENOME_FASTA``, ...).
            The assembly data report is always included.  When omitted or empty, genomic FASTA is
            explicitly requested in a fully hydrated package. File types are sent as repeated
            ``include_annotation_type`` parameters; the ``include`` parameter this
            function used to send is not part of the v2 API and was silently ignored by the server.
        client: A configured :class:`NCBIDatasetsClient` (API key, timeout); one is created if omitted.

    Returns:
        Path to the directory the package was extracted to.

    Raises:
        ValueError: If the accession or a file type is invalid.
        requests.RequestException: If the download fails (:class:`NCBIDatasetsError`, a subclass).

    Example:
        >>> # This would download the human genome (large file!)
        >>> # package_dir = download_genome_package("GCF_000001405.40", "genomes/")
        >>> # isinstance(package_dir, Path)
        >>> # True
    """
    accession = normalize_assembly_accession(accession)
    request = build_package_request(accession, include)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    zip_path = output_dir / f"{accession}.zip"

    logger.info("Downloading genome %s from NCBI Datasets", accession)
    with client_scope(client) as active:
        active.download_package(request, zip_path)
    logger.info("Downloaded genome package to %s", zip_path)

    extract_dir = output_dir / accession
    try:
        extracted = extract_package_zip(zip_path, extract_dir)
    finally:
        zip_path.unlink(missing_ok=True)
    if not extracted:
        raise NCBIDatasetsError(f"NCBI Datasets returned a package without files for {accession}")
    return extract_dir


def download_genome_package_best_effort(
    accession: str,
    output_dir: str | Path,
    include: Sequence[str] | str | None = None,
    ftp_url: Optional[str] = None,
    *,
    client: Optional[NCBIDatasetsClient] = None,
) -> Path:
    """Download genome package with a fallback strategy.

    This function tries, in order:
    1. The NCBI Datasets API (:func:`download_genome_package`)
    2. A direct FTP download, only if ``ftp_url`` is provided

    Args:
        accession: Genome accession
        output_dir: Output directory
        include: File types to include (see :func:`download_genome_package`)
        ftp_url: FTP URL for fallback download
        client: Optional configured :class:`NCBIDatasetsClient`

    Returns:
        Path to genome package directory

    Raises:
        RuntimeError: If every strategy failed; the message lists each failure.

    Example:
        >>> # Try to download with fallbacks
        >>> # package_dir = download_genome_package_best_effort("GCF_000001405.40", "genomes/")
    """
    output_dir = Path(output_dir)

    failures: list[str] = []

    try:
        # Try primary method first
        return download_genome_package(accession, output_dir, include, client=client)
    except Exception as e:
        failures.append(f"NCBI Datasets API: {e}")
        logger.warning(f"Primary download failed: {e}")

        # Try FTP fallback if provided
        if ftp_url:
            try:
                return _download_from_ftp(ftp_url, output_dir / accession)
            except Exception as e2:
                failures.append(f"FTP fallback: {e2}")
                logger.warning(f"FTP download failed: {e2}")

        detail = "; ".join(failures) or "no download strategy was attempted"
        raise RuntimeError(f"All genome download methods failed for {accession}: {detail}") from e


def is_valid_assembly_accession(accession: str) -> bool:
    """Validate NCBI assembly accession format.

    Args:
        accession: Assembly accession string

    Returns:
        True if accession has valid NCBI assembly format

    Example:
        >>> is_valid_assembly_accession("GCF_000001405.39")
        True
        >>> is_valid_assembly_accession("GCA_000001405")
        True
        >>> is_valid_assembly_accession("invalid")
        False
    """
    # NCBI assembly accession patterns
    # GCF_##########.## or GCA_##########
    pattern = r"^(GCF|GCA)_[0-9]{9,}(\.[0-9]+)?$"

    return bool(re.match(pattern, accession))


def validate_accession(accession: str) -> bool:
    """Validate genome accession format.

    Args:
        accession: Genome accession to validate

    Returns:
        True if accession format is valid

    Example:
        >>> validate_accession("GCF_000001405.39")
        True
        >>> validate_accession("invalid")
        False
    """
    if not accession:
        return False

    # Common NCBI accession patterns
    patterns = [
        r"^GC[FA]_\d{9}\.\d+$",  # GCF_000001405.39 (RefSeq)
        r"^GCA_\d{9}\.\d+$",  # GCA_000001405.39 (GenBank)
        r"^ASM\d+$",  # ASM2732 (older format)
        r"^chr\d+$",  # chr1, chr2, etc.
        r"^NC_\d{6}$",  # NC_000001 (chromosome)
        r"^NT_\d{6}$",  # NT_000001 (contig)
        r"^NW_\d{8}$",  # NW_001838827 (scaffold)
    ]

    return any(re.match(pattern, accession) for pattern in patterns)


def get_genome_metadata(accession: str, *, client: Optional[NCBIDatasetsClient] = None) -> Dict[str, Any]:
    """Retrieve genome assembly metadata from the NCBI Datasets assembly report.

    Args:
        accession: Versioned assembly accession (e.g., "GCF_003254395.2")
        client: Optional configured :class:`NCBIDatasetsClient`

    Returns:
        Dictionary with ``accession``, ``organism_name``, ``tax_id``, ``assembly_name``,
        ``assembly_level``, ``assembly_status``, ``assembly_method``, ``refseq_category``,
        ``source_database``, ``paired_accession``, ``genome_size`` (bp, an ``int``), ``gc_percent``,
        ``chromosome_count``, ``contig_count``, ``contig_n50``, ``scaffold_count``, ``scaffold_n50``,
        ``release_date`` and ``submitter``.  A field NCBI does not report is ``None``; nothing is invented.

    Raises:
        ValueError: If the accession is not a versioned GenBank/RefSeq assembly accession.
        AssemblyNotFoundError: If NCBI has no report for the accession.
        requests.RequestException: If the request fails (:class:`NCBIDatasetsError`, a subclass).
            This function used to return ``{"accession": ..., "error": ...}`` instead of raising.

    Example:
        >>> # This would query NCBI for metadata
        >>> # metadata = get_genome_metadata("GCF_003254395.2")
        >>> # metadata["organism_name"]
        >>> # 'Apis mellifera'
    """
    accession = normalize_assembly_accession(accession)
    request = build_accession_report_request(accession)
    with client_scope(client) as active:
        page = active.fetch_page(request)

    for report in page.reports:
        if report.get("accession") == accession:
            return parse_assembly_report(report)
    if page.reports:
        returned = ", ".join(str(report.get("accession")) for report in page.reports)
        raise NCBIDatasetsError(f"NCBI Datasets returned a report for {returned}, not for {accession}")
    raise AssemblyNotFoundError(f"NCBI Datasets has no assembly report for {accession}")


def list_genome_assemblies(
    organism: str | int,
    max_results: int = 10,
    *,
    reference_only: bool = False,
    assembly_source: str = "refseq",
    client: Optional[NCBIDatasetsClient] = None,
) -> List[Dict[str, Any]]:
    """Search for genome assemblies of an organism or taxon.

    Pages through ``GET /genome/taxon/{taxon}/dataset_report`` until ``max_results`` assemblies are
    collected or NCBI has no more.  Assemblies of child taxa are included, and they come in the order
    NCBI returns them, which is not sorted by release date.

    Args:
        organism: Scientific or common name (e.g., "Homo sapiens") or NCBI Taxonomy ID at any rank.
        max_results: Maximum number of assemblies to return (at least 1).
        reference_only: Only reference genome assemblies.
        assembly_source: ``"refseq"`` (GCF_, the default), ``"genbank"`` (GCA_) or ``"all"``.  Species
            that have only GenBank assemblies need ``"genbank"`` or ``"all"``.
        client: Optional configured :class:`NCBIDatasetsClient`

    Returns:
        A list of metadata dictionaries in the format of :func:`get_genome_metadata`.  An empty list
        means NCBI has no matching assembly; an unknown organism name looks the same, because NCBI
        answers both with an empty report.

    Raises:
        ValueError: If an argument is invalid.
        requests.RequestException: If the request fails (:class:`NCBIDatasetsError`, a subclass).
            This function used to return ``[]`` on failure.

    Example:
        >>> # Search for human genome assemblies
        >>> # assemblies = list_genome_assemblies("Homo sapiens", max_results=5)
        >>> # len(assemblies) <= 5
        >>> # True
    """
    if isinstance(max_results, bool) or not isinstance(max_results, int) or max_results < 1:
        raise ValueError(f"max_results must be a positive integer, got {max_results!r}")
    request = build_taxon_report_request(
        organism,
        reference_only=reference_only,
        assembly_source=assembly_source,
        page_size=min(max_results, MAX_PAGE_SIZE),
    )
    with client_scope(client) as active:
        collected = active.collect_reports(request, limit=max_results)
    reports = require_complete(collected, f"assemblies for {organism}", limit=max_results)
    return [parse_assembly_report(report) for report in reports]


def download_reference_genome(
    species: str | int,
    output_dir: str | Path,
    *,
    include: Sequence[str] | str | None = None,
    client: Optional[NCBIDatasetsClient] = None,
) -> Optional[Path]:
    """Download the reference genome of a species.

    Asks NCBI for the reference genome assemblies of the taxon (``filters.reference_only``, RefSeq copy
    preferred) and downloads the best one, ranked as in :func:`select_reference_assembly`.  The first
    assembly NCBI happens to list is not assumed to be the reference.

    Args:
        species: Species name or NCBI Taxonomy ID
        output_dir: Output directory
        include: File types to download (see :func:`download_genome_package`)
        client: Optional configured :class:`NCBIDatasetsClient`

    Returns:
        Path to the extracted genome package, or ``None`` if NCBI lists no reference genome for the
        taxon.  Many species have none, e.g. most viruses.

    Raises:
        ValueError: If an argument is invalid.
        RuntimeError: If a reference genome was found but could not be downloaded.
        requests.RequestException: If the lookup fails (:class:`NCBIDatasetsError`, a subclass).
            Lookup and download failures used to be logged and returned as ``None``.

    Example:
        >>> # Download the honey bee reference genome
        >>> # genome_path = download_reference_genome("Apis mellifera", "genomes/")
    """
    resolve_include_types(include)
    request = build_taxon_report_request(
        species,
        reference_only=True,
        exclude_paired_reports=True,
        page_size=_REFERENCE_CANDIDATES,
    )
    with client_scope(client) as active:
        collected = active.collect_reports(request, limit=_REFERENCE_CANDIDATES)
        chosen = select_reference_assembly(collected.reports)
        if chosen is None:
            logger.warning("NCBI lists no reference genome assembly for %s", species)
            return None
        if collected.total_count is not None and collected.total_count > len(collected.reports):
            logger.warning(
                "NCBI lists %d reference assemblies for %s; choosing among the first %d",
                collected.total_count,
                species,
                len(collected.reports),
            )
        accession = str(chosen["accession"])
        logger.info("Downloading reference genome %s for %s", accession, species)
        return download_genome_package_best_effort(accession, output_dir, include, client=active)


def _download_from_ftp(ftp_url: str, output_dir: Path) -> Path:
    """Download the file at ``ftp_url`` into ``output_dir``.

    Raises:
        ValueError: If the URL names a directory (empty path or trailing ``/``) instead of a file.  The
            download used to be saved as an invented ``genome_data.fasta``, or under the name of the last
            directory, which turned an FTP listing into a "genome".
    """
    parsed = urlparse(ftp_url)
    if not parsed.path or parsed.path.endswith("/"):
        raise ValueError(f"FTP URL must name a file, not a directory: {ftp_url}")
    filename = PurePosixPath(parsed.path).name

    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / filename

    logger.info(f"Downloading from FTP: {ftp_url}")
    from metainformant.core.io.download import download_with_progress

    # Use the centralized downloader so we get heartbeat + retry.
    result = download_with_progress(
        ftp_url,
        output_file,
        protocol="ftp",
        show_progress=False,
        heartbeat_interval=5,
        max_retries=3,
        chunk_size=1024 * 1024,
        timeout=300,
        resume=False,
    )
    if not result.success:
        raise RuntimeError(f"FTP download failed: {ftp_url}: {result.error or 'unknown error'}")

    return output_dir


def get_chromosome_lengths(
    accession: str,
    email: str | None = None,
    *,
    allow_anonymous: bool = False,
    names: str = "accession",
    client: Optional[NCBIDatasetsClient] = None,
) -> Dict[str, int]:
    """Get sequence lengths for a genome assembly from the NCBI Datasets sequence report.

    Reads every sequence of the assembly (chromosomes, organelles, unlocalized and unplaced scaffolds),
    following pagination, so nothing is truncated.  This replaces a nucleotide search through Entrez that
    stopped at 100 records and keyed them by a ``chr<N>`` name parsed from their titles, so records whose
    titles parse to the same name overwrote one another.

    Args:
        accession: Versioned assembly accession (e.g., 'GCF_003254395.2' for Apis mellifera)
        email: Explicit contact email, sent as the ``From`` header; otherwise ``NCBI_EMAIL``.
        allow_anonymous: Explicitly permit a request without any contact.
        names: ``"accession"`` (default) keys each sequence by its RefSeq accession, which is what the FASTA
            headers of an NCBI genome package carry.  ``"chromosome"`` keys assembled molecules by NCBI's
            chromosome name (``1``, ``X``, ``LG1``, ``MT``) and all other sequences by accession.
        client: Optional configured :class:`NCBIDatasetsClient` (used as is; ``email`` then only
            enforces the contact policy).

    Returns:
        Dictionary mapping sequence names to lengths in base pairs.

    Raises:
        ValueError: If an argument is invalid.
        NCBIContactError: If no contact is available and ``allow_anonymous`` is false.
        AssemblyNotFoundError: If NCBI lists no sequences for the accession.
        requests.RequestException: If the request fails (:class:`NCBIDatasetsError`, a subclass).
            Failures, a missing assembly and a missing Biopython used to return ``{}``.

    Example:
        >>> # Get sequence lengths for the honey bee genome (requires network)
        >>> # lengths = get_chromosome_lengths("GCF_003254395.2", names="chromosome")
        >>> # lengths["LG1"]
        >>> # 27754200
    """
    accession = normalize_assembly_accession(accession)
    if names not in SEQUENCE_NAME_MODES:
        raise ValueError(f"names must be one of {SEQUENCE_NAME_MODES}, got {names!r}")
    contact = resolve_ncbi_contact(email, allow_anonymous=allow_anonymous)
    logger.info("NCBI contact mode: %s", contact.mode)
    logger.info(f"Fetching sequence lengths for: {accession}")

    request = build_sequence_report_request(accession)
    with client_scope(client, email=contact.email) as active:
        collected = active.collect_reports(request)
    reports = require_complete(collected, f"sequences for {accession}")
    if not reports:
        raise AssemblyNotFoundError(f"NCBI Datasets lists no sequences for {accession}")

    lengths = chromosome_lengths_from_sequence_reports(reports, names=names)
    logger.info(f"Found {len(lengths)} sequence lengths")
    return lengths


def validate_genome_files(genome_dir: str | Path) -> Dict[str, Any]:
    """Validate genome assembly files.

    Searches ``genome_dir`` recursively, so a directory returned by :func:`download_genome_package`
    (files under ``ncbi_dataset/data/<accession>/``) validates; FASTA is recognised by ``.fna``, ``.fa`` or
    ``.fasta`` and GFF by ``.gff`` or ``.gff3``, each optionally gzip-compressed.

    Args:
        genome_dir: Directory containing genome files

    Returns:
        Validation results dictionary with ``valid``, ``fasta_files``, ``gff_files``, ``total_files`` (files
        at any depth) and ``issues``; ``{"valid": False, "error": ...}`` if ``genome_dir`` is not a directory.

    Example:
        >>> # Validate downloaded genome
        >>> # results = validate_genome_files("genomes/GCF_003254395.2/")
        >>> # results['valid']
        >>> # True
    """
    genome_dir = Path(genome_dir)

    if not genome_dir.exists():
        return {"valid": False, "error": "Directory does not exist"}
    if not genome_dir.is_dir():
        return {"valid": False, "error": "Path is not a directory"}

    files = [Path(root) / name for root, _dirs, names in os.walk(genome_dir) for name in names]

    def count(suffixes: tuple[str, ...]) -> int:
        stems = (path.name.lower().removesuffix(".gz") for path in files)
        return sum(1 for stem in stems if stem.endswith(suffixes))

    fasta_files = count(_FASTA_SUFFIXES)
    results: dict[str, Any] = {
        "valid": True,
        "fasta_files": fasta_files,
        "gff_files": count(_GFF_SUFFIXES),
        "total_files": len(files),
        "issues": [],
    }

    if not fasta_files:
        results["issues"].append("No FASTA files found")
        results["valid"] = False

    return results
