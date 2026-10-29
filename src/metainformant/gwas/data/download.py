"""Data download utilities for GWAS analysis.

This module provides tools for downloading reference genomes, variant databases,
and SRA sequencing data for GWAS analysis.

The GWAS Catalog association table (:func:`download_variant_data`) is split into a thin HTTP layer and
pure code that takes response text and returns values (:func:`decode_catalog_response`,
:func:`parse_catalog_study`, :func:`parse_catalog_associations` and the sans-IO :class:`AssociationsPager`),
so the response handling can be tested without a network.
"""

from __future__ import annotations

import csv
import json
import math
import re
import stat
import zipfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import IO, Any, BinaryIO, Dict, Iterable, List, Mapping, Optional, Tuple

import requests

from metainformant.core.ncbi import resolve_ncbi_contact
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def download_reference_genome(accession: str, output_dir: str | Path) -> Path:
    """Download reference genome from NCBI.

    Args:
        accession: Genome accession (e.g., "GCF_000001405.39" for GRCh38)
        output_dir: Output directory for downloaded genome

    Returns:
        Path to downloaded genome directory

    Raises:
        ValueError: If accession is invalid
        requests.RequestException: If download fails
        RuntimeError: If both the Datasets API and FTP fallback fail

    Example:
        >>> # Download human reference genome
        >>> genome_dir = download_reference_genome("GCF_000001405.39", "genomes/")
        >>> genome_dir.exists()
        True
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Validate accession format
    if not _validate_genome_accession(accession):
        raise ValueError(f"Invalid genome accession: {accession}")

    logger.info(f"Downloading reference genome {accession}")

    failures: list[str] = []

    # Try NCBI Datasets API first
    try:
        return _download_from_ncbi_datasets(accession, output_dir)
    except Exception as e:
        failures.append(f"NCBI Datasets API: {e}")
        logger.warning(f"NCBI Datasets download failed: {e}")

    # Fallback to FTP download
    try:
        return _download_from_ftp(accession, output_dir)
    except Exception as e:
        failures.append(f"NCBI FTP: {e}")
        logger.warning(f"FTP download failed: {e}")

    detail = "; ".join(failures) or "no download strategy was attempted"
    raise RuntimeError(f"All genome download methods failed for {accession}: {detail}") from None


def _download_from_ncbi_datasets(accession: str, output_dir: Path) -> Path:
    """Download genome using NCBI Datasets API."""
    api_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download"

    params = {"include": "genome,annotation,gff3,protein,cds,seq-report"}

    response = requests.get(api_url, params=params, timeout=60)
    response.raise_for_status()

    # Save as zip file
    zip_path = output_dir / f"{accession}.zip"
    with open(zip_path, "wb") as f:
        f.write(response.content)

    # Extract zip
    import glob

    extract_dir = output_dir / accession
    extract_dir.mkdir(exist_ok=True)

    _extract_zip_safely(zip_path, extract_dir)

    # Clean up zip
    zip_path.unlink()

    # Validate that we actually got a FASTA file
    if not glob.glob(str(extract_dir / "**" / "*.fna*"), recursive=True):
        import shutil

        shutil.rmtree(extract_dir)
        raise RuntimeError(f"NCBI Datasets API returned a zip without any .fna files for {accession}")

    return extract_dir


def _extract_zip_safely(zip_path: Path, extract_dir: Path) -> None:
    """Extract a downloaded ZIP without permitting path or symlink escapes.

    NCBI responses are external input.  ``ZipFile.extractall`` trusts member
    names and can write outside the requested output directory when a remote
    archive contains traversal entries.  Extraction is intentionally explicit
    so every member is checked before it is written.
    """
    destination = extract_dir.resolve()
    destination.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(zip_path, "r") as archive:
        for member in archive.infolist():
            member_name = member.filename.replace("\\", "/")
            member_path = Path(member_name)
            if member_path.is_absolute() or ".." in member_path.parts:
                raise ValueError(f"Unsafe ZIP member path: {member.filename!r}")

            target = (destination / member_path).resolve()
            if not target.is_relative_to(destination):
                raise ValueError(f"ZIP member escapes extraction directory: {member.filename!r}")

            mode = (member.external_attr >> 16) & 0o177777
            if stat.S_ISLNK(mode):
                raise ValueError(f"Symlink ZIP members are not allowed: {member.filename!r}")

            if member.is_dir() or member_name.endswith("/"):
                target.mkdir(parents=True, exist_ok=True)
                continue

            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(member, "r") as source, target.open("wb") as sink:
                _copy_zip_member(source, sink)


def _copy_zip_member(source: IO[bytes], sink: BinaryIO, chunk_size: int = 1024 * 1024) -> None:
    """Copy one validated ZIP member to its confined destination."""
    while chunk := source.read(chunk_size):
        sink.write(chunk)


def _download_from_ftp(accession: str, output_dir: Path) -> Path:
    """Download genome from NCBI FTP using ftplib.

    Constructs the NCBI FTP directory path from the accession, lists the
    remote directory to find the assembly folder, and downloads the genomic
    FASTA (.fna.gz) and GFF3 annotation (.gff.gz) files.

    Args:
        accession: Genome accession (e.g., "GCF_003254395.2")
        output_dir: Parent output directory

    Returns:
        Path to the directory containing downloaded files

    Raises:
        RuntimeError: If the accession directory is not found on the FTP server
    """
    import ftplib  # nosec B402 - retained as a documented NCBI legacy fallback; primary path is HTTPS.

    # Extract the 9-digit numeric portion after the prefix (e.g. "003254395")
    prefix = accession[:3]  # "GCF" or "GCA"
    digits = accession.split("_")[1].split(".")[0]  # "003254395"
    ftp_dir = f"/genomes/all/{prefix}/{digits[:3]}/{digits[3:6]}/{digits[6:9]}"

    genome_dir = output_dir / accession
    genome_dir.mkdir(exist_ok=True)

    ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")  # nosec B321 - NCBI fallback endpoint.
    ftp.login()

    try:
        ftp.cwd(ftp_dir)

        # The assembly directory includes a name suffix, e.g.
        # GCF_003254395.2_Amel_HAv3.1 — find it by prefix match
        dirs = ftp.nlst()
        target_dir = next((d for d in dirs if d.startswith(accession)), None)
        if not target_dir:
            raise RuntimeError(f"No directory found for {accession} on NCBI FTP at {ftp_dir}")

        ftp.cwd(target_dir)
        files = ftp.nlst()

        # Download genomic FASTA and GFF3 annotation
        found_genome = False
        for suffix in ("_genomic.fna.gz", "_genomic.gff.gz"):
            matching = [f for f in files if f.endswith(suffix) and "_from_genomic" not in f and "_rna_from" not in f]
            if matching:
                filename = matching[0]
                if suffix == "_genomic.fna.gz":
                    found_genome = True
                local_path = genome_dir / filename
                if local_path.exists():
                    logger.info(f"File already exists: {local_path}")
                    continue
                logger.info(f"Downloading {filename} from NCBI FTP...")
                with open(local_path, "wb") as out_file:
                    ftp.retrbinary(f"RETR {filename}", out_file.write)
                logger.info(f"Downloaded {filename}")
    finally:
        ftp.quit()

    if not found_genome and not any(genome_dir.glob("*_genomic.fna*")):
        raise RuntimeError(f"NCBI FTP returned no genomic FASTA for {accession}")

    return genome_dir


def _validate_genome_accession(accession: str) -> bool:
    """Validate genome accession format."""
    import re

    # NCBI RefSeq or GenBank patterns
    patterns = [
        r"^GC[FA]_\d{9}\.\d+$",  # GCF_000001405.39
        r"^GCA_\d{9}\.\d+$",  # GCA_000001405.39
    ]

    return any(re.match(pattern, accession) for pattern in patterns)


def download_variant_database(db_name: str, output_dir: str | Path) -> Path:
    """Download variant database (dbSNP, 1000 Genomes, etc.).

    Args:
        db_name: Database name ("dbsnp", "1000genomes", "hapmap")
        output_dir: Output directory

    Returns:
        Path to downloaded database

    Example:
        >>> # Download dbSNP for human
        >>> db_path = download_variant_database("dbsnp", "databases/")
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if db_name == "dbsnp":
        return _download_dbsnp(output_dir)
    elif db_name == "1000genomes":
        return _download_1000genomes(output_dir)
    else:
        raise ValueError(f"Unknown database: {db_name}")


def _download_dbsnp(output_dir: Path, chromosome: str | None = None) -> Path:
    """Download dbSNP VCF files from NCBI FTP.

    Args:
        output_dir: Output directory
        chromosome: Specific chromosome to download (e.g., '1', '22', 'X'), or None for all

    Returns:
        Path to downloaded dbSNP directory

    Note:
        Downloads from NCBI FTP. Full database is very large (~40GB compressed).
        Consider downloading specific chromosomes instead.
    """
    import ftplib  # nosec B402 - retained as a documented NCBI legacy fallback; primary path is HTTPS.

    ftp_host = "ftp.ncbi.nih.gov"
    ftp_path = "/snp/latest_release/VCF"

    db_dir = output_dir / "dbsnp"
    db_dir.mkdir(exist_ok=True)

    logger.info(f"Connecting to {ftp_host} for dbSNP download...")

    try:
        ftp = ftplib.FTP(ftp_host)  # nosec B321 - NCBI fallback endpoint.
        ftp.login()
        ftp.cwd(ftp_path)

        # List available files
        files = ftp.nlst()

        # Filter to VCF files
        vcf_files = [f for f in files if f.endswith(".vcf.gz")]

        if chromosome:
            # Download specific chromosome
            vcf_files = [f for f in vcf_files if chromosome in f]

        if not vcf_files:
            logger.warning("No dbSNP VCF files found matching criteria")
            ftp.quit()
            return db_dir

        # Download first matching file (or provide info about available files)
        for vcf_file in vcf_files[:1]:  # Limit to first file for safety
            local_path = db_dir / vcf_file
            if local_path.exists():
                logger.info(f"File already exists: {local_path}")
                continue

            logger.info(f"Downloading {vcf_file}...")
            with open(local_path, "wb") as f:
                ftp.retrbinary(f"RETR {vcf_file}", f.write)
            logger.info(f"Downloaded: {local_path}")

        ftp.quit()
        return db_dir

    except ftplib.all_errors as e:
        logger.error(f"FTP download failed: {e}")
        raise RuntimeError(f"dbSNP download failed: {e}")


def _download_1000genomes(output_dir: Path, chromosome: str | None = None) -> Path:
    """Download 1000 Genomes Project VCF data.

    Args:
        output_dir: Output directory
        chromosome: Specific chromosome to download (e.g., '1', '22'), or None for index

    The release VCFs hold every sample, so there is no per-population download: the sample panel
    (``integrated_call_samples_v3.20130502.ALL.panel``) maps samples to populations. (A ``population``
    argument was accepted and ignored.)

    Returns:
        Path to downloaded 1000 Genomes directory

    Note:
        Data from Phase 3 release (2013). For newer data, consider gnomAD.
    """
    import ftplib  # nosec B402 - retained as a documented NCBI legacy fallback; primary path is HTTPS.

    ftp_host = "ftp.1000genomes.ebi.ac.uk"
    ftp_path = "/vol1/ftp/release/20130502"

    db_dir = output_dir / "1000genomes"
    db_dir.mkdir(exist_ok=True)

    logger.info(f"Connecting to {ftp_host} for 1000 Genomes download...")

    try:
        ftp = ftplib.FTP(ftp_host)  # nosec B321 - legacy 1000 Genomes endpoint.
        ftp.login()
        ftp.cwd(ftp_path)

        # List available files
        files = ftp.nlst()

        # Filter to VCF files
        vcf_files = [f for f in files if f.endswith(".vcf.gz")]

        if chromosome:
            # Download specific chromosome
            vcf_files = [f for f in vcf_files if f"chr{chromosome}." in f or f".{chromosome}." in f]

        if not vcf_files:
            logger.warning("No 1000 Genomes VCF files found matching criteria")
            # Download the panel file instead for population info
            panel_file = "integrated_call_samples_v3.20130502.ALL.panel"
            if panel_file in files:
                local_path = db_dir / panel_file
                if not local_path.exists():
                    logger.info(f"Downloading sample panel: {panel_file}")
                    with open(local_path, "wb") as f:
                        ftp.retrbinary(f"RETR {panel_file}", f.write)
            ftp.quit()
            return db_dir

        # Download first matching file (or provide info)
        for vcf_file in vcf_files[:1]:  # Limit for safety
            local_path = db_dir / vcf_file
            if local_path.exists():
                logger.info(f"File already exists: {local_path}")
                continue

            logger.info(f"Downloading {vcf_file} (this may take a while)...")
            with open(local_path, "wb") as f:
                ftp.retrbinary(f"RETR {vcf_file}", f.write)
            logger.info(f"Downloaded: {local_path}")

        ftp.quit()
        return db_dir

    except ftplib.all_errors as e:
        logger.error(f"FTP download failed: {e}")
        raise RuntimeError(f"1000 Genomes download failed: {e}")


def download_sra_run(sra_accession: str, output_dir: str | Path, threads: int = 1) -> Path:
    """Download SRA run data.

    Args:
        sra_accession: SRA accession (e.g., "SRR123456")
        output_dir: Output directory
        threads: Number of download threads

    Returns:
        Path to downloaded data

    Example:
        >>> # Download sequencing run
        >>> data_path = download_sra_run("SRR123456", "sra_data/")
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    run_dir = output_dir / sra_accession
    run_dir.mkdir(exist_ok=True)

    # Use sra-tools or direct download
    try:
        # Try direct ENA download first (more reliable)
        return _download_sra_from_ena(sra_accession, run_dir, threads)
    except Exception as e:
        logger.warning(f"ENA download failed: {e}")

    # Fallback to SRA toolkit
    try:
        return _download_sra_with_toolkit(sra_accession, run_dir, threads)
    except Exception as e:
        logger.warning(f"SRA toolkit download failed: {e}")

    raise RuntimeError(
        f"All SRA download methods failed for {sra_accession}. "
        "Ensure either ENA FTP is accessible or SRA toolkit (fastq-dump) is installed."
    )


def _download_sra_from_ena(accession: str, output_dir: Path, threads: int) -> Path:
    """Download SRA data from ENA."""
    # ENA URL pattern (FTP). We keep the original behavior but add heartbeat + retry.
    # Query ENA API for exact FTP links to avoid guessing and 550 errors
    import requests

    from metainformant.core.io.download import download_with_progress

    api_url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport?"
        f"accession={accession}&result=read_run&fields=fastq_ftp&format=json"
    )

    try:
        resp = requests.get(api_url, timeout=15)
        resp.raise_for_status()
        metadata = resp.json()
        if not metadata or "fastq_ftp" not in metadata[0]:
            raise ValueError(f"No fastq_ftp fields returned for {accession}")

        ftp_links_raw = metadata[0]["fastq_ftp"]
        # Can be multiple links separated by semicolon (e.g. for paired end)
        patterns = [f"http://{link}" for link in ftp_links_raw.split(";") if link]

    except Exception as e:
        logger.warning(f"Failed to fetch ENA API for {accession}: {e}. Falling back to default patterns.")
        ena_url = f"http://ftp.sra.ebi.ac.uk/vol1/fastq/{accession[:6]}"
        if len(accession) > 9:
            # Simple fallback for 10 or 11 chars
            if len(accession) == 10:
                ena_url += f"/00{accession[-1]}"  # noqa: E701
            elif len(accession) == 11:
                ena_url += f"/0{accession[-2:]}"  # noqa: E701
        ena_url += f"/{accession}"

        patterns = [
            f"{ena_url}/{accession}.fastq.gz",
            f"{ena_url}/{accession}_1.fastq.gz",
            f"{ena_url}/{accession}_2.fastq.gz",
        ]

    downloaded_files: list[Path] = []
    last_error: str | None = None

    for url in patterns:
        filename = url.split("/")[-1]
        output_file = output_dir / filename
        result = download_with_progress(
            url,
            output_file,
            protocol="ftp",
            show_progress=False,
            heartbeat_interval=5,
            max_retries=3,
            chunk_size=1024 * 1024,
            timeout=60,
            resume=False,
        )
        if result.success:
            downloaded_files.append(output_file)
            logger.info(f"Downloaded {filename}")
        else:
            last_error = result.error

    if not downloaded_files:
        raise RuntimeError(f"No files downloaded from ENA for {accession}: {last_error or 'unknown error'}")

    return output_dir


def _download_sra_with_toolkit(accession: str, output_dir: Path, threads: int) -> Path:
    """Download using SRA toolkit."""
    import subprocess

    # Check if fastq-dump is available
    try:
        subprocess.run(["fastq-dump", "--version"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        raise RuntimeError("SRA toolkit (fastq-dump) not found")

    cmd = [
        "fastq-dump",
        "--split-files",  # Split paired-end reads
        "--gzip",  # Compress output
        "--outdir",
        str(output_dir),
        accession,
    ]

    logger.info(f"Running fastq-dump for {accession}")

    result = subprocess.run(cmd, capture_output=True, text=True)
    result.check_returncode()

    return output_dir


def download_sra_project(
    project_id: str,
    output_dir: str | Path,
    threads: int = 1,
    email: str | None = None,
    *,
    allow_anonymous: bool = False,
) -> List[Path]:
    """Download all runs from an SRA project.

    Args:
        project_id: SRA project ID (e.g., "PRJNA123456")
        output_dir: Output directory
        threads: Number of threads
        email: Explicit NCBI contact; otherwise ``NCBI_EMAIL``.
        allow_anonymous: Explicitly permit anonymous NCBI metadata lookup.

    Returns:
        List of paths to downloaded run directories

    Example:
        >>> # Download entire project
        >>> run_paths = download_sra_project("PRJNA123456", "project_data/")
        >>> len(run_paths) > 0
        True
    """
    # First, get list of runs in the project
    run_accessions = _get_project_runs(project_id, email=email, allow_anonymous=allow_anonymous)

    if not run_accessions:
        logger.warning(f"No runs found for project {project_id}")
        return []

    output_dir = Path(output_dir)
    project_dir = output_dir / project_id
    project_dir.mkdir(parents=True, exist_ok=True)

    downloaded_runs = []

    for run_acc in run_accessions:
        try:
            run_path = download_sra_run(run_acc, project_dir, threads)
            downloaded_runs.append(run_path)
            logger.info(f"Downloaded run {run_acc}")
        except Exception as e:
            logger.error(f"Failed to download run {run_acc}: {e}")
            continue

    return downloaded_runs


def _get_project_runs(
    project_id: str,
    email: str | None = None,
    *,
    allow_anonymous: bool = False,
) -> List[str]:
    """Get list of SRA runs for a project."""
    contact = resolve_ncbi_contact(email, allow_anonymous=allow_anonymous)
    # Query NCBI for project runs
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    # Search for runs in this project
    search_params: Dict[str, str | int] = {
        "db": "sra",
        "term": f"{project_id}[BioProject]",
        "retmax": 1000,
        "retmode": "json",
    }
    if contact.email:
        search_params["email"] = contact.email

    try:
        response = requests.get(f"{base_url}/esearch.fcgi", params=search_params, timeout=30)
        response.raise_for_status()

        data = response.json()

        if "esearchresult" not in data:
            return []

        id_list = data["esearchresult"].get("idlist", [])

        # Get run accessions from summaries
        if id_list:
            summary_params = {
                "db": "sra",
                "id": ",".join(id_list[:100]),
                "retmode": "json",
            }  # Limit for performance
            if contact.email:
                summary_params["email"] = contact.email

            summary_response = requests.get(f"{base_url}/esummary.fcgi", params=summary_params, timeout=30)
            summary_response.raise_for_status()

            summary_data = summary_response.json()

            # Import re at the module level or locally
            import re

            run_accessions = []
            if "result" in summary_data:
                for uid in id_list[:100]:
                    if uid in summary_data["result"]:
                        record = summary_data["result"][uid]
                        if not isinstance(record, dict):
                            continue

                        runs_data = record.get("runs", "")
                        if isinstance(runs_data, str) and runs_data.strip():
                            # EUtils returns an XML string inside the JSON for SRA runs
                            # e.g.: <Run acc="SRR13968413" .../>
                            matches = re.findall(r'acc="([^"]+)"', runs_data)
                            for match in matches:
                                run_accessions.append(match)
                        elif isinstance(runs_data, dict):
                            # Fallback in case they ever return structured JSON
                            runs_list = runs_data.get("run", [])
                            if isinstance(runs_list, list):
                                for run in runs_list:
                                    run_acc = run.get("@acc")
                                    if run_acc:
                                        run_accessions.append(run_acc)
                            elif isinstance(runs_list, dict):
                                run_acc = runs_list.get("@acc")
                                if run_acc:
                                    run_accessions.append(run_acc)

            return run_accessions

    except requests.RequestException as e:
        logger.error(f"Failed to get project runs: {e}")

    return []


def search_sra_for_organism(
    organism: str,
    max_results: int = 100,
    email: str | None = None,
    *,
    allow_anonymous: bool = False,
) -> List[Dict[str, Any]]:
    """Search SRA for sequencing data from a specific organism.

    Args:
        organism: Organism name
        max_results: Maximum number of results
        email: Explicit NCBI contact; otherwise ``NCBI_EMAIL``.
        allow_anonymous: Explicitly permit anonymous NCBI metadata lookup.

    Returns:
        List of SRA records

    Example:
        >>> # Find human sequencing data
        >>> results = search_sra_for_organism("Homo sapiens", max_results=10)
        >>> len(results) <= 10
        True
    """
    contact = resolve_ncbi_contact(email, allow_anonymous=allow_anonymous)
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    query = f'"{organism}"[Organism] AND "RNA-seq"[Strategy]'

    search_params: Dict[str, str | int] = {
        "db": "sra",
        "term": query,
        "retmax": min(max_results, 1000),
        "retmode": "json",
    }
    if contact.email:
        search_params["email"] = contact.email

    try:
        response = requests.get(f"{base_url}/esearch.fcgi", params=search_params, timeout=30)
        response.raise_for_status()

        search_data = response.json()

        if "esearchresult" not in search_data:
            return []

        id_list = search_data["esearchresult"].get("idlist", [])

        if not id_list:
            return []

        # Get summaries
        summary_params = {
            "db": "sra",
            "id": ",".join(id_list[:max_results]),
            "retmode": "json",
        }
        if contact.email:
            summary_params["email"] = contact.email

        summary_response = requests.get(f"{base_url}/esummary.fcgi", params=summary_params, timeout=30)
        summary_response.raise_for_status()

        summary_data = summary_response.json()

        results = []
        if "result" in summary_data:
            for uid in id_list[:max_results]:
                if uid in summary_data["result"]:
                    record = summary_data["result"][uid]
                    result = {
                        "id": uid,
                        "accession": record.get("runs", {}).get("run", {}).get("@acc"),
                        "experiment": record.get("expxml", {}).get("summary"),
                        "platform": record.get("expxml", {}).get("platform"),
                        "library_strategy": record.get("expxml", {}).get("library_strategy"),
                        "spots": record.get("runs", {}).get("run", {}).get("@spots"),
                        "bases": record.get("runs", {}).get("run", {}).get("@bases"),
                        "size_mb": record.get("runs", {}).get("run", {}).get("@size_MB"),
                        "taxon_id": record.get("taxon", {}).get("@taxid"),
                        "organism": record.get("taxon", {}).get("@scientificname"),
                    }
                    results.append(result)

        return results

    except requests.RequestException as e:
        logger.error(f"SRA search failed for {organism}: {e}")
        return []


def download_annotation(accession: str, output_dir: str | Path) -> Path:
    """Download genome annotation (GFF/GTF) files.

    Args:
        accession: Genome accession
        output_dir: Output directory

    Returns:
        Path to annotation files
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Try NCBI Datasets API for annotation
    try:
        api_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download"
        params = {"include": "gff3"}

        response = requests.get(api_url, params=params, timeout=60)
        response.raise_for_status()

        # Save annotation
        zip_path = output_dir / f"{accession}_annotation.zip"
        with open(zip_path, "wb") as f:
            f.write(response.content)

        # Extract
        extract_dir = output_dir / f"{accession}_annotation"
        extract_dir.mkdir(exist_ok=True)

        _extract_zip_safely(zip_path, extract_dir)

        zip_path.unlink()
        return extract_dir

    except Exception as e:
        logger.warning(f"NCBI Datasets annotation download failed: {e}, trying FTP fallback")

        # FTP fallback: download GFF3 annotation directly
        try:
            import ftplib  # nosec B402 - retained as a documented NCBI legacy fallback; primary path is HTTPS.

            prefix = accession[:3]
            digits = accession.split("_")[1].split(".")[0]
            ftp_dir = f"/genomes/all/{prefix}/{digits[:3]}/{digits[3:6]}/{digits[6:9]}"

            annotation_dir = output_dir / f"{accession}_annotation"
            annotation_dir.mkdir(exist_ok=True)

            ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")  # nosec B321 - NCBI fallback endpoint.
            ftp.login()
            try:
                ftp.cwd(ftp_dir)
                dirs = ftp.nlst()
                target_dir = next((d for d in dirs if d.startswith(accession)), None)
                if not target_dir:
                    raise RuntimeError(f"No directory found for {accession} on NCBI FTP")

                ftp.cwd(target_dir)
                files = ftp.nlst()

                gff_files = [f for f in files if f.endswith("_genomic.gff.gz") and "_from_genomic" not in f]
                if gff_files:
                    filename = gff_files[0]
                    local_path = annotation_dir / filename
                    logger.info(f"Downloading annotation {filename} from NCBI FTP...")
                    with open(local_path, "wb") as out_file:
                        ftp.retrbinary(f"RETR {filename}", out_file.write)
                    logger.info(f"Downloaded annotation {filename}")
                else:
                    raise RuntimeError(f"No GFF annotation file found for {accession} on NCBI FTP")
            finally:
                ftp.quit()

            return annotation_dir

        except Exception as ftp_err:
            raise RuntimeError(
                f"All annotation download methods failed for {accession}: " f"API error: {e}, FTP error: {ftp_err}"
            )


#: Base URL of the GWAS Catalog REST API, the HAL/JSON service :func:`download_variant_data` reads.
GWAS_CATALOG_API_URL = "https://www.ebi.ac.uk/gwas/rest/api"

#: Columns of the association summary table written by :func:`download_variant_data`.
SUMMARY_STATS_COLUMNS: Tuple[str, ...] = ("variant", "risk_allele", "p_value", "odds_ratio", "beta")

#: Associations requested per page unless ``page_size`` says otherwise.
DEFAULT_ASSOCIATIONS_PAGE_SIZE = 200

#: Pages one download may read before it gives up (a service that keeps announcing more pages is misbehaving).
MAX_ASSOCIATION_PAGES = 5000

_CATALOG_TIMEOUT_SECONDS = 30


@dataclass
class AssociationsPage:
    """One decoded page of the ``/studies/{accession}/associations`` response.

    Attributes:
        rows: One row per association, keyed by :data:`SUMMARY_STATS_COLUMNS`; every value is text, and a value
            the response does not carry is the empty string.
        number: Zero-based number the service gives this page, or None if the response has no ``page`` block.
        total_pages: Number of pages the service reports, or None.
        total_elements: Number of associations the study has according to the service, or None.
        has_next: Whether the response announces a further page.
    """

    rows: List[Dict[str, str]] = field(default_factory=list)
    number: Optional[int] = None
    total_pages: Optional[int] = None
    total_elements: Optional[int] = None
    has_next: bool = False


def _excerpt(text: str, limit: int = 200) -> str:
    """The first ``limit`` characters of a response body, for an error message."""
    flat = " ".join(text.split())
    return flat if len(flat) <= limit else flat[:limit] + "..."


def _json_object(payload: Any, what: str) -> Mapping[str, Any]:
    """``payload`` (JSON text, or an already decoded value) as a JSON object, or ``ValueError``."""
    decoded = payload
    if isinstance(payload, (str, bytes, bytearray)):
        try:
            decoded = json.loads(payload)
        except ValueError as exc:
            raise ValueError(f"{what} is not valid JSON: {exc}") from exc
    if not isinstance(decoded, Mapping):
        raise ValueError(f"{what} must be a JSON object, got {type(decoded).__name__}")
    return decoded


def decode_catalog_response(status_code: int, body: str, url: str) -> Dict[str, Any]:
    """Decode a GWAS Catalog REST response, accepting only HTTP 200 with a JSON object as its body.

    This is the status check of every request :func:`download_variant_data` makes. Any other status, such as a
    404 for an unknown study or a 5xx from the service, is an error: it used to be skipped, which left a
    download that looked complete but held no association table.

    Args:
        status_code: HTTP status of the response.
        body: Response body text.
        url: Requested URL, named in the error message.

    Returns:
        The decoded JSON object.

    Raises:
        requests.HTTPError: If the status is not 200.
        ValueError: If a 200 response is not a JSON object.
    """
    if status_code != 200:
        raise requests.HTTPError(
            f"GWAS Catalog request to {url} returned HTTP {status_code}, not 200: {_excerpt(body)}"
        )
    return dict(_json_object(body, f"GWAS Catalog response from {url}"))


def _fetch_catalog_json(url: str, params: Optional[Mapping[str, Any]] = None) -> Dict[str, Any]:
    """GET ``url`` from the GWAS Catalog and return its decoded JSON object (see :func:`decode_catalog_response`)."""
    response = requests.get(url, params=params, timeout=_CATALOG_TIMEOUT_SECONDS)
    return decode_catalog_response(response.status_code, response.text, url)


def _catalog_number(value: Any, field_name: str) -> float:
    """``value`` as a finite float; the GWAS Catalog sends numbers as JSON numbers, and older exports as text."""
    if isinstance(value, bool) or value is None:
        raise ValueError(f"{field_name} must be a number, got {value!r}")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must be a number, got {value!r}") from exc
    if not math.isfinite(number):
        raise ValueError(f"{field_name} must be finite, got {value!r}")
    return number


def _is_absent(value: Any) -> bool:
    """Whether a response field carries no value (missing, null or blank)."""
    return value is None or (isinstance(value, str) and not value.strip())


def _format_number(number: float) -> str:
    """Shortest text that reads back as ``number`` (``2`` for 2.0, ``0.37`` for 0.37)."""
    return format(number, ".15g")


def format_catalog_p_value(mantissa: Any, exponent: Any) -> str:
    """Text ``"<mantissa>e<exponent>"`` for the ``pvalueMantissa`` / ``pvalueExponent`` pair of an association.

    Both fields may be JSON numbers (``2`` and ``-8`` give ``"2e-8"``) or numeric text, which the earlier string
    concatenation turned into a ``TypeError`` for numbers. No p-value (both fields absent) is the empty string;
    one field without the other is malformed.

    Raises:
        ValueError: If only one of the two is present, a field is not a finite number, or the exponent is not
            a whole number.
    """
    if _is_absent(mantissa) and _is_absent(exponent):
        return ""
    if _is_absent(mantissa) or _is_absent(exponent):
        raise ValueError(f"pvalueMantissa and pvalueExponent must come as a pair, got {mantissa!r} and {exponent!r}")
    base = _catalog_number(mantissa, "pvalueMantissa")
    power = _catalog_number(exponent, "pvalueExponent")
    if not power.is_integer():
        raise ValueError(f"pvalueExponent must be a whole number, got {exponent!r}")
    return f"{_format_number(base)}e{int(power)}"


def _optional_number_text(value: Any, field_name: str) -> str:
    """``value`` formatted as text, or the empty string when the response carries none."""
    return "" if _is_absent(value) else _format_number(_catalog_number(value, field_name))


def _optional_int(value: Any, field_name: str) -> Optional[int]:
    """A non-negative whole number from the ``page`` block of a response, or None when it is absent."""
    if _is_absent(value):
        return None
    number = _catalog_number(value, field_name)
    if not number.is_integer() or number < 0:
        raise ValueError(f"{field_name} must be a non-negative whole number, got {value!r}")
    return int(number)


def _list_field(owner: Mapping[str, Any], key: str, where: str) -> List[Any]:
    """``owner[key]`` as a list: absent or null is empty, anything but a list is an error."""
    value = owner.get(key)
    if value is None:
        return []
    if not isinstance(value, list):
        raise ValueError(f"{where}.{key} must be a list, got {type(value).__name__}")
    return value


def _split_risk_allele(name: str) -> Tuple[str, str]:
    """Split a ``riskAlleleName`` such as ``"rs7329174-A"`` into ``("rs7329174", "A")``.

    The allele follows the last hyphen; a name without a hyphen is a variant with no stated allele.
    """
    text = name.strip()
    variant, separator, allele = text.rpartition("-")
    if not separator:
        return text, ""
    return variant, allele


def _association_row(index: int, association: Any) -> Dict[str, str]:
    """One summary-table row from one element of ``_embedded.associations``."""
    where = f"association {index}"
    if not isinstance(association, Mapping):
        raise ValueError(f"{where} must be a JSON object, got {type(association).__name__}")

    variants: List[str] = []
    alleles: List[str] = []
    for locus_index, locus in enumerate(_list_field(association, "loci", where)):
        if not isinstance(locus, Mapping):
            raise ValueError(f"{where}.loci[{locus_index}] must be a JSON object, got {type(locus).__name__}")
        for risk_allele in _list_field(locus, "strongestRiskAlleles", f"{where}.loci[{locus_index}]"):
            if not isinstance(risk_allele, Mapping):
                raise ValueError(f"{where} has a strongestRiskAlleles entry that is not a JSON object")
            name = risk_allele.get("riskAlleleName")
            if _is_absent(name):
                continue
            if not isinstance(name, str):
                raise ValueError(f"{where} has a riskAlleleName that is not text: {name!r}")
            variant, allele = _split_risk_allele(name)
            variants.append(variant)
            alleles.append(allele)

    try:
        p_value = format_catalog_p_value(association.get("pvalueMantissa"), association.get("pvalueExponent"))
        odds_ratio = _optional_number_text(association.get("orPerCopyNum"), "orPerCopyNum")
        beta = _optional_number_text(association.get("betaNum"), "betaNum")
    except ValueError as exc:
        raise ValueError(f"{where}: {exc}") from exc

    # A haplotype or interaction has several risk alleles but one p-value: one row, alleles joined with ';'
    return {
        "variant": ";".join(variants),
        "risk_allele": ";".join(alleles),
        "p_value": p_value,
        "odds_ratio": odds_ratio,
        "beta": beta,
    }


def parse_catalog_study(payload: Any) -> Dict[str, Optional[str]]:
    """Trait and PubMed id of a study from the ``/studies/{accession}`` response.

    Args:
        payload: The response text, or its decoded JSON object.

    Returns:
        ``{"trait": ..., "pubmed_id": ...}``; either is None when the response does not carry it (a study with
        no ``diseaseTrait`` or ``publicationInfo`` block is not an error).

    Raises:
        ValueError: If the payload is not a JSON object.
    """
    study = _json_object(payload, "GWAS Catalog study response")
    publication = study.get("publicationInfo")
    disease_trait = study.get("diseaseTrait")
    pubmed_id = publication.get("pubmedId") if isinstance(publication, Mapping) else None
    trait = disease_trait.get("traitName") if isinstance(disease_trait, Mapping) else None
    return {
        "trait": None if _is_absent(trait) else str(trait),
        "pubmed_id": None if _is_absent(pubmed_id) else str(pubmed_id),
    }


def parse_catalog_associations(payload: Any) -> AssociationsPage:
    """Decode one page of the ``/studies/{accession}/associations`` response into summary-table rows.

    This is a pure function of the response, written against the response shape ``download_variant_data`` has
    always assumed (the GWAS Catalog serves HAL JSON). No real response could be fetched while this was written,
    so the shape below is an assumption that the tests reproduce, not a recorded response::

        {
          "_embedded": {"associations": [
            {
              "pvalueMantissa": 2, "pvalueExponent": -8,      # JSON numbers (text is accepted too)
              "orPerCopyNum": 1.07, "betaNum": null,
              "loci": [{"strongestRiskAlleles": [{"riskAlleleName": "rs7329174-A"}]}]
            }
          ]},
          "_links": {"next": {"href": "..."}},                 # present while a further page exists
          "page": {"size": 200, "totalElements": 412, "totalPages": 3, "number": 0}
        }

    The ``variant`` and ``risk_allele`` columns come from ``riskAlleleName`` (``<variant>-<allele>``, split at
    the last hyphen); an association with several risk alleles (a haplotype or an interaction) keeps them in
    order, joined with ``;``. ``p_value`` is ``"<mantissa>e<exponent>"``. A field the response does not carry
    gives an empty cell, never a placeholder value; ``_embedded`` may be absent for a study with no
    associations. A page is "the last" when the ``page`` block says so, or, without a ``page`` block, when
    there is no ``_links.next``.

    Args:
        payload: The response text, or its decoded JSON object.

    Returns:
        The page's rows and its paging information.

    Raises:
        ValueError: If the payload is not a JSON object, an association has a field of the wrong type (a
            non-numeric p-value, a risk allele that is not text, ...), or a ``page`` field is not a
            non-negative whole number.
    """
    document = _json_object(payload, "GWAS Catalog associations response")

    embedded = document.get("_embedded")
    if embedded is not None and not isinstance(embedded, Mapping):
        raise ValueError(f"_embedded must be a JSON object, got {type(embedded).__name__}")
    associations = _list_field(embedded or {}, "associations", "_embedded")
    rows = [_association_row(index, association) for index, association in enumerate(associations)]

    page_block = document.get("page")
    if page_block is not None and not isinstance(page_block, Mapping):
        raise ValueError(f"page must be a JSON object, got {type(page_block).__name__}")
    number = _optional_int((page_block or {}).get("number"), "page.number")
    total_pages = _optional_int((page_block or {}).get("totalPages"), "page.totalPages")
    total_elements = _optional_int((page_block or {}).get("totalElements"), "page.totalElements")

    if number is not None and total_pages is not None:
        has_next = number + 1 < total_pages
    else:
        links = document.get("_links")
        next_link = links.get("next") if isinstance(links, Mapping) else None
        has_next = isinstance(next_link, Mapping) and not _is_absent(next_link.get("href"))

    return AssociationsPage(
        rows=rows, number=number, total_pages=total_pages, total_elements=total_elements, has_next=has_next
    )


def write_summary_stats_table(path: str | Path, rows: Iterable[Mapping[str, str]]) -> Path:
    """Write association rows as a tab-separated table headed by :data:`SUMMARY_STATS_COLUMNS`.

    The rows go through :mod:`csv`, so a tab, quote or newline inside a cell is quoted instead of adding a
    column or a row.

    Args:
        path: File to write (its directory must exist).
        rows: Rows keyed by :data:`SUMMARY_STATS_COLUMNS`, as :func:`parse_catalog_associations` returns them.

    Returns:
        ``path`` as a :class:`~pathlib.Path`.
    """
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(SUMMARY_STATS_COLUMNS), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return path


class AssociationsPager:
    """Sans-IO state machine that pages through ``/studies/{accession}/associations``.

    The pager owns the paging decisions of :func:`download_variant_data` and does no I/O: ask it for the query
    of the next request, fetch that page, hand the decoded response back::

        pager = AssociationsPager(page_size=200)
        while (query := pager.request()) is not None:
            pager.accept(parse_catalog_associations(fetch(associations_url, query)))
        pager.check_complete()
        pager.rows, pager.truncated, pager.available

    Attributes:
        rows: The association rows collected so far (never more than ``max_associations``).
        available: Number of associations the service reports for the study, or None if it never said.
        truncated: True once ``max_associations`` cut the table short of the study's associations.
        pages_read: Number of pages accepted.
    """

    def __init__(self, page_size: int = DEFAULT_ASSOCIATIONS_PAGE_SIZE, max_associations: Optional[int] = None) -> None:
        """Start a download at the first page.

        Args:
            page_size: Associations requested per page.
            max_associations: Stop after this many associations; None collects them all.

        Raises:
            ValueError: If ``page_size`` or ``max_associations`` is not a positive whole number.
        """
        if isinstance(page_size, bool) or not isinstance(page_size, int) or page_size < 1:
            raise ValueError(f"page_size must be a positive whole number, got {page_size!r}")
        if max_associations is not None and (
            isinstance(max_associations, bool) or not isinstance(max_associations, int) or max_associations < 1
        ):
            raise ValueError(f"max_associations must be a positive whole number or None, got {max_associations!r}")
        self.page_size = page_size
        self.max_associations = max_associations
        self.rows: List[Dict[str, str]] = []
        self.available: Optional[int] = None
        self.truncated = False
        self.pages_read = 0
        self._next_page: Optional[int] = 0

    def request(self) -> Optional[Dict[str, int]]:
        """Query parameters (``page`` and ``size``) of the next request, or None when no page is left.

        Raises:
            RuntimeError: If the service has already announced :data:`MAX_ASSOCIATION_PAGES` pages and still
                announces another.
        """
        if self._next_page is None:
            return None
        if self.pages_read >= MAX_ASSOCIATION_PAGES:
            raise RuntimeError(
                f"the service announces more than {MAX_ASSOCIATION_PAGES} pages of associations at "
                f"page_size={self.page_size}; read {len(self.rows)} associations and stopped"
            )
        return {"page": self._next_page, "size": self.page_size}

    def accept(self, page: AssociationsPage) -> None:
        """Take the decoded response to the request :meth:`request` last returned.

        Raises:
            ValueError: If the service answered with another page than the one requested (it ignored the
                page parameter), which would make the download repeat or skip associations.
            RuntimeError: If no page was outstanding.
        """
        requested = self._next_page
        if requested is None:
            raise RuntimeError("accept() called although every page has been read")
        if page.number is not None and page.number != requested:
            raise ValueError(f"requested association page {requested} but the service returned page {page.number}")

        self.pages_read += 1
        if page.total_elements is not None:
            self.available = page.total_elements
        self.rows.extend(page.rows)
        following = requested + 1 if page.has_next else None

        if self.max_associations is not None and len(self.rows) >= self.max_associations:
            # The table is short exactly when this cap cut rows off or the service still has pages
            self.truncated = len(self.rows) > self.max_associations or following is not None
            self.rows = self.rows[: self.max_associations]
            following = None
        self._next_page = following

    def check_complete(self) -> None:
        """Check the rows read against the study's association count, unless the table was cut on purpose.

        Raises:
            RuntimeError: If the service reported a count and it differs from the number of rows read.
        """
        if not self.truncated and self.available is not None and len(self.rows) != self.available:
            raise RuntimeError(
                f"the service reports {self.available} associations but {len(self.rows)} were read; "
                "the download is incomplete"
            )


def download_variant_data(
    study_accession: str,
    output_dir: str | Path,
    data_type: str = "vcf",
    *,
    page_size: int = DEFAULT_ASSOCIATIONS_PAGE_SIZE,
    max_associations: Optional[int] = None,
) -> Path:
    """Download variant data for a GWAS study.

    Only the association summary table of a study can be downloaded
    (``data_type="summary_stats"``): the GWAS Catalog REST API serves association
    records, not study-level genotype files. Genotype VCFs (``"vcf"``, the default,
    and ``"both"``) are not implemented, so those values raise ``NotImplementedError``
    before anything is written or requested; pass ``data_type="summary_stats"`` to
    download.

    The table has one row per association (see :data:`SUMMARY_STATS_COLUMNS`) and is read page by page until
    the service reports the last page (:class:`AssociationsPager`), so a study with more associations than one
    page holds is downloaded whole. Every request must answer HTTP 200 (:func:`decode_catalog_response`); a
    failure raises and leaves nothing on disk. When the service reports how many associations the study has,
    the rows read are checked against it. The response shape this relies on is documented at
    :func:`parse_catalog_associations`.

    Args:
        study_accession: GWAS study accession (e.g., GCST123456)
        output_dir: Output directory for downloaded files
        data_type: ``"summary_stats"``, the only data type that can be downloaded.
            ``"vcf"`` and ``"both"`` raise ``NotImplementedError``.
        page_size: Associations requested per page.
        max_associations: Keep at most this many associations. ``None`` (the default) keeps all of them. When
            the cap cuts the table short the download is flagged, not silent: a warning is logged and the
            metadata file records ``"truncated": true`` with the number of associations available.

    Returns:
        Path to downloaded data directory. It holds ``<accession>_summary_stats.tsv`` and
        ``<accession>_metadata.json`` (trait, PubMed id, the files written, ``n_associations``,
        ``n_associations_available`` and ``truncated``).

    Raises:
        ValueError: If study accession is invalid, data_type is not recognised, ``page_size`` or
            ``max_associations`` is not a positive whole number, or a response cannot be decoded
        NotImplementedError: If data_type asks for genotype VCFs ("vcf" or "both")
        requests.RequestException: If a request fails or does not answer HTTP 200
        RuntimeError: If the service reports more associations than were read, or announces more than
            ``MAX_ASSOCIATION_PAGES`` pages
    """
    if not re.fullmatch(r"GCST[0-9]+", study_accession):
        raise ValueError(f"Invalid GWAS Catalog study accession: {study_accession!r}")

    if data_type in ("vcf", "both"):
        raise NotImplementedError(
            f"data_type={data_type!r} is not implemented: the GWAS Catalog REST API publishes association "
            "summaries, not study-level genotype VCFs, so there is nothing to download. Use "
            "data_type='summary_stats' for the association table; for genotypes follow the data-availability "
            "statement of the study's publication, or use download_variant_database for dbSNP, 1000 Genomes "
            "or HapMap variants."
        )
    if data_type != "summary_stats":
        raise ValueError(
            f"Unknown data_type {data_type!r}: expected 'summary_stats' ('vcf' and 'both' are not implemented)"
        )

    pager = AssociationsPager(page_size=page_size, max_associations=max_associations)

    logger.info(f"Downloading variant data for {study_accession}")

    study_url = f"{GWAS_CATALOG_API_URL}/studies/{study_accession}"
    associations_url = f"{study_url}/associations"

    try:
        study = parse_catalog_study(_fetch_catalog_json(study_url))
        logger.info(f"Study {study_accession}: {study['trait']} (PMID: {study['pubmed_id']})")

        while (query := pager.request()) is not None:
            pager.accept(parse_catalog_associations(_fetch_catalog_json(associations_url, query)))
        pager.check_complete()
    except (requests.RequestException, ValueError, RuntimeError) as e:
        logger.error(f"Failed to download variant data for {study_accession}: {e}")
        raise

    rows = pager.rows
    if pager.truncated:
        available = pager.available if pager.available is not None else "an unknown number of"
        logger.warning(f"{study_accession}: wrote the first {len(rows)} of {available} associations (max_associations)")

    output_dir = Path(output_dir)
    study_dir = output_dir / study_accession
    study_dir.mkdir(parents=True, exist_ok=True)

    stats_file = write_summary_stats_table(study_dir / f"{study_accession}_summary_stats.tsv", rows)
    logger.info(f"Downloaded {len(rows)} associations to {stats_file}")

    metadata_file = study_dir / f"{study_accession}_metadata.json"
    with open(metadata_file, "w", encoding="utf-8") as handle:
        json.dump(
            {
                "accession": study_accession,
                "trait": study["trait"],
                "pubmed_id": study["pubmed_id"],
                "downloaded_files": [str(stats_file)],
                "n_associations": len(rows),
                "n_associations_available": pager.available,
                "truncated": pager.truncated,
                "download_date": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            },
            handle,
            indent=2,
        )

    logger.info(f"Downloaded variant data to {study_dir}")
    return study_dir


def prepare_reference_genome(output_dir: str | Path) -> Path | None:
    """Decompress and index reference genome FASTA files.

    Finds compressed FASTA files (.fna.gz), decompresses them,
    and runs bwa index and samtools faidx.

    Args:
        output_dir: Directory containing reference genome files

    Returns:
        Path to the primary uncompressed FASTA file, or None if not found
    """
    output_dir = Path(output_dir)
    import subprocess

    # Decompress
    fastas = list(output_dir.glob("**/*.fna.gz"))
    for f in fastas:
        logger.info(f"Decompressing {f}...")
        subprocess.run(["gunzip", "-f", str(f)], check=True)

    # Re-discover uncompressed FASTAs
    final_fastas = list(output_dir.glob("**/*.fna")) + list(output_dir.glob("**/*.fasta"))
    if not final_fastas:
        return None

    ref_path = final_fastas[0]

    # Indexing
    logger.info(f"Indexing reference: {ref_path}")
    subprocess.run(["bwa", "index", str(ref_path)], check=True)
    subprocess.run(["samtools", "faidx", str(ref_path)], check=True)

    return ref_path
