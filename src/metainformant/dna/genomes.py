"""Genome data retrieval and processing utilities.

This module provides tools for downloading genome assemblies, validating
accessions, and retrieving genome metadata from NCBI and other sources.
"""

from __future__ import annotations

import json
import re
import time
from pathlib import Path
from typing import Any, Dict, List, Optional
from urllib.parse import urlparse

import requests

from metainformant.core import io, logging, paths

logger = logging.get_logger(__name__)


def download_genome_package(accession: str, output_dir: str | Path,
                          include: List[str] | None = None) -> Path:
    """Download genome package from NCBI Datasets API.

    Args:
        accession: Genome accession (e.g., "GCF_000001405.39" for GRCh38)
        output_dir: Output directory for downloaded files
        include: List of file types to include (default: all)

    Returns:
        Path to downloaded genome package directory

    Raises:
        ValueError: If accession is invalid
        requests.RequestException: If download fails

    Example:
        >>> # This would download the human genome (large file!)
        >>> # package_dir = download_genome_package("GCF_000001405.39", "genomes/")
        >>> # isinstance(package_dir, Path)
        >>> # True
    """
    if not validate_accession(accession):
        raise ValueError(f"Invalid genome accession: {accession}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Use NCBI Datasets API
    api_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download"

    params = {}
    if include:
        params['include'] = ','.join(include)

    logger.info(f"Downloading genome {accession} from NCBI Datasets")

    try:
        response = requests.get(api_url, params=params, timeout=60)
        response.raise_for_status()

        # Save as zip file
        zip_path = output_dir / f"{accession}.zip"
        with open(zip_path, 'wb') as f:
            f.write(response.content)

        logger.info(f"Downloaded genome package to {zip_path}")

        # Extract zip (basic extraction)
        import zipfile
        extract_dir = output_dir / accession
        extract_dir.mkdir(exist_ok=True)

        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)

        # Remove zip file
        zip_path.unlink()

        return extract_dir

    except requests.RequestException as e:
        logger.error(f"Failed to download genome {accession}: {e}")
        raise


def download_genome_package_best_effort(accession: str, output_dir: str | Path,
                                       include: List[str] | None = None,
                                       ftp_url: Optional[str] = None) -> Path:
    """Download genome package with fallback strategies.

    This function tries multiple approaches to download genome data:
    1. NCBI Datasets API
    2. FTP download if URL provided
    3. Best-effort retrieval

    Args:
        accession: Genome accession
        output_dir: Output directory
        include: File types to include
        ftp_url: FTP URL for fallback download

    Returns:
        Path to genome package directory

    Example:
        >>> # Try to download with fallbacks
        >>> # package_dir = download_genome_package_best_effort("GCF_000001405.39", "genomes/")
    """
    output_dir = Path(output_dir)

    try:
        # Try primary method first
        return download_genome_package(accession, output_dir, include)
    except Exception as e:
        logger.warning(f"Primary download failed: {e}")

        # Try FTP fallback if provided
        if ftp_url:
            try:
                return _download_from_ftp(ftp_url, output_dir / accession)
            except Exception as e2:
                logger.warning(f"FTP download failed: {e2}")

        # Create empty directory as fallback
        fallback_dir = output_dir / accession
        fallback_dir.mkdir(parents=True, exist_ok=True)

        logger.warning(f"All download methods failed for {accession}, created empty directory")
        return fallback_dir


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
    import re

    # NCBI assembly accession patterns
    # GCF_##########.## or GCA_##########
    pattern = r'^(GCF|GCA)_[0-9]{9,}(\.[0-9]+)?$'

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
        r'^GC[FA]_\d{9}\.\d+$',  # GCF_000001405.39 (RefSeq)
        r'^GCA_\d{9}\.\d+$',     # GCA_000001405.39 (GenBank)
        r'^ASM\d+$',             # ASM2732 (older format)
        r'^chr\d+$',             # chr1, chr2, etc.
        r'^NC_\d{6}$',           # NC_000001 (chromosome)
        r'^NT_\d{6}$',           # NT_000001 (contig)
        r'^NW_\d{8}$',           # NW_001838827 (scaffold)
    ]

    return any(re.match(pattern, accession) for pattern in patterns)


def get_genome_metadata(accession: str) -> Dict[str, Any]:
    """Retrieve genome metadata from NCBI.

    Args:
        accession: Genome accession

    Returns:
        Dictionary with genome metadata

    Raises:
        ValueError: If accession is invalid
        requests.RequestException: If API request fails

    Example:
        >>> # This would query NCBI for metadata
        >>> # metadata = get_genome_metadata("GCF_000001405.39")
        >>> # "organism_name" in metadata
        >>> # True
    """
    if not validate_accession(accession):
        raise ValueError(f"Invalid genome accession: {accession}")

    # Use NCBI Datasets API for metadata
    api_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}"

    try:
        response = requests.get(api_url, timeout=30)
        response.raise_for_status()

        data = response.json()

        # Extract relevant metadata
        if 'reports' in data and data['reports']:
            report = data['reports'][0]

            metadata = {
                'accession': accession,
                'organism_name': report.get('organism', {}).get('organism_name', 'Unknown'),
                'tax_id': report.get('organism', {}).get('tax_id'),
                'assembly_level': report.get('assembly_info', {}).get('assembly_level'),
                'assembly_method': report.get('assembly_info', {}).get('assembly_method'),
                'genome_size': report.get('assembly_stats', {}).get('total_sequence_length'),
                'gc_percent': report.get('assembly_stats', {}).get('gc_percent'),
                'contig_count': report.get('assembly_stats', {}).get('contig_count'),
                'scaffold_count': report.get('assembly_stats', {}).get('scaffold_count'),
                'release_date': report.get('assembly_info', {}).get('release_date'),
                'submitter': report.get('assembly_info', {}).get('submitter'),
            }

            return metadata
        else:
            return {'accession': accession, 'error': 'No metadata found'}

    except requests.RequestException as e:
        logger.error(f"Failed to retrieve metadata for {accession}: {e}")
        return {'accession': accession, 'error': str(e)}


def list_genome_assemblies(organism: str, max_results: int = 10) -> List[Dict[str, Any]]:
    """Search for genome assemblies for an organism.

    Args:
        organism: Organism name (e.g., "Homo sapiens")
        max_results: Maximum number of results to return

    Returns:
        List of genome assembly metadata dictionaries

    Example:
        >>> # Search for human genome assemblies
        >>> assemblies = list_genome_assemblies("Homo sapiens", max_results=5)
        >>> len(assemblies) <= 5
        True
    """
    # Use NCBI Datasets API search
    api_url = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome"

    params = {
        'filters.reference_only': 'false',
        'filters.assembly_source': 'refseq',
        'taxon': organism,
        'page_size': max_results
    }

    try:
        response = requests.get(api_url, params=params, timeout=30)
        response.raise_for_status()

        data = response.json()

        assemblies = []
        if 'reports' in data:
            for report in data['reports']:
                assembly_info = {
                    'accession': report.get('accession'),
                    'organism_name': report.get('organism', {}).get('organism_name'),
                    'assembly_level': report.get('assembly_info', {}).get('assembly_level'),
                    'genome_size': report.get('assembly_stats', {}).get('total_sequence_length'),
                    'release_date': report.get('assembly_info', {}).get('release_date'),
                }
                assemblies.append(assembly_info)

        return assemblies

    except requests.RequestException as e:
        logger.error(f"Failed to search assemblies for {organism}: {e}")
        return []


def download_reference_genome(species: str, output_dir: str | Path) -> Optional[Path]:
    """Download reference genome for a species.

    Args:
        species: Species name
        output_dir: Output directory

    Returns:
        Path to downloaded genome, or None if not found

    Example:
        >>> # Download human reference genome
        >>> # genome_path = download_reference_genome("Homo sapiens", "genomes/")
    """
    # First search for assemblies
    assemblies = list_genome_assemblies(species, max_results=1)

    if not assemblies:
        logger.warning(f"No genome assemblies found for {species}")
        return None

    # Take the first (usually the reference)
    accession = assemblies[0]['accession']

    try:
        return download_genome_package_best_effort(accession, output_dir)
    except Exception as e:
        logger.error(f"Failed to download reference genome for {species}: {e}")
        return None


def _download_from_ftp(ftp_url: str, output_dir: Path) -> Path:
    """Download genome data from FTP URL."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Parse FTP URL
    parsed = urlparse(ftp_url)
    filename = Path(parsed.path).name

    if not filename:
        filename = "genome_data.fasta"

    output_file = output_dir / filename

    logger.info(f"Downloading from FTP: {ftp_url}")
    from metainformant.core.download import download_with_progress

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


def get_chromosome_lengths(accession: str) -> Dict[str, int]:
    """Get chromosome/contig lengths for a genome assembly.

    Args:
        accession: Genome accession

    Returns:
        Dictionary mapping chromosome names to lengths

    Example:
        >>> # Get chromosome lengths for human genome
        >>> # lengths = get_chromosome_lengths("GCF_000001405.39")
        >>> # isinstance(lengths, dict)
        >>> # True
    """
    # This would typically query NCBI for chromosome information
    # For now, return empty dict as placeholder
    logger.warning("get_chromosome_lengths not yet fully implemented")
    return {}


def validate_genome_files(genome_dir: str | Path) -> Dict[str, Any]:
    """Validate genome assembly files.

    Args:
        genome_dir: Directory containing genome files

    Returns:
        Validation results dictionary

    Example:
        >>> # Validate downloaded genome
        >>> # results = validate_genome_files("genomes/GCF_000001405.39/")
        >>> # results['valid']
        >>> # True
    """
    genome_dir = Path(genome_dir)

    if not genome_dir.exists():
        return {'valid': False, 'error': 'Directory does not exist'}

    # Look for common genome file types
    fasta_files = list(genome_dir.glob("*.fasta")) + list(genome_dir.glob("*.fa"))
    gff_files = list(genome_dir.glob("*.gff")) + list(genome_dir.glob("*.gff3"))

    results = {
        'valid': True,
        'fasta_files': len(fasta_files),
        'gff_files': len(gff_files),
        'total_files': len(list(genome_dir.glob("*"))),
        'issues': []
    }

    if not fasta_files:
        results['issues'].append("No FASTA files found")
        results['valid'] = False

    return results



