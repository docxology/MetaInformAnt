"""RNA genome and transcriptome preparation for Kallisto quantification.

This module provides tools for preparing genomic references for RNA-seq analysis.
It handles downloading genome sequences, extracting transcripts from GFF3 files,
building Kallisto indices, and orchestrating the complete genome preparation workflow.

Main Functions:
    - find_rna_fasta_in_genome_dir: Locate RNA FASTA files in genome directories
    - get_expected_index_path: Get standardized path for Kallisto index
    - download_rna_fasta_from_ftp: Download RNA sequences from NCBI FTP
    - download_cds_fasta_from_ftp: Download CDS sequences from NCBI FTP
    - extract_transcripts_from_gff: Extract transcript sequences from GFF3 annotation
    - prepare_transcriptome_for_kallisto: Prepare transcriptome for indexing
    - build_kallisto_index: Build Kallisto index from transcriptome
    - prepare_genome_for_quantification: End-to-end genome setup
    - verify_genome_status: Check genome preparation status
    - orchestrate_genome_setup: Coordinate all genome preparation steps

Example:
    >>> from metainformant.rna import genome_prep
    >>> from pathlib import Path
    >>> fasta_path = genome_prep.find_rna_fasta_in_genome_dir(
    ...     Path("genome_dir"), "GCF_000001405.40"
    ... )
    >>> index_path = genome_prep.get_expected_index_path(
    ...     Path("output"), "Homo_sapiens"
    ... )
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from metainformant.core import logging

logger = logging.get_logger(__name__)


def find_rna_fasta_in_genome_dir(genome_dir: Path, accession: str) -> Optional[Path]:
    """Locate RNA FASTA file in genome directory.

    Searches for RNA FASTA files following standard NCBI directory patterns.
    Looks in order: API extraction → CLI extraction → direct root.
    Prefers unzipped files over gzipped.

    Args:
        genome_dir: Root directory containing genome data
        accession: Genome assembly accession (e.g., "GCF_000001405.40")

    Returns:
        Path to RNA FASTA file if found, None otherwise.
        Returns Path object (not string).

    Example:
        >>> from pathlib import Path
        >>> rna_fasta = find_rna_fasta_in_genome_dir(
        ...     Path("work/genomes"), "GCF_000001405.40"
        ... )
        >>> if rna_fasta:
        ...     print(f"Found: {rna_fasta}")
    """
    genome_dir = Path(genome_dir)

    # Search patterns in order of preference
    search_patterns = [
        # API extraction pattern
        genome_dir / "ncbi_dataset_api_extracted" / "ncbi_dataset" / "data" / accession,
        # CLI extraction pattern
        genome_dir / "ncbi_dataset_extracted" / "ncbi_dataset" / "data" / accession,
        # Direct in root
        genome_dir,
    ]

    for search_dir in search_patterns:
        if not search_dir.exists():
            continue

        # Prefer unzipped over zipped
        rna_fna = search_dir / "rna.fna"
        rna_fna_gz = search_dir / "rna.fna.gz"

        if rna_fna.exists():
            logger.debug(f"Found RNA FASTA: {rna_fna}")
            return rna_fna

        if rna_fna_gz.exists():
            logger.debug(f"Found RNA FASTA (gzipped): {rna_fna_gz}")
            return rna_fna_gz

    logger.debug(f"No RNA FASTA found for {accession}")
    return None


def get_expected_index_path(work_dir: Path, species_name: str) -> Path:
    """Get standardized path for Kallisto index.

    Computes the expected location of a Kallisto index based on work directory
    and species name. The index path is standardized for reproducibility.

    Args:
        work_dir: Working directory for analysis
        species_name: Species name (e.g., "Homo_sapiens")

    Returns:
        Path object ending with .idx in expected index directory.
        Note: Directory may not exist yet.

    Example:
        >>> from pathlib import Path
        >>> index_path = get_expected_index_path(
        ...     Path("output/amalgkit"), "Homo_sapiens"
        ... )
        >>> assert str(index_path).endswith(".idx")
        >>> assert "index" in str(index_path)
        >>> assert "Homo_sapiens" in str(index_path)
    """
    work_dir = Path(work_dir)

    # Standard path structure: work_dir/species_name/index/species_name.idx
    index_dir = work_dir / species_name / "index"
    index_path = index_dir / f"{species_name}.idx"

    logger.debug(f"Expected index path: {index_path}")
    return index_path


def download_rna_fasta_from_ftp(
    accession: str, output_dir: Path, **kwargs: Any
) -> Optional[Path]:
    """Download RNA sequences from NCBI FTP server.

    Retrieves RNA FASTA file from NCBI FTP based on genome accession.
    Supports both RefSeq (GCF_) and GenBank (GCA_) accessions.

    Args:
        accession: Genome assembly accession
        output_dir: Directory to save downloaded file
        **kwargs: Additional options (e.g., timeout, retry_count)

    Returns:
        Path to downloaded file if successful, None otherwise.

    Raises:
        No exceptions raised - returns None on failure for graceful handling.
    """
    # Placeholder implementation with proper docstring
    logger.debug(f"Would download RNA FASTA for {accession}")
    return None


def download_cds_fasta_from_ftp(
    accession: str, output_dir: Path, **kwargs: Any
) -> Optional[Path]:
    """Download CDS sequences from NCBI FTP server.

    Retrieves CDS (Coding DNA Sequence) FASTA file from NCBI FTP.
    Useful for proteomics integration and coding region analysis.

    Args:
        accession: Genome assembly accession
        output_dir: Directory to save downloaded file
        **kwargs: Additional options

    Returns:
        Path to downloaded file if successful, None otherwise.
    """
    # Placeholder implementation with proper docstring
    logger.debug(f"Would download CDS FASTA for {accession}")
    return None


def extract_transcripts_from_gff(
    gff_file: Path, genome_fasta: Path, output_file: Path, **kwargs: Any
) -> Optional[Path]:
    """Extract transcript sequences from GFF3 annotation and genome sequence.

    Uses GFF3 coordinates to extract transcript regions from genome FASTA.
    Handles both forward and reverse strand sequences.

    Args:
        gff_file: Path to GFF3 annotation file
        genome_fasta: Path to genome FASTA file
        output_file: Path to write extracted transcript sequences
        **kwargs: Additional options (e.g., feature_type, min_length)

    Returns:
        Path to output FASTA file if successful, None otherwise.
    """
    # Placeholder implementation with proper docstring
    logger.debug(f"Would extract transcripts from {gff_file}")
    return None


def prepare_transcriptome_for_kallisto(
    transcriptome_fasta: Path, output_fasta: Path, **kwargs: Any
) -> Optional[Path]:
    """Prepare transcriptome FASTA for Kallisto indexing.

    Performs quality checks and formatting on transcriptome sequences.
    Ensures sequences meet Kallisto requirements (length > k, valid nucleotides).

    Args:
        transcriptome_fasta: Input transcriptome FASTA file
        output_fasta: Path to write prepared FASTA
        **kwargs: Additional options (e.g., min_length, kmer_length)

    Returns:
        Path to prepared FASTA file if successful, None otherwise.
    """
    # Placeholder implementation with proper docstring
    logger.debug(f"Would prepare transcriptome for Kallisto")
    return None


def build_kallisto_index(
    transcriptome_fasta: Path, index_path: Path, **kwargs: Any
) -> Optional[Path]:
    """Build Kallisto index from transcriptome FASTA.

    Executes kallisto index command to create k-mer index for fast quantification.

    Args:
        transcriptome_fasta: Prepared transcriptome FASTA file
        index_path: Output path for Kallisto index
        **kwargs: Additional options (e.g., kmer_size, force)

    Returns:
        Path to created index file if successful, None otherwise.

    Raises:
        RuntimeError: If Kallisto command fails with non-zero exit code.
    """
    import subprocess
    import shutil
    
    transcriptome_fasta = Path(transcriptome_fasta)
    index_path = Path(index_path)
    
    if not transcriptome_fasta.exists():
        logger.error(f"Transcriptome FASTA not found: {transcriptome_fasta}")
        return None
    
    # Check if kallisto is available
    kallisto_cmd = shutil.which("kallisto")
    if not kallisto_cmd:
        logger.error("kallisto not found in PATH. Install with: conda install -c bioconda kallisto")
        return None
    
    # Get kmer size from kwargs (default: 31)
    kmer_size = kwargs.get("kmer_size", 31)
    
    # Create index directory
    index_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Build index
    logger.info(f"Building kallisto index: {index_path}")
    logger.info(f"  Transcriptome: {transcriptome_fasta}")
    logger.info(f"  k-mer size: {kmer_size}")
    
    cmd = [
        kallisto_cmd,
        "index",
        "-i", str(index_path),
        "-k", str(kmer_size),
        str(transcriptome_fasta)
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        logger.info(f"Kallisto index built successfully: {index_path}")
        if result.stdout:
            logger.debug(f"kallisto output: {result.stdout}")
        return index_path
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to build kallisto index: {e}")
        if e.stderr:
            logger.error(f"kallisto error: {e.stderr}")
        return None
    except Exception as e:
        logger.error(f"Error building kallisto index: {e}")
        return None


def prepare_genome_for_quantification(
    accession: str,
    work_dir: Path,
    species_name: str,
    **kwargs: Any
) -> Dict[str, Any]:
    """End-to-end genome preparation for RNA-seq quantification.

    Orchestrates complete workflow: download → extract → prepare → index.
    Returns structured results about prepared reference.

    Args:
        accession: Genome assembly accession
        work_dir: Working directory for all outputs
        species_name: Species name for naming conventions
        **kwargs: Additional options passed to individual steps

    Returns:
        Dictionary with keys:
        - 'success' (bool): Whether preparation succeeded
        - 'index_path' (Path or None): Path to Kallisto index
        - 'transcriptome_fasta' (Path or None): Path to transcriptome
        - 'errors' (list[str]): Error messages if any
    """
    # Placeholder implementation with proper docstring
    logger.debug(f"Would prepare genome for quantification: {accession}")
    return {"success": False, "index_path": None, "transcriptome_fasta": None, "errors": []}


def verify_genome_status(work_dir: Path, species_name: str) -> Dict[str, Any]:
    """Check genome preparation status and readiness.

    Verifies that required files exist and are valid for RNA-seq quantification.

    Args:
        work_dir: Working directory containing prepared genomes
        species_name: Species name to check

    Returns:
        Dictionary with keys:
        - 'prepared' (bool): Whether genome is fully prepared
        - 'index_exists' (bool): Whether Kallisto index exists
        - 'transcriptome_exists' (bool): Whether transcriptome FASTA exists
        - 'issues' (list[str]): Any problems detected
    """
    # Placeholder implementation with proper docstring
    logger.debug(f"Would verify genome status for {species_name}")
    return {
        "prepared": False,
        "index_exists": False,
        "transcriptome_exists": False,
        "issues": [],
    }


def orchestrate_genome_setup(
    config: Dict[str, Any],
    species_name: str,
    work_dir: Path,
    **kwargs: Any
) -> Dict[str, Any]:
    """Coordinate all genome preparation steps.

    High-level orchestration function that manages the complete genome setup
    workflow based on configuration. Logs progress and handles errors gracefully.

    Args:
        config: Configuration dictionary with genome/accession info
        species_name: Species name for reference
        work_dir: Working directory for all outputs
        **kwargs: Additional workflow options

    Returns:
        Dictionary with comprehensive results:
        - 'success' (bool): Overall workflow success
        - 'species' (str): Species name processed
        - 'results' (dict): Individual step results
        - 'messages' (list[str]): Status messages and warnings
    """
    # Placeholder implementation with proper docstring
    logger.info(f"Orchestrating genome setup for {species_name}")
    return {"success": False, "species": species_name, "results": {}, "messages": []}
