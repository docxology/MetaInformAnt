"""Genome preparation utilities for RNA-seq workflows.

This module provides functions to:
- Extract transcriptome FASTA files from downloaded genome packages
- Prepare FASTA files for kallisto indexing
- Build kallisto indexes from transcriptome FASTA files
- Download RNA/CDS files directly from FTP when package download fails
- Extract transcripts from GFF files using gffread
"""

from __future__ import annotations

import gzip
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Any
from urllib.request import urlopen

from ..core.io import ensure_directory

logger = logging.getLogger(__name__)


def find_rna_fasta_in_genome_dir(genome_dir: Path, accession: str) -> Path | None:
    """Find RNA FASTA file in extracted genome directory.
    
    NCBI datasets extracts to various patterns:
    - ncbi_dataset_api_extracted/ncbi_dataset/data/{accession}/rna.fna
    - ncbi_dataset_api_extracted/ncbi_dataset/data/{accession}/rna.fna.gz
    - ncbi_dataset_extracted/ncbi_dataset/data/{accession}/rna.fna
    
    Also searches recursively for RNA FASTA files and checks for direct FTP downloads.
    
    Args:
        genome_dir: Base genome directory
        accession: NCBI assembly accession
        
    Returns:
        Path to RNA FASTA file if found, None otherwise
    """
    # Common extraction patterns
    patterns = [
        # API extraction pattern
        genome_dir / "ncbi_dataset_api_extracted" / "ncbi_dataset" / "data" / accession / "rna.fna",
        genome_dir / "ncbi_dataset_api_extracted" / "ncbi_dataset" / "data" / accession / "rna.fna.gz",
        genome_dir / "ncbi_dataset_api_extracted" / "ncbi_dataset" / "data" / accession / "transcript_sequences.fna",
        # CLI extraction pattern
        genome_dir / "ncbi_dataset_extracted" / "ncbi_dataset" / "data" / accession / "rna.fna",
        genome_dir / "ncbi_dataset_extracted" / "ncbi_dataset" / "data" / accession / "rna.fna.gz",
        # Direct file patterns (including FTP downloads)
        genome_dir / "rna.fna",
        genome_dir / "rna.fna.gz",
        genome_dir / f"{accession}_rna_from_genomic.fna",
        genome_dir / f"{accession}_rna_from_genomic.fna.gz",
    ]
    
    for pattern in patterns:
        if pattern.exists() and pattern.is_file():
            return pattern
    
    # Search recursively for RNA FASTA files with more patterns
    search_patterns = [
        "rna.fna",
        "rna.fna.gz",
        "transcript_sequences.fna",
        "*_rna_from_genomic.fna.gz",
        "*_rna_from_genomic.fna",
        "*rna*.fna*",
        "*transcript*.fna*",
    ]
    
    for pattern in search_patterns:
        matches = list(genome_dir.rglob(pattern))
        if matches:
            # Prefer unzipped files
            unzipped = [m for m in matches if not m.name.endswith(".gz")]
            if unzipped:
                return unzipped[0]
            return matches[0]
    
    return None


def download_rna_fasta_from_ftp(
    ftp_url: str,
    genome_dir: Path,
    accession: str,
    assembly_name: str | None = None,
    config: dict[str, Any] | None = None,
) -> Path | None:
    """Download RNA FASTA file directly from NCBI FTP.
    
    Args:
        ftp_url: Base FTP URL for the genome
        genome_dir: Directory to save the file
        accession: NCBI assembly accession
        assembly_name: Optional assembly name (for filename construction)
        config: Optional config dict (may contain genome.files.transcriptome_fasta)
        
    Returns:
        Path to downloaded file if successful, None otherwise
    """
    ensure_directory(genome_dir)
    
    # First try filename from config if available
    filename_patterns = []
    if config:
        genome_config = config.get("genome", {})
        files_config = genome_config.get("files", {})
        transcriptome_fasta = files_config.get("transcriptome_fasta")
        if transcriptome_fasta:
            filename_patterns.append(transcriptome_fasta)
    
    # Then try common RNA FASTA filename patterns
    if assembly_name:
        filename_patterns.extend([
            f"{accession}_{assembly_name}_rna_from_genomic.fna.gz",
            f"{accession}_{assembly_name}_rna.fna.gz",
            f"{accession}_{assembly_name}_transcript.fna.gz",
        ])
    else:
        filename_patterns.extend([
            f"{accession}_rna_from_genomic.fna.gz",
            f"{accession}_rna.fna.gz",
            f"{accession}_transcript.fna.gz",
        ])
    
    base_url = ftp_url.rstrip("/")
    
    for filename in filename_patterns:
        url = f"{base_url}/{filename}"
        output_path = genome_dir / filename
        
        # Skip if already exists
        if output_path.exists() and output_path.stat().st_size > 0:
            logger.info(f"RNA FASTA already exists: {output_path.name}")
            return output_path
        
        try:
            logger.info(f"Downloading RNA FASTA from FTP: {url}")
            req = urlopen(url, timeout=60)
            try:
                # Check HTTP status code
                status = getattr(req, 'status', None) or getattr(req, 'code', None)
                if status and status != 200:
                    logger.warning(f"HTTP {status} for {url}")
                    req.close()
                    continue
                
                # Download the file
                with open(output_path, "wb") as f:
                    shutil.copyfileobj(req, f)
                
                # Verify file was downloaded and has content
                if output_path.exists() and output_path.stat().st_size > 0:
                    logger.info(f"Downloaded RNA FASTA: {output_path.name} ({output_path.stat().st_size} bytes)")
                    return output_path
                else:
                    logger.warning(f"Downloaded file is empty or missing: {output_path}")
            finally:
                req.close()
        except Exception as e:
            logger.warning(f"Failed to download {url}: {e}")
            continue
    
    return None


def download_cds_fasta_from_ftp(
    ftp_url: str,
    genome_dir: Path,
    accession: str,
    assembly_name: str | None = None,
    config: dict[str, Any] | None = None,
) -> Path | None:
    """Download CDS FASTA file directly from NCBI FTP.
    
    Args:
        ftp_url: Base FTP URL for the genome
        genome_dir: Directory to save the file
        accession: NCBI assembly accession
        assembly_name: Optional assembly name (for filename construction)
        config: Optional config dict (may contain genome.files.cds_fasta)
        
    Returns:
        Path to downloaded file if successful, None otherwise
    """
    ensure_directory(genome_dir)
    
    # First try filename from config if available
    filename_patterns = []
    if config:
        genome_config = config.get("genome", {})
        files_config = genome_config.get("files", {})
        cds_fasta = files_config.get("cds_fasta")
        if cds_fasta:
            filename_patterns.append(cds_fasta)
    
    # Then try common CDS FASTA filename patterns
    if assembly_name:
        filename_patterns.extend([
            f"{accession}_{assembly_name}_cds_from_genomic.fna.gz",
            f"{accession}_{assembly_name}_cds.fna.gz",
        ])
    else:
        filename_patterns.extend([
            f"{accession}_cds_from_genomic.fna.gz",
            f"{accession}_cds.fna.gz",
        ])
    
    base_url = ftp_url.rstrip("/")
    
    for filename in filename_patterns:
        url = f"{base_url}/{filename}"
        output_path = genome_dir / filename
        
        # Skip if already exists
        if output_path.exists() and output_path.stat().st_size > 0:
            logger.info(f"CDS FASTA already exists: {output_path.name}")
            return output_path
        
        try:
            logger.info(f"Downloading CDS FASTA from FTP: {url}")
            req = urlopen(url, timeout=60)
            try:
                # Check HTTP status code
                status = getattr(req, 'status', None) or getattr(req, 'code', None)
                if status and status != 200:
                    logger.warning(f"HTTP {status} for {url}")
                    req.close()
                    continue
                
                # Download the file
                with open(output_path, "wb") as f:
                    shutil.copyfileobj(req, f)
                
                # Verify file was downloaded and has content
                if output_path.exists() and output_path.stat().st_size > 0:
                    logger.info(f"Downloaded CDS FASTA: {output_path.name} ({output_path.stat().st_size} bytes)")
                    return output_path
                else:
                    logger.warning(f"Downloaded file is empty or missing: {output_path}")
            finally:
                req.close()
        except Exception as e:
            logger.warning(f"Failed to download {url}: {e}")
            continue
    
    return None


def extract_transcripts_from_gff(
    gff_path: Path,
    genome_fasta_path: Path,
    output_fasta_path: Path,
) -> bool:
    """Extract transcript sequences from GFF file using gffread.
    
    Args:
        gff_path: Path to GFF/GFF3 annotation file
        genome_fasta_path: Path to genomic FASTA file
        output_fasta_path: Path to output transcript FASTA file
        
    Returns:
        True if extraction successful, False otherwise
    """
    gffread_exe = shutil.which("gffread")
    if not gffread_exe:
        logger.warning("gffread not found on PATH - cannot extract transcripts from GFF")
        return False
    
    if not gff_path.exists():
        logger.warning(f"GFF file not found: {gff_path}")
        return False
    
    if not genome_fasta_path.exists():
        logger.warning(f"Genome FASTA not found: {genome_fasta_path}")
        return False
    
    ensure_directory(output_fasta_path.parent)
    
    # Handle gzipped GFF and genome FASTA files
    import tempfile
    
    gff_input = gff_path
    genome_input = genome_fasta_path
    temp_files = []
    
    if gff_path.name.endswith(".gz"):
        # Decompress GFF for gffread compatibility
        with tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix='.gff') as tmp:
            with gzip.open(gff_path, 'rb') as f_in:
                shutil.copyfileobj(f_in, tmp)
            gff_input = Path(tmp.name)
            temp_files.append(gff_input)
    
    if genome_fasta_path.name.endswith(".gz"):
        # Decompress genome FASTA for gffread compatibility
        with tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix='.fna') as tmp:
            with gzip.open(genome_fasta_path, 'rb') as f_in:
                shutil.copyfileobj(f_in, tmp)
            genome_input = Path(tmp.name)
            temp_files.append(genome_input)
    
    cmd = [
        gffread_exe,
        "-g", str(genome_input),
        "-w", str(output_fasta_path),
        str(gff_input),
    ]
    
    try:
        logger.info(f"Extracting transcripts from GFF using gffread...")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        logger.info(f"Extracted transcripts to: {output_fasta_path}")
        
        # Clean up temporary files if created
        for tmp_file in temp_files:
            if tmp_file.exists():
                tmp_file.unlink()
        
        return output_fasta_path.exists() and output_fasta_path.stat().st_size > 0
    except subprocess.CalledProcessError as e:
        logger.error(f"gffread failed: {e.stderr}")
        # Clean up temporary files if created
        for tmp_file in temp_files:
            if tmp_file.exists():
                tmp_file.unlink()
        return False
    except Exception as e:
        logger.error(f"Error running gffread: {e}")
        # Clean up temporary files if created
        for tmp_file in temp_files:
            if tmp_file.exists():
                tmp_file.unlink()
        return False


def prepare_transcriptome_for_kallisto(
    genome_dir: Path,
    species_name: str,
    work_dir: Path,
    *,
    accession: str | None = None,
    use_cds_fallback: bool = True,
    ftp_url: str | None = None,
    assembly_name: str | None = None,
    config: dict[str, Any] | None = None,
) -> Path | None:
    """Prepare transcriptome FASTA file for kallisto indexing.
    
    Finds the RNA FASTA file in the extracted genome directory and
    copies it to the expected location for amalgkit/kallisto.
    If RNA FASTA is not found and use_cds_fallback is True, uses CDS sequences as fallback.
    If still not found, attempts to download from FTP or extract from GFF.
    
    Expected output: work_dir/fasta/{Species_Name}_rna.fasta
    
    Args:
        genome_dir: Directory containing extracted genome package
        species_name: Species name (with underscores, e.g., "Camponotus_floridanus")
        work_dir: Work directory for amalgkit workflow
        accession: Optional NCBI assembly accession (for searching)
        use_cds_fallback: If True, use CDS sequences when RNA FASTA is not available
        ftp_url: Optional FTP URL for direct file download
        assembly_name: Optional assembly name (for FTP filename construction)
        config: Optional config dict (may contain genome.files.transcriptome_fasta, etc.)
        
    Returns:
        Path to prepared FASTA file if successful, None otherwise
    """
    # Find RNA FASTA in genome directory
    if accession:
        rna_fasta = find_rna_fasta_in_genome_dir(genome_dir, accession)
    else:
        # Try to find any RNA FASTA file
        matches = list(genome_dir.rglob("rna.fna*"))
        if matches:
            rna_fasta = matches[0]
        else:
            rna_fasta = None
    
    # Fallback to CDS if RNA FASTA not found
    if (not rna_fasta or not rna_fasta.exists()) and use_cds_fallback:
        logger.info(f"RNA FASTA not found, searching for CDS sequences as fallback...")
        cds_patterns = [
            "*cds*.fna*",
            "*cds*.fa*",
            "*CDS*.fna*",
            "*CDS*.fa*",
        ]
        cds_matches = []
        for pattern in cds_patterns:
            cds_matches.extend(genome_dir.rglob(pattern))
        
        # Remove duplicates
        cds_matches = list(set(cds_matches))
        
        if cds_matches:
            # Prefer unzipped files
            unzipped = [f for f in cds_matches if not f.name.endswith('.gz')]
            if unzipped:
                rna_fasta = unzipped[0]
                logger.info(f"Using CDS file as transcriptome: {rna_fasta.name}")
            else:
                rna_fasta = cds_matches[0]
                logger.info(f"Using gzipped CDS file as transcriptome: {rna_fasta.name}")
    
    # Try FTP download if still not found
    if (not rna_fasta or not rna_fasta.exists()) and ftp_url and accession:
        logger.info(f"Attempting to download RNA FASTA from FTP...")
        downloaded_rna = download_rna_fasta_from_ftp(ftp_url, genome_dir, accession, assembly_name, config=config)
        if downloaded_rna and downloaded_rna.exists():
            rna_fasta = downloaded_rna
            # Re-search for it in case it was downloaded with a different name
            if accession:
                rna_fasta = find_rna_fasta_in_genome_dir(genome_dir, accession) or rna_fasta
    
    # Try CDS FTP download if RNA still not found
    if (not rna_fasta or not rna_fasta.exists()) and use_cds_fallback and ftp_url and accession:
        logger.info(f"Attempting to download CDS FASTA from FTP as fallback...")
        downloaded_cds = download_cds_fasta_from_ftp(ftp_url, genome_dir, accession, assembly_name, config=config)
        if downloaded_cds and downloaded_cds.exists():
            rna_fasta = downloaded_cds
            # Re-search for CDS files
            cds_patterns = ["*cds*.fna*", "*cds*.fa*", "*CDS*.fna*", "*CDS*.fa*"]
            cds_matches = []
            for pattern in cds_patterns:
                cds_matches.extend(genome_dir.rglob(pattern))
            if cds_matches:
                unzipped = [f for f in cds_matches if not f.name.endswith('.gz')]
                rna_fasta = unzipped[0] if unzipped else cds_matches[0]
    
    # Try GFF extraction as last resort
    if (not rna_fasta or not rna_fasta.exists()) and config:
        genome_config = config.get("genome", {})
        files_config = genome_config.get("files", {})
        gff_file = files_config.get("annotation_gff") or files_config.get("annotation_gff3")
        genomic_fasta_file = files_config.get("genomic_fasta")
        
        if gff_file and genomic_fasta_file and ftp_url:
            logger.info(f"Attempting to extract transcripts from GFF file...")
            gff_path = genome_dir / gff_file
            genome_fasta_path = genome_dir / genomic_fasta_file
            
            # Try to download GFF and genome FASTA if not present
            if not gff_path.exists() and ftp_url:
                base_url = ftp_url.rstrip("/")
                gff_url = f"{base_url}/{gff_file}"
                try:
                    logger.info(f"Downloading GFF from FTP: {gff_url}")
                    req = urlopen(gff_url, timeout=60)
                    try:
                        status = getattr(req, 'status', None) or getattr(req, 'code', None)
                        if status and status != 200:
                            logger.warning(f"HTTP {status} for GFF download")
                        else:
                            with open(gff_path, "wb") as f:
                                shutil.copyfileobj(req, f)
                            if gff_path.exists() and gff_path.stat().st_size > 0:
                                logger.info(f"Downloaded GFF: {gff_path.name} ({gff_path.stat().st_size} bytes)")
                            else:
                                logger.warning(f"GFF download failed or file is empty")
                    finally:
                        req.close()
                except Exception as e:
                    logger.warning(f"Failed to download GFF: {e}")
            
            if not genome_fasta_path.exists() and ftp_url:
                base_url = ftp_url.rstrip("/")
                genome_url = f"{base_url}/{genomic_fasta_file}"
                try:
                    logger.info(f"Downloading genome FASTA from FTP: {genome_url}")
                    req = urlopen(genome_url, timeout=60)
                    try:
                        status = getattr(req, 'status', None) or getattr(req, 'code', None)
                        if status and status != 200:
                            logger.warning(f"HTTP {status} for genome FASTA download")
                        else:
                            with open(genome_fasta_path, "wb") as f:
                                shutil.copyfileobj(req, f)
                            if genome_fasta_path.exists() and genome_fasta_path.stat().st_size > 0:
                                logger.info(f"Downloaded genome FASTA: {genome_fasta_path.name} ({genome_fasta_path.stat().st_size} bytes)")
                            else:
                                logger.warning(f"Genome FASTA download failed or file is empty")
                    finally:
                        req.close()
                except Exception as e:
                    logger.warning(f"Failed to download genome FASTA: {e}")
            
            # Extract transcripts if we have both files
            if gff_path.exists() and genome_fasta_path.exists():
                extracted_fasta = genome_dir / f"{accession}_transcripts_from_gff.fna"
                if extract_transcripts_from_gff(gff_path, genome_fasta_path, extracted_fasta):
                    rna_fasta = extracted_fasta
    
    if not rna_fasta or not rna_fasta.exists():
        logger.warning(f"RNA FASTA and CDS not found in {genome_dir}")
        logger.warning(f"  Tried: RNA FASTA search, CDS fallback, FTP download")
        
        # Check if GFF extraction was attempted but failed due to missing gffread
        if config:
            genome_config = config.get("genome", {})
            files_config = genome_config.get("files", {})
            gff_file = files_config.get("annotation_gff") or files_config.get("annotation_gff3")
            if gff_file:
                gffread_exe = shutil.which("gffread")
                if not gffread_exe:
                    logger.warning(f"  GFF extraction not attempted: gffread not found on PATH")
                    logger.warning(f"  Install gffread (from Cufflinks package) to enable GFF-based transcript extraction")
                else:
                    logger.warning(f"  GFF extraction attempted but failed")
        
        return None
    
    # Prepare output directory
    fasta_dir = work_dir / "fasta"
    ensure_directory(fasta_dir)
    
    # Prepare output filename
    # Replace spaces/underscores and ensure .fasta extension
    output_name = species_name.replace(" ", "_") + "_rna.fasta"
    output_path = fasta_dir / output_name
    
    # If source is gzipped, decompress; otherwise copy
    if rna_fasta.name.endswith(".gz"):
        logger.info(f"Decompressing {rna_fasta} to {output_path}")
        try:
            with gzip.open(rna_fasta, "rb") as f_in:
                with open(output_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
        except Exception as e:
            logger.error(f"Failed to decompress RNA FASTA: {e}")
            return None
    else:
        logger.info(f"Copying {rna_fasta} to {output_path}")
        try:
            shutil.copy2(rna_fasta, output_path)
        except Exception as e:
            logger.error(f"Failed to copy RNA FASTA: {e}")
            return None
    
    # Verify output file exists and has content
    if output_path.exists() and output_path.stat().st_size > 0:
        logger.info(f"Prepared transcriptome FASTA: {output_path}")
        return output_path
    else:
        logger.error(f"Output file invalid or empty: {output_path}")
        return None


def build_kallisto_index(
    fasta_path: Path,
    index_path: Path,
    *,
    kmer_size: int = 31,
    check_existing: bool = True,
) -> bool:
    """Build kallisto index from transcriptome FASTA file.
    
    Args:
        fasta_path: Path to transcriptome FASTA file
        index_path: Path where index should be created
        kmer_size: k-mer size for index (default: 31)
        check_existing: If True, skip if index already exists
        
    Returns:
        True if index was built successfully or already exists, False otherwise
    """
    # Check if index already exists
    if check_existing and index_path.exists():
        logger.info(f"Kallisto index already exists: {index_path}")
        return True
    
    # Check if kallisto is available
    kallisto_exe = shutil.which("kallisto")
    if not kallisto_exe:
        logger.error("kallisto not found on PATH")
        return False
    
    # Verify FASTA file exists
    if not fasta_path.exists():
        logger.error(f"FASTA file not found: {fasta_path}")
        return False
    
    # Ensure index directory exists
    ensure_directory(index_path.parent)
    
    # Build index
    cmd = [
        kallisto_exe,
        "index",
        "-i",
        str(index_path),
        "-k",
        str(kmer_size),
        str(fasta_path),
    ]
    
    logger.info(f"Building kallisto index: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        logger.info(f"Kallisto index built successfully: {index_path}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to build kallisto index: {e.stderr}")
        return False
    except Exception as e:
        logger.error(f"Error building kallisto index: {e}")
        return False


def get_expected_index_path(work_dir: Path, species_name: str) -> Path:
    """Get expected kallisto index path for a species.
    
    Args:
        work_dir: Work directory for amalgkit workflow
        species_name: Species name (with underscores)
        
    Returns:
        Expected path to kallisto index file
    """
    index_dir = work_dir / "index"
    index_filename = f"{species_name}_transcripts.idx"
    return index_dir / index_filename


def prepare_genome_for_quantification(
    genome_dir: Path,
    species_name: str,
    work_dir: Path,
    *,
    accession: str | None = None,
    build_index: bool = True,
    kmer_size: int = 31,
) -> dict[str, Any]:
    """Complete genome preparation for quantification.
    
    This function orchestrates the full process:
    1. Extract transcriptome FASTA from genome package
    2. Prepare FASTA in expected location
    3. Build kallisto index (if requested)
    
    Args:
        genome_dir: Directory containing extracted genome package
        species_name: Species name (with underscores)
        work_dir: Work directory for amalgkit workflow
        accession: Optional NCBI assembly accession
        build_index: If True, build kallisto index after preparing FASTA
        kmer_size: k-mer size for kallisto index
        
    Returns:
        Dictionary with preparation results:
        - success: bool
        - fasta_path: Path | None
        - index_path: Path | None
        - error: str | None
    """
    result: dict[str, Any] = {
        "success": False,
        "fasta_path": None,
        "index_path": None,
        "error": None,
    }
    
    # Step 1: Prepare transcriptome FASTA
    fasta_path = prepare_transcriptome_for_kallisto(
        genome_dir,
        species_name,
        work_dir,
        accession=accession,
    )
    
    if not fasta_path:
        result["error"] = "Failed to prepare transcriptome FASTA"
        return result
    
    result["fasta_path"] = str(fasta_path)
    
    # Step 2: Build kallisto index (if requested)
    if build_index:
        index_path = get_expected_index_path(work_dir, species_name)
        success = build_kallisto_index(
            fasta_path,
            index_path,
            kmer_size=kmer_size,
        )
        
        if success:
            result["index_path"] = str(index_path)
            result["success"] = True
        else:
            result["error"] = "Failed to build kallisto index"
    else:
        result["success"] = True
    
    return result


def verify_genome_status(
    config_path: Path,
    *,
    repo_root: Path | None = None,
) -> dict[str, Any]:
    """Verify genome and index status for a species.

    Args:
        config_path: Path to species config file
        repo_root: Repository root directory (optional)

    Returns:
        Dictionary with verification results
    """
    import json
    import yaml

    from ..core.config import load_mapping_from_file

    if repo_root is None:
        repo_root = Path.cwd()

    try:
        config = load_mapping_from_file(config_path)
    except Exception as e:
        return {
            "config_file": str(config_path),
            "error": f"Failed to load config: {e}",
            "genome_downloaded": False,
            "rna_fasta_found": False,
            "kallisto_index_found": False,
        }

    species_list = config.get("species_list", [])
    species_name = species_list[0] if species_list else "unknown"
    genome = config.get("genome", {})

    result: dict[str, Any] = {
        "config_file": str(config_path),
        "species_name": species_name,
        "genome_downloaded": False,
        "rna_fasta_found": False,
        "rna_fasta_path": None,
        "kallisto_index_found": False,
        "kallisto_index_path": None,
        "genome_accession": None,
        "work_dir": None,
    }

    if not genome:
        result["error"] = "No genome configuration found"
        return result

    accession = genome.get("accession")
    result["genome_accession"] = accession

    if not accession:
        result["error"] = "No genome accession in config"
        return result

    # Check genome download
    dest_dir_str = genome.get("dest_dir", "")
    if not dest_dir_str:
        work_dir_str = config.get("work_dir", "")
        if work_dir_str:
            work_dir_path = Path(work_dir_str).expanduser()
            if not work_dir_path.is_absolute():
                work_dir_path = repo_root / work_dir_path
            dest_dir = work_dir_path.parent / "genome"
        else:
            result["error"] = "Cannot determine genome directory"
            return result
    else:
        dest_dir = Path(dest_dir_str).expanduser()
        if not dest_dir.is_absolute():
            dest_dir = repo_root / dest_dir

    result["genome_dir"] = str(dest_dir)

    # Check download record
    download_record = dest_dir / "download_record.json"
    if download_record.exists():
        try:
            with open(download_record, "r") as f:
                record = json.load(f)
                if record.get("return_code") == 0:
                    result["genome_downloaded"] = True
        except Exception:
            pass

    # Check if extracted directory exists
    extracted_dirs = [
        dest_dir / "ncbi_dataset_api_extracted",
        dest_dir / "ncbi_dataset_extracted",
    ]
    if any(d.exists() for d in extracted_dirs):
        result["genome_downloaded"] = True

    # Find RNA FASTA
    rna_fasta = find_rna_fasta_in_genome_dir(dest_dir, accession)
    if rna_fasta:
        result["rna_fasta_found"] = True
        result["rna_fasta_path"] = str(rna_fasta)

    # Check kallisto index
    work_dir_str = config.get("work_dir", "")
    if work_dir_str:
        work_dir_path = Path(work_dir_str).expanduser()
        if not work_dir_path.is_absolute():
            work_dir_path = repo_root / work_dir_path
        result["work_dir"] = str(work_dir_path)

        index_path = get_expected_index_path(work_dir_path, species_name.replace(" ", "_"))
        if index_path.exists():
            result["kallisto_index_found"] = True
            result["kallisto_index_path"] = str(index_path)

    return result


def orchestrate_genome_setup(
    config_path: Path,
    *,
    repo_root: Path | None = None,
    skip_verify_initial: bool = False,
    skip_download: bool = False,
    skip_prepare: bool = False,
    skip_build: bool = False,
    skip_verify_final: bool = False,
    kmer_size: int = 31,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Orchestrate complete genome setup pipeline.

    Args:
        config_path: Path to species config file
        repo_root: Repository root directory (optional)
        skip_verify_initial: Skip initial verification
        skip_download: Skip genome download
        skip_prepare: Skip transcriptome preparation
        skip_build: Skip kallisto index building
        skip_verify_final: Skip final verification
        kmer_size: k-mer size for kallisto index
        dry_run: If True, only report what would be done

    Returns:
        Dictionary with orchestration results
    """
    import yaml

    from ..core.config import load_mapping_from_file

    if repo_root is None:
        repo_root = Path.cwd()

    results: dict[str, Any] = {
        "success": True,
        "steps_completed": [],
        "steps_failed": [],
        "verification_initial": None,
        "download": None,
        "prepare": None,
        "build": None,
        "verification_final": None,
    }

    config = load_mapping_from_file(config_path)
    species_list = config.get("species_list", [])
    species_name = species_list[0] if species_list else "unknown"
    genome = config.get("genome", {})

    # Step 1: Initial verification
    if not skip_verify_initial:
        logger.info("Step 1: Initial Verification")
        verification = verify_genome_status(config_path, repo_root=repo_root)
        results["verification_initial"] = verification
        if verification.get("error"):
            logger.warning(f"Initial verification warning: {verification.get('error')}")

    # Step 2: Download genome
    if not skip_download:
        logger.info("Step 2: Download Missing Genomes")
        # This would call download_genome_for_species if we had it
        # For now, we'll note that it should be called
        if not dry_run:
            logger.info("Genome download should be performed via download_genome_for_species()")
        results["download"] = {"status": "skipped" if dry_run else "pending"}

    # Step 3: Prepare transcriptomes
    if not skip_prepare:
        logger.info("Step 3: Prepare Transcriptomes")
        if not dry_run:
            # Get paths
            work_dir_str = config.get("work_dir", "")
            if work_dir_str:
                work_dir_path = Path(work_dir_str).expanduser()
                if not work_dir_path.is_absolute():
                    work_dir_path = repo_root / work_dir_path
            else:
                logger.error("Cannot determine work_dir")
                results["success"] = False
                results["steps_failed"].append("prepare")
                return results

            dest_dir_str = genome.get("dest_dir", "")
            if not dest_dir_str:
                dest_dir = work_dir_path.parent / "genome"
            else:
                dest_dir = Path(dest_dir_str).expanduser()
                if not dest_dir.is_absolute():
                    dest_dir = repo_root / dest_dir

            accession = genome.get("accession")
            fasta_path = prepare_transcriptome_for_kallisto(
                dest_dir,
                species_name.replace(" ", "_"),
                work_dir_path,
                accession=accession,
            )

            if fasta_path:
                results["prepare"] = {"success": True, "fasta_path": str(fasta_path)}
                results["steps_completed"].append("prepare")
            else:
                results["prepare"] = {"success": False, "error": "Failed to prepare transcriptome"}
                results["steps_failed"].append("prepare")
                results["success"] = False
        else:
            results["prepare"] = {"status": "dry_run"}

    # Step 4: Build kallisto index
    if not skip_build:
        logger.info("Step 4: Build Kallisto Indexes")
        if not dry_run:
            work_dir_str = config.get("work_dir", "")
            if work_dir_str:
                work_dir_path = Path(work_dir_str).expanduser()
                if not work_dir_path.is_absolute():
                    work_dir_path = repo_root / work_dir_path
            else:
                logger.error("Cannot determine work_dir")
                results["success"] = False
                results["steps_failed"].append("build")
                return results

            species_name_underscore = species_name.replace(" ", "_")
            fasta_path = work_dir_path / "fasta" / f"{species_name_underscore}_rna.fasta"
            index_path = get_expected_index_path(work_dir_path, species_name_underscore)

            if fasta_path.exists():
                success = build_kallisto_index(fasta_path, index_path, kmer_size=kmer_size)
                if success:
                    results["build"] = {"success": True, "index_path": str(index_path)}
                    results["steps_completed"].append("build")
                else:
                    results["build"] = {"success": False, "error": "Failed to build index"}
                    results["steps_failed"].append("build")
                    results["success"] = False
            else:
                results["build"] = {"success": False, "error": f"FASTA not found: {fasta_path}"}
                results["steps_failed"].append("build")
                results["success"] = False
        else:
            results["build"] = {"status": "dry_run"}

    # Step 5: Final verification
    if not skip_verify_final:
        logger.info("Step 5: Final Verification")
        verification = verify_genome_status(config_path, repo_root=repo_root)
        results["verification_final"] = verification

    return results


__all__ = [
    "find_rna_fasta_in_genome_dir",
    "download_rna_fasta_from_ftp",
    "download_cds_fasta_from_ftp",
    "extract_transcripts_from_gff",
    "prepare_transcriptome_for_kallisto",
    "build_kallisto_index",
    "get_expected_index_path",
    "prepare_genome_for_quantification",
    "verify_genome_status",
    "orchestrate_genome_setup",
]

