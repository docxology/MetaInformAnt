"""Genome preparation utilities for RNA-seq workflows.

This module provides functions to:
- Extract transcriptome FASTA files from downloaded genome packages
- Prepare FASTA files for kallisto indexing
- Build kallisto indexes from transcriptome FASTA files
"""

from __future__ import annotations

import gzip
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Any

from ..core.io import ensure_directory

logger = logging.getLogger(__name__)


def find_rna_fasta_in_genome_dir(genome_dir: Path, accession: str) -> Path | None:
    """Find RNA FASTA file in extracted genome directory.
    
    NCBI datasets extracts to various patterns:
    - ncbi_dataset_api_extracted/ncbi_dataset/data/{accession}/rna.fna
    - ncbi_dataset_api_extracted/ncbi_dataset/data/{accession}/rna.fna.gz
    - ncbi_dataset_extracted/ncbi_dataset/data/{accession}/rna.fna
    
    Also searches recursively for RNA FASTA files.
    
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
        # Direct file patterns
        genome_dir / "rna.fna",
        genome_dir / "rna.fna.gz",
    ]
    
    for pattern in patterns:
        if pattern.exists() and pattern.is_file():
            return pattern
    
    # Search recursively for RNA FASTA files
    for pattern in ["rna.fna", "rna.fna.gz", "transcript_sequences.fna", "*_rna_from_genomic.fna.gz"]:
        matches = list(genome_dir.rglob(pattern))
        if matches:
            # Prefer unzipped files
            unzipped = [m for m in matches if not m.name.endswith(".gz")]
            if unzipped:
                return unzipped[0]
            return matches[0]
    
    return None


def prepare_transcriptome_for_kallisto(
    genome_dir: Path,
    species_name: str,
    work_dir: Path,
    *,
    accession: str | None = None,
) -> Path | None:
    """Prepare transcriptome FASTA file for kallisto indexing.
    
    Finds the RNA FASTA file in the extracted genome directory and
    copies it to the expected location for amalgkit/kallisto.
    
    Expected output: work_dir/fasta/{Species_Name}_rna.fasta
    
    Args:
        genome_dir: Directory containing extracted genome package
        species_name: Species name (with underscores, e.g., "Camponotus_floridanus")
        work_dir: Work directory for amalgkit workflow
        accession: Optional NCBI assembly accession (for searching)
        
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
    
    if not rna_fasta or not rna_fasta.exists():
        logger.warning(f"RNA FASTA not found in {genome_dir}")
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


__all__ = [
    "find_rna_fasta_in_genome_dir",
    "prepare_transcriptome_for_kallisto",
    "build_kallisto_index",
    "get_expected_index_path",
    "prepare_genome_for_quantification",
]

