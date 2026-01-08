"""RNA genome and transcriptome preparation for Kallisto quantification.

This module provides tools for preparing genomic references for RNA-seq analysis.
It handles downloading genome sequences, extracting transcripts from GFF3 files,
building Kallisto indices, and orchestrating the complete genome preparation workflow.
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
import urllib.request
import urllib.error
import ssl
from pathlib import Path
from typing import Any, Dict, Optional

from metainformant.core import logging

logger = logging.get_logger(__name__)


def find_rna_fasta_in_genome_dir(genome_dir: Path, accession: str) -> Optional[Path]:
    """Locate RNA FASTA file in genome directory."""
    genome_dir = Path(genome_dir)

    # Search patterns
    search_patterns = [
        genome_dir / "ncbi_dataset_api_extracted" / "ncbi_dataset" / "data" / accession,
        genome_dir / "ncbi_dataset_extracted" / "ncbi_dataset" / "data" / accession,
        genome_dir,
    ]

    for search_dir in search_patterns:
        if not search_dir.exists():
            continue

        rna_fna = search_dir / "rna.fna"
        rna_fna_gz = search_dir / "rna.fna.gz"
        
        # Also check for files ending with _rna.fna or _rna.fna.gz
        rna_suffix_fna = list(search_dir.glob("*_rna.fna"))
        rna_suffix_fna_gz = list(search_dir.glob("*_rna.fna.gz"))
        
        if rna_fna.exists():
            return rna_fna
        if rna_fna_gz.exists():
            return rna_fna_gz
            
        if rna_suffix_fna:
            return rna_suffix_fna[0]
        if rna_suffix_fna_gz:
            return rna_suffix_fna_gz[0]

    return None

def get_expected_index_path(work_dir: Path, species_name: str) -> Path:
    """Get standardized path for Kallisto index."""
    work_dir = Path(work_dir)
    # Standard path structure: work_dir/species_name/index/species_name.idx
    # Assuming typical amalgkit structure where work_dir is the species output root?
    # Based on setup_genome.py, work_dir is typically output/amalgkit/species/work
    # Or output/amalgkit/species/
    
    # Let's assume standardized "index" folder
    index_dir = work_dir / "index"
    index_path = index_dir / f"{species_name}.idx"
    return index_path


def _download_url(url: str, output_path: Path, timeout: int = 600) -> bool:
    """Helper to download file from URL."""
    try:
        # Create context to ignore SSL errors if needed (often an issue with some legacy scientific sites)
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        
        logger.info(f"Downloading {url} to {output_path}")
        with urllib.request.urlopen(url, context=ctx, timeout=timeout) as response, open(output_path, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
        return True
    except Exception as e:
        logger.error(f"Failed to download {url}: {e}")
        if output_path.exists():
            output_path.unlink()
        return False

def download_rna_fasta_from_ftp(
    accession: str, output_dir: Path, **kwargs: Any
) -> Optional[Path]:
    """Download RNA sequences from NCBI FTP server.
    
    Expects kwargs to contain 'ftp_url' or relies on standard NCBI structure.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    ftp_url = kwargs.get('ftp_url', "")
    filename = kwargs.get('filename', f"{accession}_rna.fna.gz") # Placeholder default
    
    if not ftp_url:
        logger.warning("No FTP URL provided for download")
        return None
        
    # If ftp_url is a directory (ends with /), append filename
    if not ftp_url.endswith(".gz") and not ftp_url.endswith(".fna") and not ftp_url.endswith(".fa"):
        # We need to guess filename if not provided? 
        # Actually config provided exact maps sometimes.
        # But if the URL is "https://ftp.../GCA_.../", we usually look for "GCA_..._rna_from_genomic.fna.gz"
        
        # Standard NCBI naming: {AssemblyAccession}_{AssemblyName}_rna_from_genomic.fna.gz
        # We might need to construct it if filename not explicit.
        pass

    target_url = ftp_url
    if filename and not target_url.endswith(filename) and not target_url.endswith(".gz"):
         target_url = f"{ftp_url.rstrip('/')}/{filename}"
    
    output_path = output_dir / (filename if filename else Path(target_url).name)
    
    if output_path.exists():
        logger.info(f"File already exists: {output_path}")
        return output_path
        
    if _download_url(target_url, output_path):
        return output_path
    return None


def download_cds_fasta_from_ftp(
    accession: str, output_dir: Path, **kwargs: Any
) -> Optional[Path]:
    """Download CDS sequences."""
    # Logic is identical to RNA fasta, just potentially different default filename assumptions
    return download_rna_fasta_from_ftp(accession, output_dir, **kwargs)


def extract_transcripts_from_gff(
    gff_file: Path, genome_fasta: Path, output_file: Path, **kwargs: Any
) -> Optional[Path]:
    """Extract transcript sequences from GFF3 annotation using gffread."""
    gff_file = Path(gff_file)
    genome_fasta = Path(genome_fasta)
    output_file = Path(output_file)
    
    if not gff_file.exists() or not genome_fasta.exists():
        logger.error("GFF or Genome FASTA missing")
        return None
        
    # Check for gffread
    gffread = shutil.which("gffread")
    if not gffread:
        logger.error("gffread not found in PATH")
        return None
        
    cmd = [
        gffread,
        str(gff_file),
        "-g", str(genome_fasta),
        "-w", str(output_file)
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        if output_file.exists() and output_file.stat().st_size > 0:
            return output_file
        return None
    except subprocess.CalledProcessError as e:
        logger.error(f"gffread failed: {e}")
        return None


def prepare_transcriptome_for_kallisto(
    transcriptome_fasta: Path, output_fasta: Path, **kwargs: Any
) -> Optional[Path]:
    """Prepare transcriptome FASTA for Kallisto.
    
    Essentially check validity and unzip if needed.
    """
    transcriptome_fasta = Path(transcriptome_fasta)
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    
    if not transcriptome_fasta.exists():
        logger.error(f"Input missing: {transcriptome_fasta}")
        return None
        
    # Check if gzipped
    is_gzipped = transcriptome_fasta.name.endswith(".gz")
    
    try:
        if is_gzipped:
            with gzip.open(transcriptome_fasta, 'rt') as f_in, open(output_fasta, 'w') as f_out:
                shutil.copyfileobj(f_in, f_out)
        else:
            shutil.copyfile(transcriptome_fasta, output_fasta)
            
        return output_fasta
    except Exception as e:
        logger.error(f"Failed to prepare transcriptome: {e}")
        return None


def build_kallisto_index(
    transcriptome_fasta: Path, index_path: Path, **kwargs: Any
) -> Optional[Path]:
    """Build Kallisto index."""
    transcriptome_fasta = Path(transcriptome_fasta)
    index_path = Path(index_path)
    
    if not transcriptome_fasta.exists():
        logger.error(f"Transcriptome FASTA not found: {transcriptome_fasta}")
        return None
        
    kallisto_cmd = shutil.which("kallisto")
    if not kallisto_cmd:
        logger.error("kallisto not found in PATH")
        return None
        
    kmer_size = kwargs.get("kmer_size", 31)
    
    index_path.parent.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        kallisto_cmd,
        "index",
        "-i", str(index_path),
        "-k", str(kmer_size),
        str(transcriptome_fasta)
    ]
    
    try:
        logger.info(f"Running kallisto index on {transcriptome_fasta}")
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return index_path
    except subprocess.CalledProcessError as e:
        logger.error(f"kallisto index failed: {e.stderr}")
        return None


def verify_genome_status(work_dir: Path, species_name: str) -> Dict[str, Any]:
    """Check genome preparation status."""
    work_dir = Path(work_dir)
    
    # Try different directory structures. 
    # setup_genome.py puts things in 'genome' dir inside work root sometimes?
    # Based on config: dest_dir: output/amalgkit/pbarbatus_test5/genome
    # Config is not passed here directly, just work_dir.
    # We'll search in potential locations.
    
    # Assuming work_dir contains 'genome' folder or similar?
    # Or work_dir is output/species/work
    
    # Let's search relative to work_dir parent if needed
    base_dir = work_dir.parent # output/species
    
    genome_dir = base_dir / "genome"
    index_dir = work_dir / "index" # Based on amalgkit
    
    transcriptome_exists = False
    rna_fasta = find_rna_fasta_in_genome_dir(genome_dir, f"{species_name}*")
    # Also check typical downloaded file names
    if not rna_fasta:
         rna_fasta = next(genome_dir.glob("*_rna_from_genomic.fna*"), None)
    
    if rna_fasta and rna_fasta.exists():
        transcriptome_exists = True
        
    index_path = index_dir / f"{species_name}.idx"
    index_exists = index_path.exists()
    
    return {
        "prepared": index_exists,
        "index_exists": index_exists,
        "transcriptome_exists": transcriptome_exists,
        "issues": [],
        "genome_downloaded": transcriptome_exists # Simplify for verify script
    }


def orchestrate_genome_setup(
    config: Dict[str, Any],
    species_name: str,
    work_dir: Path,
    **kwargs: Any
) -> Dict[str, Any]:
    """Coordinate all genome setup steps."""
    work_dir = Path(work_dir)
    # The 'work_dir' from config is typically 'work/', but downloads go to 'genome/' sibling
    # Let's parse config to find exact locations if possible
    
    genome_cfg = config.get('genome', {})
    
    # 1. Determine download directory
    dest_dir_str = genome_cfg.get('dest_dir')
    if dest_dir_str:
        genome_dir = Path(dest_dir_str)
    else:
        # Fallback to sibling of work_dir
        genome_dir = work_dir.parent / "genome"
    
    genome_dir.mkdir(parents=True, exist_ok=True)
    
    # 2. Download files
    ftp_url = genome_cfg.get('ftp_url')
    if not ftp_url:
        logger.error("No FTP URL in genome config")
        return {"success": False, "messages": ["No FTP URL"]}
        
    messages = []
    
    # Check specifically for transcriptome
    files_map = genome_cfg.get('files', {})
    rna_filename = files_map.get('transcriptome_fasta')
    
    logger.info("Downloading RNA sequences...")
    rna_path = download_rna_fasta_from_ftp(
        accession=genome_cfg.get('accession', species_name),
        output_dir=genome_dir,
        ftp_url=ftp_url,
        filename=rna_filename
    )
    
    if not rna_path:
        return {"success": False, "messages": ["Failed to download RNA FASTA"]}
        
    messages.append(f"Downloaded RNA FASTA: {rna_path}")
    
    # 3. Build Index
    # We typically need to unzip first for some tools, but kallisto handles gz mostly?
    # Actually kallisto index typically accepts .gz.
    
    index_dir = work_dir / "index" # Standard amalgkit location
    index_path = get_expected_index_path(work_dir, species_name)
    
    if not kwargs.get("skip_build"):
        logger.info("Building Kallisto index...")
        idx = build_kallisto_index(rna_path, index_path)
        if idx:
            messages.append(f"Built index: {idx}")
        else:
            return {"success": False, "messages": messages + ["Index build failed"]}
            
    return {"success": True, "messages": messages, "species": species_name}
