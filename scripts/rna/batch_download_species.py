#!/usr/bin/env python3
"""
Batch download samples for multiple species in parallel with dynamic thread allocation and immediate per-sample quantification.

# ============================================================================
# CONFIGURATION
# ============================================================================
# Scope: Multi-species parallel batch download with configurable throughput
# Steps: getfastq → quant (per-sample, immediate)
# Config: Auto-discovers all config/amalgkit/amalgkit_*.yaml files
# Threads: Configurable total (default 24) distributed across species
# Batch Size: Per-sample processing (download → quant → delete immediately)
# Output: output/amalgkit/{species}/fastq/ and quant/ per species
# Dependencies: amalgkit, kallisto, wget (for ENA) or SRA Toolkit
# Virtual Env: Auto-activates if .venv exists
# Disk Management: Pre-flight checks, auto-cleanup when space low
# Reliability: Depends on download method (ENA recommended)
# ============================================================================

Features:
- Configurable total thread count (default: 24, recommended for standard systems)
- Even initial distribution across species (minimum 1 thread per species)
- Dynamic reallocation as species complete
- Automatic thread concentration on remaining species
- Immediate quantification and deletion: Each sample is quantified and FASTQs deleted as soon as download completes
- Separate thread pool for quantification operations (default: 10 threads)
- Disk space monitoring: Pre-flight checks and automatic cleanup when space is low
- Resilient error handling: Automatic retry and recovery from disk space issues

Thread Allocation:
- Downloads: 24 threads default (distributed evenly, minimum 1 per species)
- Can increase to 48+ threads on larger drives with more disk space
- Quantification: Separate pool (default: 10 threads) for immediate quant+delete operations
- Dynamic: As species complete, download threads redistribute to remaining species
- Concentration: Threads concentrate on fewer species as workflow progresses

Disk Space Management:
- Pre-flight checks: Verifies sufficient disk space before starting (default: 10GB minimum)
- Automatic cleanup: Removes partial/failed downloads when space drops below threshold (default: 5GB)
- Health checks: Monitors disk space every 10 minutes during execution
- Error recovery: Automatically attempts cleanup and retry on disk space errors

Workflow:
1. Downloads run in parallel across species with dynamic thread allocation
2. Every 60 seconds (configurable), monitor detects completed downloads:
   - Completed SRA files (not modified in 5+ min, reasonable size) → SRA→FASTQ conversion
   - Completed FASTQ files → Direct quantification
3. Pipeline for each sample: SRA→FASTQ (if needed) → Quant → Delete
4. Processes ALL species (including backlogged completed downloads)
5. No waiting for batch completion - maximum disk efficiency

Usage:
  # Default (24 threads, recommended for standard systems)
  python3 scripts/rna/batch_download_species.py
  
  # Larger drive (48 threads)
  python3 scripts/rna/batch_download_species.py --total-threads 48
  
  # Custom disk space thresholds
  python3 scripts/rna/batch_download_species.py --min-free-gb 15.0 --auto-cleanup-threshold 10.0
"""

import os
import sys
import subprocess
import time
import shutil
from pathlib import Path
from datetime import datetime
from glob import glob
from typing import Optional
from concurrent.futures import ThreadPoolExecutor, as_completed, Future

# Import setup utilities (must be before other imports)
sys.path.insert(0, str(Path(__file__).parent))
from _setup_utils import ensure_venv_activated, check_environment_or_exit

# Auto-setup and activate venv using uv
ensure_venv_activated(auto_setup=True)
check_environment_or_exit(auto_setup=True)

# Add src to path (must be before any imports)
repo_root = Path(__file__).parent.parent.parent.resolve()
src_path = repo_root / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

from metainformant.rna.workflow import load_workflow_config
from metainformant.rna.amalgkit import build_amalgkit_command, run_amalgkit
from metainformant.rna.steps import process_sample_pipeline, delete_sample_fastqs
from metainformant.rna import count_quantified_samples
from metainformant.rna.progress_tracker import get_tracker
from metainformant.core.disk import check_disk_space, get_disk_space_info
from metainformant.core.io import read_delimited, write_delimited
from metainformant.core.logging import get_logger

logger = get_logger("batch_download")


def is_sample_download_complete(sample_dir: Path) -> bool:
    """Check if a sample directory has completed download with valid FASTQ files.
    
    Args:
        sample_dir: Path to sample directory (e.g., fastq/getfastq/SRR12345 or fastq/SRR12345)
        
    Returns:
        True if sample has valid FASTQ files (non-empty .fastq.gz files, minimum 1KB each)
        
    Side effects:
        None (read-only operation)
        
    Notes:
        - Checks for .fastq.gz files
        - Validates file size (minimum 1KB to avoid empty/partial files)
    """
    if not sample_dir.exists() or not sample_dir.is_dir():
        return False
    
    # Check for FASTQ files
    fastq_files = list(sample_dir.glob("*.fastq.gz"))
    if not fastq_files:
        return False
    
    # Verify files are non-empty (at least 1KB)
    for fastq_file in fastq_files:
        if fastq_file.stat().st_size < 1024:
            return False
    
    return True


def is_sra_download_complete(sample_dir: Path) -> bool:
    """Check if a sample has completed SRA download (ready for conversion to FASTQ).
    
    Args:
        sample_dir: Path to sample directory (e.g., fastq/getfastq/SRR12345)
        
    Returns:
        True if sample has SRA file that appears complete:
        - Size > 100 MB and < 100 GB (reasonable range)
        - Not modified in last 5 minutes (download finished)
        
    Side effects:
        None (read-only operation)
        
    Notes:
        - Uses file modification time to detect completion
        - Filters out suspiciously large files (>100 GB)
    """
    if not sample_dir.exists() or not sample_dir.is_dir():
        return False
    
    # Check for SRA files
    sra_files = list(sample_dir.glob("*.sra"))
    if not sra_files:
        return False
    
    # Check if SRA file is complete
    # Criteria: file exists, has reasonable size (>100 MB), and hasn't been modified in last 5 minutes
    import time
    current_time = time.time()
    
    for sra_file in sra_files:
        size_mb = sra_file.stat().st_size / (1024 * 1024)
        mod_time = sra_file.stat().st_mtime
        age_minutes = (current_time - mod_time) / 60
        
        # File is complete if:
        # 1. Size > 100 MB (reasonable minimum)
        # 2. Size < 100 GB (not an error - suspiciously large files)
        # 3. Not modified in last 5 minutes (download finished)
        if size_mb > 100 and size_mb < 100000 and age_minutes >= 5:
            return True
    
    return False


def get_species_id_from_config(config_path: Path) -> str:
    """Extract species ID from config file.
    
    Args:
        config_path: Path to config file
        
    Returns:
        Species ID (e.g., "camponotus_floridanus" from config filename or species_list)
    """
    try:
        cfg = load_workflow_config(config_path)
        # Try to get from species_list first
        if cfg.species_list:
            # Use first species, convert to lowercase with underscores
            species_id = cfg.species_list[0].lower().replace(" ", "_")
            return species_id
    except Exception:
        pass
    
    # Fallback: extract from filename
    # e.g., "amalgkit_camponotus_floridanus.yaml" -> "camponotus_floridanus"
    stem = config_path.stem
    if stem.startswith("amalgkit_"):
        return stem.replace("amalgkit_", "")
    return stem


def process_sample_with_tracker(
    sample_id: str,
    config_path: Path,
    status: str,
    tracker=None,
    *,
    log_dir: Path | None = None,
) -> tuple[bool, str]:
    """Process a sample with progress tracking.
    
    Wraps process_sample_pipeline to add ProgressTracker integration.
    
    Args:
        sample_id: SRA accession ID
        config_path: Path to species config file
        status: "sra" or "fastq"
        tracker: ProgressTracker instance (optional)
        log_dir: Optional log directory
        
    Returns:
        Tuple of (success: bool, message: str)
    """
    if tracker:
        try:
            species_id = get_species_id_from_config(config_path)
            # Download is already complete (reported in detect_completed_downloads)
            # Now process: quant + delete
        except Exception as e:
            logger.debug(f"Failed to get species_id for tracker: {e}")
            tracker = None
    
    # Process the sample
    success, message = process_sample_pipeline(
        sample_id,
        config_path,
        status,
        log_dir=log_dir,
    )
    
    # Update tracker based on result
    if tracker:
        try:
            species_id = get_species_id_from_config(config_path)
            if success:
                # Quant and delete both completed
                tracker.on_quant_complete(species_id, sample_id)
                tracker.on_delete_complete(species_id, sample_id)
            else:
                # Quant may have failed, but check if it actually completed
                cfg = load_workflow_config(config_path)
                quant_params = dict(cfg.per_step.get("quant", {}))
                quant_dir = Path(quant_params.get("out_dir", cfg.work_dir / "quant"))
                abundance_file = quant_dir / sample_id / "abundance.tsv"
                if abundance_file.exists():
                    # Quant actually succeeded, just delete failed or something else
                    tracker.on_quant_complete(species_id, sample_id)
                    tracker.on_delete_complete(species_id, sample_id)
        except Exception as e:
            logger.debug(f"Failed to update tracker after processing {sample_id}: {e}")
    
    return success, message


def detect_completed_downloads(
    species_configs: list[tuple[str, Path]],
    processed_samples: set[str],
    quant_dirs: dict[str, Path],
    tracker=None,
) -> list[tuple[str, str, Path, str]]:
    """Detect samples that have completed download but not yet been quantified.
    
    Args:
        species_configs: List of (species_name, config_path) tuples
        processed_samples: Set of sample IDs already being processed (format: "species:sample_id")
        quant_dirs: Dictionary mapping species_name -> quant_dir Path
        
    Returns:
        List of (species_name, sample_id, sample_dir, status) tuples where status is:
        - "fastq": Has FASTQ files ready for quant
        - "sra": Has SRA file ready for conversion
        
    Side effects:
        - Updates processed_samples set to mark detected samples
        
    Notes:
        - Checks both fastq/getfastq/{sample}/ and fastq/{sample}/ directory structures
        - Skips samples already quantified or currently being processed
    """
    completed = []
    
    for species_name, config_path in species_configs:
        try:
            cfg = load_workflow_config(config_path)
            fastq_dir = Path(cfg.per_step.get("getfastq", {}).get("out_dir", cfg.work_dir / "fastq"))
            quant_dir = quant_dirs.get(species_name)
            
            if not fastq_dir.exists():
                continue
            
            # Check both structures: fastq/getfastq/{sample}/ and fastq/{sample}/
            search_dirs = []
            getfastq_subdir = fastq_dir / "getfastq"
            if getfastq_subdir.exists():
                search_dirs.append(getfastq_subdir)
            search_dirs.append(fastq_dir)
            
            for search_dir in search_dirs:
                if not search_dir.exists():
                    continue
                
                for sample_dir in search_dir.iterdir():
                    if not sample_dir.is_dir():
                        continue
                    
                    sample_id = sample_dir.name
                    sample_key = f"{species_name}:{sample_id}"
                    
                    # Skip if already processed or already quantified
                    if sample_key in processed_samples:
                        continue
                    
                    if quant_dir and (quant_dir / sample_id / "abundance.tsv").exists():
                        continue
                    
                    # Check for FASTQ files first (highest priority)
                    if is_sample_download_complete(sample_dir):
                        completed.append((species_name, sample_id, sample_dir, "fastq"))
                        processed_samples.add(sample_key)
                        # Report download complete to tracker
                        if tracker:
                            try:
                                species_id = get_species_id_from_config(config_path)
                                tracker.on_download_complete(species_id, sample_id)
                            except Exception as e:
                                logger.debug(f"Failed to update tracker for {sample_id}: {e}")
                    # Check for completed SRA files (need conversion)
                    elif is_sra_download_complete(sample_dir):
                        completed.append((species_name, sample_id, sample_dir, "sra"))
                        processed_samples.add(sample_key)
                        # Report download complete to tracker (SRA is downloaded, needs conversion)
                        if tracker:
                            try:
                                species_id = get_species_id_from_config(config_path)
                                tracker.on_download_complete(species_id, sample_id)
                            except Exception as e:
                                logger.debug(f"Failed to update tracker for {sample_id}: {e}")
        
        except Exception as e:
            logger.debug(f"Error checking downloads for {species_name}: {e}")
            continue
    
    return completed


# delete_sample_fastqs now imported from metainformant.rna.steps.getfastq

# NOTE: convert_sra_to_fastq moved to metainformant.rna.steps.getfastq.convert_sra_to_fastq
def _convert_sra_to_fastq_legacy(
    species_name: str,
    sample_id: str,
    config_path: Path,
    log_dir: Path,
) -> tuple[bool, str]:
    """Convert SRA file to FASTQ for a sample using fasterq-dump directly.
    
    Args:
        species_name: Display name for species
        sample_id: Sample ID to convert
        config_path: Path to species config file
        log_dir: Directory for log files
        
    Returns:
        Tuple of (success: bool, message: str)
    """
    import shutil
    import subprocess
    
    sample_logger = get_logger(f"convert_{species_name.replace(' ', '_')}")
    
    try:
        cfg = load_workflow_config(config_path)
        
        # Get fastq directory
        getfastq_params = dict(cfg.per_step.get("getfastq", {}))
        fastq_dir = Path(getfastq_params.get("out_dir", cfg.work_dir / "fastq"))
        
        # Find SRA file
        sample_dir = fastq_dir / "getfastq" / sample_id
        if not sample_dir.exists():
            sample_dir = fastq_dir / sample_id
        
        if not sample_dir.exists():
            return False, f"Sample directory not found: {sample_dir}"
        
        sra_files = list(sample_dir.glob("*.sra"))
        if not sra_files:
            # Check if FASTQ already exists
            fastq_files = list(sample_dir.glob("*.fastq.gz"))
            if not fastq_files:
                fastq_files = list(sample_dir.glob("*.fastq"))
            if fastq_files:
                sample_logger.info(f"✅ FASTQ files already exist for {sample_id}")
                return True, f"FASTQ files already exist for {sample_id}"
            return False, f"No SRA file found for {sample_id}"
        
        sra_file = sra_files[0]
        sample_logger.info(f"🔄 Converting SRA to FASTQ for {sample_id} (SRA: {sra_file.name})...")
        
        # Use fasterq-dump directly on the SRA file
        fasterq_dump = shutil.which("fasterq-dump")
        if not fasterq_dump:
            return False, "fasterq-dump not found in PATH"
        
        # Get thread count from config
        threads = cfg.threads or getfastq_params.get("threads", 4)
        
        # Run fasterq-dump
        # fasterq-dump needs the accession ID and will look for the SRA file
        # in the current working directory. We need to run it from the sample_dir
        # where the SRA file is located.
        cmd = [
            fasterq_dump,
            sample_id,  # Use accession ID
            "--outdir", str(sample_dir),
            "--threads", str(threads),
            "--split-files",  # Split paired-end reads
            "-p",  # Show progress
        ]
        
        log_file = log_dir / f"fasterq_dump_{species_name.replace(' ', '_')}_{sample_id}.log"
        log_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Run from the sample directory where the SRA file is located
        # fasterq-dump will find the SRA file there
        with open(log_file, "w") as f:
            result = subprocess.run(
                cmd,
                stdout=f,
                stderr=subprocess.STDOUT,
                cwd=str(sample_dir),  # Run from directory containing the SRA file
            )
        
        # Compress FASTQ files
        pigz = shutil.which("pigz") or "gzip"
        for fastq_file in sample_dir.glob(f"{sample_id}_*.fastq"):
            if not fastq_file.name.endswith('.gz'):
                subprocess.run([pigz, "-f", str(fastq_file)], check=False)
        
        # Also check for single-end FASTQ
        for fastq_file in sample_dir.glob(f"{sample_id}.fastq"):
            if not fastq_file.name.endswith('.gz'):
                subprocess.run([pigz, "-f", str(fastq_file)], check=False)
        
        # Check if FASTQ files were created
        fastq_files = list(sample_dir.glob("*.fastq.gz"))
        if not fastq_files:
            fastq_files = list(sample_dir.glob("*.fastq"))
        
        if fastq_files:
            sample_logger.info(f"✅ Converted SRA to FASTQ for {sample_id} ({len(fastq_files)} files)")
            return True, f"Converted SRA to FASTQ for {sample_id}"
        else:
            sample_logger.warning(f"⚠️  SRA conversion failed for {sample_id} (code {result.returncode})")
            return False, f"SRA conversion failed (code {result.returncode})"
    
    except Exception as e:
        sample_logger.error(f"❌ Error converting SRA for {sample_id}: {e}", exc_info=True)
        return False, str(e)


# NOTE: process_complete_sample replaced by metainformant.rna.steps.sample_pipeline.process_sample_pipeline
def _process_complete_sample_legacy(
    species_name: str,
    sample_id: str,
    config_path: Path,
    log_dir: Path,
    status: str,
) -> tuple[bool, str]:
    """Process a complete sample: convert SRA→FASTQ if needed, then quant→delete.
    
    Args:
        species_name: Display name for species
        sample_id: Sample ID to process
        config_path: Path to species config file
        log_dir: Directory for log files
        status: "fastq" (has FASTQ) or "sra" (needs conversion)
        
    Returns:
        Tuple of (success: bool, message: str)
    """
    if status == "sra":
        # First convert SRA to FASTQ
        success, message = _convert_sra_to_fastq_legacy(species_name, sample_id, config_path, log_dir)
        if not success:
            return False, f"SRA conversion failed: {message}"
        # Fall through to quant+delete
    
    # Now quant and delete (works for both "fastq" and after "sra" conversion)
    return _quantify_and_delete_sample_legacy(species_name, sample_id, config_path, log_dir)


# NOTE: quantify_and_delete_sample replaced by metainformant.rna.steps.quant.quantify_sample + delete
def _quantify_and_delete_sample_legacy(
    species_name: str,
    sample_id: str,
    config_path: Path,
    log_dir: Path,
) -> tuple[bool, str]:
    """Quantify a sample and delete its FASTQ files immediately after.
    
    Args:
        species_name: Display name for species
        sample_id: Sample ID to quantify
        config_path: Path to species config file
        log_dir: Directory for log files
        
    Returns:
        Tuple of (success: bool, message: str)
    """
    sample_logger = get_logger(f"quant_{species_name.replace(' ', '_')}")
    
    try:
        cfg = load_workflow_config(config_path)
        
        # Get metadata file
        metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
        if not metadata_file.exists():
            metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
        
        if not metadata_file.exists():
            return False, f"No metadata file found for {species_name}"
        
        # Read metadata to find this sample
        rows = list(read_delimited(metadata_file, delimiter="\t"))
        sample_rows = [row for row in rows if row.get("run") == sample_id]
        
        if not sample_rows:
            return False, f"Sample {sample_id} not found in metadata"
        
        # Get quant params
        quant_params = dict(cfg.per_step.get("quant", {}))
        quant_dir = Path(quant_params.get("out_dir", cfg.work_dir / "quant"))
        quant_dir.mkdir(parents=True, exist_ok=True)
        quant_params["out_dir"] = str(quant_dir.absolute())
        
        # Create temporary single-sample metadata
        temp_metadata = cfg.work_dir / f"metadata.quant.{sample_id}.tsv"
        write_delimited(sample_rows, temp_metadata, delimiter="\t")
        quant_params["metadata"] = str(temp_metadata.absolute())
        
        # Run quantification
        sample_logger.info(f"🔬 Quantifying {sample_id}...")
        result = run_amalgkit(
            "quant",
            quant_params,
            work_dir=None,
            log_dir=log_dir,
            step_name=f"quant_{species_name.replace(' ', '_')}_{sample_id}",
            check=False,
        )
        
        # Clean up temp metadata
        try:
            temp_metadata.unlink()
        except Exception:
            pass
        
        # Verify quantification succeeded
        abundance_file = quant_dir / sample_id / "abundance.tsv"
        if result.returncode == 0 and abundance_file.exists():
            sample_logger.info(f"✅ Quantified {sample_id}")
            
            # Delete FASTQ files immediately
            fastq_dir = Path(cfg.per_step.get("getfastq", {}).get("out_dir", cfg.work_dir / "fastq"))
            delete_sample_fastqs(sample_id, fastq_dir)
            
            return True, f"Quantified and deleted {sample_id}"
        else:
            sample_logger.warning(f"⚠️  Quantification failed for {sample_id} (code {result.returncode})")
            
            # Still delete FASTQ to free space even if quant failed
            fastq_dir = Path(cfg.per_step.get("getfastq", {}).get("out_dir", cfg.work_dir / "fastq"))
            delete_sample_fastqs(sample_id, fastq_dir)
            
            return False, f"Quantification failed (code {result.returncode})"
    
    except Exception as e:
        sample_logger.error(f"❌ Error processing {sample_id}: {e}", exc_info=True)
        return False, str(e)


# count_quantified_samples now imported from metainformant.rna


def is_species_complete(config_path: Path, completion_threshold: float = 0.95) -> bool:
    """Check if species download/quantification is complete.
    
    Args:
        config_path: Path to species config file
        completion_threshold: Fraction of samples that must be quantified (default: 0.95)
        
    Returns:
        True if species is complete (≥threshold quantified)
    """
    quantified, total = count_quantified_samples(config_path)
    if total == 0:
        return False
    return (quantified / total) >= completion_threshold


def distribute_threads(total_threads: int, num_species: int) -> list[int]:
    """Distribute threads evenly across species.
    
    Args:
        total_threads: Total number of threads to distribute
        num_species: Number of species
        
    Returns:
        List of thread counts per species (one per species, minimum 1 per species)
        
    Example:
        distribute_threads(30, 25) -> [1]*20 + [2]*5 (20 species get 1, 5 get 2)
        distribute_threads(8, 24) -> [1]*24 (all species get 1 thread minimum)
    """
    if num_species == 0:
        return []
    
    # If we have fewer threads than species, each gets at least 1
    # This shouldn't happen often, but ensures all species get at least 1 thread
    if total_threads < num_species:
        # All species get 1 thread (we'll use all available threads)
        allocation = [1] * num_species
        return allocation
    
    base_threads = total_threads // num_species
    extra_threads = total_threads % num_species
    
    # Base allocation: all species get base_threads (minimum 1)
    allocation = [max(1, base_threads)] * num_species
    
    # Distribute extra threads to first extra_threads species
    for i in range(extra_threads):
        allocation[i] += 1
    
    return allocation


def redistribute_threads(total_threads: int, active_species: list[tuple[str, Path]]) -> dict[str, int]:
    """Redistribute threads among active (incomplete) species.
    
    Args:
        total_threads: Total number of threads available
        active_species: List of (species_name, config_path) tuples for incomplete species
        
    Returns:
        Dictionary mapping species_name -> thread_count
    """
    if not active_species:
        return {}
    
    num_active = len(active_species)
    thread_allocation = distribute_threads(total_threads, num_active)
    
    return {
        species_name: threads
        for (species_name, _), threads in zip(active_species, thread_allocation)
    }


def check_process_health(process: subprocess.Popen, log_file: Path, timeout: int) -> tuple[bool, str]:
    """Check if a process is healthy (not hung).
    
    Args:
        process: The subprocess.Popen object
        log_file: Path to log file
        timeout: Timeout in seconds - if log not modified in this time, consider hung
        
    Returns:
        Tuple of (is_healthy: bool, reason: str)
    """
    if process.poll() is not None:
        return False, "Process has terminated"
    
    # Check if log file exists and has been modified recently
    if log_file.exists():
        last_modified = log_file.stat().st_mtime
        age = time.time() - last_modified
        if age > timeout:
            return False, f"Log file not modified in {age:.0f}s (hung process)"
        
        # Check log file size - if it's growing, process is active
        current_size = log_file.stat().st_size
        # If file is very small and old, might be stuck
        if current_size < 1000 and age > 300:
            return False, "Log file too small and old (likely stuck)"
    
    return True, "Process appears healthy"


def start_download_process(config_path: Path, species_name: str, threads: int, log_dir: Path, retry_num: int = 0, auto_cleanup_threshold: float = 5.0) -> tuple[subprocess.Popen, Path]:
    """Start an amalgkit getfastq process for a species.
    
    Args:
        config_path: Path to species config file
        species_name: Display name for species
        threads: Number of threads to use
        log_dir: Directory for log files
        
    Returns:
        subprocess.Popen object for the running process
    """
    try:
        cfg = load_workflow_config(config_path)
        
        # Get metadata file
        metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
        if not metadata_file.exists():
            metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
        
        if not metadata_file.exists():
            raise FileNotFoundError(f"No metadata file found for {species_name}")
        
        # Build getfastq params
        getfastq_params = dict(cfg.per_step.get("getfastq", {}))
        getfastq_params["out_dir"] = str(Path(getfastq_params.get("out_dir", cfg.work_dir / "fastq")).absolute())
        getfastq_params["metadata"] = str(metadata_file.absolute())
        getfastq_params["threads"] = threads
        
        # Enable cloud acceleration if configured
        if cfg.per_step.get("getfastq", {}).get("accelerate", True):
            getfastq_params.setdefault("aws", "yes")
            getfastq_params.setdefault("gcp", "yes")
            getfastq_params.setdefault("ncbi", "yes")
        
        # Set remove_sra to yes to save space
        getfastq_params["remove_sra"] = "yes"
        
        # Build command
        cmd = build_amalgkit_command("getfastq", getfastq_params)
        
        # Setup log file with rotation
        log_dir.mkdir(parents=True, exist_ok=True)
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        retry_suffix = f"_retry{retry_num}" if retry_num > 0 else ""
        log_file = log_dir / f"getfastq_{species_name.replace(' ', '_')}_{timestamp}{retry_suffix}.log"
        
        logger.info(f"Starting {species_name} with {threads} threads (retry {retry_num})")
        logger.info(f"  Command: {' '.join(cmd[:5])}...")
        logger.info(f"  Log: {log_file}")
        
        # Start process with enhanced logging
        log_handle = open(log_file, 'w')
        log_handle.write(f"=== Batch Download Process Started ===\n")
        log_handle.write(f"Species: {species_name}\n")
        log_handle.write(f"Threads: {threads}\n")
        log_handle.write(f"Retry: {retry_num}\n")
        log_handle.write(f"Timestamp: {datetime.now()}\n")
        log_handle.write(f"Command: {' '.join(cmd)}\n")
        log_handle.write("=" * 80 + "\n\n")
        log_handle.flush()
        
        process = subprocess.Popen(
            cmd,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            cwd=str(cfg.work_dir.parent) if cfg.work_dir else None,
        )
        
        # Write PID to log
        log_handle.write(f"\nProcess started with PID: {process.pid}\n")
        log_handle.flush()
        
        return process, log_file
        
    except Exception as e:
        # Check if it's a disk space error
        error_str = str(e).lower()
        if "no space" in error_str or "disk" in error_str or "enospc" in error_str:
            # Check disk space
            disk_info = get_disk_space_info(log_dir)
            logger.error(f"❌ Disk space error for {species_name}: {e}")
            logger.error(f"   Current free space: {disk_info['free_gb']:.1f}GB")
            if disk_info["free_gb"] < auto_cleanup_threshold:
                logger.warning("   Attempting cleanup...")
                try:
                    from scripts.rna.cleanup_partial_downloads import cleanup_species
                    result = cleanup_species(species_name, config_path, dry_run=False)
                    if result["freed_mb"] > 0:
                        logger.info(f"   Freed {result['freed_mb']/1024:.2f}GB, retrying...")
                        # Retry once after cleanup
                        return start_download_process(config_path, species_name, threads, log_dir, retry_num=retry_num, auto_cleanup_threshold=auto_cleanup_threshold)
                except Exception as cleanup_error:
                    logger.error(f"   Cleanup failed: {cleanup_error}")
        
        logger.error(f"Failed to start process for {species_name}: {e}")
        raise


def main():
    """Main entry point for batch download orchestrator.
    
    Orchestrates multi-species parallel downloads with:
    1. Auto-discovery of species configs
    2. Dynamic thread allocation across species
    3. Per-sample immediate quantification and cleanup
    4. Disk space monitoring and automatic cleanup
    5. Process health checks and recovery
    
    Configuration (command-line arguments):
        --total-threads: Total threads across all species (default: 24)
        --max-species: Maximum species to process (default: all)
        --check-interval: Seconds between completion checks (default: 300)
        --monitor-interval: Seconds between sample monitoring (default: 60)
        --quant-threads: Threads for quantification pool (default: 10)
        --min-free-gb: Minimum free disk space GB (default: 10.0)
        --auto-cleanup-threshold: GB threshold for auto-cleanup (default: 5.0)
        
    Side effects:
        - Downloads FASTQ files for all species
        - Quantifies samples immediately after download
        - Deletes FASTQ files after quantification
        - Writes logs to output/amalgkit/{species}/logs/
        
    Exit codes:
        0: Success (all species processed or completed)
        1: Failure (critical errors encountered)
    """
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Batch download samples with dynamic thread allocation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default: 24 total threads, distributed evenly across all species
  python3 scripts/rna/batch_download_species.py
  
  # 48 total threads, redistributed as species complete
  python3 scripts/rna/batch_download_species.py --total-threads 48
  
  # Process only first 10 species
  python3 scripts/rna/batch_download_species.py --max-species 10
        """
    )
    parser.add_argument(
        "--total-threads",
        type=int,
        default=24,
        help="Total number of threads to use across all species (default: 24, recommended for standard systems, can increase to 48+ on larger drives)"
    )
    parser.add_argument(
        "--max-species",
        type=int,
        default=None,
        help="Maximum number of species to process (default: all)"
    )
    parser.add_argument(
        "--check-interval",
        type=int,
        default=300,
        help="Seconds between completion checks (default: 300)"
    )
    parser.add_argument(
        "--monitor-interval",
        type=int,
        default=60,
        help="Seconds between sample download monitoring (default: 60)"
    )
    parser.add_argument(
        "--quant-threads",
        type=int,
        default=10,
        help="Number of parallel quantification operations (separate from download threads, default: 10)"
    )
    parser.add_argument(
        "--max-retries",
        type=int,
        default=3,
        help="Maximum retries for failed processes (default: 3)"
    )
    parser.add_argument(
        "--process-timeout",
        type=int,
        default=3600,
        help="Process timeout in seconds - restart if no activity (default: 3600 = 1 hour)"
    )
    parser.add_argument(
        "--health-check-interval",
        type=int,
        default=600,
        help="Health check interval in seconds (default: 600 = 10 minutes)"
    )
    parser.add_argument(
        "--min-free-gb",
        type=float,
        default=None,
        help="Minimum free disk space in GB before starting downloads (default: auto-detect based on drive size)"
    )
    parser.add_argument(
        "--auto-cleanup-threshold",
        type=float,
        default=None,
        help="Auto-cleanup partial downloads when free space drops below this GB (default: auto-detect based on drive size)"
    )
    parser.add_argument(
        "--max-batch-disk-gb",
        type=float,
        default=None,
        help="Maximum disk space per batch in GB (default: auto-calculate from batch size)"
    )
    
    args = parser.parse_args()
    
    repo_root = Path(__file__).parent.parent.parent.resolve()
    config_dir = repo_root / "config" / "amalgkit"
    
    if not config_dir.exists() or not list(config_dir.glob("amalgkit_*.yaml")):
        config_dir = repo_root / "config"
    
    # Discover all config files
    config_pattern = str(config_dir / "amalgkit_*.yaml")
    config_files = sorted(glob(config_pattern))
    
    species_configs = []
    for config_file in config_files:
        path = Path(config_file)
        stem_lower = path.stem.lower()
        if "template" in stem_lower or "test" in stem_lower:
            continue
        
        species_code = path.stem.replace("amalgkit_", "")
        display_name = species_code.replace("_", " ").title()
        species_configs.append((display_name, path))
    
    if not species_configs:
        logger.error("⚠️  No species configs found")
        return 1
    
    # Limit species if requested
    if args.max_species:
        species_configs = species_configs[:args.max_species]
    
    print("\n" + "=" * 80)
    print("DYNAMIC THREAD ALLOCATION: MULTI-SPECIES DOWNLOADS")
    print("=" * 80)
    print(f"Date: {datetime.now()}")
    
    # Check disk space before starting (use repo_root already defined)
    output_dir = repo_root / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Auto-detect disk space thresholds if not specified
    from metainformant.core.disk import detect_drive_size_category
    
    drive_category = detect_drive_size_category(output_dir)
    logger.info(f"Drive category: {drive_category}")
    
    if args.min_free_gb is None:
        # Auto-detect based on drive size
        if drive_category == "large":
            min_free_gb = 50.0
        elif drive_category == "medium":
            min_free_gb = 20.0
        else:
            min_free_gb = 10.0
        logger.info(f"Auto-detected min_free_gb: {min_free_gb} (based on {drive_category} drive)")
    else:
        min_free_gb = args.min_free_gb
        logger.info(f"Using specified min_free_gb: {min_free_gb}")
    
    if args.auto_cleanup_threshold is None:
        # Auto-detect based on drive size
        if drive_category == "large":
            auto_cleanup_threshold = 20.0
        elif drive_category == "medium":
            auto_cleanup_threshold = 10.0
        else:
            auto_cleanup_threshold = 5.0
        logger.info(f"Auto-detected auto_cleanup_threshold: {auto_cleanup_threshold} (based on {drive_category} drive)")
    else:
        auto_cleanup_threshold = args.auto_cleanup_threshold
        logger.info(f"Using specified auto_cleanup_threshold: {auto_cleanup_threshold}")
    
    disk_info = get_disk_space_info(output_dir)
    is_sufficient, disk_msg = check_disk_space(output_dir, min_free_gb=min_free_gb)
    
    print(f"\nDisk Space Check:")
    print(f"  Total: {disk_info['total_gb']:.1f}GB")
    print(f"  Used: {disk_info['used_gb']:.1f}GB ({disk_info['percent_used']})")
    print(f"  Free: {disk_info['free_gb']:.1f}GB ({disk_info['percent_free']})")
    print(f"  Status: {disk_msg}")
    sys.stdout.flush()  # Ensure output is written immediately
    
    if not is_sufficient:
        logger.warning(f"⚠️  Disk space check failed: {disk_msg}")
        logger.warning("Attempting automatic cleanup of partial downloads...")
        
        # Try to clean up partial downloads
        try:
            from scripts.rna.cleanup_partial_downloads import cleanup_species
            total_freed = 0
            for species_name, config_path in species_configs:
                if "template" in config_path.stem.lower():
                    continue
                result = cleanup_species(species_name, config_path, dry_run=False)
                total_freed += result["freed_mb"]
            
            if total_freed > 0:
                logger.info(f"✅ Freed {total_freed/1024:.2f}GB from partial downloads")
                # Re-check disk space
                is_sufficient, disk_msg = check_disk_space(output_dir, min_free_gb=min_free_gb)
                if not is_sufficient:
                    logger.error(f"❌ Still insufficient disk space after cleanup: {disk_msg}")
                    logger.error("Please free up disk space manually before continuing")
                    return 1
            else:
                logger.error("❌ No partial downloads to clean up. Please free up disk space manually")
                return 1
        except Exception as e:
            logger.error(f"❌ Failed to cleanup: {e}")
            logger.error("Please free up disk space manually before continuing")
            return 1
    
    print(f"\nThread Allocation:")
    print(f"  Downloads: {args.total_threads} threads TOTAL (distributed across all species)")
    print(f"  Quantification: {args.quant_threads} parallel operations (separate pool)")
    print(f"  Note: Each kallisto quant operation uses threads from config (typically 4-8)")
    print(f"\nSpecies discovered: {len(species_configs)}")
    print(f"Check interval: {args.check_interval} seconds")
    print(f"Monitor interval: {args.monitor_interval} seconds")
    print(f"Health check interval: {args.health_check_interval} seconds")
    print(f"Process timeout: {args.process_timeout} seconds")
    print(f"Max retries: {args.max_retries}")
    print(f"Min free space: {min_free_gb}GB")
    print(f"Auto-cleanup threshold: {auto_cleanup_threshold}GB")
    print("=" * 80)
    print("\nThread Allocation Strategy:")
    print("  • Initial: Even distribution across all species")
    print("  • Dynamic: Redistribute as species complete")
    print("  • Concentration: Threads focus on remaining species")
    print("  • Immediate: Quant+delete each sample as download completes")
    print("\nEnhanced Features:")
    print(f"  • Health checks: Detect and restart hung processes")
    print(f"  • Auto-retry: Automatically retry failed processes (up to {args.max_retries} times)")
    print(f"  • Enhanced logging: Comprehensive logs with rotation")
    print(f"  • Process monitoring: Track process activity and restart if needed")
    print(f"  • Disk space monitoring: Automatic cleanup when space is low")
    print(f"  • Resilient to disk space issues: Pre-flight checks and auto-recovery")
    print("=" * 80 + "\n")
    
    # Environment check already verified amalgkit in check_environment_or_exit()
    # Just log confirmation
    logger.info("✅ Environment check passed (amalgkit and dependencies verified)")
    print()
    
    # Filter out already-complete species
    active_species = []
    complete_species = []
    
    for species_name, config_path in species_configs:
        if is_species_complete(config_path):
            complete_species.append((species_name, config_path))
            quantified, total = count_quantified_samples(config_path)
            logger.info(f"✅ {species_name}: Already complete ({quantified}/{total} quantified)")
        else:
            active_species.append((species_name, config_path))
    
    if complete_species:
        print(f"Found {len(complete_species)} already-complete species (skipping)")
        print()
    
    if not active_species:
        logger.info("✅ All species are already complete!")
        return 0
    
    print(f"Processing {len(active_species)} active species with {args.total_threads} total threads")
    print()
    
    # Initial thread distribution
    thread_allocation = redistribute_threads(args.total_threads, active_species)
    
    # Validate allocation sums to exactly total_threads
    total_allocated = sum(thread_allocation.values())
    if total_allocated != args.total_threads:
        logger.warning(f"⚠️  Thread allocation mismatch: {total_allocated} allocated, {args.total_threads} requested")
    else:
        logger.info(f"✅ Thread allocation verified: {total_allocated} threads allocated across {len(thread_allocation)} species")
    
    print("Initial Thread Allocation (downloads):")
    for species_name, threads in sorted(thread_allocation.items()):
        print(f"  {species_name:30s}: {threads} threads")
    print(f"  {'TOTAL':30s}: {total_allocated} threads")
    print()
    
    # Setup logging directory
    repo_root = Path(__file__).parent.parent.parent.resolve()
    log_dir = repo_root / "output" / "batch_download_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # Track running processes
    running_processes: dict[str, subprocess.Popen] = {}
    process_configs: dict[str, Path] = {}
    results: dict[str, dict] = {}
    process_start_times: dict[str, float] = {}
    process_last_activity: dict[str, float] = {}
    process_retries: dict[str, int] = {}
    process_log_files: dict[str, Path] = {}
    
    # Track processed samples and quant operations
    processed_samples: set[str] = set()  # Format: "species_name:sample_id"
    quant_futures: dict[str, Future] = {}  # Format: "species_name:sample_id" -> Future
    
    # Build quant_dirs mapping for efficient lookup
    quant_dirs: dict[str, Path] = {}
    for species_name, config_path in active_species:
        try:
            cfg = load_workflow_config(config_path)
            quant_dir = Path(cfg.per_step.get("quant", {}).get("out_dir", cfg.work_dir / "quant"))
            quant_dirs[species_name] = quant_dir
        except Exception:
            pass
    
    # Create thread pool for quant operations
    quant_executor = ThreadPoolExecutor(max_workers=args.quant_threads, thread_name_prefix="quant")
    
    # Initialize progress tracker
    tracker = get_tracker()
    logger.info("Progress tracker initialized - visualization will update automatically")
    
    # Initialize tracker for all active species
    for species_name, config_path in active_species:
        try:
            species_id = get_species_id_from_config(config_path)
            cfg = load_workflow_config(config_path)
            
            # Get metadata to find all sample IDs
            metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
            if not metadata_file.exists():
                metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
            
            if metadata_file.exists():
                rows = list(read_delimited(metadata_file, delimiter="\t"))
                sample_ids = [row.get("run") for row in rows if row.get("run")]
                if sample_ids:
                    tracker.initialize_species(species_id, len(sample_ids), sample_ids)
                    logger.debug(f"Initialized tracker for {species_id} with {len(sample_ids)} samples")
        except Exception as e:
            logger.debug(f"Failed to initialize tracker for {species_name}: {e}")
    
    # Start initial processes
    for species_name, config_path in active_species:
        threads = thread_allocation[species_name]
        try:
            process, log_file = start_download_process(config_path, species_name, threads, log_dir, retry_num=0, auto_cleanup_threshold=auto_cleanup_threshold)
            running_processes[species_name] = process
            process_configs[species_name] = config_path
            process_start_times[species_name] = time.time()
            process_last_activity[species_name] = time.time()
            process_retries[species_name] = 0
            process_log_files[species_name] = log_file
            logger.info(f"Started {species_name} (PID: {process.pid}, {threads} threads, log: {log_file.name})")
        except Exception as e:
            logger.error(f"Failed to start {species_name}: {e}")
            results[species_name] = {
                "success": False,
                "message": f"Failed to start: {e}",
                "stats": {},
            }
    
    print(f"\nStarted {len(running_processes)} download processes")
    print(f"Quant thread pool: {args.quant_threads} threads")
    print(f"Monitoring downloads every {args.monitor_interval} seconds...")
    print(f"Checking completion every {args.check_interval} seconds...")
    print("=" * 80)
    print()
    
    # Monitor and redistribute
    iteration = 0
    last_monitor_check = time.time()
    last_full_check = time.time()
    last_health_check = time.time()
    
    while running_processes:
        iteration += 1
        current_time = time.time()
        
        # Add periodic status logging to ensure loop is running
        if iteration % 100 == 0:
            logger.debug(f"Main loop iteration {iteration}: {len(running_processes)} processes, {len(quant_futures)} quant futures")
        
        # Check for completed quant operations (every iteration)
        completed_quants = []
        for sample_key, future in list(quant_futures.items()):
            if future.done():
                completed_quants.append(sample_key)
                try:
                    success, message = future.result(timeout=0.1)  # Should be instant if done()
                    species_name, sample_id = sample_key.split(":", 1)
                    if success:
                        logger.info(f"✅ {species_name}/{sample_id}: {message}")
                    else:
                        logger.warning(f"⚠️  {species_name}/{sample_id}: {message}")
                except Exception as e:
                    logger.error(f"❌ Error in quant operation for {sample_key}: {e}", exc_info=True)
        
        # Remove completed futures
        for sample_key in completed_quants:
            del quant_futures[sample_key]
        
        # Monitor for completed downloads (more frequent than completion checks)
        do_monitor = (current_time - last_monitor_check) >= args.monitor_interval
        do_full_check = (current_time - last_full_check) >= args.check_interval
        do_health_check = (current_time - last_health_check) >= args.health_check_interval
        
        # Health check for hung processes and disk space
        if do_health_check:
            last_health_check = current_time
            logger.info(f"Performing health check on {len(running_processes)} processes...")
            
            # Check disk space
            disk_info = get_disk_space_info(output_dir)
            if disk_info["free_gb"] < auto_cleanup_threshold:
                logger.warning(f"⚠️  Low disk space: {disk_info['free_gb']:.1f}GB free (threshold: {auto_cleanup_threshold}GB)")
                logger.info("Attempting automatic cleanup of partial downloads...")
                try:
                    from scripts.rna.cleanup_partial_downloads import cleanup_species
                    total_freed = 0
                    for species_name_check, config_path in [(name, process_configs.get(name)) for name in running_processes.keys()]:
                        if config_path and config_path.exists():
                            result = cleanup_species(species_name_check, config_path, dry_run=False)
                            total_freed += result["freed_mb"]
                    if total_freed > 0:
                        logger.info(f"✅ Freed {total_freed/1024:.2f}GB from partial downloads")
                    else:
                        logger.warning("⚠️  No partial downloads to clean up")
                except Exception as e:
                    logger.warning(f"⚠️  Cleanup failed: {e}")
            
            unhealthy = []
            for species_name in list(running_processes.keys()):
                process = running_processes[species_name]
                log_file = process_log_files.get(species_name)
                
                if log_file:
                    is_healthy, reason = check_process_health(process, log_file, args.process_timeout)
                    if not is_healthy:
                        unhealthy.append((species_name, reason))
                        logger.warning(f"⚠️  {species_name}: {reason}")
            
            # Restart unhealthy processes
            for species_name, reason in unhealthy:
                if species_name not in running_processes:
                    continue
                
                process = running_processes[species_name]
                retry_count = process_retries.get(species_name, 0)
                
                if retry_count >= args.max_retries:
                    logger.error(f"❌ {species_name}: Max retries ({args.max_retries}) reached, giving up")
                    process.terminate()
                    del running_processes[species_name]
                    results[species_name] = {
                        "success": False,
                        "message": f"Failed after {retry_count} retries: {reason}",
                        "stats": {},
                    }
                    continue
                
                logger.warning(f"🔄 Restarting {species_name} (retry {retry_count + 1}/{args.max_retries})")
                
                # Terminate old process
                try:
                    process.terminate()
                    time.sleep(2)
                    if process.poll() is None:
                        process.kill()
                except Exception as e:
                    logger.warning(f"Error terminating {species_name}: {e}")
                
                # Restart process
                try:
                    config_path = process_configs[species_name]
                    threads = thread_allocation.get(species_name, 1)
                    new_process, new_log_file = start_download_process(
                        config_path, species_name, threads, log_dir, retry_num=retry_count + 1, auto_cleanup_threshold=auto_cleanup_threshold
                    )
                    running_processes[species_name] = new_process
                    process_start_times[species_name] = time.time()
                    process_last_activity[species_name] = time.time()
                    process_retries[species_name] = retry_count + 1
                    process_log_files[species_name] = new_log_file
                    logger.info(f"✅ Restarted {species_name} (new PID: {new_process.pid})")
                except Exception as e:
                    logger.error(f"❌ Failed to restart {species_name}: {e}")
                    del running_processes[species_name]
                    results[species_name] = {
                        "success": False,
                        "message": f"Failed to restart: {e}",
                        "stats": {},
                    }
        
        if do_monitor:
            last_monitor_check = current_time
            
            # Detect completed downloads (both FASTQ and SRA)
            completed_downloads = detect_completed_downloads(
                [(name, process_configs[name]) for name in running_processes.keys()],
                processed_samples,
                quant_dirs,
                tracker=tracker,
            )
            
            # Also check ALL species for backlogged completed downloads
            all_species_configs = [
                (name, cfg) for name, cfg in process_configs.items()
            ]
            backlogged = detect_completed_downloads(
                all_species_configs,
                processed_samples,
                quant_dirs,
                tracker=tracker,
            )
            # Add backlogged samples (avoid duplicates)
            for item in backlogged:
                if item not in completed_downloads:
                    completed_downloads.append(item)
            
            # Submit processing jobs for newly completed downloads
            for species_name, sample_id, sample_dir, status in completed_downloads:
                sample_key = f"{species_name}:{sample_id}"
                
                # Skip if already queued
                if sample_key in quant_futures:
                    continue
                
                # Submit to processing thread pool using metainformant function
                config_path = process_configs.get(species_name)
                if config_path:
                    try:
                        future = quant_executor.submit(
                            process_sample_with_tracker,
                            sample_id,
                            config_path,
                            status,
                            tracker=tracker,
                            log_dir=log_dir,
                        )
                        quant_futures[sample_key] = future
                        if status == "fastq":
                            logger.info(f"📥 {species_name}/{sample_id}: FASTQ ready, queued for quant+delete")
                        else:
                            logger.info(f"📦 {species_name}/{sample_id}: SRA complete, queued for SRA→FASTQ→quant→delete")
                    except Exception as e:
                        logger.error(f"❌ Failed to queue {species_name}/{sample_id} for processing: {e}", exc_info=True)
                else:
                    logger.warning(f"⚠️  No config path found for {species_name}, skipping {sample_id}")
        
        if do_full_check:
            last_full_check = current_time
            print(f"\n[{datetime.now().strftime('%H:%M:%S')}] Check #{iteration} (Full)")
            print("-" * 80)
            
            # Check for completed processes (only on full check)
            completed = []
            for species_name, process in list(running_processes.items()):
                if process.poll() is not None:
                    # Process finished
                    returncode = process.returncode
                    completed.append(species_name)
                    retry_count = process_retries.get(species_name, 0)
                    
                    if returncode == 0:
                        status = "✅"
                        message = "Completed successfully"
                        results[species_name] = {
                            "success": True,
                            "message": message,
                            "stats": {"return_code": returncode, "retries": retry_count},
                        }
                    else:
                        # Check if we should retry
                        if retry_count < args.max_retries:
                            logger.warning(f"⚠️  {species_name}: Process failed with code {returncode}, will retry")
                            # Restart process
                            try:
                                config_path = process_configs[species_name]
                                threads = thread_allocation.get(species_name, 1)
                                new_process, new_log_file = start_download_process(
                                    config_path, species_name, threads, log_dir, retry_num=retry_count + 1, auto_cleanup_threshold=auto_cleanup_threshold
                                )
                                running_processes[species_name] = new_process
                                process_start_times[species_name] = time.time()
                                process_last_activity[species_name] = time.time()
                                process_retries[species_name] = retry_count + 1
                                process_log_files[species_name] = new_log_file
                                logger.info(f"🔄 Restarted {species_name} after failure (retry {retry_count + 1}/{args.max_retries})")
                                continue  # Don't mark as completed
                            except Exception as e:
                                logger.error(f"❌ Failed to restart {species_name}: {e}")
                                status = "❌"
                                message = f"Failed after {retry_count} retries: {e}"
                        else:
                            status = "❌"
                            message = f"Failed after {retry_count} retries (code {returncode})"
                        
                        results[species_name] = {
                            "success": False,
                            "message": message,
                            "stats": {"return_code": returncode, "retries": retry_count},
                        }
                    
                    logger.info(f"{status} {species_name}: {message}")
                    if species_name in running_processes:
                        del running_processes[species_name]
            
            # Check for newly-complete species (even if process still running)
            # This handles cases where process is still running but work is done
            newly_complete = []
            for species_name in list(running_processes.keys()):
                config_path = process_configs[species_name]
                if is_species_complete(config_path):
                    newly_complete.append(species_name)
                    process = running_processes[species_name]
                    process.terminate()  # Stop the process since work is done
                    logger.info(f"✅ {species_name}: Work complete (terminating process)")
                    results[species_name] = {
                        "success": True,
                        "message": "Completed (detected completion)",
                        "stats": {"return_code": 0},
                    }
                    del running_processes[species_name]
            
            # Redistribute threads if any completed
            if completed or newly_complete:
                remaining_active = [
                    (name, cfg) for name, cfg in active_species
                    if name not in results or not results[name].get("success", False)
                ]
                
                if remaining_active:
                    # Get current active (running) species
                    currently_running = [
                        (name, process_configs[name])
                        for name in running_processes.keys()
                    ]
                    
                    # Calculate threads to redistribute
                    threads_in_use = sum(thread_allocation.get(name, 0) for name in running_processes.keys())
                    threads_available = args.total_threads - threads_in_use
                    
                    # Redistribute among remaining active
                    new_allocation = redistribute_threads(args.total_threads, remaining_active)
                    
                    # Update allocations for running processes
                    for species_name in list(running_processes.keys()):
                        old_threads = thread_allocation.get(species_name, 1)
                        new_threads = new_allocation.get(species_name, old_threads)
                        
                        if new_threads != old_threads:
                            logger.info(f"  {species_name}: {old_threads} → {new_threads} threads")
                            # Note: Can't change threads on running process, but will use new allocation when restarting
                            thread_allocation[species_name] = new_threads
                    
                    # Validate redistribution
                    redistributed_total = sum(new_allocation.values())
                    if redistributed_total != args.total_threads:
                        logger.warning(f"⚠️  Redistribution mismatch: {redistributed_total} allocated, {args.total_threads} requested")
                    else:
                        logger.debug(f"✅ Redistribution verified: {redistributed_total} threads")
                    
                    print(f"Redistributed threads: {len(completed + newly_complete)} completed, "
                          f"{len(running_processes)} still running ({redistributed_total} threads total)")
                else:
                    print(f"All species complete or in progress")
            
            # Status update (on full check)
            if running_processes:
                print(f"Active downloads: {len(running_processes)} species")
                print(f"Quant queue: {len(quant_futures)} samples")
                for species_name in sorted(running_processes.keys()):
                    threads = thread_allocation.get(species_name, 1)
                    process = running_processes[species_name]
                    pid = process.pid
                    status = "running" if process.poll() is None else "finished"
                    print(f"  {species_name:30s}: PID {pid}, {threads} threads, {status}")
            else:
                # Wait for remaining quant operations
                if quant_futures:
                    print(f"Downloads complete, waiting for {len(quant_futures)} quant operations...")
                else:
                    print("All processes completed")
                    break
        else:
            # Monitor-only check - just show quant queue
            if quant_futures:
                print(f"Quant queue: {len(quant_futures)} samples")
        
        # Sleep until next monitor interval
        elapsed = time.time() - current_time
        sleep_time = max(0.1, args.monitor_interval - elapsed)
        if sleep_time > 0:
            time.sleep(sleep_time)
    
    # Wait for all quant operations to complete
    if quant_futures:
        print("\nWaiting for remaining quant operations to complete...")
        for sample_key, future in quant_futures.items():
            try:
                success, message = future.result()
                species_name, sample_id = sample_key.split(":", 1)
                if success:
                    logger.info(f"✅ {species_name}/{sample_id}: {message}")
                else:
                    logger.warning(f"⚠️  {species_name}/{sample_id}: {message}")
            except Exception as e:
                logger.error(f"❌ Error in quant operation for {sample_key}: {e}")
    
    # Shutdown quant executor
    quant_executor.shutdown(wait=True)
    
    # Final summary
    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    print(f"Total species: {len(species_configs)}")
    print(f"Already complete: {len(complete_species)}")
    print(f"Processed: {len(active_species)}")
    print(f"Samples processed: {len(processed_samples)}")
    print()
    
    success_count = sum(1 for r in results.values() if r.get("success", False))
    failed_count = len(results) - success_count
    
    print(f"Successfully completed: {success_count}")
    print(f"Failed: {failed_count}")
    print()
    
    if results:
        print("Per-species results:")
        for species_name, result in sorted(results.items()):
            status = "✅" if result.get("success", False) else "❌"
            message = result.get("message", "Unknown")
            print(f"  {status} {species_name:30s}: {message}")
    
    print("=" * 80)
    
    return 0 if failed_count == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

