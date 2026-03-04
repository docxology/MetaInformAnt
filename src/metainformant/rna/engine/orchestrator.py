"""
RNA-Seq Pipeline Orchestrator.

This module provides the `StreamingPipeline` class, which manages high-throughput
RNA-seq processing workflows. It prioritizes the "streaming" approach:
1. Download FASTQ (ENA/Direct) to temp
2. Quantify (Kallisto)
3. Cleanup FASTQ immediately

This minimizes disk usage and maximizes throughput compared to the traditional
"download all, then process all" approach.
"""

import logging
import shutil
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from metainformant.core.utils.logging import get_logger
from metainformant.core.utils.watchdog import ProcessWatchdog
from metainformant.rna.retrieval.ena_downloader import ENADownloader

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Resource Profiles — species-specific tuning for workers, threads, timeouts
# ---------------------------------------------------------------------------

@dataclass
class ResourceProfile:
    """Resource allocation profile for pipeline execution.

    Presets:
        - ``high_mem``: Low concurrency, high threads per sample (complex genomes).
        - ``high_io``: High concurrency, moderate threads (large sample counts).
        - ``default``: Balanced settings suitable for most species.
        - ``minimal``: Conservative settings for constrained environments.
    """

    name: str = "default"
    workers: int = 4
    threads_per_sample: int = 2
    watchdog_timeout: int = 3600
    use_watchdog: bool = True
    description: str = ""


# Pre-defined profile presets
RESOURCE_PROFILES: Dict[str, ResourceProfile] = {
    "high_mem": ResourceProfile(
        name="high_mem",
        workers=2,
        threads_per_sample=8,
        watchdog_timeout=7200,
        use_watchdog=True,
        description="Complex genomes (e.g., Harpegnathos): low concurrency, high per-sample threads",
    ),
    "high_io": ResourceProfile(
        name="high_io",
        workers=16,
        threads_per_sample=2,
        watchdog_timeout=3600,
        use_watchdog=True,
        description="Large datasets (e.g., Apis mellifera): max concurrency, moderate threads",
    ),
    "default": ResourceProfile(
        name="default",
        workers=4,
        threads_per_sample=2,
        watchdog_timeout=3600,
        use_watchdog=True,
        description="Balanced settings for typical ant species",
    ),
    "minimal": ResourceProfile(
        name="minimal",
        workers=2,
        threads_per_sample=1,
        watchdog_timeout=1800,
        use_watchdog=False,
        description="Resource-constrained environments (CI, small VMs)",
    ),
}

# Auto-detection: species name → profile name
_SPECIES_PROFILE_MAP: Dict[str, str] = {
    "harpegnathos_saltator": "high_mem",
    "apis_mellifera": "high_io",
    "amellifera": "high_io",
}


def get_profile(name: Optional[str] = None, species: Optional[str] = None) -> ResourceProfile:
    """Return a ``ResourceProfile`` by explicit name or auto-detected from species.

    Args:
        name: Explicit profile name (``high_mem``, ``high_io``, ``default``, ``minimal``).
              Takes priority over *species* if both are provided.
        species: Species name for auto-detection.  Falls back to ``"default"``.

    Returns:
        A copy of the matching ``ResourceProfile``.
    """
    if name and name in RESOURCE_PROFILES:
        profile = RESOURCE_PROFILES[name]
    elif species:
        key = species.lower().replace(" ", "_")
        mapped = _SPECIES_PROFILE_MAP.get(key, "default")
        profile = RESOURCE_PROFILES[mapped]
    else:
        profile = RESOURCE_PROFILES["default"]

    # Return a fresh copy so callers can mutate freely
    return ResourceProfile(
        name=profile.name,
        workers=profile.workers,
        threads_per_sample=profile.threads_per_sample,
        watchdog_timeout=profile.watchdog_timeout,
        use_watchdog=profile.use_watchdog,
        description=profile.description,
    )

class StreamingPipeline:
    """
    Orchestrates streaming RNA-seq processing (Download -> Quant -> Delete).
    """

    def __init__(
        self,
        species: str,
        index_file: Path,
        work_dir: Path,
        fastq_dir: Path,
        workers: Optional[int] = None,
        threads_per_sample: Optional[int] = None,
        dry_run: bool = False,
        use_watchdog: Optional[bool] = None,
        watchdog_timeout: Optional[int] = None,
        profile: Optional[str] = None,
    ):
        """
        Initialize the pipeline.

        Args:
            species: Species name (e.g., "apis_mellifera").
            index_file: Path to Kallisto index.
            work_dir: Main Amalgkit working directory.
            fastq_dir: Temporary directory for FASTQ downloads (high I/O).
            workers: Number of parallel downloads/quantifications.
                     Defaults to profile setting if not specified.
            threads_per_sample: CPU threads for Kallisto per sample.
                     Defaults to profile setting if not specified.
            dry_run: If True, simulate actions without executing.
            use_watchdog: ProcessWatchdog to kill stalled processes (0 CPU).
                     Defaults to profile setting if not specified.
            watchdog_timeout: Seconds of inactivity before kill.
                     Defaults to profile setting if not specified.
            profile: Resource profile name ("high_mem", "high_io", "default",
                     "minimal"). If None, auto-detected from species name.
        """
        # Resolve resource profile (explicit name > species auto-detect > default)
        rp = get_profile(name=profile, species=species)
        logger.info(f"Resource profile: {rp.name} ({rp.description})")

        self.species = species
        self.index_file = Path(index_file)
        self.work_dir = Path(work_dir)
        self.fastq_dir = Path(fastq_dir)
        self.workers = workers if workers is not None else rp.workers
        self.threads = threads_per_sample if threads_per_sample is not None else rp.threads_per_sample
        self.dry_run = dry_run
        self.use_watchdog = use_watchdog if use_watchdog is not None else rp.use_watchdog
        self.watchdog_timeout = watchdog_timeout if watchdog_timeout is not None else rp.watchdog_timeout
        self.profile = rp
        
        # Computed paths
        self.quant_dir = self.work_dir / "quant"
        self.downloader = ENADownloader()

        # Ensure directories exist
        if not self.dry_run:
            self.work_dir.mkdir(parents=True, exist_ok=True)
            self.quant_dir.mkdir(parents=True, exist_ok=True)
            self.fastq_dir.mkdir(parents=True, exist_ok=True)

    def get_processed_samples(self) -> Set[str]:
        """Return set of sample IDs that have already been quantified."""
        processed = set()
        if not self.quant_dir.exists():
            return processed
        
        for sample_dir in self.quant_dir.iterdir():
            if sample_dir.is_dir() and (sample_dir / "abundance.tsv").exists():
                processed.add(sample_dir.name)
        return processed

    def _quantify_sample(self, sample_id: str, fastq_files: List[Path]) -> Tuple[bool, str]:
        """Run kallisto quant for a single sample."""
        output_dir = self.quant_dir / sample_id
        
        if self.dry_run:
            logger.info(f"[DRY RUN] Would quantify {sample_id} with {fastq_files} -> {output_dir}")
            return True, "Dry Run"

        output_dir.mkdir(parents=True, exist_ok=True)
        
        cmd = [
            "kallisto", "quant",
            "-i", str(self.index_file),
            "-o", str(output_dir),
            "-t", str(self.threads),
        ]
        
        # Handle Single vs Paired
        if len(fastq_files) >= 2:
             cmd.extend([str(fastq_files[0]), str(fastq_files[1])])
        else:
             # Estimated fragment length for single-end
             cmd.extend(["--single", "-l", "200", "-s", "30", str(fastq_files[0])])

        proc = None
        watchdog = None
        
        try:
            # Start process using Popen interaction to allow monitoring
            proc = subprocess.Popen(
                cmd, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                text=True
            )
            
            # Start Watchdog if enabled
            if self.use_watchdog:
                watchdog = ProcessWatchdog(
                    pid=proc.pid, 
                    cpu_threshold=1.0, 
                    timeout_seconds=self.watchdog_timeout,
                    on_stall=ProcessWatchdog.kill_process_tree
                )
                watchdog.start()
            
            # Wait for completion (simulates subprocess.run)
            stdout, stderr = proc.communicate(timeout=7200) # Hard timeout 2h just in case
            
            if proc.returncode != 0:
                # Check if it was killed by watchdog/signal?
                if proc.returncode == -9 or proc.returncode == -15:
                    return False, "Killed (likely by Watchdog)"
                return False, f"Kallisto failed: {stderr[:200] if stderr else 'No stderr'}"
            
            return True, "Quantified"

        except subprocess.TimeoutExpired:
            if proc:
                ProcessWatchdog.kill_process_tree(proc.pid)
            return False, "Hard Timeout (2h)"
        
        except Exception as e:
            if proc:
                ProcessWatchdog.kill_process_tree(proc.pid)
            return False, str(e)
            
        finally:
            if watchdog:
                watchdog.stop()

    def process_one(self, sample_id: str) -> Dict:
        """
        Full lifecycle for one sample: Download -> Quant -> Clean.
        Designed to be run in a parallel worker.
        """
        result = {
            "sample_id": sample_id, 
            "success": False, 
            "timing": {}, 
            "error": None
        }
        
        sample_tmp_dir = self.fastq_dir / sample_id
        t0 = time.time()

        try:
            # 1. Check if done
            if not self.dry_run and (self.quant_dir / sample_id / "abundance.tsv").exists():
                return {**result, "success": True, "status": "Skipped (Exists)"}

            # 2. Download
            if self.dry_run:
                logger.info(f"[DRY RUN] Downloading {sample_id}...")
                fastq_files = [sample_tmp_dir / f"{sample_id}_1.fastq.gz"]
                ok, msg = True, "Dry Run"
            else:
                ok, msg, fastq_files = self.downloader.download_run(sample_id, sample_tmp_dir)

            result["timing"]["download"] = time.time() - t0
            
            if not ok:
                result["error"] = f"Download failed: {msg}"
                return result

            # 3. Quantify
            t1 = time.time()
            ok, msg = self._quantify_sample(sample_id, fastq_files)
            result["timing"]["quantify"] = time.time() - t1
            
            if not ok:
                result["error"] = msg
                return result

            # 4. Cleanup
            if not self.dry_run:
                shutil.rmtree(sample_tmp_dir, ignore_errors=True)
            
            result["success"] = True

        except Exception as e:
            result["error"] = str(e)
            # Cleanup on failure
            if not self.dry_run:
                shutil.rmtree(sample_tmp_dir, ignore_errors=True)
        
        return result

    def run(self, sample_ids: List[str]):
        """Run the pipeline for a list of samples."""
        total = len(sample_ids)
        logger.info(f"Starting pipeline for {total} samples with {self.workers} workers.")
        logger.info(f"Watchdog: {'Enabled' if self.use_watchdog else 'Disabled'} (Timeout: {self.watchdog_timeout}s)")
        
        if self.dry_run:
            logger.info("[DRY RUN] Simulating processing for first 5 samples only.")
            sample_ids = sample_ids[:5]

        # Use ProcessPoolExecutor for true parallelism
        with ProcessPoolExecutor(max_workers=self.workers) as executor:
            future_to_sample = {
                executor.submit(self.process_one, sid): sid 
                for sid in sample_ids
            }
            
            completed = 0
            for future in as_completed(future_to_sample):
                completed += 1
                res = future.result()
                sid = res["sample_id"]
                
                if res["success"]:
                    logger.info(f"[{completed}/{total}] ✓ {sid}")
                else:
                    logger.warning(f"[{completed}/{total}] ✗ {sid}: {res['error']}")
