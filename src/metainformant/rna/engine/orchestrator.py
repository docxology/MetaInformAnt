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
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from metainformant.core.utils.logging import get_logger
from metainformant.rna.retrieval.ena_downloader import ENADownloader

logger = get_logger(__name__)


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
        workers: int = 4,
        threads_per_sample: int = 2,
        dry_run: bool = False,
    ):
        """
        Initialize the pipeline.

        Args:
            species: Species name (e.g., "apis_mellifera").
            index_file: Path to Kallisto index.
            work_dir: Main Amalgkit working directory.
            fastq_dir: Temporary directory for FASTQ downloads (high I/O).
            workers: Number of parallel downloads/quantifications.
            threads_per_sample: CPU threads for Kallisto per sample.
            dry_run: If True, simulate actions without executing.
        """
        self.species = species
        self.index_file = Path(index_file)
        self.work_dir = Path(work_dir)
        self.fastq_dir = Path(fastq_dir)
        self.workers = workers
        self.threads = threads_per_sample
        self.dry_run = dry_run
        
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

        try:
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
            if res.returncode != 0:
                return False, f"Kallisto failed: {res.stderr[:200]}"
            return True, "Quantified"
        except Exception as e:
            return False, str(e)

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
        
        if self.dry_run:
            logger.info("[DRY RUN] Simulating processing for first 5 samples only.")
            sample_ids = sample_ids[:5]

        # Use ProcessPoolExecutor for true parallelism
        # Note: self.process_one must be picklable. 
        # Since it's an instance method, pickle sends 'self'. 
        # This is fine as long as 'self' state is picklable (paths, ints, bools are fine).
        
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
