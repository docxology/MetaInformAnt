
"""Streaming RNA-seq Orchestrator implementation.

This module provides robust orchestration for running the amalgkit pipeline 
with ENA-first download streaming, concurrent processing, and automatic cleanup.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import logging
import os
import subprocess
import sys
import time
import urllib.request
import yaml

from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple

from metainformant.core.utils import logging as log_utils
from metainformant.rna.amalgkit import amalgkit
from metainformant.rna.amalgkit.tissue_normalizer import apply_tissue_normalization
from metainformant.rna.retrieval.ena_downloader import ENADownloader
from metainformant.rna.engine.progress_db import ProgressDB

# Determine project root if possible, or use relative paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
CONFIG_DIR = PROJECT_ROOT / "config/amalgkit"
# Use local output for logs to avoid TCC/Permission issues on external volumes
LOG_DIR = PROJECT_ROOT / "output/amalgkit"


logger = log_utils.get_logger(__name__)

class StreamingPipelineOrchestrator:
    """Orchestrator for streaming ENA download -> Quant processing."""

    def __init__(self, config_dir: Path = CONFIG_DIR, log_dir: Path = LOG_DIR,
                 db_path: Optional[Path] = None):
        self.config_dir = Path(config_dir)
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        
        # SQLite progress database
        self.db = ProgressDB(db_path) if db_path else ProgressDB()
        
        # Configure logging to file as well
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        log_file = self.log_dir / f"streaming_orchestrator_{timestamp}.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
        logger.addHandler(file_handler)

    def query_ena_fastq_urls(self, srr_id: str) -> List[str]:
        """Query ENA API for direct FASTQ download URLs."""
        url = (
            f"https://www.ebi.ac.uk/ena/portal/api/filereport"
            f"?accession={srr_id}&result=read_run"
            f"&fields=run_accession,fastq_ftp&format=tsv"
        )
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                lines = resp.read().decode().strip().split('\n')
                if len(lines) < 2:
                    return []
                ftp_field = lines[1].split('\t')[1] if '\t' in lines[1] else ''
                if not ftp_field:
                    return []
                # Ensure protocol is https
                return [f"https://{p}" if not p.startswith("http") else p for p in ftp_field.split(';') if p]
        except Exception as e:
            logger.warning(f"ENA query failed for {srr_id}: {e}")
            return []

    def download_fastq(self, srr_id: str, out_dir: Path) -> bool:
        """Download FASTQ files using ENADownloader.

        Uses the ENA Portal API to discover FASTQ URLs and downloads
        via curl with retries.

        Args:
            srr_id: SRA accession ID (e.g., SRR12345)
            out_dir: Base directory for downloads (sample subdir created automatically)

        Returns:
            True if all FASTQ files downloaded successfully, False otherwise.
        """
        sample_dir = out_dir / srr_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        downloader = ENADownloader(timeout=7200, retries=3)
        success, message, downloaded_files = downloader.download_run(srr_id, sample_dir)

        if success:
            for fq in downloaded_files:
                sz_gb = fq.stat().st_size / (1024**3)
                logger.info(f"Downloaded {fq.name}: {sz_gb:.2f} GB")
        else:
            logger.error(f"Download failed for {srr_id}: {message}")
            # Clean up partial files
            for f in sample_dir.iterdir():
                if f.is_file():
                    f.unlink(missing_ok=True)

        return success

    def quant_sample(self, config_path: Path, batch_index: int, species_name: str, threads: int) -> tuple:
        """Run amalgkit quant for a single sample using --batch.
        
        Returns:
            Tuple of (success: bool, error_msg: str or None).
        """
        with open(config_path) as f:
            cfg = yaml.safe_load(f)

        steps_cfg = cfg.get("steps", {})
        quant_cfg = steps_cfg.get("quant", {})
        # Always use the work dir for quant output so results land in the expected location
        quant_out = f"output/amalgkit/{species_name}/work"
        
        # Determine metadata path
        work_dir = Path(f"output/amalgkit/{species_name}/work")
        meta_path = None
        for mp in [work_dir / "metadata/metadata.tsv", work_dir / "metadata/metadata_selected.tsv"]:
            if mp.exists():
                meta_path = str(mp)
                break
        
        if not meta_path:
            logger.error(f"No metadata found for {species_name} batch {batch_index}")
            return False, f"No metadata found for {species_name}"

        # Build command
        cmd = [
            "amalgkit", "quant",
            "--out_dir", quant_out,
            "--metadata", meta_path,
            "--threads", str(threads),
            "--batch", str(batch_index),
        ]
        
        # Handle cleanup settings
        keep_fastq = str(quant_cfg.get("keep_fastq", "yes")).lower()
        if keep_fastq in ("no", "false"):
            cmd.extend(["--clean_fastq", "yes"])
        else:
            cmd.extend(["--clean_fastq", "no"])

        # Add genome/index paths
        genome_cfg = cfg.get("genome", {})
        index_dir = genome_cfg.get("index_dir", "")
        if index_dir and Path(index_dir).exists():
            cmd.extend(["--index_dir", index_dir])
        
        log_path = self.log_dir / f"{species_name}_quant.log"
        
        try:
            with open(log_path, 'a') as log_f:
                log_f.write(f"\n[{time.strftime('%Y-%m-%d %H:%M:%S')}] Quant batch {batch_index} START: {' '.join(cmd)}\n")
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
            
            with open(log_path, 'a') as log_f:
                log_f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Quant batch {batch_index} END (Exit {result.returncode})\n")
                if result.stdout: log_f.write(f"-- STDOUT --\n{result.stdout}\n")
                if result.stderr: log_f.write(f"-- STDERR --\n{result.stderr}\n")
            
            if result.returncode == 0:
                return True, None
            
            # Extract meaningful error from amalgkit output
            error_msg = "Quantification Failed"
            for line in (result.stdout + result.stderr).splitlines():
                line_lower = line.strip().lower()
                if any(kw in line_lower for kw in ["error", "exiting", "no sample", "not found", "failed", "exception"]):
                    error_msg = line.strip()[:120]
                    break
            return False, error_msg
        except subprocess.TimeoutExpired:
            logger.error(f"Quant timeout batch {batch_index} after 2h")
            with open(log_path, 'a') as log_f:
                log_f.write(f"TIMEOUT after 7200s\n")
            return False, "Quant timeout (>2h)"
        except Exception as e:
            logger.error(f"Quant exception batch {batch_index}: {e}")
            with open(log_path, 'a') as log_f:
                log_f.write(f"EXCEPTION: {e}\n")
            return False, f"Exception: {e}"

    def is_quantified(self, species_name: str, srr_id: str) -> bool:
        """Check if sample is already quantified (DB-first, filesystem fallback)."""
        db_state = self.db.get_state(species_name, srr_id)
        if db_state == "quantified":
            return True
        # Filesystem fallback for samples not yet in DB
        quant_dir = Path(f"output/amalgkit/{species_name}/work/quant/{srr_id}")
        return any(quant_dir.glob("*_abundance.tsv")) if quant_dir.exists() else False

    def process_single_sample(self, srr_id: str, batch_idx: int, fastq_dir: Path, 
                            config_path: Path, species_name: str, threads: int) -> Dict[str, Any]:
        """Processing unit for a single sample."""
        result = {
            "srr": srr_id, "batch": batch_idx, 
            "downloaded": False, "quantified": False, 
            "skipped": False, "error": None
        }

        if self.is_quantified(species_name, srr_id):
            result.update({"downloaded": True, "quantified": True, "skipped": True})
            return result
        
        # Download
        self.db.set_state(species_name, srr_id, "downloading")
        if not self.download_fastq(srr_id, fastq_dir):
            error_msg = "Download Failed (all sources: ENA FTP/HTTP, NCBI)"
            self.db.set_state(species_name, srr_id, "failed", error=error_msg)
            result["error"] = error_msg
            return result
        self.db.set_state(species_name, srr_id, "downloaded")
        result["downloaded"] = True
        
        # Quantify
        self.db.set_state(species_name, srr_id, "quantifying")
        success, quant_error = self.quant_sample(config_path, batch_idx, species_name, threads)
        if success:
            self.db.set_state(species_name, srr_id, "quantified")
        else:
            error_msg = quant_error or "Quantification Failed"
            self.db.set_state(species_name, srr_id, "failed", error=error_msg)
        result["quantified"] = success
        if not success:
            result["error"] = quant_error or "Quantification Failed"
        
        return result

    def verify_genome_index(self, config_path: Path, species_name: str) -> bool:
        """Verify Kallisto index exists."""
        try:
            with open(config_path) as f:
                cfg = yaml.safe_load(f)
        except Exception as e:
            logger.error(f"Failed to load config {config_path}: {e}")
            return False
            
        # Check multiple possible index locations
        # 1. From quant config (has the correct capitalized path)
        quant_index_dir = cfg.get('steps', {}).get('quant', {}).get('index_dir', '')
        # 2. From genome config dest_dir + /index
        genome_dest = cfg.get('genome', {}).get('dest_dir', '')
        genome_index_dir = f"{genome_dest}/index" if genome_dest else ''
        # 3. Legacy genome.index_dir
        legacy_index_dir = cfg.get('genome', {}).get('index_dir', '')
        
        search_dirs = [
            quant_index_dir,
            genome_index_dir,
            legacy_index_dir,
            f"output/amalgkit/{species_name}/work/index",
            f"output/amalgkit/shared/genome/{species_name}/index",
        ]
        
        logger.info(f"Verifying index for {species_name}...")
        for d in search_dirs:
            if not d:
                continue
                
            path_obj = Path(d)
            abs_path = path_obj.resolve() if path_obj.exists() else path_obj.absolute()
            
            logger.info(f"  Checking {d} -> {abs_path}")
            
            if not path_obj.exists():
                logger.warning(f"  Path does not exist: {d}")
                continue
                
            # Check permissions by attempting listdir
            try:
                files = list(path_obj.glob("*.idx"))
                if files:
                    logger.info(f"  Genome index found in {d}: {[f.name for f in files]}")
                    return True
                else:
                    logger.warning(f"  Path exists but no .idx files found in {d}")
                    try:
                        contents = os.listdir(d) 
                        logger.info(f"  Directory contents: {contents[:5]}...")
                    except Exception as e:
                        logger.error(f"  Failed to list directory {d}: {e}")
            except Exception as e:
                logger.error(f"  Permission/Error accessing {d}: {e}")
                
        logger.error(f"No genome index found for {species_name}")
        return False

    def run_tissue_normalization(self, metadata_path: Path):
        """Run tissue normalization on the metadata file."""
        mapping_path = self.config_dir / "tissue_mapping.yaml"
        patches_path = self.config_dir / "tissue_patches.yaml"
        
        if not mapping_path.exists():
            logger.warning(f"Tissue mapping not found at {mapping_path}, skipping normalization.")
            return

        logger.info(f"Normalizing tissues in {metadata_path}")
        try:
            import pandas as pd
            df = pd.read_csv(metadata_path, sep="\t", low_memory=False)
            
            # Apply normalization
            df_norm = apply_tissue_normalization(
                df, 
                mapping_path=mapping_path,
                patches_path=patches_path if patches_path.exists() else None,
                output_column="tissue_normalized" # Amalgkit curate looks for 'tissue' usually, but let's be safe
            )
            
            # For curation to work seamlessly without code changes in amalgkit, 
            # we might want to update the 'tissue' column itself or ensure curate uses 'tissue_normalized'
            # The 'curate' step often uses the 'tissue' column. 
            # Strategy: Overwrite 'tissue' with normalized values, keep original in 'tissue_original'
            if "tissue_normalized" in df_norm.columns:
                df_norm["tissue_original"] = df_norm["tissue"]
                df_norm["tissue"] = df_norm["tissue_normalized"]
                # Drop temporary column if desired, or keep it
            
            df_norm.to_csv(metadata_path, sep="\t", index=False)
            logger.info("Tissue normalization complete and saved.")
            
        except Exception as e:
            logger.error(f"Tissue normalization failed: {e}")

    def process_species(self, config_name: str, max_gb: float, workers: int, threads: int) -> bool:
        """Main processing loop for a species."""
        config_path = self.config_dir / config_name
        if not config_path.exists():
            logger.error(f"Config not found: {config_path}")
            return False
            
        species_name = config_name.replace("amalgkit_", "").replace(".yaml", "")
        logger.info(f"=== Processing {species_name} ===")
        
        # Autonomous Preprocessing Detection
        work_dir = Path(f"output/amalgkit/{species_name}/work")
        metadata_path = work_dir / "metadata/metadata.tsv"
        if not metadata_path.exists():
            metadata_path = work_dir / "metadata/metadata_selected.tsv"
            
        needs_prep = False
        if not metadata_path.exists():
            logger.info(f"Metadata missing for {species_name}. Queuing autonomous preprocessing...")
            needs_prep = True
        elif not self.verify_genome_index(config_path, species_name):
            logger.info(f"Genome index missing for {species_name}. Queuing autonomous preprocessing...")
            needs_prep = True
            
        if needs_prep:
            logger.info(f"Running preprocessing stages (config -> select -> metadata -> index) for {species_name}...")
            prep_cmd = [
                "python3", "scripts/rna/run_workflow.py",
                "--config", str(config_path),
                "--no-progress",
                "--steps", "config", "select", "metadata", "index"
            ]
            try:
                prep_result = subprocess.run(prep_cmd, capture_output=True, text=True, timeout=7200)
                if prep_result.returncode != 0:
                    logger.error(f"Preprocessing failed for {species_name}:\n{prep_result.stderr}\n{prep_result.stdout}")
                    return False
                logger.info(f"Preprocessing complete for {species_name}.")
            except subprocess.TimeoutExpired:
                logger.error(f"Preprocessing timed out after 2h for {species_name}.")
                return False
            except Exception as e:
                logger.error(f"Preprocessing exception for {species_name}: {e}")
                return False
                
            # Post-prep explicit verification
            if not self.verify_genome_index(config_path, species_name):
                logger.error(f"Genome index still missing for {species_name} after successful prep spawn.")
                return False

            metadata_path = work_dir / "metadata/metadata.tsv"
            if not metadata_path.exists():
                metadata_path = work_dir / "metadata/metadata_selected.tsv"
                
            if not metadata_path.exists():
                logger.error(f"No metadata generated for {species_name} after successful prep spawn.")
                return False

        # Apply Normalization Step
        self.run_tissue_normalization(metadata_path)

        # Load samples
        import pandas as pd
        df = pd.read_csv(metadata_path, sep="\t")
        
        # Filter logic (simplified from run_all_species.py)
        # Ensure total_bases is numeric
        df["total_bases"] = pd.to_numeric(df["total_bases"], errors="coerce").fillna(0)
        
        # Filter by size
        max_bases = max_gb * 1e9
        filtered = df[df["total_bases"] <= max_bases].copy()
        filtered = filtered.sort_values("total_bases")
        
        logger.info(f"Samples: {len(df)} total -> {len(filtered)} filtered (<= {max_gb} GB)")
        
        if filtered.empty:
            logger.info("No samples to process.")
            return True

        # Register samples in progress DB and reconcile existing results
        srr_col = "run" if "run" in filtered.columns else "run_accession"
        all_srr_ids = filtered[srr_col].tolist()
        self.db.init_species(species_name, all_srr_ids)
        
        quant_dir = Path(f"output/amalgkit/{species_name}/work/quant")
        reconciled = self.db.reconcile(species_name, quant_dir)
        if reconciled:
            logger.info(f"Reconciled {reconciled} already-quantified samples from filesystem")

        # Mark all filtered samples as sampled so amalgkit quant processes them
        filtered["is_sampled"] = "yes"
        
        # Write sorted metadata for batch processing
        filtered.to_csv(metadata_path, sep="\t", index=False)
        
        # ThreadPool Execution
        fastq_dir = Path(f"output/amalgkit/{species_name}/work/getfastq")
        threads_per_worker = max(1, threads // workers)
        
        srr_col = "run" if "run" in filtered.columns else "run_accession"
        
        quantified_count = 0
        
        # Per-sample timeout: 2h for download + quant combined
        SAMPLE_TIMEOUT = 7200
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            future_to_srr = {}
            for i, (_, row) in enumerate(filtered.iterrows()):
                srr = row[srr_col]
                batch_idx = i + 1
                f = executor.submit(
                    self.process_single_sample, 
                    srr, batch_idx, fastq_dir, config_path, species_name, threads_per_worker
                )
                future_to_srr[f] = srr
                
            for future in concurrent.futures.as_completed(future_to_srr):
                srr = future_to_srr[future]
                try:
                    res = future.result(timeout=SAMPLE_TIMEOUT)
                    if res["quantified"]:
                        quantified_count += 1
                        status = "Skipped (Done)" if res.get("skipped") else "Done"
                        logger.info(f"sample {res['srr']}: {status}")
                    elif res["error"]:
                        logger.warning(f"sample {res['srr']}: Failed ({res['error']})")
                except concurrent.futures.TimeoutError:
                    logger.error(f"sample {srr}: TIMEOUT after {SAMPLE_TIMEOUT}s — skipping")
                    future.cancel()
                except Exception as e:
                    logger.error(f"sample {srr}: Unexpected error — {e}")

        # Downstream Steps
        if quantified_count > 0:
            logger.info("Running downstream steps (merge, curate, sanity)...")
            workflow_log = self.log_dir / f"{species_name}_workflow.log"
            
            try:
                with open(workflow_log, "w") as f:
                    cmd = ["python3", "scripts/rna/run_workflow.py", 
                           "--config", str(config_path), "--no-progress", 
                           "--steps", "merge", "curate", "sanity"]
                    
                    result = subprocess.run(
                        cmd, stdout=f, stderr=subprocess.STDOUT,
                        timeout=1800,  # 30-minute timeout for downstream steps
                    )
                
                if result.returncode == 0:
                    logger.info("  ✓ Downstream steps complete!")
                else:
                    logger.error(f"  ⚠ Downstream steps had errors. See: {workflow_log}")

            except subprocess.TimeoutExpired:
                logger.error(
                    f"  ⚠ Downstream steps timed out after 30 minutes for {species_name}. "
                    f"See: {workflow_log}"
                )
            except Exception as e:
                logger.error(f"Failed to run downstream steps: {e}")
            
        return True

    def run_all(self, species_list: List[str], max_gb: float, workers: int, threads: int):
        """Run pipeline for all listed species configurations."""
        results = {}
        for config_name in species_list:
            try:
                success = self.process_species(config_name, max_gb, workers, threads)
                results[config_name] = "Success" if success else "Failed"
            except Exception as e:
                logger.error(f"Fatal error processing {config_name}: {e}")
                results[config_name] = f"Error: {e}"
        
        logger.info("Final Results:")
        for k, v in results.items():
            logger.info(f"{k}: {v}")

