
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

# Determine project root if possible, or use relative paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
CONFIG_DIR = PROJECT_ROOT / "config/amalgkit"
# Use local output for logs to avoid TCC/Permission issues on external volumes
LOG_DIR = PROJECT_ROOT / "output/amalgkit"
print(f"DEBUG: LOG_DIR is {LOG_DIR}")

logger = log_utils.get_logger(__name__)

class StreamingPipelineOrchestrator:
    """Orchestrator for streaming ENA download -> Quant processing."""

    def __init__(self, config_dir: Path = CONFIG_DIR, log_dir: Path = LOG_DIR):
        self.config_dir = Path(config_dir)
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        
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
        """Download FASTQ files directly from ENA using curl."""
        urls = self.query_ena_fastq_urls(srr_id)
        if not urls:
            return False

        sample_dir = out_dir / srr_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        # Check if already downloaded
        existing_fq = list(sample_dir.glob("*.fastq.gz"))
        if len(existing_fq) >= len(urls) and len(existing_fq) > 0:
            return True

        success = True
        for url in urls:
            fname = url.split('/')[-1]
            fpath = sample_dir / fname
            if fpath.exists() and fpath.stat().st_size > 0:
                continue
            
            try:
                # Use curl for robust downloading
                cmd = ["curl", "-L", "-f", "-o", str(fpath), "--retry", "3",
                       "--retry-delay", "5", "-s", "--show-error", url]
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
                
                if result.returncode != 0:
                    logger.error(f"Download failed for {fname}: {result.stderr}")
                    fpath.unlink(missing_ok=True)
                    success = False
                else:
                    sz_gb = fpath.stat().st_size / (1024**3)
                    logger.info(f"Downloaded {fname}: {sz_gb:.2f} GB")
            except Exception as e:
                logger.error(f"Download exception for {fname}: {e}")
                fpath.unlink(missing_ok=True)
                success = False

        return success

    def quant_sample(self, config_path: Path, batch_index: int, species_name: str, threads: int) -> bool:
        """Run amalgkit quant for a single sample using --batch."""
        with open(config_path) as f:
            cfg = yaml.safe_load(f)

        steps_cfg = cfg.get("steps", {})
        quant_cfg = steps_cfg.get("quant", {})
        quant_out = quant_cfg.get("out_dir", f"blue/amalgkit/{species_name}/work")
        
        # Determine metadata path
        work_dir = Path(f"blue/amalgkit/{species_name}/work")
        meta_path = None
        for mp in [work_dir / "metadata/metadata.tsv", work_dir / "metadata/metadata_selected.tsv"]:
            if mp.exists():
                meta_path = str(mp)
                break
        
        if not meta_path:
            logger.error(f"No metadata found for {species_name} batch {batch_index}")
            return False

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
        
        # Amalgkit usually handles fasta_dir via index building, but can pass if needed
        # fasta_dir = genome_cfg.get("fasta_dir", "")
        
        log_path = self.log_dir / f"{species_name}_quant.log"
        
        try:
            with open(log_path, 'a') as log_f:
                log_f.write(f"\n[{time.strftime('%Y-%m-%d %H:%M:%S')}] Quant batch {batch_index} START: {' '.join(cmd)}\n")
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
            
            with open(log_path, 'a') as log_f:
                log_f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Quant batch {batch_index} END (Exit {result.returncode})\n")
                if result.stdout: log_f.write(f"-- STDOUT --\n{result.stdout}\n")
                if result.stderr: log_f.write(f"-- STDERR --\n{result.stderr}\n")
            
            return result.returncode == 0
        except Exception as e:
            logger.error(f"Quant exception batch {batch_index}: {e}")
            with open(log_path, 'a') as log_f:
                log_f.write(f"EXCEPTION: {e}\n")
            return False

    def is_quantified(self, species_name: str, srr_id: str) -> bool:
        """Check if sample is already quantified."""
        quant_dir = Path(f"blue/amalgkit/{species_name}/work/quant/{srr_id}")
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
        if not self.download_fastq(srr_id, fastq_dir):
            result["error"] = "ENA Download Failed"
            return result
        result["downloaded"] = True
        
        # Quantify
        success = self.quant_sample(config_path, batch_idx, species_name, threads)
        result["quantified"] = success
        if not success:
            result["error"] = "Quantification Failed"
        
        return result

    def verify_genome_index(self, config_path: Path, species_name: str) -> bool:
        """Verify Kallisto index exists."""
        try:
            with open(config_path) as f:
                cfg = yaml.safe_load(f)
        except Exception as e:
            logger.error(f"Failed to load config {config_path}: {e}")
            return False
            
        index_dir = cfg.get('genome', {}).get('index_dir', '')
        search_dirs = [
            index_dir,
            f"blue/amalgkit/{species_name}/work/index",
            f"blue/amalgkit/shared/genome/{species_name}/index",
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
        
        # Verify Index
        if not self.verify_genome_index(config_path, species_name):
            return False

        # Metadata Handling
        work_dir = Path(f"blue/amalgkit/{species_name}/work")
        metadata_path = work_dir / "metadata/metadata.tsv"
        
        if not metadata_path.exists():
            # Try selected metadata
            metadata_path = work_dir / "metadata/metadata_selected.tsv"
        
        if not metadata_path.exists():
            logger.error(f"No metadata found for {species_name}")
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

        # Write sorted metadata for batch processing
        filtered.to_csv(metadata_path, sep="\t", index=False)
        
        # ThreadPool Execution
        fastq_dir = Path(f"blue/amalgkit/{species_name}/fastq/getfastq")
        threads_per_worker = max(1, threads // workers)
        
        srr_col = "run" if "run" in filtered.columns else "run_accession"
        
        quantified_count = 0
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            futures = []
            for i, (_, row) in enumerate(filtered.iterrows()):
                srr = row[srr_col]
                batch_idx = i + 1
                futures.append(executor.submit(
                    self.process_single_sample, 
                    srr, batch_idx, fastq_dir, config_path, species_name, threads_per_worker
                ))
                
            for future in concurrent.futures.as_completed(futures):
                res = future.result()
                if res["quantified"]:
                    quantified_count += 1
                    status = "Skipped (Done)" if res.get("skipped") else "Done"
                    logger.info(f"sample {res['srr']}: {status}")
                elif res["error"]:
                    logger.warning(f"sample {res['srr']}: Failed ({res['error']})")

        # Downstream Steps
        if quantified_count > 0:
            logger.info("Running downstream steps (merge, curate, sanity)...")
            workflow_log = self.log_dir / f"{species_name}_workflow.log"
            
            try:
                with open(workflow_log, "w") as f:
                    cmd = ["python3", "scripts/rna/run_workflow.py", 
                           "--config", str(config_path), "--no-progress", 
                           "--steps", "merge", "curate", "sanity"]
                    
                    result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
                
                if result.returncode == 0:
                    logger.info("  ✓ Downstream steps complete!")
                else:
                    logger.error(f"  ⚠ Downstream steps had errors. See: {workflow_log}")

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

