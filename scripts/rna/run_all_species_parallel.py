#!/usr/bin/env python3
"""Cross-species parallel pipeline orchestrator.

**ARCHITECTURE NOTE:** 
Historically, `amalgkit` processed samples chunk-by-chunk and species-by-species. If 
a species had 24 samples running on 16 threads, the final 8 samples would run alone, 
leaving 8 CPU threads idle, blocking the entire pipeline from moving to the next species.

This orchestrator shatters that boundary. It collects samples from ALL species into a 
single, massive `work_queue`. It submits the entire queue to a generic `ThreadPoolExecutor`.
As soon as any sample finishes (getfastq → quant → cleanup), that thread slot immediately 
picks up the next pending sample — even if it belongs to a completely different species.
This guarantees 100% CPU utilization across the entire RNA-seq dataset.

When the final sample for a *specific species* completes, the orchestrator detects the 
species is fully quantified and asynchronously triggers its post-processing (merge, cstmm).

Usage:
    .venv/bin/python scripts/rna/run_all_species_parallel.py --max-workers 16
    .venv/bin/python scripts/rna/run_all_species_parallel.py --species acromyrmex_echinatior amellifera
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import os
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

REPO_ROOT = Path(__file__).parent.parent.parent.resolve()
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(Path(__file__).parent))

from _setup_utils import check_environment_or_exit, ensure_venv_activated
from metainformant.core.utils.optional_deps import enable_optional_warnings, suppress_optional_warnings

suppress_optional_warnings()
ensure_venv_activated(auto_setup=True)
check_environment_or_exit(auto_setup=True)
enable_optional_warnings()

from metainformant.core.utils.logging import get_logger
from metainformant.rna.amalgkit.amalgkit import run_amalgkit
from metainformant.rna.engine.workflow import (
    apply_step_defaults,
    load_workflow_config,
    plan_workflow,
    sanitize_params_for_cli,
)
from metainformant.rna.engine.workflow_cleanup import check_disk_space, check_disk_space_or_fail, cleanup_fastqs, get_quantified_samples
from metainformant.rna.engine.workflow_core import AmalgkitWorkflowConfig, WorkflowStepResult
from metainformant.rna.engine.orchestration import run_workflow_for_species

logger = get_logger("xspecies_parallel")


class SpeciesContext:
    """Per-species config, step functions (as amalgkit CLI calls), and completion tracking."""

    def __init__(self, config_path: Path, max_workers: int):
        self.config_path = config_path
        self.species_name = config_path.stem.replace("amalgkit_", "")
        self.max_workers = max_workers

        # Load config
        self.config = load_workflow_config(config_path)
        apply_step_defaults(self.config)

        # Plan steps and separate chunk vs post-process
        planned = plan_workflow(self.config)
        self.chunk_step_params: Dict[str, Dict] = {}  # step_name -> params
        self.post_steps: List[Tuple[str, Dict]] = []

        for name, params in planned:
            if name == "integrate":
                continue
            elif name in ("getfastq", "quant"):
                self.chunk_step_params[name] = params
            elif name in ("metadata", "config", "select"):
                continue  # Already done in setup
            else:
                self.post_steps.append((name, params))

        # Load metadata
        selected_meta = self.config.work_dir / "metadata" / "metadata_selected.tsv"
        fallback_meta = self.config.work_dir / "metadata" / "metadata.tsv"
        meta_path = selected_meta if selected_meta.exists() else fallback_meta

        self.fieldnames = None
        self.samples = []
        if meta_path.exists():
            with open(meta_path, "r") as f:
                reader = csv.DictReader(f, delimiter="\t")
                self.fieldnames = reader.fieldnames
                self.samples = list(reader)

        # Symlink getfastq -> fastq/getfastq
        fastq_src_dir = self.config.work_dir.parent / "fastq" / "getfastq"
        fastq_link_dir = self.config.work_dir / "getfastq"
        if fastq_src_dir.exists() and not fastq_link_dir.exists():
            try:
                fastq_link_dir.symlink_to(fastq_src_dir.resolve(), target_is_directory=True)
            except Exception:
                pass

        # Thread-safe tracking
        self._lock = threading.Lock()
        
        # Pre-calculate already quantified samples
        self.quantified_samples = get_quantified_samples(self.config)

        self.total_samples = len(self.samples)
        self.completed_samples = 0
        self.failed_samples = 0

    def mark_done(self, sample_id: str, success: bool) -> bool:
        """Mark a sample as done. Returns True if this was the last sample for the species."""
        with self._lock:
            if success:
                self.completed_samples += 1
            else:
                self.failed_samples += 1
            total_done = self.completed_samples + self.failed_samples
            logger.info(
                f"  [{self.species_name}] {sample_id} {'✓' if success else '✗'} "
                f"({total_done}/{self.total_samples})"
            )
            return total_done >= self.total_samples


def process_sample(ctx: SpeciesContext, sample_row: Dict[str, str], idx: int) -> bool:
    """Process one sample: getfastq → quant → cleanup."""
    sample_id = sample_row.get("run", str(idx))

    # Disk space guardrail: skip sample if < 20 GB free
    ok, free_gb = check_disk_space(ctx.config.work_dir, min_free_gb=20.0)
    if not ok:
        logger.warning(
            f"  [{ctx.species_name}/{sample_id}] Skipping — only {free_gb:.1f} GB free (need 20 GB)"
        )
        return False

    # Write single-sample metadata
    single_meta = ctx.config.work_dir / "metadata" / f"metadata_xsp_{sample_id}.tsv"
    with open(single_meta, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=ctx.fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(sample_row)

    # Run each chunk step (getfastq, quant)
    for step_name in ["getfastq", "quant"]:
        if step_name not in ctx.chunk_step_params:
            continue

        params = ctx.chunk_step_params[step_name].copy()
        params["metadata"] = str(single_meta)

        # Throttle threads
        base_threads = int(params.get("threads", 1))
        params["threads"] = max(1, base_threads // ctx.max_workers)
        if "jobs" in params:
            params["jobs"] = params["threads"]

        try:
            if step_name == "getfastq":
                check_disk_space_or_fail(ctx.config.work_dir, min_free_gb=20.0, step_name=f"getfastq/{sample_id}")
                
                # Size guardrail: skip oversized samples BEFORE any download attempt.
                total_bases_str = str(sample_row.get("total_bases", "")).strip()
                try:
                    if not total_bases_str or total_bases_str.lower() in ["na", "none", "nan"]:
                        total_bases = 0
                    else:
                        total_bases = int(float(total_bases_str))
                except (ValueError, TypeError):
                    logger.warning(f"  [{ctx.species_name}/{sample_id}] Invalid total_bases value: '{total_bases_str}', defaulting to 0 for limit check.")
                    total_bases = 0
                    
                if total_bases > 4_000_000_000:
                    logger.warning(
                        f"  [{ctx.species_name}/{sample_id}] Sample too large "
                        f"({total_bases:,} bases > 4B limit). AUTO-SKIPPING to prioritize smaller samples first."
                    )
                    return False

                # Hybrid Architecture: Try custom ENA exact downloader first
                ena_script = REPO_ROOT / "scripts" / "rna" / "download_ena.py"
                getfastq_out = ctx.config.work_dir / "fastq"
                if "out_dir" in params:
                    getfastq_out = Path(params["out_dir"])
                
                ena_cmd = [sys.executable, str(ena_script), sample_id, str(getfastq_out)]
                ena_result = subprocess.run(ena_cmd, capture_output=True, text=True)
                
                if ena_result.returncode == 0:
                    logger.info(f"  [{ctx.species_name}/{sample_id}] ENA Direct Download successful!")
                    
                    # Stage flat ENA files into the subdirectory structure amalgkit quant expects.
                    # amalgkit quant looks for FASTQs at {quant_out_dir}/getfastq/{SRR_ID}/ (quant.py line 260)
                    # but ENA downloads put them flat in {getfastq_out}/SRR_1.fastq.gz
                    quant_params = ctx.chunk_step_params.get("quant", {})
                    quant_out_dir = Path(quant_params.get("out_dir", str(ctx.config.work_dir)))
                    staging_dir = quant_out_dir / "getfastq" / sample_id
                    staging_dir.mkdir(parents=True, exist_ok=True)
                    
                    for pattern in [f"{sample_id}_*.fastq.gz", f"{sample_id}.fastq.gz"]:
                        for fq_file in getfastq_out.glob(pattern):
                            link_target = staging_dir / fq_file.name
                            if not link_target.exists():
                                try:
                                    link_target.symlink_to(fq_file.resolve())
                                    logger.debug(f"  Staged {fq_file.name} → {staging_dir}")
                                except OSError:
                                    pass
                    
                    continue  # Skip amalgkit getfastq — ENA already got the files
                else:
                    logger.info(f"  [{ctx.species_name}/{sample_id}] ENA failed or not found. Falling back to AWS/NCBI fasterq-dump...")

            safe_params = sanitize_params_for_cli(step_name, params)
            result = run_amalgkit(step_name, safe_params)

            if result.returncode != 0:
                logger.error(f"  [{ctx.species_name}/{sample_id}] {step_name} failed (code {result.returncode})")
                return False
                
            # If amalgkit getfastq fallback succeeded, stage its flat output files for quant
            if step_name == "getfastq":
                quant_params = ctx.chunk_step_params.get("quant", {})
                quant_out_dir = Path(quant_params.get("out_dir", str(ctx.config.work_dir)))
                staging_dir = quant_out_dir / "getfastq" / sample_id
                staging_dir.mkdir(parents=True, exist_ok=True)
                
                fallback_out = ctx.config.work_dir / "fastq"
                if "out_dir" in params:
                    fallback_out = Path(params["out_dir"])
                    
                for pattern in [f"{sample_id}_*.fastq.gz", f"{sample_id}.fastq.gz"]:
                    for fq_file in fallback_out.glob(pattern):
                        link_target = staging_dir / fq_file.name
                        if not link_target.exists():
                            try:
                                link_target.symlink_to(fq_file.resolve())
                                logger.debug(f"  Staged fallback {fq_file.name} → {staging_dir}")
                            except OSError:
                                pass
            
            # Additional validation loop: verify downloaded/staged fastq size
            if step_name == "getfastq":
                quant_params = ctx.chunk_step_params.get("quant", {})
                quant_out_dir = Path(quant_params.get("out_dir", str(ctx.config.work_dir)))
                staging_dir = quant_out_dir / "getfastq" / sample_id
                
                fastq_files = list(staging_dir.glob("*.fastq*"))
                if not fastq_files:
                    logger.error(f"  [{ctx.species_name}/{sample_id}] getfastq succeeded but no FASTQ files were found in staging!")
                    return False
                
                for fq in fastq_files:
                    if fq.stat().st_size < 100:  # Less than 100 bytes is effectively empty/corrupt
                        logger.error(f"  [{ctx.species_name}/{sample_id}] Corrupt/Empty FASTQ downloaded or staged: {fq.name} ({fq.stat().st_size} bytes)")
                        return False


        except Exception as e:
            logger.error(f"  [{ctx.species_name}/{sample_id}] {step_name} exception: {e}")
            return False

    # Cleanup FASTQs on success
    logger.info(f"  [{ctx.species_name}/{sample_id}] ✓ Cleaning FASTQs")
    cleanup_fastqs(ctx.config, [sample_id])
    return True


def run_post_processing(ctx: SpeciesContext):
    """Run merge/curate/sanity after all samples for a species are done."""
    if not ctx.post_steps:
        return
    logger.info(f"━━━ [{ctx.species_name}] All samples done → post-processing ━━━")
    for step_name, params in ctx.post_steps:
        try:
            result = run_amalgkit(step_name, params)
            status = "✓" if result.returncode == 0 else f"✗ (code {result.returncode})"
            logger.info(f"  [{ctx.species_name}] {step_name} {status}")
        except Exception as e:
            logger.error(f"  [{ctx.species_name}] {step_name} exception: {e}")


def discover_configs(species_filter: Optional[List[str]] = None) -> List[Path]:
    config_dir = REPO_ROOT / "config" / "amalgkit"
    configs = []
    for p in sorted(config_dir.glob("amalgkit_*.yaml")):
        stem = p.stem.lower()
        if "template" in stem or "test" in stem or "cross_species" in stem:
            continue
        if species_filter:
            species_name = p.stem.replace("amalgkit_", "")
            if species_name not in species_filter:
                continue
        configs.append(p)
    return configs


def main():
    parser = argparse.ArgumentParser(description="Cross-species parallel pipeline")
    parser.add_argument("--max-workers", type=int, default=12, help="Max concurrent samples")
    parser.add_argument("--species", nargs="+", help="Only these species (default: all)")
    parser.add_argument("--dry-run", action="store_true", help="Show plan only")
    args = parser.parse_args()

    configs = discover_configs(args.species)
    if not configs:
        logger.error("No species configs found!")
        return 1

    logger.info(f"{'='*60}")
    logger.info(f" Cross-Species Parallel Pipeline")
    logger.info(f" Species: {len(configs)}, Max workers: {args.max_workers}")
    logger.info(f"{'='*60}")

    # Phase 1: Build species contexts from configs that already have metadata
    # (Skip running metadata/select setup — that's slow NCBI queries. Only process ready species.)
    species_contexts: List[SpeciesContext] = []
    for cfg_path in configs:
        species_name = cfg_path.stem.replace("amalgkit_", "")
        try:
            ctx = SpeciesContext(cfg_path, args.max_workers)
            if ctx.total_samples == 0:
                logger.info(f"  [{species_name}] No metadata/samples yet — skipping")
                continue

            species_contexts.append(ctx)
            logger.info(f"  [{species_name}] {ctx.total_samples} samples ready")
        except Exception as e:
            logger.error(f"  [{species_name}] Error loading config: {e}")
            continue

    # Build unified work queue
    work_queue = []
    completed_contexts = []
    for ctx in species_contexts:
        pending_samples = []
        for idx, sample in enumerate(ctx.samples):
            sample_id = sample.get("run", str(idx))
            if sample_id in ctx.quantified_samples:
                ctx.completed_samples += 1
            else:
                pending_samples.append((ctx, sample, idx))
                
        if not pending_samples and ctx.total_samples > 0:
            logger.info(f"  [{ctx.species_name}] All {ctx.total_samples} samples already quantified! Queueing for immediate post-processing.")
            completed_contexts.append(ctx)
        else:
            work_queue.extend(pending_samples)

    logger.info(f"\n Total: {len(work_queue)} pending samples across {len(species_contexts)} species")

    if args.dry_run:
        for ctx in species_contexts:
            logger.info(f"  {ctx.species_name}: {ctx.total_samples} samples")
        return 0

    # Phase 2: Continuous pool — all samples, all species, one pool
    logger.info(f"\n Starting continuous pool ({args.max_workers} workers)...")
    
    completed_species = {ctx.species_name for ctx in completed_contexts}
    post_lock = threading.Lock()

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        # Immediately submit post-processing for species that were already completely 100% quantified
        for ctx in completed_contexts:
            executor.submit(run_post_processing, ctx)
            
        future_to_info = {}
        for ctx, sample, idx in work_queue:
            fut = executor.submit(process_sample, ctx, sample, idx)
            future_to_info[fut] = (ctx, sample.get("run", str(idx)))

        for fut in concurrent.futures.as_completed(future_to_info):
            ctx, sample_id = future_to_info[fut]
            try:
                success = fut.result()
                is_last = ctx.mark_done(sample_id, success)

                if is_last:
                    with post_lock:
                        if ctx.species_name not in completed_species:
                            completed_species.add(ctx.species_name)
                            executor.submit(run_post_processing, ctx)
            except Exception as exc:
                logger.error(f"  [{ctx.species_name}/{sample_id}] Exception: {exc}")
                ctx.mark_done(sample_id, False)

    # Summary
    logger.info(f"\n{'='*60}")
    logger.info(f" FINAL SUMMARY")
    logger.info(f"{'='*60}")
    for ctx in species_contexts:
        s = "✓" if ctx.failed_samples == 0 else "✗"
        logger.info(f"  {s} {ctx.species_name}: {ctx.completed_samples}/{ctx.total_samples} OK, {ctx.failed_samples} fail")

    return 1 if any(c.failed_samples > 0 for c in species_contexts) else 0


if __name__ == "__main__":
    sys.exit(main())
