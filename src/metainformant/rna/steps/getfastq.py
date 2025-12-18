"""Step runner for `amalgkit getfastq` (download raw FASTQ files)."""

from __future__ import annotations

import os
import shutil
import subprocess
import time
import urllib.request
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ...core.errors import ValidationError
from ...core.io import read_delimited
from ...core.logging import get_logger
from ..amalgkit import getfastq as _getfastq
from .download_progress import DownloadProgressMonitor


def _inject_robust_defaults(raw_params: Mapping[str, Any] | None) -> dict[str, Any]:
    """Return a copy of params with robust download defaults injected when absent.

    Heuristics:
    - Prefer parallel-fastq-dump when available; otherwise disable it to let
      amalgkit use sra-tools directly.
    - Provide explicit executables for PFD and prefetch when discoverable.
    - Use NCBI_EMAIL from environment when not specified.
    - Enable printing of helper tool outputs for better diagnostics.
    """
    params: dict[str, Any] = {}
    if raw_params:
        params.update({str(k): v for k, v in raw_params.items()})

    # Entrez email from environment if missing
    if not params.get("entrez_email"):
        email = os.environ.get("NCBI_EMAIL")
        if email:
            params["entrez_email"] = email

    # Tool discovery
    pfd_path = shutil.which("parallel-fastq-dump")
    prefetch_path = shutil.which("prefetch")

    # Normalize pfd parameter (handle "yes"/"no" strings and booleans)
    # Determine if PFD should be enabled based on user preference and tool availability
    pfd_value = params.get("pfd")
    if "pfd" not in params:
        # Not specified - enable if tool is available
        pfd_enabled = pfd_path is not None
        params["pfd"] = pfd_enabled
    else:
        # User specified a value - normalize it
        if isinstance(pfd_value, str):
            pfd_enabled = pfd_value.lower() in ("yes", "true", "1")
            params["pfd"] = "yes" if pfd_enabled else "no"
        else:
            pfd_enabled = bool(pfd_value)
            params["pfd"] = pfd_enabled

    # Only set pfd_exe if PFD is actually enabled AND tool is available
    # If PFD is disabled, explicitly remove pfd_exe to prevent amalgkit from checking
    if pfd_enabled and pfd_path:
        if not params.get("pfd_exe"):
            params["pfd_exe"] = pfd_path
    else:
        # PFD is disabled or not available - explicitly remove pfd_exe to prevent amalgkit from checking
        params.pop("pfd_exe", None)
    
    if prefetch_path and not params.get("prefetch_exe"):
        params["prefetch_exe"] = prefetch_path

    # Improve diagnostics by default unless user opted out
    params.setdefault("pfd_print", True)
    params.setdefault("fastp_print", True)
    # Avoid fastp bottlenecks by default for initial fetch; downstream QC is optional
    params.setdefault("fastp", False)

    # Enable optional accelerated sources (ENA/AWS) unless caller opts out
    params.setdefault("accelerate", True)
    return params


def run(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
):
    """Run `amalgkit getfastq` with robust defaults, verification, and targeted retries.

    Strategy:
    - Bulk attempt (metadata or id) via `amalgkit getfastq` using robust sra-tools path.
    - Verify per-SRR outputs under out_dir/getfastq/<SRR>.
    - Targeted retries for missing SRRs: first via `amalgkit getfastq --id <SRR>`,
      then fall back to direct sra-tools: prefetch + fasterq-dump (+ gzip/pigz).
    """
    effective_params = _inject_robust_defaults(params)
    out_dir = Path(str(effective_params.get("out_dir", work_dir or "."))).expanduser()
    
    # Extract accelerate flag and apply to amalgkit params
    accelerate_enabled = bool(effective_params.pop("accelerate", False))
    if accelerate_enabled:
        # Enable cloud acceleration
        effective_params.setdefault("aws", "yes")
        effective_params.setdefault("gcp", "yes")
        effective_params.setdefault("ncbi", "yes")
    
    # Remove any max_size/min_size params - these are not amalgkit parameters
    # (they're internal to prefetch, which amalgkit manages automatically)
    effective_params.pop("max_size", None)
    effective_params.pop("min_size", None)
    
    # Initialize progress monitor if enabled
    progress_monitor: DownloadProgressMonitor | None = None
    show_progress = effective_params.get("show_progress", True)
    if show_progress:
        update_interval = float(effective_params.get("progress_update_interval", 2.0))
        use_bars = effective_params.get("progress_style", "bar") == "bar"
        progress_monitor = DownloadProgressMonitor(
            out_dir=out_dir,
            update_interval=update_interval,
            use_progress_bars=use_bars,
            show_summary=not use_bars,
        )
        progress_monitor.start_monitoring()
        get_logger(__name__).info("📊 Progress tracking enabled for getfastq step")
    
    # 1) Bulk attempt (metadata or id)
    logger = get_logger(__name__)
    meta_path_used = effective_params.get("metadata", "inferred")
    if meta_path_used and meta_path_used != "inferred":
        logger.info(f"📋 Using metadata file: {meta_path_used}")
        # Verify it has run column
        try:
            rows = list(read_delimited(str(meta_path_used), delimiter="\t"))
            if rows:
                sample_count = len(rows)
                has_run = 'run' in rows[0]
                logger.info(f"   ✓ Metadata file has {sample_count} samples, 'run' column: {has_run}")
        except Exception as e:
            logger.warning(f"   ⚠ Could not verify metadata file: {e}")
    
    logger.info("🚀 Starting bulk download via amalgkit getfastq...")
    bulk_result = _getfastq(
        effective_params,
        work_dir=work_dir,
        log_dir=log_dir,
        step_name="getfastq",
        check=False,
    )
    logger.info(f"   Bulk download completed with return code: {bulk_result.returncode}")

    # 2) Determine SRR list to verify
    srr_list: list[str] = []
    id_val = effective_params.get("id")
    if id_val:
        if isinstance(id_val, (list, tuple)):
            srr_list = [str(srr).strip() for srr in id_val]
        else:
            srr_list = [str(id_val).strip()]
    else:
        meta_path = effective_params.get("metadata")
        if not meta_path or meta_path == "inferred":
            # Try to find metadata file - use species-specific work_dir from params if available
            # First priority: use work_dir parameter if provided
            if work_dir:
                actual_work_dir = Path(work_dir)
            else:
                # Infer from out_dir structure: output/amalgkit/<species>/fastq -> output/amalgkit/<species>/work
                actual_work_dir = Path(out_dir).parent / "work"
            
            # Try multiple metadata file locations in order of preference
            # Skip pivot tables - they lack Run IDs and are not suitable for getfastq
            candidate_paths = [
                actual_work_dir / "metadata" / "metadata.filtered.tissue.tsv",
                actual_work_dir / "metadata" / "metadata.tsv",
            ]
            
            # Find first existing metadata file
            meta_path = None
            for candidate in candidate_paths:
                if candidate.exists():
                    meta_path = str(candidate)
                    break
            
            # If still not found, use the default metadata.tsv path (will fail with clear error)
            if not meta_path:
                meta_path = str(actual_work_dir / "metadata" / "metadata.tsv")
            
            # Validate selected metadata file has 'run' column
            if meta_path:
                logger = get_logger(__name__)
                logger.info(f"📋 Validating metadata file: {meta_path}")
                try:
                    rows = list(read_delimited(str(meta_path), delimiter="\t"))
                    if rows and 'run' not in rows[0]:
                        # This is likely a pivot table - skip it and try fallback
                        logger.warning(
                            f"   ⚠ Metadata file lacks 'run' column (likely pivot table). "
                            f"Skipping and trying fallback."
                        )
                        # Try fallback to metadata.tsv
                        fallback = actual_work_dir / "metadata" / "metadata.tsv"
                        if fallback.exists() and fallback != Path(meta_path):
                            meta_path = str(fallback)
                            logger.info(f"   ✓ Using fallback metadata file: {fallback}")
                        else:
                            raise ValidationError(
                                f"Metadata file {meta_path} lacks 'run' column and no fallback available"
                            )
                    elif rows:
                        logger.info(f"   ✓ Metadata file validated: {len(rows)} samples, has 'run' column")
                except ValueError:
                    # Re-raise validation errors
                    raise
                except Exception as e:
                    logger.warning(f"   ⚠ Could not validate metadata file {meta_path}: {e}")
            
            # Debug: print what we're looking for
            get_logger(__name__).debug(f"Looking for metadata at: {meta_path}")
            get_logger(__name__).debug(f"Work dir: {actual_work_dir}")
            get_logger(__name__).debug(f"Out dir: {out_dir}")
        try:
            rows = list(read_delimited(str(meta_path), delimiter="\t"))
            if rows:
                header = rows[0].keys()
                run_key = "run" if "run" in header else next((k for k in header if k.lower() == "run"), None)
                if run_key:
                    for row in rows:
                        val = (row.get(run_key) or "").strip()
                        if val:
                            srr_list.append(val)
        except Exception as e:
            # If metadata cannot be read, proceed without verification
            get_logger(__name__).warning(f"Could not read metadata file {meta_path}: {e}")
            # Continue with empty metadata

    srr_list = sorted(set(srr_list))

    def _has_fastq(srr: str) -> bool:
        d = out_dir / "getfastq" / srr
        return any(
            (d / f).exists()
            for f in (
                f"{srr}_1.fastq.gz",
                f"{srr}_2.fastq.gz",
                f"{srr}_1.fastq",
                f"{srr}_2.fastq",
                f"{srr}.fastq.gz",
                f"{srr}.fastq",
            )
        )

    def _has_sra(srr: str) -> bool:
        """Check if SRA file exists for a sample."""
        srr_dir = out_dir / "getfastq" / srr
        return (srr_dir / f"{srr}.sra").exists()

    # 2.5) Check for samples with SRA but no FASTQ (conversion failure)
    logger = get_logger(__name__)
    if srr_list:
        logger.info(f"🔍 Verifying {len(srr_list)} samples for FASTQ files...")
    
    conversion_needed: list[str] = []
    if srr_list:
        for srr in srr_list:
            if _has_sra(srr) and not _has_fastq(srr):
                conversion_needed.append(srr)
                logger.warning(
                    f"   ⚠ Sample {srr} has SRA file but no FASTQ files. "
                    f"This indicates SRA-to-FASTQ conversion failed. Attempting automatic conversion."
                )
    
    # Attempt automatic conversion for samples with SRA but no FASTQ
    if conversion_needed:
        logger.info(f"🔄 Attempting automatic SRA-to-FASTQ conversion for {len(conversion_needed)} samples")
        threads = effective_params.get("threads", 6)
        for idx, srr in enumerate(conversion_needed, 1):
            logger.info(f"   [{idx}/{len(conversion_needed)}] Converting {srr}...")
            srr_dir = out_dir / "getfastq" / srr
            sra_file = srr_dir / f"{srr}.sra"
            if sra_file.exists():
                success, message, fastq_files = convert_sra_to_fastq(
                    srr,
                    sra_file,
                    srr_dir,
                    threads=threads,
                    log_dir=log_dir,
                )
                if success:
                    logger.info(f"   ✓ Successfully converted {srr}: {len(fastq_files)} FASTQ files")
                else:
                    logger.error(f"   ✗ Failed to convert {srr}: {message}")
            else:
                logger.warning(f"   ⚠ SRA file not found for {srr} (expected at {sra_file})")

    # 3) Targeted retries for missing SRRs
    missing = [s for s in srr_list if not _has_fastq(s)] if srr_list else []
    
    if missing:
        logger.info(f"🔄 Retrying download for {len(missing)} missing samples: {missing[:5]}{'...' if len(missing) > 5 else ''}")
    
    # Register missing samples with progress monitor for tracking
    if progress_monitor and missing:
        for idx, srr in enumerate(missing, 1):
            progress_monitor.register_thread(idx, srr)
    
    for idx, srr in enumerate(missing, 1):
        # Attempt 0: Accelerated sources (ENA HTTP FASTQ or AWS ODP S3 .sra)
        if accelerate_enabled:
            try:
                srr_dir = Path(out_dir) / "getfastq" / srr
                srr_dir.mkdir(parents=True, exist_ok=True)
                # AWS ODP direct .sra
                aws_url = f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{srr}/{srr}"
                try:
                    with urllib.request.urlopen(aws_url) as resp:
                        if resp.status == 200:
                            # stream to file
                            sra_path = srr_dir / f"{srr}.sra"
                            with open(sra_path, "wb") as out_f:
                                shutil.copyfileobj(resp, out_f)
                except Exception as e:
                    # Log error but continue - file copy issues shouldn't block entire download
                    get_logger(__name__).debug(f"Could not copy file for {srr}: {e}")
                # If we now have a local .sra, try fasterq-dump quickly
                if (srr_dir / f"{srr}.sra").exists():
                    fasterq_bin = shutil.which("fasterq-dump")
                    if fasterq_bin:
                        subprocess.run(
                            [
                                fasterq_bin,
                                "--threads",
                                str(effective_params.get("threads", 6)),
                                "--split-files",
                                "--size-check", "off",  # Disable disk size check to prevent "disk-limit exceeded" errors
                                "-O",
                                str(srr_dir),
                                str(srr_dir / f"{srr}.sra"),
                            ],
                            check=False,
                        )
                        # compress
                        pigz = shutil.which("pigz")
                        if (srr_dir / f"{srr}_1.fastq").exists():
                            subprocess.run([pigz or "gzip", "-f", str(srr_dir / f"{srr}_1.fastq")], check=False)
                        if (srr_dir / f"{srr}_2.fastq").exists():
                            subprocess.run([pigz or "gzip", "-f", str(srr_dir / f"{srr}_2.fastq")], check=False)
                        if _has_fastq(srr):
                            if progress_monitor:
                                progress_monitor.unregister_thread(idx, success=True)
                            continue
            except Exception:
                pass
        # Attempt A: amalgkit getfastq --id <SRR> using robust flags
        per_params = dict(effective_params)
        per_params["id"] = srr
        per_params["redo"] = True
        ak_res = _getfastq(
            per_params,
            work_dir=work_dir,
            log_dir=log_dir,
            step_name="getfastq",
            check=False,
        )
        if _has_fastq(srr):
            if progress_monitor:
                progress_monitor.unregister_thread(idx, success=True)
            continue
        # Attempt B: direct sra-tools fallback: prefetch + fasterq-dump (+ gzip)
        prefetch_bin = shutil.which("prefetch")
        fasterq_bin = shutil.which("fasterq-dump")
        if prefetch_bin and fasterq_bin:
            srr_dir = out_dir / "getfastq" / srr
            srr_dir.mkdir(parents=True, exist_ok=True)
            try:
                subprocess.run([prefetch_bin, "--output-directory", str(srr_dir), srr], check=False)
                sra_path = srr_dir / f"{srr}.sra"
                if sra_path.exists():
                    subprocess.run(
                        [
                            fasterq_bin,
                            "--threads",
                            str(effective_params.get("threads", 6)),
                            "--split-files",
                            "-O",
                            str(srr_dir),
                            str(sra_path),
                        ],
                        check=False,
                    )
                    # gzip if present
                    pigz = shutil.which("pigz")
                    if (srr_dir / f"{srr}_1.fastq").exists():
                        if pigz:
                            subprocess.run([pigz, "-f", str(srr_dir / f"{srr}_1.fastq")], check=False)
                        else:
                            subprocess.run(["gzip", "-f", str(srr_dir / f"{srr}_1.fastq")], check=False)
                    if (srr_dir / f"{srr}_2.fastq").exists():
                        if pigz:
                            subprocess.run([pigz, "-f", str(srr_dir / f"{srr}_2.fastq")], check=False)
                        else:
                            subprocess.run(["gzip", "-f", str(srr_dir / f"{srr}_2.fastq")], check=False)
            except Exception:
                pass
        
        # Unregister from progress monitor
        if progress_monitor:
            success = _has_fastq(srr)
            progress_monitor.unregister_thread(idx, success=success)

    # Stop progress monitoring
    if progress_monitor:
        progress_monitor.stop_monitoring()

    # 4) Final status: success if all SRRs (if known) have FASTQs
    final_missing = [s for s in srr_list if not _has_fastq(s)] if srr_list else []
    
    # Detailed diagnostics for failures
    if final_missing:
        logger = get_logger(__name__)
        logger.error(f"Failed to obtain FASTQ files for {len(final_missing)} samples: {final_missing[:10]}{'...' if len(final_missing) > 10 else ''}")
        
        # Distinguish between download failures and conversion failures
        download_failed: list[str] = []
        conversion_failed: list[str] = []
        for srr in final_missing:
            if _has_sra(srr):
                conversion_failed.append(srr)
                logger.error(
                    f"  {srr}: SRA file exists but FASTQ conversion failed. "
                    f"This indicates a conversion problem, not a download problem."
                )
            else:
                download_failed.append(srr)
                logger.error(
                    f"  {srr}: No SRA file found. This indicates a download problem."
                )
        
        if conversion_failed:
            logger.error(
                f"Conversion failures ({len(conversion_failed)} samples): "
                f"These samples were downloaded but failed to convert to FASTQ format. "
                f"Check SRA file integrity and sra-tools availability."
            )
        if download_failed:
            logger.error(
                f"Download failures ({len(download_failed)} samples): "
                f"These samples could not be downloaded from SRA. "
                f"Check network connectivity and SRA accession validity."
            )
    
    # Determine return code: prioritize missing samples over bulk result
    # If we have missing samples, return error code 1
    # Otherwise, use bulk_result return code
    rc = 0 if not final_missing and bulk_result.returncode == 0 else (1 if final_missing else bulk_result.returncode)
    
    # Final summary
    if srr_list:
        successful = len(srr_list) - len(final_missing)
        logger.info("=" * 80)
        logger.info(f"📊 getfastq Step Summary:")
        logger.info(f"   Total samples: {len(srr_list)}")
        logger.info(f"   ✓ Successful: {successful}")
        if final_missing:
            logger.info(f"   ✗ Failed: {len(final_missing)}")
        if conversion_needed:
            logger.info(f"   🔄 Auto-converted: {len([s for s in conversion_needed if _has_fastq(s)])}")
        logger.info(f"   Return code: {rc}")
        logger.info("=" * 80)

    # 5) Optional cleanup: remove raw FASTQ files if cleanup_raw is enabled
    cleanup_raw = effective_params.get("cleanup_raw", False)
    if cleanup_raw and rc == 0 and not final_missing:
        get_logger(__name__).info("Cleaning up raw FASTQ files after successful download")
        try:
            for srr in srr_list:
                srr_dir = out_dir / "getfastq" / srr
                if srr_dir.exists():
                    # Remove raw FASTQ files, keep gzipped versions if they exist
                    for fastq_file in srr_dir.glob("*.fastq"):
                        if not fastq_file.name.endswith('.gz'):
                            fastq_file.unlink()
                            get_logger(__name__).debug(f"Removed raw FASTQ file: {fastq_file}")
        except Exception as e:
            get_logger(__name__).warning(f"Failed to cleanup raw FASTQ files: {e}")

    result = subprocess.CompletedProcess(["amalgkit", "getfastq"], rc, stdout="", stderr="")
    if check and rc != 0:
        raise subprocess.CalledProcessError(rc, result.args)
    return result


def convert_sra_to_fastq(
    sample_id: str,
    sra_file: Path,
    output_dir: Path,
    *,
    threads: int = 4,
    log_dir: Path | None = None,
) -> tuple[bool, str, list[Path]]:
    """Convert a local SRA file to FASTQ format.
    
    Prefers parallel-fastq-dump (works better with local files) and falls back
    to fasterq-dump if needed. Automatically compresses output FASTQ files.
    
    Args:
        sample_id: SRA accession ID (e.g., "SRR1234567")
        sra_file: Path to the SRA file
        output_dir: Directory where FASTQ files should be written
        threads: Number of threads for conversion
        log_dir: Optional directory for log files
        
    Returns:
        Tuple of (success: bool, message: str, fastq_files: list[Path])
        fastq_files contains paths to created FASTQ files (may be empty if failed)
    """
    logger = get_logger(__name__)
    
    if not sra_file.exists():
        return False, f"SRA file not found: {sra_file}", []
    
    # Check if FASTQ files already exist
    fastq_files_existing = list(output_dir.glob(f"{sample_id}_*.fastq.gz"))
    if not fastq_files_existing:
        fastq_files_existing = list(output_dir.glob(f"{sample_id}_*.fastq"))
    if not fastq_files_existing:
        fastq_files_existing = list(output_dir.glob(f"{sample_id}.fastq.gz"))
    if not fastq_files_existing:
        fastq_files_existing = list(output_dir.glob(f"{sample_id}.fastq"))
    
    if fastq_files_existing:
        logger.info(f"FASTQ files already exist for {sample_id}")
        return True, f"FASTQ files already exist for {sample_id}", fastq_files_existing
    
    sra_size_gb = sra_file.stat().st_size / 1e9
    logger.info(f"Converting SRA to FASTQ for {sample_id} (SRA: {sra_file.name}, {sra_size_gb:.2f} GB)...")
    
    conversion_start = time.time()
    
    # Try parallel-fastq-dump first (works better with local files)
    parallel_fastq_dump = shutil.which("parallel-fastq-dump")
    if parallel_fastq_dump:
        cmd = [
            parallel_fastq_dump,
            "-s", str(sra_file),  # SRA file path
            "-O", str(output_dir),  # Output directory
            "--threads", str(threads),
            "--gzip",  # Compress directly
        ]
        
        log_file = None
        if log_dir:
            log_file = log_dir / f"parallel_fastq_dump_{sample_id}.log"
            log_file.parent.mkdir(parents=True, exist_ok=True)
        
        try:
            # Set timeout: 2 hours per GB of SRA file (minimum 30 minutes, maximum 6 hours)
            sra_size_gb = sra_file.stat().st_size / 1e9
            timeout_seconds = max(1800, min(21600, int(sra_size_gb * 7200)))  # 30min to 6 hours
            
            with open(log_file, "w") if log_file else open(os.devnull, "w") as f:
                result = subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    cwd=str(output_dir),
                    check=False,
                    timeout=timeout_seconds,
                )
            
            # Check for output files
            fastq_files = list(output_dir.glob(f"{sample_id}_*.fastq.gz"))
            if not fastq_files:
                fastq_files = list(output_dir.glob(f"{sample_id}.fastq.gz"))
            
            if fastq_files:
                elapsed = time.time() - conversion_start
                total_size_mb = sum(f.stat().st_size for f in fastq_files) / 1e6
                logger.info(
                    f"Converted SRA to FASTQ for {sample_id} using parallel-fastq-dump "
                    f"({len(fastq_files)} files, {total_size_mb:.1f} MB, {elapsed:.1f}s)"
                )
                return True, f"Converted SRA to FASTQ for {sample_id}", fastq_files
            
            # If parallel-fastq-dump failed, fall through to fasterq-dump
            logger.warning(f"parallel-fastq-dump failed (code {result.returncode}), trying fasterq-dump")
        except subprocess.TimeoutExpired:
            elapsed = time.time() - conversion_start
            logger.warning(
                f"parallel-fastq-dump timed out for {sample_id} after {elapsed/60:.1f} minutes, trying fasterq-dump"
            )
        except Exception as e:
            logger.warning(f"parallel-fastq-dump error: {e}, trying fasterq-dump")
    
    # Fallback to fasterq-dump
    # Find the real fasterq-dump binary (not the wrapper)
    # Check common system locations first
    fasterq_dump = None
    for path in ["/usr/bin/fasterq-dump", "/usr/local/bin/fasterq-dump"]:
        if shutil.which(path) or Path(path).exists():
            fasterq_dump = path
            break
    
    # If not found in system locations, use which but filter out wrapper
    if not fasterq_dump:
        which_result = shutil.which("fasterq-dump")
        if which_result and "temp" not in which_result:
            fasterq_dump = which_result
        elif which_result:
            # Wrapper found, try to find real binary
            # The wrapper calls /usr/bin/fasterq-dump, so try that
            if Path("/usr/bin/fasterq-dump").exists():
                fasterq_dump = "/usr/bin/fasterq-dump"
    
    if not fasterq_dump:
        return False, "Neither parallel-fastq-dump nor fasterq-dump found in PATH", []
    
    # fasterq-dump needs the accession ID and looks for SRA file in current directory
    # or we can pass the file path directly (per its usage: fasterq-dump <path>)
    # We'll run from the directory containing the SRA file
    sra_dir = sra_file.parent
    
    cmd = [
        fasterq_dump,
        sample_id,  # Use accession ID
        "--outdir", str(output_dir),
        "--threads", str(threads),
        "--split-files",  # Split paired-end reads
        "--size-check", "off",  # Disable disk size check to prevent "disk-limit exceeded" errors
        "-p",  # Show progress
    ]
    
    log_file = None
    if log_dir:
        log_file = log_dir / f"fasterq_dump_{sample_id}.log"
        log_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Set timeout: 2 hours per GB of SRA file (minimum 30 minutes, maximum 6 hours)
    timeout_seconds = max(1800, min(21600, int(sra_size_gb * 7200)))  # 30min to 6 hours
    
    try:
        with open(log_file, "w") if log_file else open(os.devnull, "w") as f:
            result = subprocess.run(
                cmd,
                stdout=f,
                stderr=subprocess.STDOUT,
                cwd=str(sra_dir),  # Run from directory containing SRA file
                check=False,
                timeout=timeout_seconds,
            )
        
        # Compress FASTQ files (fasterq-dump doesn't compress automatically)
        pigz = shutil.which("pigz") or "gzip"
        for fastq_file in output_dir.glob(f"{sample_id}_*.fastq"):
            if not fastq_file.name.endswith('.gz'):
                subprocess.run([pigz, "-f", str(fastq_file)], check=False, capture_output=True)
        
        # Also check for single-end FASTQ
        for fastq_file in output_dir.glob(f"{sample_id}.fastq"):
            if not fastq_file.name.endswith('.gz'):
                subprocess.run([pigz, "-f", str(fastq_file)], check=False, capture_output=True)
        
        # Check if FASTQ files were created
        fastq_files = list(output_dir.glob("*.fastq.gz"))
        if not fastq_files:
            fastq_files = list(output_dir.glob("*.fastq"))
        
        if fastq_files:
            elapsed = time.time() - conversion_start
            total_size_mb = sum(f.stat().st_size for f in fastq_files) / 1e6
            logger.info(
                f"Converted SRA to FASTQ for {sample_id} using fasterq-dump "
                f"({len(fastq_files)} files, {total_size_mb:.1f} MB, {elapsed:.1f}s)"
            )
            return True, f"Converted SRA to FASTQ for {sample_id}", fastq_files
        else:
            elapsed = time.time() - conversion_start
            logger.warning(
                f"SRA conversion failed for {sample_id} (code {result.returncode}, {elapsed:.1f}s elapsed)"
            )
            return False, f"SRA conversion failed (code {result.returncode})", []
    except subprocess.TimeoutExpired:
        elapsed = time.time() - conversion_start
        logger.error(
            f"SRA conversion timed out for {sample_id} after {elapsed/60:.1f} minutes "
            f"(timeout was {timeout_seconds/60:.1f} minutes). Process may be stuck."
        )
        return False, f"SRA conversion timed out after {elapsed/60:.1f} minutes", []
    except Exception as e:
        elapsed = time.time() - conversion_start
        logger.error(
            f"Error converting SRA for {sample_id} after {elapsed:.1f}s: {e}", exc_info=True
        )
        return False, str(e), []


def delete_sample_fastqs(sample_id: str, fastq_dir: Path) -> None:
    """Delete FASTQ files for a specific sample.
    
    Searches for FASTQ files in both getfastq subdirectory and direct structure,
    and removes them to free disk space.
    
    Args:
        sample_id: SRA accession ID (e.g., "SRR1234567")
        fastq_dir: Base directory containing FASTQ files
    """
    logger = get_logger(__name__)
    
    # Check both structures: fastq/getfastq/{sample}/ and fastq/{sample}/
    sample_dirs = []
    getfastq_subdir = fastq_dir / "getfastq" / sample_id
    if getfastq_subdir.exists():
        sample_dirs.append(getfastq_subdir)
    direct_dir = fastq_dir / sample_id
    if direct_dir.exists():
        sample_dirs.append(direct_dir)
    
    # Delete sample directories
    for sample_dir in sample_dirs:
        if sample_dir.exists() and sample_dir.is_dir():
            try:
                shutil.rmtree(sample_dir)
                logger.debug(f"Deleted FASTQ directory: {sample_dir}")
            except Exception as e:
                logger.warning(f"Failed to delete {sample_dir}: {e}")
    
    # Also check for loose FASTQ files
    for pattern in [f"{sample_id}_*.fastq*", f"{sample_id}.fastq*"]:
        for fastq_file in fastq_dir.rglob(pattern):
            try:
                fastq_file.unlink()
                logger.debug(f"Deleted loose FASTQ: {fastq_file.name}")
            except Exception as e:
                logger.warning(f"Failed to delete {fastq_file}: {e}")
