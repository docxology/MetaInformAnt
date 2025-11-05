"""Step runner for `amalgkit getfastq` (download raw FASTQ files)."""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import urllib.request
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ...core.io import read_delimited
from ..amalgkit import getfastq as _getfastq


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

    # Enable parallel-fastq-dump for faster downloads when available
    if "pfd" not in params:
        params["pfd"] = True if pfd_path else False  # Use PFD when available
    if pfd_path and not params.get("pfd_exe"):
        params["pfd_exe"] = pfd_path
    if prefetch_path and not params.get("prefetch_exe"):
        params["prefetch_exe"] = prefetch_path

    # If user explicitly enabled PFD, provide its path when available
    if bool(params.get("pfd")) and pfd_path and not params.get("pfd_exe"):
        params["pfd_exe"] = pfd_path

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
    
    # 1) Bulk attempt (metadata or id)
    bulk_result = _getfastq(
        effective_params,
        work_dir=work_dir,
        log_dir=log_dir,
        step_name="getfastq",
        check=False,
    )

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
            candidate_paths = [
                actual_work_dir / "metadata" / "pivot_qualified.tsv",
                actual_work_dir / "metadata" / "pivot_selected.tsv",
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
            
            # Debug: print what we're looking for
            import logging
            logging.getLogger(__name__).debug(f"Looking for metadata at: {meta_path}")
            logging.getLogger(__name__).debug(f"Work dir: {actual_work_dir}")
            logging.getLogger(__name__).debug(f"Out dir: {out_dir}")
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
            import logging
            logging.getLogger(__name__).warning(f"Could not read metadata file {meta_path}: {e}")
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

    # 3) Targeted retries for missing SRRs
    missing = [s for s in srr_list if not _has_fastq(s)] if srr_list else []
    for srr in missing:
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
                    import logging
                    logging.getLogger(__name__).debug(f"Could not copy file for {srr}: {e}")
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

    # 4) Final status: success if all SRRs (if known) have FASTQs
    final_missing = [s for s in srr_list if not _has_fastq(s)] if srr_list else []
    rc = 0 if not final_missing and bulk_result.returncode == 0 else (1 if final_missing else bulk_result.returncode)

    # 5) Optional cleanup: remove raw FASTQ files if cleanup_raw is enabled
    cleanup_raw = effective_params.get("cleanup_raw", False)
    if cleanup_raw and rc == 0 and not final_missing:
        logging.getLogger(__name__).info("Cleaning up raw FASTQ files after successful download")
        try:
            for srr in srr_list:
                srr_dir = out_dir / "getfastq" / srr
                if srr_dir.exists():
                    # Remove raw FASTQ files, keep gzipped versions if they exist
                    for fastq_file in srr_dir.glob("*.fastq"):
                        if not fastq_file.name.endswith('.gz'):
                            fastq_file.unlink()
                            logging.getLogger(__name__).debug(f"Removed raw FASTQ file: {fastq_file}")
        except Exception as e:
            logging.getLogger(__name__).warning(f"Failed to cleanup raw FASTQ files: {e}")

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
    logger = logging.getLogger(__name__)
    
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
    
    logger.info(f"Converting SRA to FASTQ for {sample_id} (SRA: {sra_file.name})...")
    
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
            with open(log_file, "w") if log_file else open(os.devnull, "w") as f:
                result = subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    cwd=str(output_dir),
                    check=False,
                )
            
            # Check for output files
            fastq_files = list(output_dir.glob(f"{sample_id}_*.fastq.gz"))
            if not fastq_files:
                fastq_files = list(output_dir.glob(f"{sample_id}.fastq.gz"))
            
            if fastq_files:
                logger.info(f"Converted SRA to FASTQ for {sample_id} using parallel-fastq-dump ({len(fastq_files)} files)")
                return True, f"Converted SRA to FASTQ for {sample_id}", fastq_files
            
            # If parallel-fastq-dump failed, fall through to fasterq-dump
            logger.warning(f"parallel-fastq-dump failed (code {result.returncode}), trying fasterq-dump")
        except Exception as e:
            logger.warning(f"parallel-fastq-dump error: {e}, trying fasterq-dump")
    
    # Fallback to fasterq-dump
    fasterq_dump = shutil.which("fasterq-dump")
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
        "-p",  # Show progress
    ]
    
    log_file = None
    if log_dir:
        log_file = log_dir / f"fasterq_dump_{sample_id}.log"
        log_file.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        with open(log_file, "w") if log_file else open(os.devnull, "w") as f:
            result = subprocess.run(
                cmd,
                stdout=f,
                stderr=subprocess.STDOUT,
                cwd=str(sra_dir),  # Run from directory containing SRA file
                check=False,
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
            logger.info(f"Converted SRA to FASTQ for {sample_id} using fasterq-dump ({len(fastq_files)} files)")
            return True, f"Converted SRA to FASTQ for {sample_id}", fastq_files
        else:
            logger.warning(f"SRA conversion failed for {sample_id} (code {result.returncode})")
            return False, f"SRA conversion failed (code {result.returncode})", []
    except Exception as e:
        logger.error(f"Error converting SRA for {sample_id}: {e}", exc_info=True)
        return False, str(e), []


def delete_sample_fastqs(sample_id: str, fastq_dir: Path) -> None:
    """Delete FASTQ files for a specific sample.
    
    Searches for FASTQ files in both getfastq subdirectory and direct structure,
    and removes them to free disk space.
    
    Args:
        sample_id: SRA accession ID (e.g., "SRR1234567")
        fastq_dir: Base directory containing FASTQ files
    """
    logger = logging.getLogger(__name__)
    
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
