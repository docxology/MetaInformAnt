"""Step runner for `amalgkit getfastq` (download raw FASTQ files)."""

from __future__ import annotations

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

    # Prefer sra-tools path for maximum robustness; only use PFD when explicitly requested
    if "pfd" not in params:
        params["pfd"] = False  # default to no PFD due to observed instability
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
    id_val = str(effective_params.get("id", "")).strip()
    if id_val:
        srr_list = [id_val]
    else:
        meta_path = effective_params.get("metadata")
        if not meta_path or meta_path == "inferred":
            meta_path = str(Path(out_dir).parent / "metadata" / "metadata.tsv")
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
            logging.getLogger(__name__).warning(f"Could not read metadata file {metadata_file}: {e}")
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
        if bool(effective_params.get("accelerate")):
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
    result = subprocess.CompletedProcess(["amalgkit", "getfastq"], rc, stdout="", stderr="")
    if check and rc != 0:
        raise subprocess.CalledProcessError(rc, result.args)
    return result
