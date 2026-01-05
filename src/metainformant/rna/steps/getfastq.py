"""FASTQ download step for RNA-seq workflow.

This step downloads FASTQ files from SRA or other sources for selected samples.
"""

from __future__ import annotations

import time
import json
from pathlib import Path
from typing import Any, Dict
from urllib.parse import urlparse

from metainformant.core import logging
from metainformant.rna.steps import StepResult

from metainformant.core.download import _utc_iso

logger = logging.get_logger(__name__)


def run_step(step_params: Dict[str, Any]) -> StepResult:
    """Execute the FASTQ download step.

    Args:
        step_params: Parameters for the download step
            - work_dir: Working directory
            - selected_samples: Path to selected samples file
            - download_method: Download method ("sra", "ena", "auto")
            - threads: Number of threads for download
            - max_concurrent: Maximum concurrent downloads

    Returns:
        StepResult with download results
    """
    start_time = time.time()

    try:
        work_dir = Path(step_params.get('work_dir', '.'))
        selected_samples_file = step_params.get('selected_samples',
                                              work_dir / "selection" / "selected_samples.json")
        download_method = step_params.get('download_method', 'ena')
        threads = step_params.get('threads', 1)
        max_concurrent = step_params.get('max_concurrent', min(threads, 4))

        # Load selected samples
        from metainformant.core import io
        if isinstance(selected_samples_file, str):
            selected_samples_file = Path(selected_samples_file)

        if not selected_samples_file.exists():
            raise FileNotFoundError(f"Selected samples file not found: {selected_samples_file}")

        selected_samples = io.load_json(selected_samples_file)

        # Create FASTQ directory
        fastq_dir = work_dir / "fastq"
        fastq_dir.mkdir(parents=True, exist_ok=True)

        download_results = {}

        # Download samples by species
        for species, samples in selected_samples.items():
            species_dir = fastq_dir / species.replace(' ', '_')
            species_dir.mkdir(exist_ok=True)

            logger.info(f"Downloading FASTQ files for {species} ({len(samples)} samples)")

            species_results = download_species_fastq(
                samples, species_dir, download_method, threads, max_concurrent
            )

            download_results[species] = species_results

        # Save download results
        results_file = fastq_dir / "download_results.json"
        io.dump_json(download_results, results_file)

        # Calculate summary
        total_downloaded = sum(len(results) for results in download_results.values())
        successful_downloads = sum(
            1 for species_results in download_results.values()
            for result in species_results.values()
            if result.get('success', False)
        )

        execution_time = time.time() - start_time

        return StepResult(
            success=True,
            outputs={
                'fastq_dir': str(fastq_dir),
                'download_results': str(results_file),
                'total_samples': total_downloaded,
                'successful_downloads': successful_downloads,
                'download_method': download_method
            },
            metadata={
                'threads': threads,
                'max_concurrent': max_concurrent,
                'execution_time': execution_time,
                'success_rate': successful_downloads / total_downloaded if total_downloaded > 0 else 0.0
            },
            execution_time=execution_time
        )

    except Exception as e:
        execution_time = time.time() - start_time
        logger.error(f"FASTQ download step failed: {e}")

        return StepResult(
            success=False,
            outputs={},
            metadata={},
            error_message=str(e),
            execution_time=execution_time
        )


def download_species_fastq(samples: Dict[str, Any], output_dir: Path,
                          method: str, threads: int, max_concurrent: int) -> Dict[str, Any]:
    """Download FASTQ files for samples from one species.

    Args:
        samples: Sample information dictionary
        output_dir: Output directory for FASTQ files
        method: Download method
        threads: Number of threads
        max_concurrent: Maximum concurrent downloads

    Returns:
        Dictionary with download results per sample
    """
    results = {}

    # Process samples (in a real implementation, would use parallel processing)
    for sample_id, sample_info in samples.items():
        accession = sample_info.get('original_info', {}).get('accession')

        if not accession:
            results[sample_id] = {
                'success': False,
                'error': 'No accession available'
            }
            continue

        logger.info(f"Downloading FASTQ for {accession}")

        try:
            result = download_sample_fastq(accession, output_dir, method, threads)
            results[sample_id] = result

            if result['success']:
                logger.info(f"Successfully downloaded {accession}")
            else:
                logger.warning(f"Failed to download {accession}: {result.get('error', 'Unknown error')}")

        except Exception as e:
            logger.error(f"Error downloading {accession}: {e}")
            results[sample_id] = {
                'success': False,
                'error': str(e)
            }

    return results


def download_sample_fastq(accession: str, output_dir: Path, method: str, threads: int) -> Dict[str, Any]:
    """Download FASTQ files for a single sample.

    Args:
        accession: Sample accession
        output_dir: Output directory
        method: Download method
        threads: Number of threads

    Returns:
        Download result dictionary
    """
    if method == 'ena':
        return download_from_ena(accession, output_dir, threads)
    elif method == 'sra':
        return download_from_sra(accession, output_dir, threads)
    elif method == 'auto':
        # Try ENA first, then SRA
        result = download_from_ena(accession, output_dir, threads)
        if not result['success']:
            result = download_from_sra(accession, output_dir, threads)
        return result
    else:
        return {
            'success': False,
            'error': f'Unknown download method: {method}'
        }


def download_from_ena(accession: str, output_dir: Path, threads: int) -> Dict[str, Any]:
    """Download FASTQ from ENA (European Nucleotide Archive).

    Args:
        accession: Sample accession
        output_dir: Output directory
        threads: Number of threads

    Returns:
        Download result
    """
    try:
        import requests

        from metainformant.core.download import download_with_progress

        output_dir.mkdir(parents=True, exist_ok=True)

        # ENA filereport API returns fastq URLs and sizes.
        # Docs: https://www.ebi.ac.uk/ena/portal/api/
        api_url = (
            "https://www.ebi.ac.uk/ena/portal/api/filereport"
            f"?accession={accession}&result=read_run&fields=fastq_http,fastq_bytes&format=tsv"
        )
        resp = requests.get(api_url, timeout=60)
        resp.raise_for_status()
        lines = [ln for ln in resp.text.splitlines() if ln.strip()]
        if len(lines) < 2:
            raise RuntimeError(f"ENA filereport returned no rows for {accession}")

        header = lines[0].split("\t")
        row = lines[1].split("\t")
        record = dict(zip(header, row, strict=False))

        fastq_http = (record.get("fastq_http") or "").strip()
        if not fastq_http:
            raise RuntimeError(f"ENA has no fastq_http for {accession}")

        urls = [u for u in fastq_http.split(";") if u]
        if not urls:
            raise RuntimeError(f"ENA has no FASTQ URLs for {accession}")

        downloaded: list[str] = []
        for url in urls:
            filename = Path(urlparse(url).path).name or f"{accession}.fastq.gz"
            dest = output_dir / filename
            result = download_with_progress(
                url,
                dest,
                protocol="https",
                show_progress=True,
                heartbeat_interval=5,
                max_retries=3,
                chunk_size=1024 * 1024,
                timeout=300,
                resume=True,
            )
            if not result.success:
                raise RuntimeError(f"Failed downloading {url}: {result.error or 'unknown error'}")
            downloaded.append(str(dest))

        return {"success": True, "method": "ena", "files": downloaded, "accession": accession}

    except Exception as e:
        return {
            'success': False,
            'method': 'ena',
            'error': str(e),
            'accession': accession
        }


def download_from_sra(accession: str, output_dir: Path, threads: int) -> Dict[str, Any]:
    """Download FASTQ from SRA (Sequence Read Archive).

    Args:
        accession: Sample accession
        output_dir: Output directory
        threads: Number of threads

    Returns:
        Download result
    """
    try:
        import shutil
        import subprocess

        output_dir.mkdir(parents=True, exist_ok=True)

        fasterq = shutil.which("fasterq-dump")
        if not fasterq:
            raise RuntimeError("SRA toolkit not available: missing fasterq-dump on PATH")

        # Heartbeat during subprocess: track total bytes of FASTQ outputs in output_dir.
        hb_path = (output_dir / ".downloads")
        hb_path.mkdir(parents=True, exist_ok=True)
        hb_file = hb_path / f"{accession}.sra_toolkit.heartbeat.json"

        def _write_hb(status: str, errors: list[str] | None = None) -> None:
            files = list(output_dir.glob(f"{accession}*.fastq")) + list(output_dir.glob(f"{accession}*.fastq.gz"))
            bytes_done = sum(p.stat().st_size for p in files if p.exists())
            payload = {
                "url": accession,
                "destination": str(output_dir),
                "started_at": _utc_iso(),
                "last_update": _utc_iso(),
                "bytes_downloaded": bytes_done,
                "total_bytes": None,
                "progress_percent": None,
                "speed_mbps": None,
                "eta_seconds": None,
                "status": status,
                "errors": errors or [],
            }
            with open(hb_file, "w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=2, sort_keys=True)
                fh.write("\n")

        cmd = [
            fasterq,
            "--split-files",
            "--threads",
            str(max(1, int(threads))),
            "--outdir",
            str(output_dir),
            accession,
        ]
        logger.info(f"Running: {' '.join(cmd)}")
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        start = time.time()
        last = 0.0
        stderr_tail: list[str] = []
        while True:
            rc = proc.poll()
            now = time.time()
            if now - last >= 5:
                _write_hb("downloading")
                last = now
            if rc is not None:
                break
            # Avoid busy loop
            time.sleep(1)

        out, err = proc.communicate(timeout=30)
        if err:
            stderr_tail = err.splitlines()[-20:]
        if proc.returncode != 0:
            _write_hb("failed", errors=stderr_tail)
            raise RuntimeError(f"fasterq-dump failed for {accession}: {stderr_tail[-1] if stderr_tail else 'unknown'}")

        _write_hb("completed")

        # Collect produced FASTQs
        files = sorted([str(p) for p in output_dir.glob(f"{accession}*.fastq*")])
        if not files:
            raise RuntimeError(f"fasterq-dump produced no FASTQ files for {accession}")

        return {"success": True, "method": "sra", "files": files, "accession": accession}

    except Exception as e:
        return {
            'success': False,
            'method': 'sra',
            'error': str(e),
            'accession': accession
        }



run_getfastq = run_step




