"""FASTQ download step for RNA-seq workflow.

This step downloads FASTQ files from SRA or other sources for selected samples.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Dict

from metainformant.core import logging
from metainformant.rna.steps import StepResult

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
    import requests

    try:
        # ENA API for FASTQ download
        # This is a simplified implementation - real implementation would need proper ENA API calls

        # For now, simulate download
        logger.info(f"Simulating ENA download for {accession}")

        # Create dummy FASTQ files for demonstration
        fastq1 = output_dir / f"{accession}_1.fastq.gz"
        fastq2 = output_dir / f"{accession}_2.fastq.gz"

        # In real implementation, would download actual files
        # For now, just create empty files to simulate success
        fastq1.touch()
        fastq2.touch()

        return {
            'success': True,
            'method': 'ena',
            'files': [str(fastq1), str(fastq2)],
            'accession': accession
        }

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
    import subprocess

    try:
        # Use fasterq-dump or similar tool
        # This is a simplified implementation

        logger.info(f"Simulating SRA download for {accession}")

        # Create dummy FASTQ files
        fastq1 = output_dir / f"{accession}_1.fastq.gz"
        fastq2 = output_dir / f"{accession}_2.fastq.gz"

        fastq1.touch()
        fastq2.touch()

        return {
            'success': True,
            'method': 'sra',
            'files': [str(fastq1), str(fastq2)],
            'accession': accession
        }

    except Exception as e:
        return {
            'success': False,
            'method': 'sra',
            'error': str(e),
            'accession': accession
        }
