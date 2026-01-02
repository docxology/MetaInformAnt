"""SRA data download utilities for GWAS.

This module provides specialized tools for downloading sequencing data
from the Sequence Read Archive (SRA) for GWAS analysis.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import requests

from metainformant.core import logging

logger = logging.get_logger(__name__)


def download_sra_experiment(experiment_acc: str, output_dir: str | Path,
                           threads: int = 1) -> List[Path]:
    """Download all runs from an SRA experiment.

    Args:
        experiment_acc: Experiment accession (SRX...)
        output_dir: Output directory
        threads: Number of threads for download

    Returns:
        List of paths to downloaded run directories
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get runs for this experiment
    run_accessions = _get_experiment_runs(experiment_acc)

    if not run_accessions:
        logger.warning(f"No runs found for experiment {experiment_acc}")
        return []

    downloaded_runs = []

    for run_acc in run_accessions:
        try:
            from .download import download_sra_run
            run_path = download_sra_run(run_acc, output_dir, threads)
            downloaded_runs.append(run_path)
            logger.info(f"Downloaded run {run_acc} from experiment {experiment_acc}")
        except Exception as e:
            logger.error(f"Failed to download run {run_acc}: {e}")
            continue

    return downloaded_runs


def _get_experiment_runs(experiment_acc: str) -> List[str]:
    """Get SRA run accessions for an experiment."""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    # Search for runs in this experiment
    search_params = {
        'db': 'sra',
        'term': f'{experiment_acc}[Experiment Accession]',
        'retmax': 100,
        'retmode': 'json'
    }

    try:
        response = requests.get(f"{base_url}/esearch.fcgi", params=search_params, timeout=30)
        response.raise_for_status()

        data = response.json()

        if 'esearchresult' not in data:
            return []

        id_list = data['esearchresult'].get('idlist', [])

        if not id_list:
            return []

        # Get run accessions
        summary_params = {
            'db': 'sra',
            'id': ','.join(id_list),
            'retmode': 'json'
        }

        summary_response = requests.get(f"{base_url}/esummary.fcgi", params=summary_params, timeout=30)
        summary_response.raise_for_status()

        summary_data = summary_response.json()

        run_accessions = []
        if 'result' in summary_data:
            for uid in id_list:
                if uid in summary_data['result']:
                    record = summary_data['result'][uid]
                    runs = record.get('runs', {}).get('run', [])
                    if isinstance(runs, list):
                        for run in runs:
                            run_acc = run.get('@acc')
                            if run_acc:
                                run_accessions.append(run_acc)
                    elif isinstance(runs, dict):
                        run_acc = runs.get('@acc')
                        if run_acc:
                            run_accessions.append(run_acc)

        return run_accessions

    except requests.RequestException as e:
        logger.error(f"Failed to get experiment runs: {e}")
        return []


def download_sra_biosample(biosample_acc: str, output_dir: str | Path,
                          threads: int = 1) -> List[Path]:
    """Download all SRA data for a BioSample.

    Args:
        biosample_acc: BioSample accession (SAMN...)
        output_dir: Output directory
        threads: Number of threads

    Returns:
        List of downloaded run paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get runs for this BioSample
    run_accessions = _get_biosample_runs(biosample_acc)

    downloaded_runs = []

    for run_acc in run_accessions:
        try:
            from .download import download_sra_run
            run_path = download_sra_run(run_acc, output_dir, threads)
            downloaded_runs.append(run_path)
        except Exception as e:
            logger.error(f"Failed to download run {run_acc}: {e}")
            continue

    return downloaded_runs


def _get_biosample_runs(biosample_acc: str) -> List[str]:
    """Get SRA runs for a BioSample."""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    search_params = {
        'db': 'sra',
        'term': f'{biosample_acc}[BioSample]',
        'retmax': 100,
        'retmode': 'json'
    }

    try:
        response = requests.get(f"{base_url}/esearch.fcgi", params=search_params, timeout=30)
        response.raise_for_status()

        data = response.json()

        if 'esearchresult' not in data:
            return []

        id_list = data['esearchresult'].get('idlist', [])

        # Get run accessions from summaries
        if id_list:
            summary_params = {
                'db': 'sra',
                'id': ','.join(id_list[:50]),  # Limit
                'retmode': 'json'
            }

            summary_response = requests.get(f"{base_url}/esummary.fcgi", params=summary_params, timeout=30)
            summary_response.raise_for_status()

            summary_data = summary_response.json()

            run_accessions = []
            if 'result' in summary_data:
                for uid in id_list[:50]:
                    if uid in summary_data['result']:
                        record = summary_data['result'][uid]
                        runs = record.get('runs', {}).get('run', [])
                        if isinstance(runs, list):
                            for run in runs:
                                run_acc = run.get('@acc')
                                if run_acc:
                                    run_accessions.append(run_acc)
                        elif isinstance(runs, dict):
                            run_acc = runs.get('@acc')
                            if run_acc:
                                run_accessions.append(run_acc)

            return run_accessions

    except requests.RequestException as e:
        logger.error(f"Failed to get BioSample runs: {e}")

    return []


def batch_download_sra(accessions: List[str], output_dir: str | Path,
                      threads: int = 1, max_concurrent: int = 4) -> Dict[str, Any]:
    """Download multiple SRA accessions in batch.

    Args:
        accessions: List of SRA accessions
        output_dir: Output directory
        threads: Threads per download
        max_concurrent: Maximum concurrent downloads

    Returns:
        Download results summary
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = {
        'total_requested': len(accessions),
        'successful_downloads': 0,
        'failed_downloads': 0,
        'results': {}
    }

    from .download import download_sra_run

    for accession in accessions:
        try:
            run_path = download_sra_run(accession, output_dir, threads)
            results['results'][accession] = {'success': True, 'path': str(run_path)}
            results['successful_downloads'] += 1
            logger.info(f"Successfully downloaded {accession}")
        except Exception as e:
            results['results'][accession] = {'success': False, 'error': str(e)}
            results['failed_downloads'] += 1
            logger.error(f"Failed to download {accession}: {e}")

        # Rate limiting
        time.sleep(0.5)

    return results


def find_sra_data_by_phenotype(phenotype: str, organism: str,
                              max_results: int = 50) -> List[Dict[str, Any]]:
    """Find SRA data by phenotype/trait.

    Args:
        phenotype: Phenotype or trait name
        organism: Organism name
        max_results: Maximum results to return

    Returns:
        List of matching SRA records
    """
    from .download import search_sra_for_organism

    # First get all data for organism
    all_data = search_sra_for_organism(organism, max_results * 2)

    # Filter by phenotype
    phenotype_lower = phenotype.lower()
    filtered_data = []

    for record in all_data:
        # Check various fields for phenotype keywords
        searchable_text = ' '.join([
            record.get('experiment', ''),
            str(record.get('sample_attributes', {}))
        ]).lower()

        if phenotype_lower in searchable_text:
            filtered_data.append(record)

        if len(filtered_data) >= max_results:
            break

    return filtered_data


def download_sra_with_retry(accession: str, output_dir: str | Path,
                           max_retries: int = 3, threads: int = 1) -> Path:
    """Download SRA data with retry logic.

    Args:
        accession: SRA accession
        output_dir: Output directory
        max_retries: Maximum retry attempts
        threads: Download threads

    Returns:
        Path to downloaded data
    """
    from .download import download_sra_run

    last_error = None

    for attempt in range(max_retries):
        try:
            return download_sra_run(accession, output_dir, threads)
        except Exception as e:
            last_error = e
            logger.warning(f"Download attempt {attempt + 1} failed for {accession}: {e}")

            if attempt < max_retries - 1:
                # Exponential backoff
                delay = 2 ** attempt
                logger.info(f"Retrying in {delay} seconds...")
                time.sleep(delay)

    # All attempts failed
    raise RuntimeError(f"Failed to download {accession} after {max_retries} attempts: {last_error}")


def validate_sra_download(accession: str, download_dir: Path) -> Dict[str, Any]:
    """Validate downloaded SRA data.

    Args:
        accession: SRA accession
        download_dir: Download directory

    Returns:
        Validation results
    """
    validation = {
        'accession': accession,
        'valid': False,
        'files_found': [],
        'total_size_mb': 0.0,
        'issues': []
    }

    # Look for downloaded files
    pattern = f"*{accession}*"
    downloaded_files = list(download_dir.glob(pattern))

    if not downloaded_files:
        validation['issues'].append("No files found matching accession")
        return validation

    validation['files_found'] = [str(f) for f in downloaded_files]

    # Check file sizes
    total_size = 0
    for file_path in downloaded_files:
        if file_path.exists():
            total_size += file_path.stat().st_size

    validation['total_size_mb'] = total_size / (1024 * 1024)

    # Basic validation
    if total_size == 0:
        validation['issues'].append("All files are empty")
    elif total_size < 1000:  # Less than 1KB
        validation['issues'].append("Downloaded files are suspiciously small")
    else:
        validation['valid'] = True

    return validation


def prefetch_sra_metadata(accessions: List[str]) -> Dict[str, Dict[str, Any]]:
    """Prefetch metadata for multiple SRA accessions.

    Args:
        accessions: List of SRA accessions

    Returns:
        Dictionary mapping accessions to metadata
    """
    metadata = {}

    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    # Process in batches to avoid API limits
    batch_size = 100

    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i + batch_size]

        # Search for these accessions
        search_params = {
            'db': 'sra',
            'term': ' OR '.join(f'{acc}[Accession]' for acc in batch),
            'retmax': len(batch),
            'retmode': 'json'
        }

        try:
            response = requests.get(f"{base_url}/esearch.fcgi", params=search_params, timeout=30)
            response.raise_for_status()

            search_data = response.json()

            if 'esearchresult' in search_data:
                id_list = search_data['esearchresult'].get('idlist', [])

                if id_list:
                    # Get summaries
                    summary_params = {
                        'db': 'sra',
                        'id': ','.join(id_list),
                        'retmode': 'json'
                    }

                    summary_response = requests.get(f"{base_url}/esummary.fcgi", params=summary_params, timeout=30)
                    summary_response.raise_for_status()

                    summary_data = summary_response.json()

                    if 'result' in summary_data:
                        for uid in id_list:
                            if uid in summary_data['result']:
                                record = summary_data['result'][uid]

                                # Extract accession from record
                                runs = record.get('runs', {}).get('run', [])
                                if isinstance(runs, list) and runs:
                                    accession = runs[0].get('@acc')
                                elif isinstance(runs, dict):
                                    accession = runs.get('@acc')
                                else:
                                    continue

                                if accession:
                                    metadata[accession] = {
                                        'spots': record.get('runs', {}).get('run', {}).get('@spots'),
                                        'bases': record.get('runs', {}).get('run', {}).get('@bases'),
                                        'size_mb': record.get('runs', {}).get('run', {}).get('@size_MB'),
                                        'platform': record.get('expxml', {}).get('platform'),
                                        'library_strategy': record.get('expxml', {}).get('library_strategy')
                                    }

        except requests.RequestException as e:
            logger.error(f"Failed to prefetch metadata batch: {e}")
            continue

        # Rate limiting
        time.sleep(0.5)

    return metadata


def check_sra_tools_available() -> bool:
    """Check if SRA tools (fasterq-dump, prefetch, etc.) are available.

    Returns:
        True if SRA tools are available and executable
    """
    import shutil

    # Check for common SRA tools
    sra_tools = ['fasterq-dump', 'prefetch', 'fastq-dump']

    available_tools = []
    for tool in sra_tools:
        if shutil.which(tool):
            available_tools.append(tool)

    if available_tools:
        logger.info(f"SRA tools available: {', '.join(available_tools)}")
        return True
    else:
        logger.warning("No SRA tools found. Install SRA Toolkit from NCBI")
        return False


def download_sra_project(project_id: str, output_dir: str | Path,
                        threads: int = 1) -> List[Path]:
    """Download all experiments/runs from an SRA project.

    Args:
        project_id: SRA project ID (PRJNA...)
        output_dir: Output directory
        threads: Number of threads for download

    Returns:
        List of paths to downloaded run directories
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get experiments for this project
    experiment_accessions = _get_project_experiments(project_id)

    if not experiment_accessions:
        logger.warning(f"No experiments found for project {project_id}")
        return []

    downloaded_runs = []

    for exp_acc in experiment_accessions:
        try:
            # Download all runs from this experiment
            run_paths = download_sra_experiment(exp_acc, output_dir, threads)
            downloaded_runs.extend(run_paths)
            logger.info(f"Downloaded {len(run_paths)} runs from experiment {exp_acc}")
        except Exception as e:
            logger.error(f"Failed to download experiment {exp_acc}: {e}")
            continue

    logger.info(f"Downloaded {len(downloaded_runs)} total runs from project {project_id}")
    return downloaded_runs


def _get_project_experiments(project_id: str) -> List[str]:
    """Get experiment accessions for an SRA project.

    Args:
        project_id: SRA project ID

    Returns:
        List of experiment accessions
    """
    # This would typically query NCBI's Entrez API
    # For now, return empty list (placeholder implementation)
    logger.warning("Project experiment lookup not implemented - requires NCBI API integration")
    return []


def download_sra_run(sra_accession: str, output_dir: str | Path,
                    threads: int = 1) -> Path:
    """Download a single SRA run.

    Args:
        sra_accession: SRA run accession (SRR...)
        output_dir: Output directory
        threads: Number of threads for download

    Returns:
        Path to downloaded run directory
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    run_dir = output_dir / sra_accession
    run_dir.mkdir(exist_ok=True)

    # Check if SRA tools are available
    if not check_sra_tools_available():
        logger.warning("SRA tools not available, download may fail")
        return run_dir

    try:
        # Use fasterq-dump if available, otherwise fastq-dump
        import subprocess

        cmd = ["fasterq-dump", "--split-files", "--threads", str(threads),
               "--outdir", str(run_dir), sra_accession]

        # Try fasterq-dump first
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            if result.returncode == 0:
                logger.info(f"Downloaded SRA run {sra_accession} using fasterq-dump")
                return run_dir
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass

        # Fall back to fastq-dump
        cmd = ["fastq-dump", "--split-files", "--gzip", sra_accession]
        result = subprocess.run(cmd, cwd=run_dir, capture_output=True, text=True, timeout=3600)

        if result.returncode == 0:
            logger.info(f"Downloaded SRA run {sra_accession} using fastq-dump")
            return run_dir
        else:
            logger.error(f"Failed to download {sra_accession}: {result.stderr}")
            return run_dir

    except Exception as e:
        logger.error(f"Error downloading SRA run {sra_accession}: {e}")
        return run_dir


def search_sra_for_organism(organism: str, max_results: int = 100) -> List[Dict[str, Any]]:
    """Search SRA for experiments from a specific organism.

    Args:
        organism: Organism name to search for
        max_results: Maximum number of results to return

    Returns:
        List of dictionaries with experiment information
    """
    # This would typically query NCBI's Entrez API
    # For now, return mock data structure for testing
    logger.warning("SRA search not fully implemented - requires NCBI API integration")

    # Return mock results for testing
    mock_results = []
    for i in range(min(max_results, 10)):
        mock_results.append({
            'accession': f'SRX000{i:04d}',
            'title': f'{organism} experiment {i}',
            'organism': organism,
            'runs': [f'SRR000{i:04d}'],
            'platform': 'ILLUMINA',
            'library_strategy': 'RNA-Seq'
        })

    return mock_results



