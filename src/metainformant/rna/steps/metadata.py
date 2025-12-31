"""Metadata retrieval step for RNA-seq workflow.

This step retrieves metadata for SRA samples from NCBI databases,
including sample information, experimental details, and quality metrics.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Dict

import requests

from metainformant.core import io, logging
from metainformant.rna.steps import StepResult

logger = logging.get_logger(__name__)


def run_step(step_params: Dict[str, Any]) -> StepResult:
    """Execute the metadata retrieval step.

    Args:
        step_params: Parameters for the metadata step
            - species_list: List of species to retrieve metadata for
            - work_dir: Working directory for outputs
            - threads: Number of threads to use
            - api_key: NCBI API key (optional)
            - email: Email for NCBI requests (required)

    Returns:
        StepResult with metadata retrieval results
    """
    start_time = time.time()

    try:
        species_list = step_params.get('species_list', [])
        work_dir = Path(step_params.get('work_dir', '.'))
        threads = step_params.get('threads', 1)
        api_key = step_params.get('api_key')
        email = step_params.get('email')

        if not email:
            raise ValueError("Email is required for NCBI requests")

        if not species_list:
            raise ValueError("Species list cannot be empty")

        # Create output directory
        metadata_dir = work_dir / "metadata"
        metadata_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Retrieving metadata for {len(species_list)} species")

        all_metadata = {}

        for species in species_list:
            logger.info(f"Retrieving metadata for {species}")

            try:
                # Search for SRA samples for this species
                species_metadata = search_sra_metadata(species, api_key, email)

                if species_metadata:
                    # Save species-specific metadata
                    species_file = metadata_dir / f"{species.replace(' ', '_')}_metadata.json"
                    io.dump_json(species_metadata, species_file)

                    all_metadata[species] = species_metadata
                    logger.info(f"Retrieved {len(species_metadata)} samples for {species}")
                else:
                    logger.warning(f"No metadata found for {species}")

            except Exception as e:
                logger.error(f"Failed to retrieve metadata for {species}: {e}")
                continue

        # Save combined metadata
        combined_file = work_dir / "all_species_metadata.json"
        io.dump_json(all_metadata, combined_file)

        # Create sample list for downstream steps
        sample_list = create_sample_list(all_metadata)
        sample_file = work_dir / "sample_list.txt"

        with open(sample_file, 'w') as f:
            for sample in sample_list:
                f.write(f"{sample}\n")

        execution_time = time.time() - start_time

        return StepResult(
            success=True,
            outputs={
                'metadata_dir': str(metadata_dir),
                'combined_metadata': str(combined_file),
                'sample_list': str(sample_file),
                'species_count': len(all_metadata),
                'total_samples': sum(len(samples) for samples in all_metadata.values())
            },
            metadata={
                'species_processed': list(all_metadata.keys()),
                'samples_per_species': {species: len(samples) for species, samples in all_metadata.items()},
                'execution_time': execution_time
            },
            execution_time=execution_time
        )

    except Exception as e:
        execution_time = time.time() - start_time
        logger.error(f"Metadata step failed: {e}")

        return StepResult(
            success=False,
            outputs={},
            metadata={},
            error_message=str(e),
            execution_time=execution_time
        )


def search_sra_metadata(species: str, api_key: Optional[str] = None,
                       email: Optional[str] = None, max_samples: int = 1000) -> Dict[str, Any]:
    """Search NCBI SRA for metadata of a specific species.

    Args:
        species: Species name (e.g., "Drosophila melanogaster")
        api_key: NCBI API key
        email: Email for NCBI requests
        max_samples: Maximum number of samples to retrieve

    Returns:
        Dictionary mapping sample accessions to metadata
    """
    # Use NCBI Entrez API to search SRA
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    # Search query
    query = f'"{species}"[Organism] AND "RNA-seq"[Strategy]'

    # First, search for IDs
    search_params = {
        'db': 'sra',
        'term': query,
        'retmax': min(max_samples, 10000),
        'retmode': 'json',
        'usehistory': 'y'
    }

    if api_key:
        search_params['api_key'] = api_key
    if email:
        search_params['email'] = email

    try:
        response = requests.get(f"{base_url}/esearch.fcgi", params=search_params, timeout=30)
        response.raise_for_status()

        search_data = response.json()

        if 'esearchresult' not in search_data:
            return {}

        id_list = search_data['esearchresult'].get('idlist', [])
        webenv = search_data['esearchresult'].get('webenv')
        query_key = search_data['esearchresult'].get('querykey')

        if not id_list:
            return {}

        # Limit to max_samples
        id_list = id_list[:max_samples]

        # Fetch summaries
        summary_params = {
            'db': 'sra',
            'id': ','.join(id_list),
            'retmode': 'json'
        }

        if webenv and query_key:
            summary_params.update({'WebEnv': webenv, 'query_key': query_key})

        if api_key:
            summary_params['api_key'] = api_key
        if email:
            summary_params['email'] = email

        # Respect NCBI rate limits
        time.sleep(0.4)

        summary_response = requests.get(f"{base_url}/esummary.fcgi", params=summary_params, timeout=30)
        summary_response.raise_for_status()

        summary_data = summary_response.json()

        metadata = {}
        if 'result' in summary_data:
            for uid in id_list:
                if uid in summary_data['result']:
                    sample_data = summary_data['result'][uid]

                    # Extract relevant metadata
                    sample_info = {
                        'accession': sample_data.get('runs', {}).get('run', {}).get('@acc', ''),
                        'experiment': sample_data.get('expxml', {}).get('summary', ''),
                        'platform': sample_data.get('expxml', {}).get('platform', ''),
                        'library_strategy': sample_data.get('expxml', {}).get('library_strategy', ''),
                        'library_layout': sample_data.get('expxml', {}).get('library_layout', ''),
                        'taxon_id': sample_data.get('taxon', {}).get('@taxid', ''),
                        'scientific_name': sample_data.get('taxon', {}).get('@scientificname', ''),
                        'sample_attributes': sample_data.get('sample', {}).get('attributes', {}),
                        'size_MB': sample_data.get('runs', {}).get('run', {}).get('@size_MB', 0),
                        'spots': sample_data.get('runs', {}).get('run', {}).get('@spots', 0),
                        'bases': sample_data.get('runs', {}).get('run', {}).get('@bases', 0)
                    }

                    metadata[uid] = sample_info

        return metadata

    except requests.RequestException as e:
        logger.error(f"SRA metadata search failed for {species}: {e}")
        return {}


def create_sample_list(metadata: Dict[str, Dict[str, Any]]) -> list[str]:
    """Create a list of sample accessions from metadata.

    Args:
        metadata: Metadata dictionary from search_sra_metadata

    Returns:
        List of sample accessions
    """
    sample_list = []

    for species_samples in metadata.values():
        for sample_info in species_samples.values():
            accession = sample_info.get('accession', '')
            if accession and accession not in sample_list:
                sample_list.append(accession)

    return sample_list


def validate_metadata(metadata: Dict[str, Any]) -> Dict[str, Any]:
    """Validate retrieved metadata for completeness and consistency.

    Args:
        metadata: Metadata dictionary to validate

    Returns:
        Validation results
    """
    validation_results = {
        'total_samples': 0,
        'valid_samples': 0,
        'missing_accessions': 0,
        'missing_platform': 0,
        'missing_library_strategy': 0,
        'issues': []
    }

    for species, samples in metadata.items():
        for sample_id, sample_info in samples.items():
            validation_results['total_samples'] += 1

            # Check for required fields
            if not sample_info.get('accession'):
                validation_results['missing_accessions'] += 1
                validation_results['issues'].append(f"Sample {sample_id} missing accession")

            if not sample_info.get('platform'):
                validation_results['missing_platform'] += 1

            if not sample_info.get('library_strategy'):
                validation_results['missing_library_strategy'] += 1

            # Count valid samples
            if (sample_info.get('accession') and
                sample_info.get('platform') and
                sample_info.get('library_strategy')):
                validation_results['valid_samples'] += 1

    return validation_results


