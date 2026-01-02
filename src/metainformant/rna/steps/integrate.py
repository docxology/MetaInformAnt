"""Integration step for RNA-seq workflow.

This step integrates local FASTQ files with metadata retrieved from NCBI SRA,
harmonizing sample information and preparing for downstream analysis.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Dict, List

from metainformant.core import io, logging, paths
from metainformant.rna.steps import StepResult

logger = logging.get_logger(__name__)


def run_step(step_params: Dict[str, Any]) -> StepResult:
    """Execute the integration step.

    Args:
        step_params: Parameters for the integration step
            - work_dir: Working directory
            - metadata_file: Path to metadata JSON file
            - local_fastq_dir: Directory containing local FASTQ files (optional)
            - species_list: List of species to process
            - threads: Number of threads to use

    Returns:
        StepResult with integration results
    """
    start_time = time.time()

    try:
        work_dir = Path(step_params.get('work_dir', '.'))
        metadata_file = step_params.get('metadata_file', work_dir / "all_species_metadata.json")
        local_fastq_dir = step_params.get('local_fastq_dir')
        species_list = step_params.get('species_list', [])
        threads = step_params.get('threads', 1)

        # Load metadata
        if isinstance(metadata_file, str):
            metadata_file = Path(metadata_file)

        if not metadata_file.exists():
            raise FileNotFoundError(f"Metadata file not found: {metadata_file}")

        metadata = io.load_json(metadata_file)
        logger.info(f"Loaded metadata for {len(metadata)} species")

        # Create integration directory
        integration_dir = work_dir / "integration"
        integration_dir.mkdir(parents=True, exist_ok=True)

        integrated_data = {}

        for species in species_list:
            if species not in metadata:
                logger.warning(f"Species {species} not found in metadata")
                continue

            logger.info(f"Integrating data for {species}")

            species_samples = metadata[species]
            integrated_samples = integrate_species_data(
                species, species_samples, local_fastq_dir, integration_dir
            )

            integrated_data[species] = integrated_samples

            logger.info(f"Integrated {len(integrated_samples)} samples for {species}")

        # Save integrated data
        integrated_file = integration_dir / "integrated_metadata.json"
        io.dump_json(integrated_data, integrated_file)

        # Create sample information summary
        summary = create_integration_summary(integrated_data)
        summary_file = integration_dir / "integration_summary.json"
        io.dump_json(summary, summary_file)

        execution_time = time.time() - start_time

        return StepResult(
            success=True,
            outputs={
                'integration_dir': str(integration_dir),
                'integrated_metadata': str(integrated_file),
                'summary': str(summary_file),
                'species_integrated': len(integrated_data),
                'total_samples': sum(len(samples) for samples in integrated_data.values())
            },
            metadata={
                'integration_summary': summary,
                'execution_time': execution_time
            },
            execution_time=execution_time
        )

    except Exception as e:
        execution_time = time.time() - start_time
        logger.error(f"Integration step failed: {e}")

        return StepResult(
            success=False,
            outputs={},
            metadata={},
            error_message=str(e),
            execution_time=execution_time
        )


def integrate_species_data(species: str, species_samples: Dict[str, Any],
                          local_fastq_dir: str | None,
                          integration_dir: Path) -> Dict[str, Any]:
    """Integrate sample data for a specific species.

    Args:
        species: Species name
        species_samples: Sample metadata for the species
        local_fastq_dir: Directory with local FASTQ files
        integration_dir: Output directory for integration results

    Returns:
        Dictionary with integrated sample information
    """
    integrated_samples = {}

    # Look for local FASTQ files if directory provided
    local_files = {}
    if local_fastq_dir:
        local_fastq_path = Path(local_fastq_dir)
        if local_fastq_path.exists():
            # Find FASTQ files
            fastq_files = list(local_fastq_path.glob("*.fastq")) + \
                         list(local_fastq_path.glob("*.fastq.gz")) + \
                         list(local_fastq_path.glob("*.fq")) + \
                         list(local_fastq_path.glob("*.fq.gz"))

            for fastq_file in fastq_files:
                # Try to match filename to accession
                filename = fastq_file.stem
                if filename.endswith('.fastq') or filename.endswith('.fq'):
                    filename = filename.rsplit('.', 1)[0]

                local_files[filename] = str(fastq_file)

    # Integrate each sample
    for sample_id, sample_info in species_samples.items():
        accession = sample_info.get('accession', '')

        integrated_sample = {
            'original_metadata': sample_info,
            'accession': accession,
            'species': species,
            'local_fastq_files': [],
            'sra_available': bool(accession),
            'integration_status': 'pending'
        }

        # Check for local FASTQ files
        if accession in local_files:
            integrated_sample['local_fastq_files'].append(local_files[accession])
            integrated_sample['integration_status'] = 'local_available'
        elif local_files:
            # Look for partial matches
            matching_files = [f for name, f in local_files.items() if accession in name]
            if matching_files:
                integrated_sample['local_fastq_files'].extend(matching_files)
                integrated_sample['integration_status'] = 'local_available'

        # Validate sample information
        integrated_sample['validation'] = validate_sample_info(sample_info)

        integrated_samples[sample_id] = integrated_sample

    return integrated_samples


def validate_sample_info(sample_info: Dict[str, Any]) -> Dict[str, Any]:
    """Validate sample information for completeness.

    Args:
        sample_info: Sample metadata dictionary

    Returns:
        Validation results
    """
    validation = {
        'is_valid': True,
        'issues': [],
        'completeness_score': 0.0
    }

    required_fields = ['accession', 'platform', 'library_strategy']
    optional_fields = ['spots', 'bases', 'size_MB', 'taxon_id', 'scientific_name']

    total_fields = len(required_fields) + len(optional_fields)
    present_fields = 0

    # Check required fields
    for field in required_fields:
        if not sample_info.get(field):
            validation['issues'].append(f"Missing required field: {field}")
            validation['is_valid'] = False
        else:
            present_fields += 1

    # Check optional fields
    for field in optional_fields:
        if sample_info.get(field):
            present_fields += 1

    validation['completeness_score'] = present_fields / total_fields

    # Additional validation
    if sample_info.get('library_strategy') != 'RNA-seq':
        validation['issues'].append("Library strategy is not RNA-seq")
        validation['is_valid'] = False

    return validation


def create_integration_summary(integrated_data: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """Create a summary of integration results.

    Args:
        integrated_data: Integrated sample data

    Returns:
        Summary statistics
    """
    summary = {
        'total_species': len(integrated_data),
        'total_samples': 0,
        'samples_with_local_files': 0,
        'samples_sra_only': 0,
        'validation_summary': {
            'valid_samples': 0,
            'invalid_samples': 0,
            'avg_completeness': 0.0
        },
        'platform_distribution': {},
        'library_strategy_distribution': {},
        'species_breakdown': {}
    }

    total_completeness = 0.0

    for species, samples in integrated_data.items():
        species_summary = {
            'sample_count': len(samples),
            'local_files': 0,
            'sra_available': 0,
            'valid_samples': 0
        }

        for sample_info in samples.values():
            summary['total_samples'] += 1

            if sample_info.get('local_fastq_files'):
                summary['samples_with_local_files'] += 1
                species_summary['local_files'] += 1
            elif sample_info.get('sra_available'):
                summary['samples_sra_only'] += 1
                species_summary['sra_available'] += 1

            # Validation
            validation = sample_info.get('validation', {})
            if validation.get('is_valid', False):
                summary['validation_summary']['valid_samples'] += 1
                species_summary['valid_samples'] += 1
            else:
                summary['validation_summary']['invalid_samples'] += 1

            total_completeness += validation.get('completeness_score', 0.0)

            # Platform distribution
            platform = sample_info.get('original_metadata', {}).get('platform', 'Unknown')
            summary['platform_distribution'][platform] = \
                summary['platform_distribution'].get(platform, 0) + 1

            # Library strategy distribution
            strategy = sample_info.get('original_metadata', {}).get('library_strategy', 'Unknown')
            summary['library_strategy_distribution'][strategy] = \
                summary['library_strategy_distribution'].get(strategy, 0) + 1

        summary['species_breakdown'][species] = species_summary

    # Calculate averages
    if summary['total_samples'] > 0:
        summary['validation_summary']['avg_completeness'] = total_completeness / summary['total_samples']

    return summary


def merge_metadata_sources(primary_metadata: Dict[str, Any],
                          secondary_metadata: Dict[str, Any]) -> Dict[str, Any]:
    """Merge metadata from multiple sources.

    Args:
        primary_metadata: Primary metadata source
        secondary_metadata: Secondary metadata source

    Returns:
        Merged metadata
    """
    merged = primary_metadata.copy()

    for species, samples in secondary_metadata.items():
        if species not in merged:
            merged[species] = {}

        for sample_id, sample_info in samples.items():
            if sample_id not in merged[species]:
                merged[species][sample_id] = sample_info
            else:
                # Merge sample information
                existing = merged[species][sample_id]
                for key, value in sample_info.items():
                    if key not in existing or not existing[key]:
                        existing[key] = value

    return merged



run_integrate = run_step
