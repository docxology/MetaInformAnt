"""Configuration step for RNA-seq workflow.

This step generates configuration files for downstream analysis steps,
setting up parameters for quantification, normalization, and quality control.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Dict

from metainformant.core import io, logging
from metainformant.rna.steps import StepResult

logger = logging.get_logger(__name__)


def run_step(step_params: Dict[str, Any]) -> StepResult:
    """Execute the configuration step.

    Args:
        step_params: Parameters for the config step
            - work_dir: Working directory
            - integrated_metadata: Path to integrated metadata file
            - species_list: List of species to configure
            - quantification_method: Method for quantification ("kallisto", "salmon", "star")
            - normalization_method: Method for normalization ("tmm", "cstmm")
            - threads: Number of threads to use

    Returns:
        StepResult with configuration results
    """
    start_time = time.time()

    try:
        work_dir = Path(step_params.get('work_dir', '.'))
        integrated_metadata_file = step_params.get('integrated_metadata',
                                                 work_dir / "integration" / "integrated_metadata.json")
        species_list = step_params.get('species_list', [])
        quantification_method = step_params.get('quantification_method', 'kallisto')
        normalization_method = step_params.get('normalization_method', 'cstmm')
        threads = step_params.get('threads', 1)

        # Load integrated metadata
        if isinstance(integrated_metadata_file, str):
            integrated_metadata_file = Path(integrated_metadata_file)

        if not integrated_metadata_file.exists():
            raise FileNotFoundError(f"Integrated metadata file not found: {integrated_metadata_file}")

        integrated_metadata = io.load_json(integrated_metadata_file)

        # Create config directory
        config_dir = work_dir / "config"
        config_dir.mkdir(parents=True, exist_ok=True)

        configurations = {}

        for species in species_list:
            if species not in integrated_metadata:
                logger.warning(f"Species {species} not found in integrated metadata")
                continue

            logger.info(f"Generating configuration for {species}")

            species_config = generate_species_config(
                species,
                integrated_metadata[species],
                quantification_method,
                normalization_method,
                threads
            )

            # Save species-specific config
            config_file = config_dir / f"{species.replace(' ', '_')}_config.json"
            io.dump_json(species_config, config_file)

            configurations[species] = species_config

        # Generate master configuration
        master_config = generate_master_config(configurations, work_dir, threads)
        master_config_file = config_dir / "master_config.json"
        io.dump_json(master_config, master_config_file)

        execution_time = time.time() - start_time

        return StepResult(
            success=True,
            outputs={
                'config_dir': str(config_dir),
                'master_config': str(master_config_file),
                'species_configs': {species: str(config_dir / f"{species.replace(' ', '_')}_config.json")
                                  for species in configurations.keys()},
                'configured_species': len(configurations)
            },
            metadata={
                'quantification_method': quantification_method,
                'normalization_method': normalization_method,
                'threads': threads,
                'execution_time': execution_time
            },
            execution_time=execution_time
        )

    except Exception as e:
        execution_time = time.time() - start_time
        logger.error(f"Config step failed: {e}")

        return StepResult(
            success=False,
            outputs={},
            metadata={},
            error_message=str(e),
            execution_time=execution_time
        )


def generate_species_config(species: str, species_metadata: Dict[str, Any],
                          quantification_method: str, normalization_method: str,
                          threads: int) -> Dict[str, Any]:
    """Generate configuration for a specific species.

    Args:
        species: Species name
        species_metadata: Metadata for species samples
        quantification_method: Quantification method
        normalization_method: Normalization method
        threads: Number of threads

    Returns:
        Species-specific configuration dictionary
    """
    config = {
        'species': species,
        'quantification': {
            'method': quantification_method,
            'threads': threads,
            'parameters': get_quantification_params(quantification_method)
        },
        'normalization': {
            'method': normalization_method,
            'parameters': get_normalization_params(normalization_method)
        },
        'samples': {},
        'reference_genome': get_reference_genome_config(species),
        'output_directories': {
            'quantification': f"quant/{species.replace(' ', '_')}",
            'normalized': f"normalized/{species.replace(' ', '_')}",
            'logs': f"logs/{species.replace(' ', '_')}"
        }
    }

    # Configure each sample
    for sample_id, sample_info in species_metadata.items():
        sample_config = {
            'accession': sample_info.get('accession', ''),
            'has_local_fastq': bool(sample_info.get('local_fastq_files')),
            'sra_available': sample_info.get('sra_available', False),
            'platform': sample_info.get('original_metadata', {}).get('platform'),
            'library_layout': sample_info.get('original_metadata', {}).get('library_layout'),
            'spots': sample_info.get('original_metadata', {}).get('spots', 0),
            'bases': sample_info.get('original_metadata', {}).get('bases', 0)
        }

        config['samples'][sample_id] = sample_config

    return config


def generate_master_config(species_configs: Dict[str, Dict[str, Any]],
                          work_dir: Path, threads: int) -> Dict[str, Any]:
    """Generate master configuration for all species.

    Args:
        species_configs: Individual species configurations
        work_dir: Working directory
        threads: Total threads available

    Returns:
        Master configuration dictionary
    """
    master_config = {
        'workflow_type': 'rna_seq_multi_species',
        'work_dir': str(work_dir),
        'total_species': len(species_configs),
        'total_samples': sum(len(config['samples']) for config in species_configs.values()),
        'thread_distribution': distribute_threads(species_configs, threads),
        'species_configs': list(species_configs.keys()),
        'pipeline_stages': [
            'metadata_retrieval',
            'data_download',
            'quantification',
            'normalization',
            'correlation_analysis',
            'quality_control'
        ],
        'output_structure': {
            'metadata': 'metadata/',
            'quantification': 'quant/',
            'normalized': 'normalized/',
            'correlations': 'correlations/',
            'logs': 'logs/',
            'reports': 'reports/'
        }
    }

    return master_config


def get_quantification_params(method: str) -> Dict[str, Any]:
    """Get default parameters for quantification method.

    Args:
        method: Quantification method name

    Returns:
        Parameter dictionary
    """
    params = {
        'kallisto': {
            'kmer_size': 31,
            'bootstrap_samples': 100,
            'fragment_length': 200,
            'fragment_sd': 20
        },
        'salmon': {
            'kmer_size': 31,
            'validate_mappings': True,
            'range_factorization_bins': 4,
            'mapping_cache_memory_limit': 5e8
        },
        'star': {
            'out_filter_multimap_nmax': 20,
            'out_filter_mismatch_nmax': 999,
            'out_filter_mismatch_nover_read_len': 0.04,
            'quant_mode': 'TranscriptomeSAM'
        }
    }

    return params.get(method, {})


def get_normalization_params(method: str) -> Dict[str, Any]:
    """Get default parameters for normalization method.

    Args:
        method: Normalization method name

    Returns:
        Parameter dictionary
    """
    params = {
        'tmm': {
            'trim_m': 0.3,
            'trim_a': 0.05,
            'bcv': 0.1
        },
        'cstmm': {
            'reference_species': None,
            'log_ratio_trim': 0.3,
            'sum_trim': 0.05,
            'do_weighting': True,
            'acutoff': -1e10
        },
        'quantile': {
            'ties': False
        }
    }

    return params.get(method, {})


def get_reference_genome_config(species: str) -> Dict[str, Any]:
    """Get reference genome configuration for a species.

    Args:
        species: Species name

    Returns:
        Reference genome configuration
    """
    # Default reference genome configurations
    default_configs = {
        'Drosophila melanogaster': {
            'source': 'NCBI',
            'accession': 'GCF_000001215.4',
            'annotation_source': 'Ensembl',
            'transcriptome_available': True
        },
        'Homo sapiens': {
            'source': 'NCBI',
            'accession': 'GCF_000001405.39',
            'annotation_source': 'GENCODE',
            'transcriptome_available': True
        },
        'Mus musculus': {
            'source': 'NCBI',
            'accession': 'GCF_000001635.27',
            'annotation_source': 'Ensembl',
            'transcriptome_available': True
        }
    }

    return default_configs.get(species, {
        'source': 'unknown',
        'accession': None,
        'annotation_source': 'unknown',
        'transcriptome_available': False
    })


def distribute_threads(species_configs: Dict[str, Dict[str, Any]], total_threads: int) -> Dict[str, int]:
    """Distribute threads across species based on sample count.

    Args:
        species_configs: Species configuration dictionaries
        total_threads: Total threads available

    Returns:
        Dictionary mapping species to thread count
    """
    if not species_configs:
        return {}

    # Calculate sample counts per species
    sample_counts = {species: len(config['samples']) for species, config in species_configs.items()}
    total_samples = sum(sample_counts.values())

    if total_samples == 0:
        # Equal distribution
        threads_per_species = total_threads // len(species_configs)
        return {species: threads_per_species for species in species_configs.keys()}

    # Distribute based on sample count
    thread_distribution = {}
    remaining_threads = total_threads

    for species, sample_count in sample_counts.items():
        # Proportional allocation
        species_threads = int(total_threads * (sample_count / total_samples))
        species_threads = max(1, min(species_threads, remaining_threads))  # At least 1, at most remaining

        thread_distribution[species] = species_threads
        remaining_threads -= species_threads

    # Distribute remaining threads to first species
    if remaining_threads > 0 and thread_distribution:
        first_species = next(iter(thread_distribution.keys()))
        thread_distribution[first_species] += remaining_threads

    return thread_distribution


def validate_configuration(config: Dict[str, Any]) -> Dict[str, Any]:
    """Validate a configuration dictionary.

    Args:
        config: Configuration to validate

    Returns:
        Validation results
    """
    validation = {
        'is_valid': True,
        'issues': [],
        'warnings': []
    }

    # Check required fields
    required_fields = ['species', 'quantification', 'normalization', 'samples']
    for field in required_fields:
        if field not in config:
            validation['issues'].append(f"Missing required field: {field}")
            validation['is_valid'] = False

    # Validate quantification method
    if 'quantification' in config:
        quant_config = config['quantification']
        valid_methods = ['kallisto', 'salmon', 'star', 'rsem']
        if quant_config.get('method') not in valid_methods:
            validation['issues'].append(f"Invalid quantification method: {quant_config.get('method')}")
            validation['is_valid'] = False

    # Validate normalization method
    if 'normalization' in config:
        norm_config = config['normalization']
        valid_methods = ['tmm', 'cstmm', 'quantile', 'upperquartile']
        if norm_config.get('method') not in valid_methods:
            validation['issues'].append(f"Invalid normalization method: {norm_config.get('method')}")
            validation['is_valid'] = False

    # Check sample configurations
    if 'samples' in config:
        samples = config['samples']
        if not samples:
            validation['warnings'].append("No samples configured")

        for sample_id, sample_config in samples.items():
            if not sample_config.get('accession'):
                validation['issues'].append(f"Sample {sample_id} missing accession")
                validation['is_valid'] = False

    return validation


