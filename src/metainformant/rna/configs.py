"""Configuration utilities for RNA-seq workflow.

This module provides utilities for loading, validating, and managing
RNA-seq workflow configurations.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from metainformant.core import io, logging

logger = logging.get_logger(__name__)


def load_workflow_config(config_path: str | Path) -> Dict[str, Any]:
    """Load RNA-seq workflow configuration from file.

    Args:
        config_path: Path to configuration file

    Returns:
        Configuration dictionary

    Raises:
        FileNotFoundError: If config file doesn't exist
        ValueError: If config is invalid
    """
    config_path = Path(config_path)

    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    config = io.load_json(config_path)

    # Validate configuration
    validate_config(config)

    logger.info(f"Loaded configuration from {config_path}")
    return config


def validate_config(config: Dict[str, Any]) -> None:
    """Validate RNA-seq workflow configuration.

    Args:
        config: Configuration dictionary to validate

    Raises:
        ValueError: If configuration is invalid
    """
    required_fields = ['species_list', 'work_dir']

    for field in required_fields:
        if field not in config:
            raise ValueError(f"Missing required configuration field: {field}")

    if not config['species_list']:
        raise ValueError("Species list cannot be empty")

    # Validate species list
    if not isinstance(config['species_list'], list):
        raise ValueError("Species list must be a list")

    logger.info("Configuration validation passed")


def create_default_config(species_list: list[str], work_dir: str | Path) -> Dict[str, Any]:
    """Create a default RNA-seq workflow configuration.

    Args:
        species_list: List of species for analysis
        work_dir: Working directory path

    Returns:
        Default configuration dictionary
    """
    config = {
        'species_list': species_list,
        'work_dir': str(work_dir),
        'threads': 4,
        'memory_gb': 16,
        'quantification_method': 'kallisto',
        'normalization_method': 'cstmm',
        'reference_genome_dir': str(Path(work_dir) / 'reference_genomes'),
        'output_dirs': {
            'metadata': 'metadata/',
            'config': 'config/',
            'selection': 'selection/',
            'fastq': 'fastq/',
            'quant': 'quant/',
            'merged': 'merged/',
            'cstmm': 'cstmm/',
            'curated': 'curated/',
            'csca': 'csca/',
            'logs': 'logs/',
            'reports': 'reports/'
        },
        'quality_thresholds': {
            'min_reads': 1000000,
            'min_mapping_rate': 0.8,
            'max_duplicates': 0.5
        },
        'amalgkit_params': {
            'metadata': {},
            'integrate': {},
            'config': {},
            'select': {},
            'getfastq': {'ENA': True},
            'quant': {'kallisto': True},
            'merge': {},
            'cstmm': {},
            'curate': {},
            'csca': {},
            'sanity': {}
        }
    }

    return config


def save_config(config: Dict[str, Any], output_path: str | Path) -> None:
    """Save configuration to file.

    Args:
        output_path: Output file path
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    io.dump_json(config, output_path)
    logger.info(f"Saved configuration to {output_path}")


def update_config_from_env(config: Dict[str, Any]) -> Dict[str, Any]:
    """Update configuration with environment variables.

    Args:
        config: Base configuration dictionary

    Returns:
        Updated configuration
    """
    import os

    # Update threads from environment
    if 'AK_THREADS' in os.environ:
        try:
            config['threads'] = int(os.environ['AK_THREADS'])
        except ValueError:
            logger.warning("Invalid AK_THREADS value, using default")

    # Update work directory from environment
    if 'AK_WORK_DIR' in os.environ:
        config['work_dir'] = os.environ['AK_WORK_DIR']

    return config


def get_config_template() -> Dict[str, Any]:
    """Get a configuration template with documentation.

    Returns:
        Configuration template with comments
    """
    template = {
        '_comment': 'RNA-seq workflow configuration template',
        'species_list': ['species1', 'species2'],  # List of species to analyze
        'work_dir': '/path/to/work/directory',     # Working directory
        'threads': 8,                             # Number of threads to use
        'quantification_method': 'kallisto',       # Quantification method
        'normalization_method': 'cstmm',          # Normalization method
        'reference_genome_dir': '/path/to/genomes', # Reference genome directory
        'quality_thresholds': {
            'min_reads': 1000000,     # Minimum reads per sample
            'min_mapping_rate': 0.8,  # Minimum mapping rate
            'max_duplicates': 0.5     # Maximum duplicate rate
        }
    }

    return template


class RNAPipelineConfig:
    """Configuration class for RNA-seq pipeline.

    Provides a structured configuration object for RNA-seq analysis workflows.
    """

    def __init__(self, work_dir: str | Path, threads: int = 8,
                 species_list: Optional[list[str]] = None,
                 quantification_method: str = "kallisto",
                 normalization_method: str = "cstmm",
                 reference_genome_dir: Optional[str | Path] = None):
        """Initialize RNA pipeline configuration.

        Args:
            work_dir: Working directory for the pipeline
            threads: Number of threads to use
            species_list: List of species to analyze
            quantification_method: Method for quantification ('kallisto', 'salmon', etc.)
            normalization_method: Method for normalization ('cstmm', 'tmm', etc.)
            reference_genome_dir: Directory containing reference genomes
        """
        self.work_dir = Path(work_dir)
        self.threads = threads
        self.species_list = species_list or []
        self.quantification_method = quantification_method
        self.normalization_method = normalization_method
        self.reference_genome_dir = Path(reference_genome_dir) if reference_genome_dir else None

        # Ensure work directory exists
        self.work_dir.mkdir(parents=True, exist_ok=True)

    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return {
            'work_dir': str(self.work_dir),
            'threads': self.threads,
            'species_list': self.species_list,
            'quantification_method': self.quantification_method,
            'normalization_method': self.normalization_method,
            'reference_genome_dir': str(self.reference_genome_dir) if self.reference_genome_dir else None,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> RNAPipelineConfig:
        """Create configuration from dictionary."""
        return cls(**data)

    def __str__(self) -> str:
        return f"RNAPipelineConfig(work_dir={self.work_dir}, threads={self.threads}, species={len(self.species_list)})"



