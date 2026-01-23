"""Configuration management for GWAS workflows.

This module provides tools for creating, validating, and managing
GWAS analysis configurations.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

from metainformant.core import io, logging

logger = logging.get_logger(__name__)


def create_gwas_config_template(output_path: str | Path) -> None:
    """Create a GWAS configuration template file.

    Args:
        output_path: Path to save the template
    """
    template = {
        # Study information
        "study_name": "My GWAS Study",
        "description": "Genome-wide association study description",
        # Input data
        "vcf_file": "/path/to/genotypes.vcf.gz",
        "phenotype_file": "/path/to/phenotypes.txt",
        "covariate_file": "/optional/path/to/covariates.txt",
        # Reference data
        "reference_genome": "GCF_000001405.39",  # GRCh38
        "population": "EUR",  # Population for PCA
        # Analysis parameters
        "maf_threshold": 0.01,  # Minor allele frequency filter
        "missing_threshold": 0.05,  # Missing data threshold
        "hwe_p_threshold": 1e-6,  # Hardy-Weinberg equilibrium p-value
        "significance_threshold": 5e-8,  # GWAS significance threshold
        # Quality control
        "qc_filters": {
            "remove_related_samples": True,
            "remove_heterozygosity_outliers": True,
            "remove_missingness_outliers": True,
        },
        # Population structure
        "pca_components": 10,  # Number of PCA components
        "kinship_method": "vanraden",  # Kinship matrix method
        # Association testing
        "association_method": "linear",  # 'linear' or 'logistic'
        "multiple_testing_correction": "bonferroni",  # 'bonferroni', 'fdr', 'genomic_control'
        # Visualization
        "create_plots": True,
        "plot_types": ["manhattan", "qq", "pca", "regional"],
        "output_format": "png",  # 'png', 'pdf', 'svg'
        # Performance
        "threads": 4,
        "memory_gb": 16,
        # Output
        "output_dir": "output/gwas/",
        "intermediate_files": True,  # Keep intermediate files
        # Advanced options
        "advanced": {
            "imputation_method": None,  # For missing data imputation
            "conditional_analysis": False,  # Conditional association testing
            "fine_mapping": False,  # Fine-mapping analysis
            "meta_analysis": False,  # Meta-analysis across studies
        },
    }

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    io.dump_json(template, output_path)
    logger.info(f"Created GWAS config template: {output_path}")


def validate_config_parameters(config: Dict[str, Any]) -> List[str]:
    """Validate GWAS configuration parameters.

    Args:
        config: Configuration dictionary

    Returns:
        List of validation error messages (empty if valid)
    """
    errors = []

    # Required parameters
    required_params = ["vcf_file", "phenotype_file", "output_dir"]
    for param in required_params:
        if param not in config:
            errors.append(f"Missing required parameter: {param}")

    # File paths
    if "vcf_file" in config:
        vcf_path = Path(config["vcf_file"])
        if not vcf_path.exists():
            errors.append(f"VCF file does not exist: {vcf_path}")

    if "phenotype_file" in config:
        pheno_path = Path(config["phenotype_file"])
        if not pheno_path.exists():
            errors.append(f"Phenotype file does not exist: {pheno_path}")

    # Numerical parameters
    if "maf_threshold" in config:
        maf = config["maf_threshold"]
        if not (0 < maf < 0.5):
            errors.append(f"MAF threshold must be between 0 and 0.5: {maf}")

    if "missing_threshold" in config:
        missing = config["missing_threshold"]
        if not (0 <= missing <= 1):
            errors.append(f"Missing threshold must be between 0 and 1: {missing}")

    if "significance_threshold" in config:
        sig = config["significance_threshold"]
        if sig <= 0:
            errors.append(f"Significance threshold must be positive: {sig}")

    # Choice parameters
    valid_methods = ["linear", "logistic"]
    if "association_method" in config:
        method = config["association_method"]
        if method not in valid_methods:
            errors.append(f"Invalid association method '{method}'. Must be one of: {valid_methods}")

    valid_corrections = ["bonferroni", "fdr", "genomic_control"]
    if "multiple_testing_correction" in config:
        correction = config["multiple_testing_correction"]
        if correction not in valid_corrections:
            errors.append(f"Invalid correction method '{correction}'. Must be one of: {valid_corrections}")

    valid_kinship = ["vanraden", "centered", "normalized"]
    if "kinship_method" in config:
        kinship = config["kinship_method"]
        if kinship not in valid_kinship:
            errors.append(f"Invalid kinship method '{kinship}'. Must be one of: {valid_kinship}")

    # Threads and memory
    if "threads" in config:
        threads = config["threads"]
        if not isinstance(threads, int) or threads < 1:
            errors.append(f"Threads must be a positive integer: {threads}")

    if "memory_gb" in config:
        memory = config["memory_gb"]
        if not isinstance(memory, (int, float)) or memory <= 0:
            errors.append(f"Memory must be a positive number: {memory}")

    return errors


def merge_config_defaults(config: Dict[str, Any]) -> Dict[str, Any]:
    """Merge configuration with default values.

    Args:
        config: User configuration

    Returns:
        Configuration with defaults filled in
    """
    defaults = {
        "maf_threshold": 0.01,
        "missing_threshold": 0.05,
        "hwe_p_threshold": 1e-6,
        "significance_threshold": 5e-8,
        "pca_components": 10,
        "kinship_method": "vanraden",
        "association_method": "linear",
        "multiple_testing_correction": "bonferroni",
        "create_plots": True,
        "plot_types": ["manhattan", "qq", "pca"],
        "output_format": "png",
        "threads": 4,
        "memory_gb": 16,
        "intermediate_files": True,
        "population": "EUR",
    }

    # Deep merge
    merged = defaults.copy()
    merged.update(config)

    return merged


def load_gwas_config(config_path: str | Path) -> Dict[str, Any]:
    """Load and validate GWAS configuration.

    Args:
        config_path: Path to configuration file

    Returns:
        Validated configuration dictionary

    Raises:
        FileNotFoundError: If config file doesn't exist
        ValueError: If configuration is invalid
    """
    config_path = Path(config_path)

    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    if config_path.suffix.lower() in [".yaml", ".yml"]:
        config = io.load_yaml(config_path)
    elif config_path.suffix.lower() == ".toml":
        config = io.load_toml(config_path)
    else:
        config = io.load_json(config_path)

    # Merge with defaults
    config = merge_config_defaults(config)

    # Validate
    errors = validate_config_parameters(config)
    if errors:
        raise ValueError(f"Configuration validation failed:\n" + "\n".join(f"- {error}" for error in errors))

    logger.info(f"Loaded and validated GWAS configuration from {config_path}")
    return config


def save_gwas_config(config: Dict[str, Any], output_path: str | Path) -> None:
    """Save GWAS configuration to file.

    Args:
        config: Configuration dictionary
        output_path: Output file path
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    io.dump_json(config, output_path)
    logger.info(f"Saved GWAS configuration to {output_path}")


def create_config_from_vcf(vcf_path: str | Path, phenotype_path: str | Path, output_dir: str | Path) -> Dict[str, Any]:
    """Create a basic GWAS configuration from VCF and phenotype files.

    Args:
        vcf_path: Path to VCF file
        phenotype_path: Path to phenotype file
        output_dir: Output directory

    Returns:
        Basic GWAS configuration
    """
    vcf_path = Path(vcf_path)
    phenotype_path = Path(phenotype_path)
    output_dir = Path(output_dir)

    config = {
        "vcf_file": str(vcf_path),
        "phenotype_file": str(phenotype_path),
        "output_dir": str(output_dir),
        # Basic parameters
        "study_name": f"GWAS_{vcf_path.stem}",
        "maf_threshold": 0.01,
        "missing_threshold": 0.05,
        "significance_threshold": 5e-8,
        "threads": 4,
        # Enable all standard analyses
        "create_plots": True,
        "plot_types": ["manhattan", "qq", "pca", "regional"],
        "intermediate_files": True,
    }

    return config


def update_config_for_runtime(config: Dict[str, Any]) -> Dict[str, Any]:
    """Update configuration with runtime-specific settings.

    Args:
        config: Base configuration

    Returns:
        Updated configuration
    """
    import os

    # Update from environment variables
    env_mappings = {"GWAS_THREADS": "threads", "GWAS_MEMORY": "memory_gb", "GWAS_OUTPUT_DIR": "output_dir"}

    for env_var, config_key in env_mappings.items():
        if env_var in os.environ:
            try:
                if config_key in ["threads", "memory_gb"]:
                    config[config_key] = int(os.environ[env_var])
                else:
                    config[config_key] = os.environ[env_var]
                logger.info(f"Updated {config_key} from environment: {config[config_key]}")
            except ValueError:
                logger.warning(f"Invalid value for {env_var}: {os.environ[env_var]}")

    return config


def get_config_summary(config: Dict[str, Any]) -> Dict[str, Any]:
    """Get a summary of configuration settings.

    Args:
        config: Configuration dictionary

    Returns:
        Configuration summary
    """
    summary = {
        "input_files": {
            "vcf": config.get("vcf_file"),
            "phenotype": config.get("phenotype_file"),
            "covariates": config.get("covariate_file"),
        },
        "quality_filters": {
            "maf_threshold": config.get("maf_threshold"),
            "missing_threshold": config.get("missing_threshold"),
            "hwe_threshold": config.get("hwe_p_threshold"),
        },
        "analysis_parameters": {
            "association_method": config.get("association_method"),
            "correction_method": config.get("multiple_testing_correction"),
            "significance_threshold": config.get("significance_threshold"),
            "pca_components": config.get("pca_components"),
        },
        "performance": {"threads": config.get("threads"), "memory_gb": config.get("memory_gb")},
        "output": {
            "directory": config.get("output_dir"),
            "plots": config.get("create_plots"),
            "plot_types": config.get("plot_types"),
            "format": config.get("output_format"),
        },
    }

    return summary


def validate_input_files(config: Dict[str, Any]) -> List[str]:
    """Validate that input files exist and are readable.

    Args:
        config: Configuration dictionary

    Returns:
        List of file validation errors
    """
    errors = []

    # Check VCF file
    if "vcf_file" in config:
        vcf_path = Path(config["vcf_file"])
        if not vcf_path.exists():
            errors.append(f"VCF file does not exist: {vcf_path}")
        elif not vcf_path.is_file():
            errors.append(f"VCF path is not a file: {vcf_path}")

    # Check phenotype file
    if "phenotype_file" in config:
        pheno_path = Path(config["phenotype_file"])
        if not pheno_path.exists():
            errors.append(f"Phenotype file does not exist: {pheno_path}")
        elif not pheno_path.is_file():
            errors.append(f"Phenotype path is not a file: {pheno_path}")

    # Check covariate file (optional)
    if "covariate_file" in config:
        cov_path = Path(config["covariate_file"])
        if cov_path and not cov_path.exists():
            errors.append(f"Covariate file does not exist: {cov_path}")

    # Check output directory
    if "output_dir" in config:
        output_dir = Path(config["output_dir"])
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            errors.append(f"Cannot create output directory {output_dir}: {e}")

    return errors


def estimate_runtime(config: Dict[str, Any]) -> Dict[str, Any]:
    """Estimate GWAS runtime based on configuration and data size.

    Args:
        config: Configuration dictionary

    Returns:
        Runtime estimate information
    """
    estimate = {
        "estimated_hours": 0.0,
        "estimated_cores": config.get("threads", 4),
        "estimated_memory_gb": config.get("memory_gb", 16),
        "limiting_factors": [],
    }

    # Basic estimation (simplified)
    # In practice, this would analyze actual file sizes and complexity

    # Quality control time
    estimate["estimated_hours"] += 0.5  # 30 minutes

    # PCA time
    pca_comp = config.get("pca_components", 10)
    estimate["estimated_hours"] += pca_comp * 0.1  # 6 minutes per component

    # Association testing time
    estimate["estimated_hours"] += 2.0  # 2 hours for typical GWAS

    # Multiple testing correction
    estimate["estimated_hours"] += 0.2  # 12 minutes

    # Plotting time
    if config.get("create_plots", True):
        estimate["estimated_hours"] += 1.0  # 1 hour for comprehensive plots

    # Scale by thread count (Amdahl's law approximation)
    threads = config.get("threads", 4)
    if threads > 1:
        serial_fraction = 0.2  # Assume 20% serial time
        speedup = 1 / (serial_fraction + (1 - serial_fraction) / threads)
        estimate["estimated_hours"] /= speedup

    return estimate
