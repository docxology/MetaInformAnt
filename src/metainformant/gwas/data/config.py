"""Configuration management for GWAS workflows.

This module provides tools for creating, validating, and managing
GWAS analysis configurations.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

from metainformant.core.utils import logging
from metainformant.core import io

logger = logging.get_logger(__name__)


class AttrDict(dict):
    """Dictionary that allows attribute-style access.

    Nested dicts are automatically converted to AttrDict on access,
    allowing chained attribute access like ``cfg.genome.accession``.
    Missing keys return ``None`` instead of raising ``KeyError``.
    """

    def __getattr__(self, name: str) -> Any:
        try:
            value = self[name]
            # Convert nested dicts to AttrDict on access
            if isinstance(value, dict) and not isinstance(value, AttrDict):
                value = AttrDict(value)
                self[name] = value
            return value
        except KeyError:
            return None

    def __setattr__(self, name: str, value: Any) -> None:
        self[name] = value

    def __delattr__(self, name: str) -> None:
        try:
            del self[name]
        except KeyError:
            raise AttributeError(name)


def _to_attrdict(obj: Any) -> Any:
    """Recursively convert dicts to AttrDict."""
    if isinstance(obj, dict) and not isinstance(obj, AttrDict):
        return AttrDict({k: _to_attrdict(v) for k, v in obj.items()})
    elif isinstance(obj, list):
        return [_to_attrdict(item) for item in obj]
    return obj


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

    Accepts both flat format (``vcf_file``, ``phenotype_file``, ``output_dir``)
    and nested YAML format (``variants.vcf_files``, ``samples.phenotype_file``,
    ``work_dir`` / ``output.results_dir``).

    File existence checks are skipped since paths may not exist yet at config
    validation time.

    Args:
        config: Configuration dictionary (flat or nested)

    Returns:
        List of validation error messages (empty if valid)
    """
    errors = []

    # ---- Required parameters (accept both flat and nested alternatives) ----

    # VCF requirement: vcf_file OR variants.vcf_files OR variants.calling OR variants section
    has_vcf = "vcf_file" in config
    if not has_vcf:
        variants = config.get("variants", {})
        if isinstance(variants, dict):
            vcf_files = variants.get("vcf_files")
            if vcf_files and isinstance(vcf_files, list) and len(vcf_files) > 0:
                has_vcf = True
            # Also accept variant calling config as an alternative source
            calling = variants.get("calling", {})
            if isinstance(calling, dict) and calling:
                has_vcf = True
            # Accept if variants section exists at all (files may be produced later)
            if variants:
                has_vcf = True
    if not has_vcf:
        # Also accept vcf_path (used by workflow module)
        if "vcf_path" in config:
            has_vcf = True
    if not has_vcf:
        errors.append("Missing required parameter: vcf_file (or variants.vcf_files)")

    # Phenotype requirement: phenotype_file OR samples.phenotype_file
    has_pheno = "phenotype_file" in config
    if not has_pheno:
        samples = config.get("samples", {})
        if isinstance(samples, dict) and samples.get("phenotype_file"):
            has_pheno = True
    if not has_pheno:
        if "phenotype_path" in config:
            has_pheno = True
    if not has_pheno:
        errors.append("Missing required parameter: phenotype_file (or samples.phenotype_file)")

    # Output dir requirement: output_dir OR work_dir OR output.results_dir
    has_output = "output_dir" in config or "work_dir" in config
    if not has_output:
        output = config.get("output", {})
        if isinstance(output, dict) and output.get("results_dir"):
            has_output = True
    if not has_output:
        errors.append("Missing required parameter: output_dir (or work_dir or output.results_dir)")

    # NOTE: File existence checks are intentionally skipped here.
    # Paths may not exist yet when loading config before workflow execution.
    # Use validate_input_files() for runtime file existence validation.

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
    valid_methods = ["linear", "logistic", "mixed"]
    method = config.get("association_method")
    assoc_dict = config.get("association")
    if not method and isinstance(assoc_dict, dict):
        method = assoc_dict.get("model")
    if method and method not in valid_methods:
        errors.append(f"Invalid association method '{method}'. Must be one of: {valid_methods}")

    valid_corrections = ["bonferroni", "fdr", "genomic_control"]
    correction = config.get("multiple_testing_correction")
    corr_dict = config.get("correction")
    if not correction and isinstance(corr_dict, dict):
        correction = corr_dict.get("method")
    if correction and correction not in valid_corrections:
        errors.append(f"Invalid correction method '{correction}'. Must be one of: {valid_corrections}")

    valid_kinship = ["vanraden", "ibs", "astle", "yang"]
    kinship = config.get("kinship_method")
    struct_dict = config.get("structure")
    if not kinship and isinstance(struct_dict, dict):
        kinship = struct_dict.get("kinship_method")
    if kinship and kinship not in valid_kinship:
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


def load_gwas_config(config_path: str | Path) -> AttrDict:
    """Load and validate GWAS configuration.

    Returns an :class:`AttrDict` so that nested sections can be accessed via
    attribute syntax (e.g. ``cfg.genome.accession``).

    Args:
        config_path: Path to configuration file

    Returns:
        Validated configuration as AttrDict (supports both dict and attribute access)

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

    # Convert to AttrDict recursively for attribute-style access
    return _to_attrdict(config)


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


def estimate_runtime(
    config: Dict[str, Any],
    *,
    n_samples: Optional[int] = None,
    n_variants: Optional[int] = None,
) -> Dict[str, Any]:
    """Estimate GWAS runtime based on configuration and data dimensions.

    Uses data-driven scaling models calibrated from empirical GWAS runs.
    When ``n_samples`` and ``n_variants`` are provided, estimates are scaled
    using known computational complexity for each pipeline step (e.g.
    O(n²·m) for kinship, O(n·m) for association testing).

    When dimensions are not provided, falls back to empirical baselines
    calibrated from an Apis mellifera 188-sample, 3000-variant pilot run
    (~111s for association testing at ~4000 tests/sec).

    Args:
        config: Configuration dictionary.
        n_samples: Number of samples in the dataset (optional).
        n_variants: Number of variants in the dataset (optional).

    Returns:
        Runtime estimate dictionary with keys:
            - estimated_hours: Total estimated hours
            - estimated_seconds: Total estimated seconds
            - per_step_seconds: Per-step breakdown in seconds
            - estimated_cores: Thread count
            - estimated_memory_gb: Memory estimate
            - limiting_factors: List of bottleneck descriptions
            - data_driven: Whether real data dimensions were used
    """
    threads = config.get("threads", 4)
    pca_comp = config.get("pca_components", 10)

    # --- Empirical baselines (seconds) calibrated from Apis mellifera pilot ---
    # Pilot: 188 samples, 3000 variants, 111s total for 293,450 tests
    PILOT_N = 188
    PILOT_M = 3000

    # Per-step pilot timings (seconds) — empirical
    pilot_steps = {
        "parse_vcf":          5.0,    # O(n·m)
        "qc_filters":         5.0,    # O(n·m)
        "ld_pruning":         5.0,    # O(w²·m)
        "population_structure": 15.0, # O(n²·m) dominated by kinship
        "association_testing": 111.0, # O(n·m) — 4000 tests/sec
        "multiple_testing":   2.0,    # O(m)
        "summary_stats":      2.0,    # O(m)
        "annotation":         3.0,    # O(m)
    }

    if config.get("create_plots", True):
        pilot_steps["visualization"] = 10.0  # O(m)

    # Complexity model per step
    step_models = {
        "parse_vcf":          "n_m",
        "qc_filters":         "n_m",
        "ld_pruning":         "m",
        "population_structure": "n2_m",
        "association_testing": "n_m",
        "multiple_testing":   "m",
        "summary_stats":      "m",
        "annotation":         "m",
        "visualization":      "m",
    }

    data_driven = n_samples is not None and n_variants is not None
    target_n = n_samples if n_samples else PILOT_N
    target_m = n_variants if n_variants else PILOT_M

    per_step: Dict[str, float] = {}
    limiting_factors: List[str] = []

    for step_name, pilot_secs in pilot_steps.items():
        model = step_models.get(step_name, "n_m")
        factor = _compute_scaling_factor(model, PILOT_N, PILOT_M, target_n, target_m, pca_comp)
        per_step[step_name] = pilot_secs * factor

    # PCA component scaling for population_structure
    if pca_comp > 10:
        per_step["population_structure"] *= (pca_comp / 10) ** 2

    total_seconds = sum(per_step.values())

    # Apply Amdahl's law for thread parallelism
    if threads > 1:
        serial_fraction = 0.2  # ~20% serial (I/O, setup)
        speedup = 1 / (serial_fraction + (1 - serial_fraction) / threads)
        total_seconds /= speedup
        per_step = {k: v / speedup for k, v in per_step.items()}

    # Identify limiting factors
    if data_driven:
        if target_n > 1000:
            limiting_factors.append(
                f"Kinship matrix O(n²): {target_n}² samples dominates structure step"
            )
        if target_m > 1_000_000:
            limiting_factors.append(
                f"Association testing O(n·m): {target_n}×{target_m} = {target_n * target_m:,} tests"
            )

    # Memory estimate: GRM is n×n floats (8 bytes each)
    grm_memory_gb = (target_n ** 2 * 8) / (1024 ** 3) if data_driven else 0.5
    estimated_memory = max(config.get("memory_gb", 16), grm_memory_gb * 2)

    return {
        "estimated_hours": total_seconds / 3600,
        "estimated_seconds": total_seconds,
        "per_step_seconds": per_step,
        "estimated_cores": threads,
        "estimated_memory_gb": round(estimated_memory, 1),
        "limiting_factors": limiting_factors,
        "data_driven": data_driven,
    }


def _compute_scaling_factor(
    model: str,
    pilot_n: int,
    pilot_m: int,
    target_n: int,
    target_m: int,
    k: int = 10,
) -> float:
    """Compute scaling factor for a given complexity model.

    Args:
        model: Complexity model ("n_m", "n2_m", "m", "m_k2").
        pilot_n: Pilot sample count.
        pilot_m: Pilot variant count.
        target_n: Target sample count.
        target_m: Target variant count.
        k: PCA components or region size.

    Returns:
        Multiplicative scaling factor.
    """
    if pilot_n <= 0 or pilot_m <= 0:
        return 1.0

    if model == "n_m":
        return (target_n * target_m) / (pilot_n * pilot_m)
    elif model == "n2_m":
        return (target_n ** 2 * target_m) / (pilot_n ** 2 * pilot_m)
    elif model == "m":
        return target_m / pilot_m
    elif model == "m_k2":
        return (target_m * k ** 2) / (pilot_m * 100)  # pilot k=10
    else:
        return (target_n * target_m) / (pilot_n * pilot_m)
