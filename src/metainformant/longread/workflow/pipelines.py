"""Pre-defined pipeline configurations for long-read analysis workflows.

Provides factory functions for generating pipeline configuration dictionaries
for QC, assembly, methylation, and structural variant calling pipelines.
Supports loading from YAML/JSON files and validation of configuration values.

Each configuration dictionary follows a consistent structure with:
- pipeline_name: Identifier for the pipeline
- steps: Ordered list of step configurations
- parameters: Pipeline-wide parameter defaults
- output: Output directory and format preferences
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Valid pipeline names
VALID_PIPELINE_NAMES = frozenset({"qc", "assembly", "methylation", "sv", "full"})

# Parameter bounds for validation
_PARAM_BOUNDS: dict[str, dict[str, tuple[float | int, float | int]]] = {
    "qc": {
        "min_length": (0, 10_000_000),
        "min_quality": (0.0, 60.0),
        "max_search_length": (10, 10_000),
        "min_adapter_identity": (0.5, 1.0),
    },
    "assembly": {
        "min_overlap": (100, 1_000_000),
        "k": (5, 31),
        "w": (3, 50),
        "polish_iterations": (0, 20),
        "min_minimizer_matches": (1, 100),
        "max_overhang": (100, 100_000),
    },
    "methylation": {
        "min_coverage": (1, 1000),
        "significance_threshold": (0.0, 1.0),
        "min_difference": (0.0, 1.0),
        "methylation_threshold": (0.0, 1.0),
    },
    "sv": {
        "min_sv_size": (10, 10_000_000),
        "min_support": (1, 1000),
        "min_mapping_quality": (0, 60),
        "max_cluster_distance": (10, 100_000),
    },
}


def get_qc_pipeline_config(
    min_length: int = 1000,
    min_quality: float = 7.0,
    trim_adapters: bool = True,
    max_search_length: int = 200,
    min_adapter_identity: float = 0.75,
) -> dict[str, Any]:
    """Get QC pipeline configuration.

    The QC pipeline performs: read loading, length filtering, quality filtering,
    adapter trimming, metric computation, and QC plot generation.

    Args:
        min_length: Minimum read length to keep (bp).
        min_quality: Minimum mean Phred quality score.
        trim_adapters: Whether to perform adapter trimming.
        max_search_length: Maximum bases to search for adapters at each end.
        min_adapter_identity: Minimum identity for adapter matching.

    Returns:
        Configuration dictionary for the QC pipeline.
    """
    return {
        "pipeline_name": "qc",
        "parameters": {
            "min_length": min_length,
            "min_quality": min_quality,
            "trim_adapters": trim_adapters,
            "max_search_length": max_search_length,
            "min_adapter_identity": min_adapter_identity,
        },
        "steps": [
            {
                "name": "read_input",
                "description": "Load reads from input source",
                "depends_on": [],
            },
            {
                "name": "filter_length",
                "description": "Filter reads by minimum length",
                "depends_on": ["read_input"],
                "params": {"min_length": min_length},
            },
            {
                "name": "filter_quality",
                "description": "Filter reads by mean quality score",
                "depends_on": ["filter_length"],
                "params": {"min_q": min_quality},
            },
            {
                "name": "trim_adapters",
                "description": "Detect and trim adapter sequences",
                "depends_on": ["filter_quality"],
                "enabled": trim_adapters,
                "params": {
                    "max_search_length": max_search_length,
                    "min_identity": min_adapter_identity,
                },
            },
            {
                "name": "compute_metrics",
                "description": "Calculate read length and quality statistics",
                "depends_on": ["trim_adapters"] if trim_adapters else ["filter_quality"],
            },
            {
                "name": "generate_qc_plots",
                "description": "Generate QC visualization plots",
                "depends_on": ["compute_metrics"],
            },
            {
                "name": "export_report",
                "description": "Export QC report to file",
                "depends_on": ["compute_metrics"],
            },
        ],
        "output": {
            "format": "json",
            "include_plots": True,
        },
    }


def get_assembly_pipeline_config(
    min_overlap: int = 2000,
    k: int = 15,
    w: int = 10,
    polish_iterations: int = 2,
    min_minimizer_matches: int = 3,
    max_overhang: int = 1000,
    min_read_length: int = 1000,
    min_read_quality: float = 7.0,
) -> dict[str, Any]:
    """Get assembly pipeline configuration.

    The assembly pipeline performs: read filtering, overlap detection using
    minimizer sketching, overlap graph construction, consensus generation,
    iterative polishing, and assembly statistics computation.

    Args:
        min_overlap: Minimum overlap length between reads (bp).
        k: K-mer size for minimizer computation.
        w: Window size for minimizer computation.
        polish_iterations: Number of consensus polishing iterations.
        min_minimizer_matches: Minimum shared minimizers for candidate pairs.
        max_overhang: Maximum unaligned overhang at overlap ends.
        min_read_length: Minimum read length for assembly input.
        min_read_quality: Minimum mean quality for assembly input.

    Returns:
        Configuration dictionary for the assembly pipeline.
    """
    return {
        "pipeline_name": "assembly",
        "parameters": {
            "min_overlap": min_overlap,
            "k": k,
            "w": w,
            "polish_iterations": polish_iterations,
            "min_minimizer_matches": min_minimizer_matches,
            "max_overhang": max_overhang,
            "min_read_length": min_read_length,
            "min_read_quality": min_read_quality,
        },
        "steps": [
            {
                "name": "filter_reads",
                "description": "Filter reads by length and quality for assembly",
                "depends_on": [],
                "params": {
                    "min_length": min_read_length,
                    "min_q": min_read_quality,
                },
            },
            {
                "name": "find_overlaps",
                "description": "Find pairwise overlaps using minimizer sketching",
                "depends_on": ["filter_reads"],
                "params": {
                    "min_overlap": min_overlap,
                    "k": k,
                    "w": w,
                    "min_minimizer_matches": min_minimizer_matches,
                    "max_overhang": max_overhang,
                },
            },
            {
                "name": "build_graph",
                "description": "Build overlap graph and remove contained reads",
                "depends_on": ["find_overlaps"],
            },
            {
                "name": "generate_consensus",
                "description": "Generate consensus sequence from overlapping reads",
                "depends_on": ["build_graph"],
            },
            {
                "name": "polish",
                "description": "Polish consensus sequence by iterative re-alignment",
                "depends_on": ["generate_consensus"],
                "params": {"iterations": polish_iterations},
            },
            {
                "name": "calculate_stats",
                "description": "Calculate assembly quality statistics",
                "depends_on": ["polish"],
            },
            {
                "name": "assembly_plots",
                "description": "Generate assembly quality plots",
                "depends_on": ["calculate_stats"],
            },
        ],
        "output": {
            "format": "json",
            "include_plots": True,
        },
    }


def get_methylation_pipeline_config(
    modification_types: list[str] | None = None,
    min_coverage: int = 5,
    significance_threshold: float = 0.05,
    min_difference: float = 0.2,
    methylation_threshold: float = 0.5,
    regions: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Get methylation pipeline configuration.

    The methylation pipeline performs: read loading, 5mC calling, 6mA calling,
    per-region aggregation, differential methylation analysis, and visualization.

    Args:
        modification_types: Types of modifications to detect. Defaults to
            ["5mC", "6mA"].
        min_coverage: Minimum read coverage for aggregation.
        significance_threshold: P-value threshold for differential analysis.
        min_difference: Minimum absolute methylation difference to report.
        methylation_threshold: Probability threshold for calling a site modified.
        regions: Optional list of genomic region dictionaries for aggregation.

    Returns:
        Configuration dictionary for the methylation pipeline.
    """
    if modification_types is None:
        modification_types = ["5mC", "6mA"]

    return {
        "pipeline_name": "methylation",
        "parameters": {
            "modification_types": modification_types,
            "min_coverage": min_coverage,
            "significance_threshold": significance_threshold,
            "min_difference": min_difference,
            "methylation_threshold": methylation_threshold,
            "regions": regions or [],
        },
        "steps": [
            {
                "name": "load_reads",
                "description": "Load reads from input source",
                "depends_on": [],
            },
            {
                "name": "call_5mc",
                "description": "Call 5-methylcytosine modifications",
                "depends_on": ["load_reads"],
                "enabled": "5mC" in modification_types,
                "params": {"threshold": methylation_threshold},
            },
            {
                "name": "call_6ma",
                "description": "Call N6-methyladenine modifications",
                "depends_on": ["load_reads"],
                "enabled": "6mA" in modification_types,
                "params": {"threshold": methylation_threshold},
            },
            {
                "name": "aggregate_regions",
                "description": "Aggregate methylation by genomic region",
                "depends_on": ["call_5mc", "call_6ma"],
                "params": {"min_coverage": min_coverage},
            },
            {
                "name": "differential_analysis",
                "description": "Differential methylation analysis between groups",
                "depends_on": ["aggregate_regions"],
                "params": {
                    "min_coverage": min_coverage,
                    "min_difference": min_difference,
                    "alpha": significance_threshold,
                },
            },
            {
                "name": "methylation_plots",
                "description": "Generate methylation visualization plots",
                "depends_on": ["aggregate_regions"],
            },
            {
                "name": "methylation_report",
                "description": "Export methylation analysis report",
                "depends_on": ["differential_analysis", "methylation_plots"],
            },
        ],
        "output": {
            "format": "json",
            "include_plots": True,
        },
    }


def get_sv_pipeline_config(
    min_sv_size: int = 50,
    min_support: int = 3,
    sv_types: list[str] | None = None,
    min_mapping_quality: int = 20,
    max_cluster_distance: int = 500,
) -> dict[str, Any]:
    """Get SV calling pipeline configuration.

    The SV pipeline performs: alignment loading, structural variant detection
    (deletions, insertions, inversions, duplications, translocations), variant
    phasing, annotation, and visualization.

    Args:
        min_sv_size: Minimum structural variant size (bp).
        min_support: Minimum supporting read count.
        sv_types: Types of SVs to detect. Defaults to all types.
        min_mapping_quality: Minimum mapping quality for supporting alignments.
        max_cluster_distance: Maximum distance between breakpoints to cluster.

    Returns:
        Configuration dictionary for the SV pipeline.
    """
    if sv_types is None:
        sv_types = ["DEL", "INS", "INV", "DUP", "BND"]

    return {
        "pipeline_name": "sv",
        "parameters": {
            "min_sv_size": min_sv_size,
            "min_support": min_support,
            "sv_types": sv_types,
            "min_mapping_quality": min_mapping_quality,
            "max_cluster_distance": max_cluster_distance,
        },
        "steps": [
            {
                "name": "load_alignments",
                "description": "Load long-read alignments from BAM",
                "depends_on": [],
            },
            {
                "name": "detect_svs",
                "description": "Detect all structural variants from alignments",
                "depends_on": ["load_alignments"],
                "params": {
                    "min_size": min_sv_size,
                    "min_support": min_support,
                    "min_mapping_quality": min_mapping_quality,
                },
            },
            {
                "name": "detect_insertions",
                "description": "Specialized insertion detection",
                "depends_on": ["load_alignments"],
                "enabled": "INS" in sv_types,
                "params": {"min_size": min_sv_size, "min_support": min_support},
            },
            {
                "name": "detect_inversions",
                "description": "Specialized inversion detection",
                "depends_on": ["load_alignments"],
                "enabled": "INV" in sv_types,
                "params": {"min_size": min_sv_size, "min_support": min_support},
            },
            {
                "name": "phase_variants",
                "description": "Phase structural variants to haplotypes",
                "depends_on": ["detect_svs"],
            },
            {
                "name": "sv_summary",
                "description": "Generate SV type and size distribution summary",
                "depends_on": ["detect_svs", "detect_insertions", "detect_inversions"],
            },
            {
                "name": "sv_plots",
                "description": "Generate SV visualization plots",
                "depends_on": ["sv_summary"],
            },
        ],
        "output": {
            "format": "json",
            "include_plots": True,
        },
    }


def load_pipeline_config(config_path: Path | str) -> dict[str, Any]:
    """Load pipeline configuration from a YAML or JSON file.

    Supports YAML (.yaml, .yml) and JSON (.json) formats. The file should
    contain a top-level dictionary with at minimum a 'pipeline_name' key.

    Args:
        config_path: Path to the configuration file.

    Returns:
        Pipeline configuration dictionary.

    Raises:
        FileNotFoundError: If the config file does not exist.
        ValueError: If the file format is unsupported or content is invalid.
    """
    from metainformant.core.utils.config import load_mapping_from_file

    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Pipeline config file not found: {config_path}")

    raw_config = load_mapping_from_file(config_path)

    if not isinstance(raw_config, dict):
        raise ValueError(f"Pipeline config must be a dictionary, got {type(raw_config).__name__}")

    # Apply LR_ environment variable overrides if the config utility supports it
    try:
        from metainformant.core.utils.config import apply_env_overrides
        raw_config = apply_env_overrides(raw_config, prefix="LR")
    except (ImportError, AttributeError):
        pass

    logger.info("Loaded pipeline config from %s: pipeline=%s", config_path, raw_config.get("pipeline_name", "unknown"))
    return raw_config


def validate_pipeline_config(config: dict[str, Any], pipeline_name: str) -> list[str]:
    """Validate a pipeline configuration dictionary.

    Checks for required fields, valid parameter ranges, and step dependency
    consistency. Returns a list of error messages (empty if valid).

    Args:
        config: Pipeline configuration dictionary.
        pipeline_name: Expected pipeline name for validation context.

    Returns:
        List of validation error strings. Empty list means valid.
    """
    errors: list[str] = []

    # Check required top-level keys
    if "pipeline_name" not in config:
        errors.append("Missing required key 'pipeline_name'")
    elif config["pipeline_name"] != pipeline_name:
        errors.append(
            f"Pipeline name mismatch: expected '{pipeline_name}', got '{config['pipeline_name']}'"
        )

    if pipeline_name not in VALID_PIPELINE_NAMES:
        errors.append(f"Unknown pipeline name '{pipeline_name}'. Valid: {sorted(VALID_PIPELINE_NAMES)}")

    # Validate parameters
    params = config.get("parameters", {})
    bounds = _PARAM_BOUNDS.get(pipeline_name, {})

    for param_name, (lower, upper) in bounds.items():
        if param_name in params:
            value = params[param_name]
            if isinstance(value, (int, float)):
                if value < lower or value > upper:
                    errors.append(
                        f"Parameter '{param_name}' value {value} out of range [{lower}, {upper}]"
                    )

    # Validate step dependencies
    steps = config.get("steps", [])
    step_names = {step.get("name", "") for step in steps if isinstance(step, dict)}

    for step in steps:
        if not isinstance(step, dict):
            errors.append(f"Step must be a dictionary, got {type(step).__name__}")
            continue

        if "name" not in step:
            errors.append("Step missing required 'name' field")
            continue

        depends_on = step.get("depends_on", [])
        for dep in depends_on:
            if dep not in step_names:
                errors.append(
                    f"Step '{step['name']}' depends on unknown step '{dep}'"
                )

    if errors:
        logger.warning("Pipeline config validation found %d errors for '%s'", len(errors), pipeline_name)
    else:
        logger.info("Pipeline config validation passed for '%s'", pipeline_name)

    return errors
