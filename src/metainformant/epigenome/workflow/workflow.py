"""Epigenome workflow orchestration and analysis pipelines.

This module provides integrated workflows for analyzing DNA methylation,
ChIP-seq, and ATAC-seq data, combining multiple epigenomic assays
into comprehensive analysis pipelines.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Any, Optional, Iterator, Tuple, Set
from dataclasses import dataclass, field

from metainformant.core import logging, errors, validation, io, config, workflow, paths

logger = logging.get_logger(__name__)


@dataclass
class EpigenomeConfig:
    """Configuration for epigenome analysis workflows."""

    # Input/output paths
    input_dir: str | Path = "data/epigenome"
    output_dir: str | Path = "output/epigenome"
    temp_dir: str | Path = "temp/epigenome"

    # Analysis parameters
    methylation_threshold: float = 0.8  # Minimum methylation level for calling methylated sites
    chipseq_qvalue_threshold: float = 0.05  # FDR threshold for ChIP-seq peaks
    atacseq_qvalue_threshold: float = 0.05  # FDR threshold for ATAC-seq peaks

    # Peak filtering
    min_peak_length: int = 50  # Minimum peak length in bp
    max_peak_length: int = 10000  # Maximum peak length in bp

    # Coverage thresholds
    min_methylation_coverage: int = 5  # Minimum reads for methylation calling
    min_chip_coverage: int = 10  # Minimum coverage for ChIP-seq analysis
    min_atac_coverage: int = 5  # Minimum coverage for ATAC-seq analysis

    # Integration parameters
    max_distance_for_integration: int = 1000  # Max distance for associating peaks with methylation
    correlation_threshold: float = 0.3  # Minimum correlation for epigenetic associations

    # Workflow options
    run_methylation: bool = True
    run_chipseq: bool = True
    run_atacseq: bool = True
    integrate_results: bool = True
    generate_reports: bool = True

    # Performance settings
    n_threads: int = 4
    memory_limit_gb: float = 8.0
    chunk_size: int = 1000000  # Process data in chunks

    def __post_init__(self):
        """Validate configuration after initialization."""
        validation.validate_range(self.methylation_threshold, 0.0, 1.0, "methylation_threshold")
        validation.validate_range(self.chipseq_qvalue_threshold, 0.0, 1.0, "chipseq_qvalue_threshold")
        validation.validate_range(self.atacseq_qvalue_threshold, 0.0, 1.0, "atacseq_qvalue_threshold")
        validation.validate_range(self.correlation_threshold, -1.0, 1.0, "correlation_threshold")
        validation.validate_range(self.n_threads, 1, 64, "n_threads")
        validation.validate_range(self.memory_limit_gb, 0.1, 1024.0, "memory_limit_gb")


def load_epigenome_config(config_path: str | Path | None = None) -> EpigenomeConfig:
    """Load epigenome configuration from file or create default.

    Args:
        config_path: Path to configuration file (optional)

    Returns:
        EpigenomeConfig object
    """
    if config_path is None:
        # Try to find default config
        default_paths = [
            "config/epigenome.yaml",
            "config/epigenome.yml",
            "config/epigenome.json",
        ]
        for path in default_paths:
            if paths.validate_path_exists(path, raise_error=False):
                config_path = path
                break

    if config_path:
        try:
            config_data = config.load_config_file(config_path)
            return EpigenomeConfig(**config_data)
        except Exception as e:
            logger.warning(f"Failed to load config from {config_path}: {e}")
            logger.info("Using default configuration")

    return EpigenomeConfig()


def run_methylation_workflow(
    input_dir: str | Path, output_dir: str | Path, config: Optional[EpigenomeConfig] = None
) -> Dict[str, Any]:
    """Run DNA methylation analysis workflow.

    Args:
        input_dir: Directory containing methylation data files
        output_dir: Directory for output results
        config: Workflow configuration

    Returns:
        Dictionary with workflow results and statistics
    """
    if config is None:
        config = EpigenomeConfig()

    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    validation.validate_path_is_dir(input_dir, "input_dir")
    paths.ensure_directory(output_dir)

    logger.info(f"Starting methylation workflow: {input_dir} -> {output_dir}")

    results = {
        "workflow_type": "methylation",
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "config": config.__dict__,
        "processed_files": [],
        "statistics": {},
        "errors": [],
    }

    try:
        # Find methylation data files
        methylation_files = _find_methylation_files(input_dir)
        logger.info(f"Found {len(methylation_files)} methylation data files")

        if not methylation_files:
            raise errors.WorkflowError(f"No methylation data files found in {input_dir}")

        # Process each methylation file
        all_sites = {}
        all_stats = []

        for file_path in methylation_files:
            try:
                logger.info(f"Processing methylation file: {file_path}")

                # Load methylation data
                from .methylation import load_methylation_bedgraph, load_methylation_cov

                file_ext = file_path.suffix.lower()
                if file_ext == ".bedgraph" or file_ext == ".bg":
                    sites = load_methylation_bedgraph(file_path, min_coverage=config.min_methylation_coverage)
                elif file_ext == ".cov":
                    sites = load_methylation_cov(file_path, min_coverage=config.min_methylation_coverage)
                else:
                    logger.warning(f"Unsupported file format: {file_path}")
                    continue

                # Calculate statistics
                from .methylation import calculate_methylation_statistics

                stats = calculate_methylation_statistics(sites)

                # Store results
                sample_name = file_path.stem
                all_sites[sample_name] = sites
                all_stats.append(
                    {
                        "sample": sample_name,
                        "file": str(file_path),
                        "statistics": stats,
                    }
                )

                # Save individual sample results
                sample_output = output_dir / f"{sample_name}_methylation.json"
                io.dump_json(
                    {
                        "sample": sample_name,
                        "statistics": stats,
                        "sites_count": sum(len(sites_list) for sites_list in sites.values()),
                        "chromosomes": list(sites.keys()),
                    },
                    sample_output,
                )

                results["processed_files"].append(str(file_path))

            except Exception as e:
                error_msg = f"Failed to process {file_path}: {e}"
                logger.error(error_msg)
                results["errors"].append(error_msg)

        # Generate combined results
        results["statistics"] = {
            "total_samples": len(all_sites),
            "total_sites": sum(
                sum(len(sites_list) for sites_list in sample_sites.values()) for sample_sites in all_sites.values()
            ),
            "sample_statistics": all_stats,
        }

        # Save combined statistics
        stats_output = output_dir / "methylation_summary.json"
        io.dump_json(results["statistics"], stats_output)

        # Generate report
        if config.generate_reports:
            from .methylation import generate_methylation_report

            report_path = output_dir / "methylation_report.txt"
            report = generate_methylation_report(all_sites, output_path=report_path)
            results["report_path"] = str(report_path)

        logger.info(f"Methylation workflow completed. Processed {len(results['processed_files'])} files")

    except Exception as e:
        error_msg = f"Methylation workflow failed: {e}"
        logger.error(error_msg)
        results["errors"].append(error_msg)
        raise errors.WorkflowError(error_msg) from e

    return results


def run_chipseq_workflow(
    input_dir: str | Path, output_dir: str | Path, config: Optional[EpigenomeConfig] = None
) -> Dict[str, Any]:
    """Run ChIP-seq analysis workflow.

    Args:
        input_dir: Directory containing ChIP-seq data files
        output_dir: Directory for output results
        config: Workflow configuration

    Returns:
        Dictionary with workflow results and statistics
    """
    if config is None:
        config = EpigenomeConfig()

    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    validation.validate_path_is_dir(input_dir, "input_dir")
    paths.ensure_directory(output_dir)

    logger.info(f"Starting ChIP-seq workflow: {input_dir} -> {output_dir}")

    results = {
        "workflow_type": "chipseq",
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "config": config.__dict__,
        "processed_files": [],
        "statistics": {},
        "errors": [],
    }

    try:
        # Find ChIP-seq peak files
        peak_files = _find_chipseq_files(input_dir)
        logger.info(f"Found {len(peak_files)} ChIP-seq peak files")

        if not peak_files:
            raise errors.WorkflowError(f"No ChIP-seq peak files found in {input_dir}")

        # Process each peak file
        all_peaks = []
        all_stats = []

        for file_path in peak_files:
            try:
                logger.info(f"Processing ChIP-seq file: {file_path}")

                # Load peak data
                from .chipseq import load_chip_peaks, filter_peaks_by_score

                format_type = _detect_peak_format(file_path)
                peaks = load_chip_peaks(file_path, format=format_type)

                # Filter peaks
                filtered_peaks = filter_peaks_by_score(peaks, min_score=config.chipseq_qvalue_threshold)

                # Filter by length
                filtered_peaks = [
                    p for p in filtered_peaks if config.min_peak_length <= p.length <= config.max_peak_length
                ]

                # Calculate statistics
                from .chipseq import calculate_peak_statistics

                stats = calculate_peak_statistics(filtered_peaks)

                # Store results
                sample_name = file_path.stem
                all_peaks.extend(filtered_peaks)
                all_stats.append(
                    {
                        "sample": sample_name,
                        "file": str(file_path),
                        "statistics": stats,
                        "peaks_count": len(filtered_peaks),
                    }
                )

                # Save individual sample results
                sample_output = output_dir / f"{sample_name}_chipseq.json"
                io.dump_json(
                    {
                        "sample": sample_name,
                        "statistics": stats,
                        "peaks_count": len(filtered_peaks),
                    },
                    sample_output,
                )

                # Save filtered peaks
                peak_output = output_dir / f"{sample_name}_peaks.narrowPeak"
                from .chipseq import save_chip_peaks

                save_chip_peaks(filtered_peaks, peak_output)

                results["processed_files"].append(str(file_path))

            except Exception as e:
                error_msg = f"Failed to process {file_path}: {e}"
                logger.error(error_msg)
                results["errors"].append(error_msg)

        # Generate combined results
        results["statistics"] = {
            "total_samples": len(all_stats),
            "total_peaks": sum(stat["peaks_count"] for stat in all_stats),
            "sample_statistics": all_stats,
        }

        # Save combined statistics
        stats_output = output_dir / "chipseq_summary.json"
        io.dump_json(results["statistics"], stats_output)

        # Generate report
        if config.generate_reports:
            from .chipseq import generate_chip_report

            report_path = output_dir / "chipseq_report.txt"
            report = generate_chip_report(all_peaks, output_path=report_path)
            results["report_path"] = str(report_path)

        logger.info(f"ChIP-seq workflow completed. Processed {len(results['processed_files'])} files")

    except Exception as e:
        error_msg = f"ChIP-seq workflow failed: {e}"
        logger.error(error_msg)
        results["errors"].append(error_msg)
        raise errors.WorkflowError(error_msg) from e

    return results


def run_atacseq_workflow(
    input_dir: str | Path, output_dir: str | Path, config: Optional[EpigenomeConfig] = None
) -> Dict[str, Any]:
    """Run ATAC-seq analysis workflow.

    Args:
        input_dir: Directory containing ATAC-seq data files
        output_dir: Directory for output results
        config: Workflow configuration

    Returns:
        Dictionary with workflow results and statistics
    """
    if config is None:
        config = EpigenomeConfig()

    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    validation.validate_path_is_dir(input_dir, "input_dir")
    paths.ensure_directory(output_dir)

    logger.info(f"Starting ATAC-seq workflow: {input_dir} -> {output_dir}")

    results = {
        "workflow_type": "atacseq",
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "config": config.__dict__,
        "processed_files": [],
        "statistics": {},
        "errors": [],
    }

    try:
        # Find ATAC-seq peak files
        peak_files = _find_atacseq_files(input_dir)
        logger.info(f"Found {len(peak_files)} ATAC-seq peak files")

        if not peak_files:
            raise errors.WorkflowError(f"No ATAC-seq peak files found in {input_dir}")

        # Process each peak file
        all_peaks = []
        all_stats = []

        for file_path in peak_files:
            try:
                logger.info(f"Processing ATAC-seq file: {file_path}")

                # Load peak data
                from .atacseq import load_atac_peaks, filter_peaks_by_score

                format_type = _detect_peak_format(file_path)
                peaks = load_atac_peaks(file_path, format=format_type)

                # Filter peaks (using accessibility score as threshold)
                filtered_peaks = [p for p in peaks if p.accessibility_score >= config.atacseq_qvalue_threshold]

                # Filter by length
                filtered_peaks = [
                    p for p in filtered_peaks if config.min_peak_length <= p.length <= config.max_peak_length
                ]

                # Calculate statistics
                from .atacseq import calculate_atac_statistics

                stats = calculate_atac_statistics(filtered_peaks)

                # Store results
                sample_name = file_path.stem
                all_peaks.extend(filtered_peaks)
                all_stats.append(
                    {
                        "sample": sample_name,
                        "file": str(file_path),
                        "statistics": stats,
                        "peaks_count": len(filtered_peaks),
                    }
                )

                # Save individual sample results
                sample_output = output_dir / f"{sample_name}_atacseq.json"
                io.dump_json(
                    {
                        "sample": sample_name,
                        "statistics": stats,
                        "peaks_count": len(filtered_peaks),
                    },
                    sample_output,
                )

                # Save filtered peaks
                peak_output = output_dir / f"{sample_name}_peaks.narrowPeak"
                from .atacseq import save_atac_peaks

                save_atac_peaks(filtered_peaks, peak_output)

                results["processed_files"].append(str(file_path))

            except Exception as e:
                error_msg = f"Failed to process {file_path}: {e}"
                logger.error(error_msg)
                results["errors"].append(error_msg)

        # Generate combined results
        results["statistics"] = {
            "total_samples": len(all_stats),
            "total_peaks": sum(stat["peaks_count"] for stat in all_stats),
            "sample_statistics": all_stats,
        }

        # Save combined statistics
        stats_output = output_dir / "atacseq_summary.json"
        io.dump_json(results["statistics"], stats_output)

        # Generate report
        if config.generate_reports:
            from .atacseq import generate_atac_report

            report_path = output_dir / "atacseq_report.txt"
            report = generate_atac_report(all_peaks, output_path=report_path)
            results["report_path"] = str(report_path)

        logger.info(f"ATAC-seq workflow completed. Processed {len(results['processed_files'])} files")

    except Exception as e:
        error_msg = f"ATAC-seq workflow failed: {e}"
        logger.error(error_msg)
        results["errors"].append(error_msg)
        raise errors.WorkflowError(error_msg) from e

    return results


def integrate_epigenome_results(
    methylation_results: Dict[str, Any],
    chipseq_results: Dict[str, Any],
    atacseq_results: Dict[str, Any],
    output_dir: str | Path,
    config: Optional[EpigenomeConfig] = None,
) -> Dict[str, Any]:
    """Integrate results from multiple epigenomic assays.

    Args:
        methylation_results: Results from methylation workflow
        chipseq_results: Results from ChIP-seq workflow
        atacseq_results: Results from ATAC-seq workflow
        output_dir: Directory for integrated results
        config: Workflow configuration

    Returns:
        Dictionary with integrated analysis results
    """
    if config is None:
        config = EpigenomeConfig()

    output_dir = Path(output_dir)
    paths.ensure_directory(output_dir)

    logger.info(f"Integrating epigenome results to {output_dir}")

    integrated_results = {
        "integration_type": "epigenome_multi_assay",
        "output_dir": str(output_dir),
        "config": config.__dict__,
        "methylation_summary": methylation_results.get("statistics", {}),
        "chipseq_summary": chipseq_results.get("statistics", {}),
        "atacseq_summary": atacseq_results.get("statistics", {}),
        "integrated_analyses": {},
        "errors": [],
    }

    try:
        # Cross-assay correlations and associations
        analyses = {}

        # 1. Methylation-ChIP-seq associations
        if (
            methylation_results.get("statistics", {}).get("total_samples", 0) > 0
            and chipseq_results.get("statistics", {}).get("total_samples", 0) > 0
        ):
            try:
                methylation_chip_associations = _analyze_methylation_chip_associations(
                    methylation_results, chipseq_results, config
                )
                analyses["methylation_chipseq"] = methylation_chip_associations
            except Exception as e:
                logger.warning(f"Methylation-ChIP-seq association analysis failed: {e}")
                integrated_results["errors"].append(f"methylation_chipseq: {e}")

        # 2. Methylation-ATAC-seq associations
        if (
            methylation_results.get("statistics", {}).get("total_samples", 0) > 0
            and atacseq_results.get("statistics", {}).get("total_samples", 0) > 0
        ):
            try:
                methylation_atac_associations = _analyze_methylation_atac_associations(
                    methylation_results, atacseq_results, config
                )
                analyses["methylation_atacseq"] = methylation_atac_associations
            except Exception as e:
                logger.warning(f"Methylation-ATAC-seq association analysis failed: {e}")
                integrated_results["errors"].append(f"methylation_atacseq: {e}")

        # 3. ChIP-seq-ATAC-seq associations
        if (
            chipseq_results.get("statistics", {}).get("total_samples", 0) > 0
            and atacseq_results.get("statistics", {}).get("total_samples", 0) > 0
        ):
            try:
                chip_atac_associations = _analyze_chip_atac_associations(chipseq_results, atacseq_results, config)
                analyses["chipseq_atacseq"] = chip_atac_associations
            except Exception as e:
                logger.warning(f"ChIP-seq-ATAC-seq association analysis failed: {e}")
                integrated_results["errors"].append(f"chipseq_atacseq: {e}")

        integrated_results["integrated_analyses"] = analyses

        # Save integrated results
        output_file = output_dir / "integrated_epigenome_results.json"
        io.dump_json(integrated_results, output_file)

        # Generate integration report
        if config.generate_reports:
            report_path = output_dir / "epigenome_integration_report.txt"
            report = _generate_integration_report(integrated_results, report_path)
            integrated_results["report_path"] = str(report_path)

        logger.info(f"Epigenome integration completed. Performed {len(analyses)} analyses")

    except Exception as e:
        error_msg = f"Epigenome integration failed: {e}"
        logger.error(error_msg)
        integrated_results["errors"].append(error_msg)
        raise errors.WorkflowError(error_msg) from e

    return integrated_results


def _find_methylation_files(input_dir: Path) -> List[Path]:
    """Find methylation data files in input directory."""
    extensions = ["*.bedgraph", "*.bg", "*.cov", "*.bedGraph", "*.BEDGRAPH"]
    files = []
    for ext in extensions:
        files.extend(input_dir.glob(f"**/{ext}"))
    return sorted(list(set(files)))


def _find_chipseq_files(input_dir: Path) -> List[Path]:
    """Find ChIP-seq peak files in input directory."""
    extensions = ["*.narrowPeak", "*.broadPeak", "*.bed", "*.narrowpeak", "*.broadpeak"]
    files = []
    for ext in extensions:
        files.extend(input_dir.glob(f"**/{ext}"))
    return sorted(list(set(files)))


def _find_atacseq_files(input_dir: Path) -> List[Path]:
    """Find ATAC-seq peak files in input directory."""
    # ATAC-seq often uses same formats as ChIP-seq
    return _find_chipseq_files(input_dir)


def _detect_peak_format(file_path: Path) -> str:
    """Detect peak file format from filename or content."""
    filename = file_path.name.lower()

    if "narrow" in filename:
        return "narrowpeak"
    elif "broad" in filename:
        return "broadpeak"
    else:
        # Default to narrowPeak format
        return "narrowpeak"


def _analyze_methylation_chip_associations(
    methylation_results: Dict[str, Any], chipseq_results: Dict[str, Any], config: EpigenomeConfig
) -> Dict[str, Any]:
    """Analyze associations between methylation and ChIP-seq data.

    Performs coordinate-based intersection analysis to identify overlapping
    epigenetic marks between methylation sites and ChIP-seq peaks.
    """
    associations = {
        "analysis_type": "methylation_chipseq_association",
        "max_distance": config.max_distance_for_integration,
        "correlation_threshold": config.correlation_threshold,
        "findings": [],
        "statistics": {},
    }

    # Extract methylation sites
    meth_sites = methylation_results.get("sites", [])
    chip_peaks = chipseq_results.get("peaks", [])

    if not meth_sites:
        associations["findings"].append({
            "type": "warning",
            "description": "No methylation sites provided for analysis",
        })
        return associations

    if not chip_peaks:
        associations["findings"].append({
            "type": "warning",
            "description": "No ChIP-seq peaks provided for analysis",
        })
        return associations

    # Perform intersection analysis
    max_dist = config.max_distance_for_integration
    overlapping = 0
    proximal = 0
    distal = 0

    for meth in meth_sites:
        meth_chrom = meth.get("chrom", meth.get("chromosome"))
        meth_pos = meth.get("position", meth.get("pos"))
        if meth_chrom is None or meth_pos is None:
            continue

        min_distance = float("inf")
        for peak in chip_peaks:
            peak_chrom = peak.get("chrom", peak.get("chromosome"))
            if peak_chrom != meth_chrom:
                continue

            peak_start = peak.get("start", peak.get("position"))
            peak_end = peak.get("end", peak_start + 1 if peak_start else None)
            if peak_start is None:
                continue

            # Calculate distance
            if peak_start <= meth_pos <= peak_end:
                distance = 0
            elif meth_pos < peak_start:
                distance = peak_start - meth_pos
            else:
                distance = meth_pos - peak_end

            min_distance = min(min_distance, distance)

        # Categorize based on distance
        if min_distance == 0:
            overlapping += 1
        elif min_distance <= max_dist:
            proximal += 1
        elif min_distance < float("inf"):
            distal += 1

    total_sites = len(meth_sites)
    associations["statistics"] = {
        "total_methylation_sites": total_sites,
        "total_chip_peaks": len(chip_peaks),
        "overlapping": overlapping,
        "proximal": proximal,
        "distal": distal,
        "overlap_rate": overlapping / total_sites if total_sites > 0 else 0,
        "proximal_rate": proximal / total_sites if total_sites > 0 else 0,
    }

    # Generate findings based on statistics
    if overlapping > 0:
        associations["findings"].append({
            "type": "overlap",
            "description": f"{overlapping} methylation sites ({associations['statistics']['overlap_rate']:.1%}) overlap with ChIP-seq peaks",
            "count": overlapping,
        })

    if proximal > 0:
        associations["findings"].append({
            "type": "proximal",
            "description": f"{proximal} methylation sites are within {max_dist}bp of ChIP-seq peaks",
            "count": proximal,
        })

    return associations


def _analyze_methylation_atac_associations(
    methylation_results: Dict[str, Any], atacseq_results: Dict[str, Any], config: EpigenomeConfig
) -> Dict[str, Any]:
    """Analyze associations between methylation and ATAC-seq data.

    Analyzes chromatin accessibility around methylated sites to identify
    relationships between methylation and open chromatin regions.
    """
    associations = {
        "analysis_type": "methylation_atacseq_association",
        "max_distance": config.max_distance_for_integration,
        "correlation_threshold": config.correlation_threshold,
        "findings": [],
        "statistics": {},
    }

    # Extract methylation sites and ATAC-seq peaks
    meth_sites = methylation_results.get("sites", [])
    atac_peaks = atacseq_results.get("peaks", atacseq_results.get("accessible_regions", []))

    if not meth_sites:
        associations["findings"].append({
            "type": "warning",
            "description": "No methylation sites provided for analysis",
        })
        return associations

    if not atac_peaks:
        associations["findings"].append({
            "type": "warning",
            "description": "No ATAC-seq peaks/accessible regions provided for analysis",
        })
        return associations

    # Analyze methylation in accessible vs closed chromatin
    max_dist = config.max_distance_for_integration
    in_accessible = 0
    near_accessible = 0
    in_closed = 0

    meth_values_accessible = []
    meth_values_closed = []

    for meth in meth_sites:
        meth_chrom = meth.get("chrom", meth.get("chromosome"))
        meth_pos = meth.get("position", meth.get("pos"))
        meth_level = meth.get("methylation_level", meth.get("beta", None))
        if meth_chrom is None or meth_pos is None:
            continue

        # Check if in accessible region
        is_accessible = False
        min_distance = float("inf")

        for peak in atac_peaks:
            peak_chrom = peak.get("chrom", peak.get("chromosome"))
            if peak_chrom != meth_chrom:
                continue

            peak_start = peak.get("start", peak.get("position"))
            peak_end = peak.get("end", peak_start + 1 if peak_start else None)
            if peak_start is None:
                continue

            if peak_start <= meth_pos <= peak_end:
                is_accessible = True
                min_distance = 0
                break
            else:
                if meth_pos < peak_start:
                    distance = peak_start - meth_pos
                else:
                    distance = meth_pos - peak_end
                min_distance = min(min_distance, distance)

        if is_accessible:
            in_accessible += 1
            if meth_level is not None:
                meth_values_accessible.append(meth_level)
        elif min_distance <= max_dist:
            near_accessible += 1
            if meth_level is not None:
                meth_values_accessible.append(meth_level)
        else:
            in_closed += 1
            if meth_level is not None:
                meth_values_closed.append(meth_level)

    total_sites = len(meth_sites)
    associations["statistics"] = {
        "total_methylation_sites": total_sites,
        "total_accessible_regions": len(atac_peaks),
        "sites_in_accessible": in_accessible,
        "sites_near_accessible": near_accessible,
        "sites_in_closed": in_closed,
        "accessible_rate": (in_accessible + near_accessible) / total_sites if total_sites > 0 else 0,
    }

    # Calculate mean methylation levels if available
    if meth_values_accessible:
        mean_accessible = sum(meth_values_accessible) / len(meth_values_accessible)
        associations["statistics"]["mean_methylation_accessible"] = mean_accessible
    if meth_values_closed:
        mean_closed = sum(meth_values_closed) / len(meth_values_closed)
        associations["statistics"]["mean_methylation_closed"] = mean_closed

    # Generate findings
    if in_accessible > 0:
        associations["findings"].append({
            "type": "accessible_overlap",
            "description": f"{in_accessible} methylation sites ({in_accessible/total_sites:.1%}) located within accessible chromatin",
            "count": in_accessible,
        })

    if meth_values_accessible and meth_values_closed:
        mean_acc = sum(meth_values_accessible) / len(meth_values_accessible)
        mean_clo = sum(meth_values_closed) / len(meth_values_closed)
        if mean_acc < mean_clo:
            associations["findings"].append({
                "type": "methylation_pattern",
                "description": f"Lower methylation in accessible regions (mean: {mean_acc:.3f}) vs closed chromatin (mean: {mean_clo:.3f})",
                "accessible_mean": mean_acc,
                "closed_mean": mean_clo,
            })

    return associations


def _analyze_chip_atac_associations(
    chipseq_results: Dict[str, Any], atacseq_results: Dict[str, Any], config: EpigenomeConfig
) -> Dict[str, Any]:
    """Analyze associations between ChIP-seq and ATAC-seq data."""
    associations = {
        "analysis_type": "chipseq_atacseq_association",
        "max_distance": config.max_distance_for_integration,
        "correlation_threshold": config.correlation_threshold,
        "findings": [],
    }

    # Placeholder for actual association analysis
    associations["findings"].append(
        {
            "type": "placeholder",
            "description": "ChIP-seq-ATAC-seq association analysis framework ready",
            "note": "Actual implementation would analyze TF binding sites in open chromatin regions",
        }
    )

    return associations


def _generate_integration_report(results: Dict[str, Any], output_path: Path) -> str:
    """Generate integration analysis report."""
    report_lines = []
    report_lines.append("=" * 70)
    report_lines.append("EPIGENOME INTEGRATION ANALYSIS REPORT")
    report_lines.append("=" * 70)
    report_lines.append("")

    # Summary statistics
    methylation_samples = results.get("methylation_summary", {}).get("total_samples", 0)
    chipseq_samples = results.get("chipseq_summary", {}).get("total_samples", 0)
    atacseq_samples = results.get("atacseq_summary", {}).get("total_samples", 0)

    report_lines.append("Dataset Summary:")
    report_lines.append(f"  Methylation samples: {methylation_samples}")
    report_lines.append(f"  ChIP-seq samples: {chipseq_samples}")
    report_lines.append(f"  ATAC-seq samples: {atacseq_samples}")
    report_lines.append("")

    # Integrated analyses
    analyses = results.get("integrated_analyses", {})
    if analyses:
        report_lines.append("Integrated Analyses:")
        for analysis_name, analysis_data in analyses.items():
            report_lines.append(
                f"  {analysis_name.replace('_', ' ').title()}: {len(analysis_data.get('findings', []))} findings"
            )
        report_lines.append("")

    # Errors
    errors = results.get("errors", [])
    if errors:
        report_lines.append("Errors/Warnings:")
        for error in errors:
            report_lines.append(f"  - {error}")
        report_lines.append("")

    report = "\n".join(report_lines)

    with open(output_path, "w") as f:
        f.write(report)

    logger.info(f"Integration report saved to {output_path}")
    return report
