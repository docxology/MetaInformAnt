"""GWAS workflow orchestration."""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import Any

from ..core.io import dump_json, ensure_directory, read_delimited
from .association import run_gwas
from .calling import call_variants_bcftools, call_variants_gatk
from .config import GWASWorkflowConfig
from .correction import bonferroni_correction, fdr_correction, genomic_control
from .download import download_reference_genome, download_variant_data
from .quality import apply_qc_filters, parse_vcf_full
from .structure import estimate_population_structure
from .visualization import manhattan_plot, qq_plot, regional_plot

logger = logging.getLogger(__name__)


def execute_gwas_workflow(config: GWASWorkflowConfig, *, check: bool = False) -> dict[str, Any]:
    """Execute the full GWAS workflow in order.

    Workflow steps:
    1. Configuration loading
    2. Genome preparation (if needed)
    3. Variant acquisition (download, call, or use existing)
    4. Quality control
    5. Population structure analysis
    6. Phenotype loading
    7. Association testing
    8. Multiple testing correction
    9. Visualization
    10. Results export

    Args:
        config: GWAS workflow configuration
        check: If True, only validate configuration without execution

    Returns:
        Dictionary with workflow results and metadata
    """
    logger.info("execute_gwas_workflow: Starting GWAS workflow execution")
    ensure_directory(config.work_dir)

    if check:
        logger.info("execute_gwas_workflow: Configuration check mode - validation only")
        return {"status": "validated", "config": config}

    start_time = datetime.utcnow()
    results: dict[str, Any] = {
        "start_time": start_time.isoformat(),
        "config": {
            "work_dir": str(config.work_dir),
            "threads": config.threads,
        },
        "steps": [],
    }

    try:
        # Step 1: Genome preparation (if needed)
        if config.genome:
            logger.info("execute_gwas_workflow: Preparing reference genome")
            genome_result = _prepare_genome(config)
            results["steps"].append({"step": "genome_prepare", "result": genome_result})

        # Step 2: Variant acquisition
        logger.info("execute_gwas_workflow: Acquiring variant data")
        variant_result = _acquire_variants(config)
        results["steps"].append({"step": "variant_acquisition", "result": variant_result})

        # Get VCF path for downstream steps
        vcf_path = None
        if variant_result.get("status") == "success":
            if variant_result.get("method") == "existing_vcf":
                vcf_files = variant_result.get("vcf_files", [])
                if vcf_files:
                    vcf_path = Path(vcf_files[0])
            elif variant_result.get("output_vcf"):
                vcf_path = Path(variant_result["output_vcf"])

        if not vcf_path or not vcf_path.exists():
            raise ValueError(f"VCF file not available: {vcf_path}")

        # Step 3: Quality control
        logger.info("execute_gwas_workflow: Applying quality control filters")
        qc_result = _apply_quality_control(config, vcf_path)
        results["steps"].append({"step": "quality_control", "result": qc_result})

        # Use QC'd VCF if available, otherwise original
        qc_vcf_path = qc_result.get("output_vcf")
        if qc_vcf_path and Path(qc_vcf_path).exists():
            vcf_path = Path(qc_vcf_path)

        # Step 4: Population structure
        logger.info("execute_gwas_workflow: Analyzing population structure")
        structure_result = _analyze_population_structure(config, vcf_path)
        results["steps"].append({"step": "population_structure", "result": structure_result})

        # Step 5: Association testing
        logger.info("execute_gwas_workflow: Running association tests")
        association_result = _run_association_tests(config, vcf_path)
        results["steps"].append({"step": "association_testing", "result": association_result})
        
        # Check if association failed critically - fail workflow
        if association_result.get("status") == "failed":
            error_msg = association_result.get("error", "")
            # Fail workflow if phenotype file is missing or other critical errors
            if "not found" in error_msg.lower() or "cannot read" in error_msg.lower() or "phenotype" in error_msg.lower():
                raise ValueError(f"Association testing failed: {error_msg}")

        # Step 6: Multiple testing correction
        logger.info("execute_gwas_workflow: Applying multiple testing correction")
        correction_result = _apply_corrections(config, association_result)
        results["steps"].append({"step": "multiple_testing", "result": correction_result})

        # Step 7: Visualization
        logger.info("execute_gwas_workflow: Generating visualizations")
        viz_result = _generate_visualizations(config, association_result)
        results["steps"].append({"step": "visualization", "result": viz_result})

        # Step 8: Results export
        logger.info("execute_gwas_workflow: Exporting results")
        export_result = _export_results(config, results)
        results["steps"].append({"step": "results_export", "result": export_result})

        end_time = datetime.utcnow()
        duration = (end_time - start_time).total_seconds()
        results["end_time"] = end_time.isoformat()
        results["duration_seconds"] = duration
        results["status"] = "completed"

        logger.info(f"execute_gwas_workflow: GWAS workflow completed in {duration:.2f} seconds")

    except Exception as exc:
        end_time = datetime.utcnow()
        duration = (end_time - start_time).total_seconds()
        results["end_time"] = end_time.isoformat()
        results["duration_seconds"] = duration
        results["status"] = "failed"
        results["error"] = str(exc)
        logger.error(f"execute_gwas_workflow: Workflow failed after {duration:.2f} seconds: {exc}", exc_info=True)
        # Don't raise - return failed status instead

    # Write results summary
    results_path = config.work_dir / "workflow_results.json"
    dump_json(results, results_path, indent=2)
    logger.info(f"execute_gwas_workflow: Results written to {results_path}")

    return results


def _prepare_genome(config: GWASWorkflowConfig) -> dict[str, Any]:
    """Prepare reference genome (download if needed)."""
    if not config.genome:
        return {"status": "skipped", "message": "No genome configuration provided"}

    accession = config.genome.get("accession")
    if not accession:
        return {"status": "skipped", "message": "No genome accession provided"}

    dest_dir = config.genome.get("dest_dir", config.work_dir / "genome")
    include = config.genome.get("include", ["genome", "gff3"])
    ftp_url = config.genome.get("ftp_url")

    result = download_reference_genome(
        accession=accession,
        dest_dir=dest_dir,
        include=include,
        ftp_url=ftp_url,
    )

    return result


def _acquire_variants(config: GWASWorkflowConfig) -> dict[str, Any]:
    """Acquire variant data (download, call, or use existing)."""
    variants_cfg = config.variants or {}

    # Option 1: Pre-existing VCF files
    if vcf_files := variants_cfg.get("vcf_files"):
        if isinstance(vcf_files, list) and vcf_files:
            logger.info(f"Using {len(vcf_files)} pre-existing VCF file(s)")
            return {
                "status": "success",
                "method": "existing_vcf",
                "vcf_files": vcf_files,
                "count": len(vcf_files),
            }

    # Option 2: Download from public database
    if download_cfg := variants_cfg.get("download"):
        source = download_cfg.get("source")
        accession = download_cfg.get("accession")
        region = download_cfg.get("region")
        dest_dir = variants_cfg.get("dest_dir", config.work_dir / "variants")

        if source:
            result = download_variant_data(
                source=source,
                accession=accession,
                region=region,
                dest_dir=dest_dir,
            )
            return result

    # Option 3: Call variants from BAM/CRAM
    if calling_cfg := variants_cfg.get("calling"):
        bam_files = calling_cfg.get("bam_files", [])
        reference = calling_cfg.get("reference")
        method = calling_cfg.get("method", "bcftools")
        dest_dir = variants_cfg.get("dest_dir", config.work_dir / "variants")
        threads = config.threads

        if not bam_files:
            return {
                "status": "failed",
                "error": "No BAM files specified for variant calling",
            }

        if not reference:
            return {
                "status": "failed",
                "error": "Reference genome not specified for variant calling",
            }

        # Prepare output VCF path
        output_vcf = Path(dest_dir) / "called_variants.vcf.gz"
        ensure_directory(output_vcf.parent)

        logger.info(f"Calling variants from {len(bam_files)} BAM file(s) using {method}")

        # Call variants based on method
        if method == "bcftools":
            result = call_variants_bcftools(
                bam_files=bam_files,
                reference=reference,
                output_vcf=output_vcf,
                threads=threads,
            )
        elif method == "gatk":
            result = call_variants_gatk(
                bam_files=bam_files,
                reference=reference,
                output_vcf=output_vcf,
                threads=threads,
            )
        else:
            return {
                "status": "failed",
                "error": f"Unsupported variant calling method: {method}",
            }

        if result.get("status") == "success":
            result["method"] = "calling"
            result["calling_method"] = method
        return result

    return {"status": "skipped", "message": "No variant source specified"}


def _apply_quality_control(config: GWASWorkflowConfig, vcf_path: str | Path) -> dict[str, Any]:
    """Apply quality control filters to variants."""
    qc_config = config.qc or {}
    variants_dir = config.work_dir / "variants"
    output_vcf = variants_dir / "variants_qc_filtered.vcf"

    result = apply_qc_filters(
        vcf_path=vcf_path,
        config=qc_config,
        output_vcf=output_vcf,
    )

    return result


def _analyze_population_structure(config: GWASWorkflowConfig, vcf_path: str | Path) -> dict[str, Any]:
    """Analyze population structure (PCA, kinship)."""
    structure_config = config.structure or {}
    structure_dir = config.work_dir / "structure"

    result = estimate_population_structure(
        vcf_path=vcf_path,
        config=structure_config,
        output_dir=structure_dir,
    )

    return result


def _run_association_tests(config: GWASWorkflowConfig, vcf_path: str | Path) -> dict[str, Any]:
    """Run association tests."""
    association_config = config.association or {}
    samples_config = config.samples or {}

    phenotype_path = samples_config.get("phenotype_file")
    if not phenotype_path:
        return {"status": "failed", "error": "Phenotype file not specified"}

    results_dir = config.work_dir / "results"
    result = run_gwas(
        vcf_path=vcf_path,
        phenotype_path=phenotype_path,
        config=association_config,
        output_dir=results_dir,
    )

    return result


def _apply_corrections(config: GWASWorkflowConfig, association_result: dict[str, Any]) -> dict[str, Any]:
    """Apply multiple testing corrections."""
    correction_config = config.correction or {}
    method = correction_config.get("method", "bonferroni")
    alpha = correction_config.get("alpha", 0.05)

    results = association_result.get("results", [])
    if not results:
        return {"status": "failed", "error": "No association results to correct"}

    pvalues = [r.get("p_value", 1.0) for r in results]

    correction_results: dict[str, Any] = {}

    if method == "bonferroni":
        bonf_result = bonferroni_correction(pvalues, alpha=alpha)
        correction_results["bonferroni"] = bonf_result
    elif method == "fdr":
        fdr_result = fdr_correction(pvalues, alpha=alpha)
        correction_results["fdr"] = fdr_result
    elif method == "both":
        bonf_result = bonferroni_correction(pvalues, alpha=alpha)
        fdr_result = fdr_correction(pvalues, alpha=alpha)
        correction_results["bonferroni"] = bonf_result
        correction_results["fdr"] = fdr_result

    # Genomic control
    gc_result = genomic_control(pvalues=pvalues)
    correction_results["genomic_control"] = gc_result

    return {
        "status": "success",
        "method": method,
        "corrections": correction_results,
    }


def _generate_visualizations(config: GWASWorkflowConfig, association_result: dict[str, Any]) -> dict[str, Any]:
    """Generate visualizations (Manhattan plots, Q-Q plots, and optionally comprehensive suite).
    
    If comprehensive=True in config.output, generates all available plot types.
    Otherwise, generates standard Manhattan and QQ plots.
    """
    output_config = config.output or {}
    plots_dir = output_config.get("plots_dir", config.work_dir / "plots")
    ensure_directory(plots_dir)

    results = association_result.get("results", [])
    if not results:
        return {"status": "failed", "error": "No results to visualize"}

    # Check if comprehensive visualization is requested
    comprehensive = output_config.get("comprehensive_plots", False)
    
    if comprehensive:
        logger.info("_generate_visualizations: Generating comprehensive visualization suite")
        try:
            from .visualization_comprehensive import generate_all_plots
            
            # Get paths to required files
            results_file = association_result.get("output_file")
            if not results_file or not Path(results_file).exists():
                logger.warning("_generate_visualizations: Results file not found for comprehensive plots")
                comprehensive = False  # Fall back to basic plots
            else:
                # Try to get PCA and kinship files
                pca_file = config.work_dir / "structure" / "pca_components.tsv"
                kinship_file = config.work_dir / "structure" / "kinship_matrix.tsv"
                
                comprehensive_result = generate_all_plots(
                    association_results=Path(results_file),
                    output_dir=Path(plots_dir),
                    pca_file=pca_file if pca_file.exists() else None,
                    kinship_file=kinship_file if kinship_file.exists() else None,
                    significance_threshold=output_config.get("significance_threshold", 5e-8),
                )
                
                return {
                    "status": "success",
                    "plots_dir": str(plots_dir),
                    "mode": "comprehensive",
                    "visualizations": comprehensive_result,
                }
        except Exception as e:
            logger.error(f"_generate_visualizations: Comprehensive plots failed: {e}")
            comprehensive = False  # Fall back to basic plots
    
    # Generate standard plots (basic or fallback)
    logger.info("_generate_visualizations: Generating standard visualizations")
    viz_results: dict[str, Any] = {}

    # Manhattan plot
    manhattan_path = Path(plots_dir) / "manhattan.png"
    manhattan_result = manhattan_plot(
        results=results,
        output_path=manhattan_path,
        title="GWAS Manhattan Plot",
    )
    viz_results["manhattan"] = manhattan_result

    # Q-Q plot
    qq_path = Path(plots_dir) / "qq_plot.png"
    pvalues = [r.get("p_value", 1.0) for r in results]
    qq_result = qq_plot(
        pvalues=pvalues,
        output_path=qq_path,
        title="GWAS Q-Q Plot",
    )
    viz_results["qq"] = qq_result

    return {
        "status": "success",
        "plots_dir": str(plots_dir),
        "mode": "standard",
        "visualizations": viz_results,
    }


def _export_results(config: GWASWorkflowConfig, workflow_results: dict[str, Any]) -> dict[str, Any]:
    """Export results to files."""
    output_config = config.output or {}
    results_dir = output_config.get("results_dir", config.work_dir / "results")
    ensure_directory(results_dir)

    # Results are already written by association module
    # Write summary JSON
    summary_path = Path(results_dir) / "summary.json"
    dump_json(workflow_results, summary_path, indent=2)

    return {
        "status": "success",
        "results_dir": str(results_dir),
        "summary_path": str(summary_path),
    }

