"""GWAS workflow execution functions.

This module provides the main workflow execution functions:
- execute_gwas_workflow: Execute the complete GWAS workflow from config
- run_gwas: Run complete GWAS with VCF/phenotype paths and config
- run_multi_trait_gwas: Run GWAS for multiple traits sharing common steps
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from metainformant.core.utils import logging

from .workflow_config import (
    _extract_genotype_matrix,
    _load_phenotypes,
    _load_phenotypes_by_id,
    _normalize_config,
    validate_gwas_config,
)

logger = logging.get_logger(__name__)

# Import required functions from GWAS modules
from metainformant.gwas.analysis.annotation import annotate_variants_with_genes
from metainformant.gwas.analysis.association import association_test_linear
from metainformant.gwas.analysis.correction import bonferroni_correction, fdr_correction
from metainformant.gwas.analysis.heritability import estimate_heritability, partition_heritability_by_chromosome
from metainformant.gwas.analysis.ld_pruning import ld_prune
from metainformant.gwas.analysis.mixed_model import association_test_mixed, run_mixed_model_gwas
from metainformant.gwas.analysis.quality import apply_qc_filters, check_haplodiploidy, parse_vcf_full
from metainformant.gwas.analysis.structure import compute_kinship_matrix, compute_pca
from metainformant.gwas.data.metadata import get_population_labels, load_sample_metadata, validate_metadata
from metainformant.gwas.visualization.config import apply_style, style_from_config
from metainformant.gwas.visualization.general import generate_all_plots
from metainformant.gwas.visualization.interactive.finemapping import compute_credible_set


def execute_gwas_workflow(config: Dict[str, Any], *, check: bool = False) -> Dict[str, Any]:
    """Execute the complete GWAS workflow.

    Args:
        config: GWAS configuration dictionary
        check: If True, only validate configuration without executing

    Returns:
        Workflow results dictionary
    """
    logger.info("Starting GWAS workflow execution")

    # Normalize nested YAML config to flat format
    config = _normalize_config(config)

    if check:
        logger.info("Running in check mode - validating configuration only")
        # Validate configuration
        errors = validate_gwas_config(config)
        if errors:
            return {"status": "invalid", "success": False, "errors": errors, "config_valid": False, "config": config}
        return {"status": "validated", "success": True, "errors": [], "config_valid": True, "config": config}

    # Execute workflow steps
    results: Dict[str, Any] = {
        "success": True,
        "status": "running",
        "steps_completed": [],
        "steps": [],
        "errors": [],
        "outputs": {},
    }

    try:
        # Step 1: Load and validate data
        logger.info("Step 1: Loading and validating data")
        vcf_data = parse_vcf_full(config["vcf_path"])
        results["steps_completed"].append("data_loading")
        results["steps"].append("data_loading")

        # Step 2: Apply QC filters
        logger.info("Step 2: Applying quality control filters")
        qc_config = config.get("quality_control", {})
        qc_result = apply_qc_filters(vcf_data, qc_config)
        filtered_data = qc_result.get("filtered_data", vcf_data)
        results["steps_completed"].append("quality_control")
        results["steps"].append("quality_control")

        # Step 3: Population structure analysis
        logger.info("Step 3: Analyzing population structure")
        genotype_matrix = _extract_genotype_matrix(filtered_data)
        pca_result = compute_pca(genotype_matrix)
        kinship_matrix = compute_kinship_matrix(genotype_matrix)
        results["steps_completed"].append("population_structure")
        results["steps"].append("population_structure")

        # Step 4: Association testing
        logger.info("Step 4: Performing association testing")
        phenotypes = _load_phenotypes(config["phenotype_path"])

        association_results = []
        for i, phenotype in enumerate(phenotypes):
            result = association_test_linear(
                filtered_data["genotypes"][:, i],  # Use first trait for now
                phenotype,
                covariates=pca_result[0][:, :10],  # Use first 10 PCs as covariates
            )
            association_results.append(result)

        results["steps_completed"].append("association_testing")
        results["steps"].append("association_testing")
        results["outputs"]["association_results"] = association_results

        # Step 5: Multiple testing correction
        logger.info("Step 5: Applying multiple testing correction")
        p_values = [r["p_value"] for r in association_results]
        corrected_results = fdr_correction(p_values)
        results["steps_completed"].append("multiple_testing_correction")
        results["steps"].append("multiple_testing_correction")

        # Step 6: Generate plots
        logger.info("Step 6: Generating visualization plots")
        plot_results = generate_all_plots(
            association_results,
            config.get("output_dir", config.get("work_dir", ".")),
            pca_file=Path(config.get("work_dir", ".")) / "pca_results.npy",
            kinship_file=Path(config.get("work_dir", ".")) / "kinship_matrix.npy",
            vcf_file=Path(config["vcf_path"]),
        )
        results["steps_completed"].append("visualization")
        results["steps"].append("visualization")
        results["outputs"]["plots"] = plot_results

        results["status"] = "completed"
        logger.info("GWAS workflow completed successfully")

    except Exception as e:
        logger.error(f"GWAS workflow failed: {e}")
        results["success"] = False
        results["status"] = "failed"
        results["error"] = str(e)
        results["errors"].append(str(e))

    return results


def run_gwas(
    vcf_path: Union[str, Path],
    phenotype_path: Union[str, Path],
    config: Dict[str, Any],
    output_dir: Union[str, Path] | None = None,
) -> Dict[str, Any]:
    """Run complete GWAS workflow.

    Supports linear, logistic, and mixed model association testing.
    Integrates LD pruning, haplodiploidy checking, summary statistics output,
    and SNP-to-gene annotation.

    Args:
        vcf_path: Path to VCF file
        phenotype_path: Path to phenotype file
        config: GWAS configuration dictionary (flat or nested YAML format)
        output_dir: Output directory (optional)

    Returns:
        Dictionary with GWAS results and metadata

    Raises:
        ValueError: If required parameters are missing
        FileNotFoundError: If input files don't exist
    """
    import tempfile
    from pathlib import Path

    logger.info("Starting GWAS workflow")

    # Normalize nested YAML config to flat format
    config = _normalize_config(config)

    # Validate inputs
    vcf_path = Path(vcf_path)
    phenotype_path = Path(phenotype_path)

    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    if not phenotype_path.exists():
        raise FileNotFoundError(f"Phenotype file not found: {phenotype_path}")

    # Create output directory
    if output_dir is None:
        temp_dir = tempfile.mkdtemp()
        output_dir = Path(temp_dir)
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    results: Dict[str, Any] = {"config": config, "output_dir": str(output_dir), "steps_completed": [], "results": {}}

    try:
        # Step 1: Load and parse VCF
        logger.info("Step 1: Parsing VCF file")
        vcf_data = parse_vcf_full(vcf_path)
        results["steps_completed"].append("parse_vcf")
        results["results"]["vcf_summary"] = {
            "num_variants": len(vcf_data.get("variants", [])),
            "num_samples": len(vcf_data.get("samples", [])),
        }

        # Step 1b: Sample subsetting (optional)
        sample_list = config.get("sample_list")
        sample_subset = config.get("sample_subset")
        if sample_list or sample_subset:
            from metainformant.gwas.analysis.quality import subset_vcf_data

            meta_for_subset = None
            if sample_subset:
                metadata_path_sub = config.get("metadata_file") or (
                    config.get("samples", {}).get("metadata_file") if isinstance(config.get("samples"), dict) else None
                )
                if metadata_path_sub and Path(metadata_path_sub).exists():
                    meta_result_sub = load_sample_metadata(metadata_path_sub)
                    if meta_result_sub.get("status") == "success":
                        meta_for_subset = meta_result_sub.get("metadata")

            vcf_data = subset_vcf_data(
                vcf_data,
                sample_list_file=sample_list,
                subset_config=sample_subset,
                metadata=meta_for_subset,
            )
            results["steps_completed"].append("sample_subsetting")
            results["results"]["vcf_summary"]["num_samples_after_subset"] = len(vcf_data.get("samples", []))
            logger.info(f"After subsetting: {len(vcf_data.get('samples', []))} samples")

        # Step 2: Apply QC filters
        logger.info("Step 2: Applying QC filters")
        qc_config = config.get("quality_control", {})
        qc_result = apply_qc_filters(vcf_data, qc_config)
        filtered_data = qc_result.get("filtered_data", vcf_data)
        results["steps_completed"].append("qc_filters")
        results["results"]["qc_summary"] = {
            "variants_before_qc": qc_result.get("num_variants_before", 0),
            "variants_after_qc": qc_result.get("num_variants_after", len(filtered_data.get("variants", []))),
            "samples_after_qc": len(filtered_data.get("samples", [])),
        }

        # Step 2b: Haplodiploidy check (optional, for haplodiploid species)
        haplodiploidy_config = config.get("haplodiploidy", {})
        if isinstance(haplodiploidy_config, dict) and haplodiploidy_config.get("enabled", False):
            logger.info("Step 2b: Checking haplodiploidy")
            het_threshold = haplodiploidy_config.get("het_threshold", 0.05)
            haplo_result = check_haplodiploidy(filtered_data, het_threshold=het_threshold)
            results["results"]["haplodiploidy"] = haplo_result
            results["steps_completed"].append("haplodiploidy_check")

            # Optionally exclude haploid samples
            if haplodiploidy_config.get("exclude_haploid", False) and haplo_result["haploid_samples"]:
                diploid_indices = haplo_result["diploid_samples"]
                genotypes_by_sample = filtered_data.get("genotypes", [])
                if genotypes_by_sample and diploid_indices:
                    filtered_data = dict(filtered_data)
                    filtered_data["genotypes"] = [genotypes_by_sample[i] for i in diploid_indices]
                    old_samples = filtered_data.get("samples", [])
                    if old_samples:
                        filtered_data["samples"] = [old_samples[i] for i in diploid_indices]
                    logger.info(f"Excluded {haplo_result['n_haploid']} haploid samples")

        # Get genotypes (samples x variants format)
        genotypes_by_sample = filtered_data.get("genotypes", [])

        # Transpose to variants x samples for analysis functions
        if genotypes_by_sample and genotypes_by_sample[0]:
            n_samples = len(genotypes_by_sample)
            n_variants = len(genotypes_by_sample[0])
            genotypes_by_variant = [[genotypes_by_sample[s][v] for s in range(n_samples)] for v in range(n_variants)]
        else:
            genotypes_by_variant = []
            n_samples = 0
            n_variants = 0

        # Step 3: LD pruning (optional, for PCA)
        ld_config = config.get("ld_pruning", {})
        ld_pruned_indices = None
        if isinstance(ld_config, dict) and ld_config.get("enabled", False) and genotypes_by_variant:
            logger.info("Step 3: LD pruning before PCA")
            ld_pruned_indices = ld_prune(
                genotypes_by_variant,
                window_size=ld_config.get("window_size", 50),
                step_size=ld_config.get("step_size", 5),
                r2_threshold=ld_config.get("r2_threshold", 0.2),
            )
            ld_pruned_genotypes = [genotypes_by_variant[i] for i in ld_pruned_indices]
            results["steps_completed"].append("ld_pruning")
            results["results"]["ld_pruning"] = {
                "variants_before": n_variants,
                "variants_after": len(ld_pruned_indices),
            }
        else:
            ld_pruned_genotypes = genotypes_by_variant

        # Step 4: Population structure analysis (PCA on LD-pruned, kinship on all)
        logger.info("Step 4: Computing population structure")
        kinship_result = None
        if genotypes_by_sample:
            # PCA on LD-pruned genotypes (transposed back to samples x variants)
            if ld_pruned_genotypes:
                pca_genotypes_by_sample = [
                    [ld_pruned_genotypes[v][s] for v in range(len(ld_pruned_genotypes))] for s in range(n_samples)
                ]
            else:
                pca_genotypes_by_sample = genotypes_by_sample

            pca_result = compute_pca(pca_genotypes_by_sample, n_components=min(10, n_samples))

            # Kinship on all genotypes (not LD-pruned)
            kinship_result = compute_kinship_matrix(genotypes_by_sample)

            results["steps_completed"].append("population_structure")
            results["results"]["pca"] = pca_result
            results["results"]["kinship"] = kinship_result

        # Step 5: Load phenotypes
        logger.info("Step 5: Loading phenotypes")
        trait_name = config.get("trait", config.get("trait_name"))
        phenotypes = _load_phenotypes(phenotype_path, trait=trait_name)
        results["steps_completed"].append("load_phenotypes")
        results["results"]["phenotype_summary"] = {
            "num_samples": len(phenotypes),
            "trait_name": trait_name or "unknown",
        }

        # Step 5b: Load sample metadata (optional)
        metadata = None
        metadata_path = (
            config.get("metadata_file") or config.get("samples", {}).get("metadata_file")
            if isinstance(config.get("samples"), dict)
            else None
        )
        if metadata_path and Path(metadata_path).exists():
            logger.info("Step 5b: Loading sample metadata")
            try:
                meta_result = load_sample_metadata(metadata_path)
                if meta_result.get("status") == "success":
                    metadata = meta_result.get("metadata", {})
                    sample_ids = filtered_data.get("samples", [])
                    if sample_ids:
                        validation = validate_metadata(metadata, sample_ids)
                        if validation.get("missing_samples"):
                            logger.warning(f"Metadata missing for {len(validation['missing_samples'])} samples")
                    results["steps_completed"].append("load_metadata")
                    results["results"]["metadata_summary"] = {
                        "n_samples_with_metadata": meta_result.get("n_samples", 0),
                        "columns": meta_result.get("columns", []),
                    }
            except Exception as e:
                logger.warning(f"Metadata loading failed: {e}")

        # Align phenotypes to genotyped samples by ID when possible
        sample_ids = filtered_data.get("samples", [])
        pheno_by_id = _load_phenotypes_by_id(phenotype_path, trait=trait_name)
        if pheno_by_id and sample_ids:
            # ID-based alignment: match phenotype to each genotyped sample
            aligned_phenotypes = []
            for sid in sample_ids:
                if sid in pheno_by_id:
                    aligned_phenotypes.append(pheno_by_id[sid])
            if aligned_phenotypes:
                phenotypes = aligned_phenotypes
                logger.info(
                    f"ID-based phenotype alignment: {len(aligned_phenotypes)}/{len(sample_ids)} samples matched"
                )

        # Fallback: positional alignment
        min_samples = min(len(phenotypes), n_samples) if genotypes_by_variant else 0
        if min_samples > 0 and len(phenotypes) != n_samples:
            logger.warning(
                f"Sample count mismatch: {n_samples} genotyped, {len(phenotypes)} phenotyped. "
                f"Using first {min_samples} samples."
            )
            phenotypes = phenotypes[:min_samples]

        # Step 6: Association testing
        model = config.get("model", "linear")
        logger.info(f"Step 6: Running {model} association tests")
        assoc_results = []

        if genotypes_by_variant and phenotypes:
            if model == "mixed":
                # Mixed model GWAS
                if kinship_result and kinship_result.get("status") == "success":
                    kinship_matrix = kinship_result["kinship_matrix"]
                    assoc_results = run_mixed_model_gwas(
                        genotype_matrix=genotypes_by_variant,
                        phenotypes=phenotypes[:min_samples],
                        kinship_matrix=kinship_matrix,
                        variant_info=filtered_data.get("variants"),
                    )
                else:
                    logger.warning("Kinship matrix unavailable, falling back to linear model")
                    model = "linear"

            if model in ("linear", "logistic"):
                # Linear/logistic model: test each variant
                from metainformant.gwas.analysis.association import association_test_logistic

                for i, genotype in enumerate(genotypes_by_variant):
                    try:
                        geno_trimmed = genotype[:min_samples]
                        pheno_trimmed = phenotypes[:min_samples]

                        if model == "logistic":
                            result = association_test_logistic(geno_trimmed, [int(p) for p in pheno_trimmed])
                        else:
                            result = association_test_linear(geno_trimmed, pheno_trimmed)

                        result["variant_id"] = f"variant_{i}"
                        if filtered_data.get("variants") and i < len(filtered_data["variants"]):
                            vinfo = filtered_data["variants"][i]
                            result["chrom"] = vinfo.get("chrom", "")
                            result["pos"] = vinfo.get("pos", 0)
                        assoc_results.append(result)
                    except Exception as e:
                        logger.warning(f"Association test failed for variant {i}: {e}")

        results["steps_completed"].append("association_testing")
        results["results"]["association_results"] = assoc_results

        # Compute lambda GC for inflation assessment
        if assoc_results:
            import math

            valid_p = [
                r.get("p_value", 1.0)
                for r in assoc_results
                if r.get("p_value") is not None and 0 < r.get("p_value", 1.0) <= 1.0
            ]
            if valid_p:
                valid_p_sorted = sorted(valid_p)
                median_p = valid_p_sorted[len(valid_p_sorted) // 2]
                # Lambda GC = median(chi2) / 0.456 where chi2 = qchisq(1-p, 1)
                # Approximation: -2 * log(median_p) / 1.386 (median of chi2(1))
                if median_p > 0:
                    lambda_gc = -2.0 * math.log(median_p) / 1.386
                    results["results"]["lambda_gc"] = round(lambda_gc, 4)
                    if lambda_gc > 1.1:
                        logger.warning(f"Genomic inflation detected: lambda_GC = {lambda_gc:.3f} (>1.1)")
                    elif lambda_gc < 0.9:
                        logger.warning(f"Genomic deflation detected: lambda_GC = {lambda_gc:.3f} (<0.9)")
                    else:
                        logger.info(f"Lambda GC = {lambda_gc:.3f} (within expected range)")

        # Step 7: Multiple testing correction
        logger.info("Step 7: Applying multiple testing correction")
        if assoc_results:
            p_values = [r.get("p_value", 1.0) for r in assoc_results]
            fdr_result = fdr_correction(p_values)

            if isinstance(fdr_result, dict):
                q_values = fdr_result.get("adjusted_p_values", [])
            else:
                _, q_values = fdr_result

            for i, q_val in enumerate(q_values):
                if i < len(assoc_results):
                    assoc_results[i]["q_value"] = q_val

            bonf_result = bonferroni_correction(p_values)
            if isinstance(bonf_result, dict):
                for i, is_sig in enumerate(bonf_result.get("significant", [])):
                    if i < len(assoc_results):
                        assoc_results[i]["bonferroni_significant"] = is_sig

        results["steps_completed"].append("multiple_testing_correction")

        # Step 7b: Heritability estimation (optional)
        if kinship_result and kinship_result.get("status") == "success" and phenotypes:
            logger.info("Step 7b: Estimating SNP heritability")
            try:
                km = kinship_result["kinship_matrix"]
                km_size = len(km) if isinstance(km, list) else (km.shape[0] if hasattr(km, "shape") else 0)
                pheno_for_h2 = phenotypes[:km_size] if len(phenotypes) > km_size else phenotypes
                logger.info(
                    f"Heritability inputs: kinship {km_size}x{km_size}, "
                    f"phenotypes {len(pheno_for_h2)} (from {len(phenotypes)} total)"
                )
                h2_result = estimate_heritability(km, pheno_for_h2)
                if h2_result.get("status") == "success":
                    results["steps_completed"].append("heritability")
                    results["results"]["heritability"] = {
                        "h2": h2_result.get("h2"),
                        "h2_se": h2_result.get("h2_se"),
                        "sigma_g": h2_result.get("sigma_g"),
                        "sigma_e": h2_result.get("sigma_e"),
                    }
                    logger.info(
                        f"SNP heritability: h2 = {h2_result.get('h2', 0):.3f} (SE = {h2_result.get('h2_se', 0):.3f})"
                    )
                else:
                    logger.warning(f"Heritability estimation returned non-success: {h2_result}")
            except Exception as e:
                logger.warning(f"Heritability estimation failed: {e}", exc_info=True)

        # Step 8: Write summary statistics
        logger.info("Step 8: Writing summary statistics")
        if assoc_results and filtered_data.get("variants"):
            variant_info = filtered_data["variants"]
            try:
                stats_path = output_dir / "summary_statistics.tsv"
                write_summary_statistics(assoc_results, variant_info, stats_path)
                results["steps_completed"].append("summary_statistics")
                results["results"]["summary_stats_path"] = str(stats_path)

                # Write significant hits
                sig_path = output_dir / "significant_hits.tsv"
                threshold = config.get("significance_threshold", 5e-8)
                write_significant_hits(assoc_results, variant_info, sig_path, threshold=threshold)

                # Write JSON summary
                summary_path = output_dir / "results_summary.json"
                summary = create_results_summary(assoc_results, summary_path, threshold=threshold)
                results["results"]["summary"] = summary
            except Exception as e:
                logger.warning(f"Summary statistics writing failed: {e}")

        # Step 9: SNP-to-gene annotation (optional)
        annotation_config = config.get("annotation", {})
        if isinstance(annotation_config, dict) and annotation_config.get("enabled", False):
            gff_path = annotation_config.get("gff3_file")
            if gff_path and Path(gff_path).exists():
                logger.info("Step 9: Annotating variants with genes")
                try:
                    window_kb = annotation_config.get("window_kb", 50)
                    annotate_variants_with_genes(assoc_results, gff_path, window_kb=window_kb)
                    results["steps_completed"].append("annotation")
                except Exception as e:
                    logger.warning(f"Variant annotation failed: {e}")
            else:
                logger.info("Step 9: Skipping annotation (GFF3 file not found)")

        # Step 9b: Fine-mapping credible sets (optional)
        finemapping_config = config.get("finemapping", {})
        if isinstance(finemapping_config, dict) and finemapping_config.get("enabled", False) and assoc_results:
            logger.info("Step 9b: Computing fine-mapping credible sets")
            try:
                credible_level = finemapping_config.get("credible_level", 0.95)
                cs_result = compute_credible_set(assoc_results, credible_level=credible_level)
                if cs_result.get("status") == "success":
                    results["steps_completed"].append("fine_mapping")
                    results["results"]["fine_mapping"] = {
                        "credible_set_size": cs_result.get("credible_set_size"),
                        "credible_level": credible_level,
                    }
            except Exception as e:
                logger.warning(f"Fine-mapping failed: {e}")

        # Apply visualization style from config
        try:
            style = style_from_config(config)
            apply_style(style)
        except Exception:
            pass  # Non-critical, use defaults

        # Step 10: Generate plots
        logger.info("Step 10: Generating visualization plots")
        try:
            from metainformant.gwas.visualization.interactive.suite import generate_all_plots as gen_plots

            plot_results = gen_plots(
                association_results=assoc_results,
                output_dir=output_dir,
                significance_threshold=config.get("significance_threshold", 5e-8),
            )
            results["steps_completed"].append("visualization")
            results["results"]["plots"] = plot_results
        except Exception as e:
            logger.warning(f"Plot generation failed: {e}")

        # Step 11: Enhanced visualizations (using new modules)
        logger.info("Step 11: Generating enhanced visualizations")
        enhanced_panels = []
        try:
            from metainformant.gwas.visualization.interactive.composite import (
                gwas_summary_panel,
                population_structure_panel,
            )

            # GWAS summary panel
            pca_data_for_viz = results["results"].get("pca")
            kinship_for_viz = kinship_result.get("kinship_matrix") if kinship_result else None
            logger.info(
                f"Enhanced viz inputs: assoc_results={len(assoc_results)}, "
                f"pca={'yes' if pca_data_for_viz else 'no'}, "
                f"kinship={'yes' if kinship_for_viz is not None else 'no'}, "
                f"metadata={'yes' if metadata else 'no'}"
            )
            try:
                summary_panel = gwas_summary_panel(
                    assoc_results,
                    pca_data=pca_data_for_viz,
                    kinship_matrix=kinship_for_viz,
                    output_file=output_dir / "gwas_summary_panel.png",
                )
                enhanced_panels.append("gwas_summary_panel")
                logger.info(f"GWAS summary panel: {summary_panel.get('status')}")
            except Exception as e:
                logger.warning(f"GWAS summary panel failed: {e}", exc_info=True)

            # Population structure panel (if PCA and kinship available)
            if pca_data_for_viz and kinship_for_viz is not None:
                try:
                    pop_panel = population_structure_panel(
                        pca_data_for_viz,
                        kinship_for_viz,
                        metadata=metadata,
                        output_file=output_dir / "population_structure_panel.png",
                    )
                    enhanced_panels.append("population_structure_panel")
                    logger.info(f"Population structure panel: {pop_panel.get('status')}")
                except Exception as e:
                    logger.warning(f"Population structure panel failed: {e}", exc_info=True)

            # Heritability bar chart (if per-chromosome available)
            h2_data = results["results"].get("heritability")
            if h2_data and h2_data.get("per_chromosome"):
                try:
                    from metainformant.gwas.analysis.heritability import heritability_bar_chart

                    heritability_bar_chart(h2_data, output_file=output_dir / "heritability_bar.png")
                    enhanced_panels.append("heritability_bar")
                    logger.info("Heritability bar chart generated")
                except Exception as e:
                    logger.warning(f"Heritability bar chart failed: {e}", exc_info=True)

            if enhanced_panels:
                results["steps_completed"].append("enhanced_visualization")
                results["results"]["enhanced_panels"] = enhanced_panels
                logger.info(f"Enhanced visualization: {len(enhanced_panels)} panels generated")
            else:
                logger.warning("No enhanced visualization panels were generated")
        except ImportError as e:
            logger.warning(f"Enhanced visualization import failed: {e}")
        except Exception as e:
            logger.warning(f"Enhanced visualization failed: {e}", exc_info=True)

        results["status"] = "success"
        logger.info("GWAS workflow completed successfully")

    except Exception as e:
        results["status"] = "failed"
        results["error"] = str(e)
        logger.error(f"GWAS workflow failed: {e}")

    return results


def run_multi_trait_gwas(
    vcf_path: Union[str, Path],
    phenotype_path: Union[str, Path],
    traits: List[str],
    config: Dict[str, Any],
    output_dir: Union[str, Path] | None = None,
) -> Dict[str, Any]:
    """Run GWAS for multiple traits, sharing VCF parse/QC/PCA/kinship.

    VCF parsed once -> QC once -> PCA once -> kinship once ->
    for each trait: load phenotype column -> association -> correction -> output.

    Args:
        vcf_path: Path to VCF file.
        phenotype_path: Path to phenotype file (TSV with multiple trait columns).
        traits: List of trait column names to analyze.
        config: GWAS configuration dictionary.
        output_dir: Base output directory; each trait gets a subdirectory.

    Returns:
        Combined results dict with per-trait results.
    """
    import tempfile

    logger.info(f"Starting multi-trait GWAS: {len(traits)} traits")

    config = _normalize_config(config)
    vcf_path = Path(vcf_path)
    phenotype_path = Path(phenotype_path)

    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    if not phenotype_path.exists():
        raise FileNotFoundError(f"Phenotype file not found: {phenotype_path}")

    if output_dir is None:
        output_dir = Path(tempfile.mkdtemp())
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    results: Dict[str, Any] = {
        "status": "running",
        "traits": traits,
        "output_dir": str(output_dir),
        "steps_completed": [],
        "trait_results": {},
    }

    try:
        # Shared step 1: Parse VCF
        logger.info("Multi-trait: Parsing VCF (shared)")
        vcf_data = parse_vcf_full(vcf_path)
        results["steps_completed"].append("parse_vcf")

        # Shared step 2: QC
        logger.info("Multi-trait: Applying QC (shared)")
        qc_config = config.get("quality_control", {})
        qc_result = apply_qc_filters(vcf_data, qc_config)
        filtered_data = qc_result.get("filtered_data", vcf_data)
        results["steps_completed"].append("qc_filters")

        # Shared step 2b: Haplodiploidy
        haplodiploidy_config = config.get("haplodiploidy", {})
        if isinstance(haplodiploidy_config, dict) and haplodiploidy_config.get("enabled", False):
            het_threshold = haplodiploidy_config.get("het_threshold", 0.05)
            haplo_result = check_haplodiploidy(filtered_data, het_threshold=het_threshold)
            if haplodiploidy_config.get("exclude_haploid", False) and haplo_result["haploid_samples"]:
                diploid_indices = haplo_result["diploid_samples"]
                genotypes_by_sample = filtered_data.get("genotypes", [])
                if genotypes_by_sample and diploid_indices:
                    filtered_data = dict(filtered_data)
                    filtered_data["genotypes"] = [genotypes_by_sample[i] for i in diploid_indices]
                    old_samples = filtered_data.get("samples", [])
                    if old_samples:
                        filtered_data["samples"] = [old_samples[i] for i in diploid_indices]
            results["steps_completed"].append("haplodiploidy_check")

        genotypes_by_sample = filtered_data.get("genotypes", [])
        n_samples = len(genotypes_by_sample)
        n_variants = len(genotypes_by_sample[0]) if genotypes_by_sample else 0
        if genotypes_by_sample and genotypes_by_sample[0]:
            genotypes_by_variant = [[genotypes_by_sample[s][v] for s in range(n_samples)] for v in range(n_variants)]
        else:
            genotypes_by_variant = []

        # Shared step 3: LD pruning
        ld_config = config.get("ld_pruning", {})
        if isinstance(ld_config, dict) and ld_config.get("enabled", False) and genotypes_by_variant:
            ld_pruned_indices = ld_prune(
                genotypes_by_variant,
                window_size=ld_config.get("window_size", 50),
                step_size=ld_config.get("step_size", 5),
                r2_threshold=ld_config.get("r2_threshold", 0.2),
            )
            ld_pruned_genotypes = [genotypes_by_variant[i] for i in ld_pruned_indices]
            results["steps_completed"].append("ld_pruning")
        else:
            ld_pruned_genotypes = genotypes_by_variant

        # Shared step 4: PCA + kinship
        kinship_result = None
        pca_result = None
        if genotypes_by_sample:
            if ld_pruned_genotypes:
                pca_geno = [
                    [ld_pruned_genotypes[v][s] for v in range(len(ld_pruned_genotypes))] for s in range(n_samples)
                ]
            else:
                pca_geno = genotypes_by_sample
            pca_result = compute_pca(pca_geno, n_components=min(10, n_samples))
            kinship_result = compute_kinship_matrix(genotypes_by_sample)
            results["steps_completed"].append("population_structure")

        # Per-trait loop
        for trait_name in traits:
            logger.info(f"Multi-trait: Running trait '{trait_name}'")
            trait_dir = output_dir / trait_name
            trait_dir.mkdir(parents=True, exist_ok=True)

            trait_result: Dict[str, Any] = {"trait": trait_name, "steps": []}

            try:
                # Load phenotypes for this trait
                pheno_by_id = _load_phenotypes_by_id(phenotype_path, trait=trait_name)
                sample_ids = filtered_data.get("samples", [])
                phenotypes = [pheno_by_id[sid] for sid in sample_ids if sid in pheno_by_id]
                min_s = min(len(phenotypes), n_samples) if genotypes_by_variant else 0
                phenotypes = phenotypes[:min_s]

                # Association testing
                model = config.get("model", "linear")
                assoc_results = []
                if genotypes_by_variant and phenotypes:
                    if model == "mixed" and kinship_result and kinship_result.get("status") == "success":
                        assoc_results = run_mixed_model_gwas(
                            genotype_matrix=genotypes_by_variant,
                            phenotypes=phenotypes[:min_s],
                            kinship_matrix=kinship_result["kinship_matrix"],
                            variant_info=filtered_data.get("variants"),
                        )
                    else:
                        for i, genotype in enumerate(genotypes_by_variant):
                            try:
                                result = association_test_linear(genotype[:min_s], phenotypes[:min_s])
                                result["variant_id"] = f"variant_{i}"
                                if filtered_data.get("variants") and i < len(filtered_data["variants"]):
                                    vinfo = filtered_data["variants"][i]
                                    result["chrom"] = vinfo.get("chrom", "")
                                    result["pos"] = vinfo.get("pos", 0)
                                assoc_results.append(result)
                            except Exception:
                                continue

                trait_result["steps"].append("association_testing")
                trait_result["n_tests"] = len(assoc_results)

                # Multiple testing correction
                if assoc_results:
                    p_values = [r.get("p_value", 1.0) for r in assoc_results]
                    fdr_result = fdr_correction(p_values)
                    if isinstance(fdr_result, dict):
                        q_values = fdr_result.get("adjusted_p_values", [])
                    else:
                        _, q_values = fdr_result
                    for i, q_val in enumerate(q_values):
                        if i < len(assoc_results):
                            assoc_results[i]["q_value"] = q_val
                    trait_result["steps"].append("correction")

                # Summary stats
                if assoc_results and filtered_data.get("variants"):
                    stats_path = trait_dir / "summary_statistics.tsv"
                    write_summary_statistics(assoc_results, filtered_data["variants"], stats_path)
                    trait_result["steps"].append("summary_statistics")

                trait_result["status"] = "success"
                trait_result["association_results"] = assoc_results

            except Exception as e:
                trait_result["status"] = "failed"
                trait_result["error"] = str(e)
                logger.warning(f"Trait '{trait_name}' failed: {e}")

            results["trait_results"][trait_name] = trait_result

        results["status"] = "success"
        # Provide a combined results dict compatible with visualization
        if results["trait_results"]:
            first_trait = list(results["trait_results"].values())[0]
            results["results"] = {
                "association_results": first_trait.get("association_results", []),
                "pca": pca_result,
                "kinship": kinship_result,
                "vcf_summary": {
                    "num_variants": n_variants,
                    "num_samples": n_samples,
                },
            }

    except Exception as e:
        results["status"] = "failed"
        results["error"] = str(e)
        logger.error(f"Multi-trait GWAS failed: {e}")

    return results
