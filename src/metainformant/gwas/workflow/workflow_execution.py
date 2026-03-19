"""GWAS workflow execution functions.

This module provides the main workflow execution functions:
- execute_gwas_workflow: Execute the complete GWAS workflow from config
- run_gwas: Run complete GWAS with VCF/phenotype paths and config
- run_multi_trait_gwas: Run GWAS for multiple traits sharing common steps
"""

from __future__ import annotations

import json
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
from metainformant.gwas.analysis.correction import bonferroni_correction, fdr_correction, genomic_control
from metainformant.gwas.analysis.heritability import estimate_heritability, heritability_bar_chart, partition_heritability_by_chromosome
from metainformant.gwas.analysis.ld_pruning import ld_prune
from metainformant.gwas.analysis.mixed_model import association_test_mixed, run_mixed_model_gwas
from metainformant.gwas.analysis.quality import apply_qc_filters, check_haplodiploidy, parse_vcf_full
from metainformant.gwas.analysis.structure import compute_kinship_matrix, compute_pca
from metainformant.gwas.analysis.summary_stats import (
    create_results_summary,
    write_significant_hits,
    write_summary_statistics,
)
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

        # Step 2.5: LD pruning (remove correlated variants before PCA)
        logger.info("Step 2.5: LD pruning")
        genotype_matrix = _extract_genotype_matrix(filtered_data)
        variants_info = filtered_data.get("variants", [])

        ld_config = config.get("ld_pruning", {})
        r2_thresh = ld_config.get("r2_threshold", 0.2)
        window = ld_config.get("window_size", 50)
        step = ld_config.get("step_size", 5)

        if genotype_matrix and len(genotype_matrix) > 1:
            kept_indices = ld_prune(
                genotype_matrix,
                r2_threshold=r2_thresh,
                window_size=window,
                step_size=step,
            )
            genotype_matrix = [genotype_matrix[i] for i in kept_indices]
            variants_info = [variants_info[i] for i in kept_indices] if variants_info else []
            logger.info(f"LD pruning: {len(kept_indices)} variants retained for structure analysis")
        else:
            logger.info("LD pruning: skipped (insufficient variants)")

        results["steps_completed"].append("ld_pruning")
        results["steps"].append("ld_pruning")

        # Step 3: Population structure analysis
        logger.info("Step 3: Analyzing population structure")
        n_pcs_config = config.get("structure", {}).get("pca_components", 10)
        kinship_method = config.get("structure", {}).get("kinship_method", "vanraden")

        # genotype_matrix is variant-major (variants × samples) from _extract_genotype_matrix.
        # compute_pca and compute_kinship_matrix expect sample-major (samples × variants).
        n_variants_gm = len(genotype_matrix)
        n_samples_gm = len(genotype_matrix[0]) if genotype_matrix else 0
        gm_sample_major = [
            [genotype_matrix[v][s] for v in range(n_variants_gm)]
            for s in range(n_samples_gm)
        ]
        logger.info(f"Transposed genotype matrix: {n_samples_gm} samples × {n_variants_gm} variants (sample-major)")

        pca_result = compute_pca(gm_sample_major, n_components=n_pcs_config)
        kinship_result = compute_kinship_matrix(gm_sample_major, method=kinship_method)

        # Extract the actual kinship matrix from the result dict
        # compute_kinship_matrix returns {"status": ..., "kinship_matrix": [[...]], ...}
        raw_kinship = kinship_result.get("kinship_matrix", [])
        if hasattr(raw_kinship, "tolist"):
            kinship_list = raw_kinship.tolist()
        elif isinstance(raw_kinship, list):
            kinship_list = raw_kinship
        else:
            kinship_list = []

        # Store full kinship metadata for output
        results["outputs"]["kinship_metadata"] = {
            k: v for k, v in kinship_result.items() if k != "kinship_matrix"
        }

        results["steps_completed"].append("population_structure")
        results["steps"].append("population_structure")

        # Step 4: Association testing (model-aware dispatch)
        model = config.get("model", "linear")
        logger.info(f"Step 4: Performing association testing (model={model})")
        trait_name = config.get("trait", config.get("trait_name", config.get("default_trait")))
        traits = _load_phenotypes(config["phenotype_path"], trait=trait_name)

        # PCA results is a dict with 'pcs' list
        pcs = pca_result.get("pcs", [])

        association_results = []
        if genotype_matrix and traits:
            n_samples = len(genotype_matrix[0])
            actual_traits = traits[:n_samples]
            actual_pcs = pcs[:n_samples]

            # Prepare covariates from PCA
            n_avail_pcs = len(actual_pcs[0]) if actual_pcs and actual_pcs[0] else 0
            n_pcs_to_use = min(n_avail_pcs, 5, max(0, n_samples - 2))
            if n_pcs_to_use > 0:
                variant_covariates = [
                    [actual_pcs[s][pc] for s in range(len(actual_pcs))]
                    for pc in range(n_pcs_to_use)
                ]
            else:
                variant_covariates = None

            # Import logistic if needed
            if model == "logistic":
                from metainformant.gwas.analysis.association import association_test_logistic

            for i, genotype in enumerate(genotype_matrix):
                try:
                    if model == "mixed":
                        result = association_test_mixed(
                            genotype,
                            actual_traits,
                            kinship_matrix=kinship_list,
                            covariates=variant_covariates,
                        )
                    elif model == "logistic":
                        result = association_test_logistic(
                            genotype,
                            [int(t) for t in actual_traits],
                            covariates=variant_covariates,
                        )
                    else:  # "linear" or default
                        result = association_test_linear(
                            genotype,
                            actual_traits,
                            covariates=variant_covariates,
                        )
                except Exception as exc:
                    # Fall back to linear on model-specific failures
                    logger.debug(f"Variant {i}: {model} model failed ({exc}), falling back to linear")
                    result = association_test_linear(
                        genotype,
                        actual_traits,
                        covariates=variant_covariates,
                    )

                # Attach variant info to result for plotter
                if i < len(variants_info):
                    vinfo = variants_info[i]
                    result["chrom"] = vinfo.get("chrom", "1")
                    result["pos"] = vinfo.get("pos", 0)
                    result["variant_id"] = vinfo.get("id", f"var_{i}")

                # Compute MAF from genotype vector
                valid_genos = [g for g in genotype if g >= 0]
                if valid_genos:
                    alt_freq = sum(valid_genos) / (2.0 * len(valid_genos))
                    result["maf"] = min(alt_freq, 1.0 - alt_freq)
                else:
                    result["maf"] = 0.0

                association_results.append(result)

            logger.info(f"Association testing complete: {len(association_results)} variants tested with {model} model")

        results["steps_completed"].append("association_testing")
        results["steps"].append("association_testing")
        results["outputs"]["association_results"] = association_results

        # Step 5: Multiple testing correction & genomic control
        logger.info("Step 5: Multiple testing correction")
        p_values = [
            r.get("p_value", r.get("pval", 1.0))
            for r in association_results
            if r.get("p_value", r.get("pval")) is not None
        ]

        if p_values:
            # Bonferroni
            bonf_result = bonferroni_correction(p_values)
            n_bonf = bonf_result.get("n_significant", 0) if isinstance(bonf_result, dict) else 0
            logger.info(f"Bonferroni correction: {n_bonf} variants significant")

            # FDR (Benjamini-Hochberg)
            fdr_result = fdr_correction(p_values)
            n_fdr = fdr_result.get("n_significant", 0) if isinstance(fdr_result, dict) else 0
            logger.info(f"FDR correction (BH): {n_fdr} variants significant")

            # Genomic control lambda
            gc_result = genomic_control(p_values=p_values)
            lambda_gc = gc_result.get("lambda_gc", 1.0) if isinstance(gc_result, dict) else 1.0
            logger.info(f"Genomic control: λ_GC = {lambda_gc:.4f}")

            # Attach correction results to output
            results["outputs"]["bonferroni"] = bonf_result if isinstance(bonf_result, dict) else {}
            results["outputs"]["fdr"] = fdr_result if isinstance(fdr_result, dict) else {}
            results["outputs"]["genomic_control"] = {"lambda_gc": lambda_gc}

            # Mark significant variants in association results
            bonf_significant = bonf_result.get("significant", []) if isinstance(bonf_result, dict) else []
            fdr_significant = fdr_result.get("significant", []) if isinstance(fdr_result, dict) else []
            fdr_adjusted = fdr_result.get("adjusted_p_values", []) if isinstance(fdr_result, dict) else []
            for idx, ar in enumerate(association_results):
                if idx < len(bonf_significant):
                    ar["bonferroni_significant"] = bonf_significant[idx]
                if idx < len(fdr_significant):
                    ar["fdr_significant"] = fdr_significant[idx]
                if idx < len(fdr_adjusted):
                    ar["fdr_p_value"] = fdr_adjusted[idx]

        results["steps_completed"].append("multiple_testing_correction")
        results["steps"].append("multiple_testing_correction")

        # Step 5.5: Heritability estimation (REML)
        logger.info("Step 5.5: Estimating SNP heritability (REML)")
        if kinship_list and traits:
            actual_traits_h2 = traits[:len(kinship_list)]
            h2_result = estimate_heritability(kinship_list, actual_traits_h2)
            if h2_result.get("status") == "success":
                logger.info(
                    f"Heritability: h² = {h2_result['h2']:.4f} ± {h2_result['h2_se']:.4f} "
                    f"(σ²_g={h2_result['sigma_g']:.4f}, σ²_e={h2_result['sigma_e']:.4f})"
                )
            else:
                logger.warning(f"Heritability estimation: {h2_result.get('message', 'failed')}")
            results["outputs"]["heritability"] = h2_result

        # Step 5.5a: Per-chromosome heritability partitioning
        logger.info("Step 5.5a: Partitioning heritability by chromosome")
        if genotype_matrix and variants_info and traits:
            try:
                # Build per-chromosome kinship matrices (sample-major)
                chrom_indices: dict[str, list[int]] = {}
                for idx, vinfo in enumerate(variants_info):
                    chrom = vinfo.get("chrom", "unknown")
                    chrom_indices.setdefault(chrom, []).append(idx)

                n_samples_h2 = len(genotype_matrix[0]) if genotype_matrix else 0
                per_chrom_kinship: dict[int, list[list[float]]] = {}
                for chrom_id, (chrom_name, indices) in enumerate(sorted(chrom_indices.items())):
                    if len(indices) < 2:
                        continue
                    # Extract variant-major genotypes for this chromosome
                    chrom_geno = [genotype_matrix[i] for i in indices]
                    # Transpose to sample-major for kinship
                    n_vars_chr = len(chrom_geno)
                    chrom_sm = [
                        [chrom_geno[v][s] for v in range(n_vars_chr)]
                        for s in range(n_samples_h2)
                    ]
                    k_result = compute_kinship_matrix(chrom_sm, method=config.get("structure", {}).get("kinship_method", "vanraden"))
                    if k_result.get("status") == "success":
                        per_chrom_kinship[chrom_id] = k_result["kinship_matrix"]

                if per_chrom_kinship:
                    actual_traits_part = traits[:n_samples_h2]
                    partition_result = partition_heritability_by_chromosome(per_chrom_kinship, actual_traits_part)
                    if partition_result.get("status") == "success":
                        logger.info(
                            f"Per-chromosome h²: total={partition_result['total_h2']:.4f} across "
                            f"{partition_result['n_chromosomes']} chromosomes"
                        )
                        results["outputs"]["heritability_partition"] = partition_result

                        # Generate heritability bar chart
                        results_dir_h2 = Path(config.get("results_dir", config.get("output_dir", ".")))
                        results_dir_h2.mkdir(parents=True, exist_ok=True)
                        bar_result = heritability_bar_chart(
                            partition_result,
                            output_file=results_dir_h2 / "heritability_bar_chart.png",
                        )
                        if bar_result.get("status") == "success":
                            logger.info(f"Heritability bar chart saved to {bar_result.get('output_path')}")
                    else:
                        logger.warning(f"Heritability partition: {partition_result.get('message', 'failed')}")
            except Exception as e:
                logger.warning(f"Heritability partitioning failed: {e}", exc_info=True)

        results["steps_completed"].append("heritability_estimation")
        results["steps"].append("heritability_estimation")

        # Step 5.5b: Write summary statistics to TSV
        logger.info("Step 5.5b: Writing summary statistics")
        out_dir = Path(config.get("output_dir", config.get("work_dir", ".")))
        out_dir.mkdir(parents=True, exist_ok=True)
        results_dir = Path(config.get("results_dir", config.get("output_dir", config.get("work_dir", "."))))
        results_dir.mkdir(parents=True, exist_ok=True)

        if association_results and variants_info:
            try:
                summary_path = results_dir / "summary_statistics.tsv"
                write_summary_statistics(association_results, variants_info, summary_path)
                logger.info(f"Summary statistics written to {summary_path}")

                sig_threshold = config.get("significance_threshold", 5e-8)
                sig_path = results_dir / "significant_hits.tsv"
                write_significant_hits(association_results, variants_info, sig_path, threshold=sig_threshold)
                logger.info(f"Significant hits written to {sig_path}")
            except Exception as exc:
                logger.warning(f"Summary stats output failed: {exc}")

        results["steps_completed"].append("summary_statistics")
        results["steps"].append("summary_statistics")

        # Step 6: Generate visualization plots
        logger.info("Step 6: Generating visualization plots")
        plot_results = {}
        if association_results:
            # Save PCA intermediate data
            pca_file = out_dir / "pca_results.json"
            kinship_file = out_dir / "kinship_results.json"

            pcs_list = pca_result.get("pcs", [])
            explained_var = pca_result.get("explained_variance_ratio", [])
            if pcs_list:
                n_samples_pca = len(pcs_list)
                n_pcs_pca = len(pcs_list[0])
                transposed_pcs = [[pcs_list[s][p] for s in range(n_samples_pca)] for p in range(n_pcs_pca)]
                pca_plot_data = {
                    "components": transposed_pcs,
                    "variance": explained_var,
                    "explained_variance": explained_var,
                    "loadings": [],
                }
                with open(pca_file, "w") as f:
                    json.dump(pca_plot_data, f)
            else:
                pca_file = None

            if kinship_list:
                with open(kinship_file, "w") as f:
                    json.dump(kinship_list, f)
            else:
                kinship_file = None

            # 6a: Core plots (Manhattan, QQ, PCA, Kinship)
            plot_results = generate_all_plots(
                association_results=association_results,
                output_dir=results_dir,
                pca_file=pca_file,
                kinship_file=kinship_file,
                vcf_file=Path(config["vcf_path"]),
            )

            # 6b: Effect size distribution plot
            try:
                from metainformant.gwas.visualization.general import effect_size_plot
                fig = effect_size_plot(association_results, output_path=results_dir / "effect_size_plot.png")
                if fig is not None:
                    plot_results["effect_size"] = str(results_dir / "effect_size_plot.png")
                    import matplotlib.pyplot as plt
                    plt.close(fig)
            except Exception as exc:
                logger.warning(f"Effect size plot failed: {exc}")

            # 6c: MAF spectrum plot
            try:
                maf_values = [r.get("maf", 0.0) for r in association_results if r.get("maf") is not None]
                if maf_values and any(m > 0 for m in maf_values):
                    import matplotlib.pyplot as plt
                    fig, ax = plt.subplots(figsize=(8, 5))
                    ax.hist(maf_values, bins=30, edgecolor="black", alpha=0.7, color="teal")
                    ax.set_xlabel("Minor Allele Frequency")
                    ax.set_ylabel("Number of Variants")
                    ax.set_title(f"MAF Spectrum (n={len(maf_values)} variants)")
                    ax.grid(True, alpha=0.3)
                    plt.tight_layout()
                    maf_path = results_dir / "maf_spectrum_plot.png"
                    fig.savefig(maf_path, dpi=300, bbox_inches="tight")
                    plot_results["maf_spectrum"] = str(maf_path)
                    logger.info(f"Saved MAF spectrum plot to {maf_path}")
                    plt.close(fig)
            except Exception as exc:
                logger.warning(f"MAF spectrum plot failed: {exc}")

            # 6d: Regional plot for top hit
            try:
                from metainformant.gwas.visualization.general import regional_plot
                if association_results:
                    top_hit = min(association_results, key=lambda r: r.get("p_value", 1.0))
                    top_chrom = str(top_hit.get("chrom", "1"))
                    top_pos = top_hit.get("pos", 0)
                    region_start = max(0, top_pos - 500000)
                    region_end = top_pos + 500000
                    fig = regional_plot(
                        association_results, top_chrom, region_start, region_end,
                        output_path=results_dir / "regional_plot.png",
                    )
                    if fig is not None:
                        plot_results["regional"] = str(results_dir / "regional_plot.png")
                        import matplotlib.pyplot as plt
                        plt.close(fig)
            except Exception as exc:
                logger.warning(f"Regional plot failed: {exc}")

            logger.info(f"Generated {len(plot_results)} plots: {list(plot_results.keys())}")

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
