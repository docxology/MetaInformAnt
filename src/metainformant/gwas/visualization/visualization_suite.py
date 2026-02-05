"""Comprehensive GWAS visualization suite.

This module provides functions to generate complete sets of GWAS plots.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def _get_p_value_column(df: Any) -> str:
    """Find the p-value column name in a DataFrame.

    Checks common column name variants in priority order.

    Args:
        df: pandas DataFrame with GWAS results

    Returns:
        The column name containing p-values

    Raises:
        KeyError: If no p-value column is found
    """
    candidates = ["p_value", "P", "pvalue", "p_val", "PVALUE", "P_VALUE"]
    for col in candidates:
        if col in df.columns:
            return col
    raise KeyError(f"No p-value column found. Available columns: {list(df.columns)}")


def _get_chrom_column(df: Any) -> str:
    """Find the chromosome column name in a DataFrame."""
    candidates = ["CHROM", "CHR", "chrom", "chr", "chromosome", "CHROMOSOME"]
    for col in candidates:
        if col in df.columns:
            return col
    raise KeyError(f"No chromosome column found. Available columns: {list(df.columns)}")


def _get_pos_column(df: Any) -> str:
    """Find the position column name in a DataFrame."""
    candidates = ["POS", "BP", "pos", "bp", "position", "POSITION"]
    for col in candidates:
        if col in df.columns:
            return col
    raise KeyError(f"No position column found. Available columns: {list(df.columns)}")


def generate_all_plots(
    association_results: Union[Path, List[Any]],
    output_dir: Path,
    pca_file: Optional[Path] = None,
    kinship_file: Optional[Path] = None,
    vcf_file: Optional[Path] = None,
    significance_threshold: float = 5e-8,
    metadata: Optional[Dict[str, Dict[str, Any]]] = None,
    kinship_matrix: Optional[Any] = None,
    genotype_matrix: Optional[List[List[int]]] = None,
    pca_data: Optional[Dict[str, Any]] = None,
    finemapping_results: Optional[List[Dict[str, Any]]] = None,
    h2_data: Optional[Dict[str, Any]] = None,
    phenotypes_dict: Optional[Dict[str, List[float]]] = None,
    phenotype_values: Optional[List[float]] = None,
) -> Dict[str, Any]:
    """Generate a comprehensive set of GWAS plots.

    Args:
        association_results: Path to association results file, or a list of results
        output_dir: Output directory for plots
        pca_file: Optional PCA file for population structure plots
        kinship_file: Optional kinship matrix file
        vcf_file: Optional VCF file for variant plots
        significance_threshold: P-value threshold for significance
        metadata: Optional sample metadata dict {sample_id: {latitude, longitude,
            population, ...}} for geographic and population-aware plots.
        kinship_matrix: Optional kinship matrix (ndarray or list of lists) for
            dendrogram and clustermap visualizations.
        genotype_matrix: Optional genotype matrix (variants x samples) for
            LD decay, genotype-phenotype, and top-hits plots.
        pca_data: Optional PCA data dict with 'pcs', 'explained_variance_ratio',
            and optionally 'sample_ids' for interactive and composite PCA plots.
        finemapping_results: Optional list of fine-mapping association results
            for credible set and PIP plots.
        h2_data: Optional heritability data dict from partition_heritability_by_chromosome
            for heritability bar chart.
        phenotypes_dict: Optional mapping of trait_name -> list of values for
            multi-trait correlation matrix.
        phenotype_values: Optional list of phenotype values (one per sample) for
            genotype-phenotype and PCA-correlation plots.

    Returns:
        Dictionary with paths to generated plots, summary statistics,
        'status' ('completed' or 'failed'), and 'num_plots_generated' count.

    Example:
        >>> results = generate_all_plots(results_file, output_dir)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results: Dict[str, Any] = {"plots": {}, "statistics": {}, "files_created": []}

    # Handle empty list input gracefully
    if isinstance(association_results, list):
        if len(association_results) == 0:
            results["status"] = "completed"
            results["num_plots_generated"] = 0
            return results
        # If list of dicts, convert to DataFrame
        try:
            import pandas as pd

            df = pd.DataFrame(association_results)
        except Exception as e:
            logger.error(f"Failed to convert list to DataFrame: {e}")
            results["status"] = "failed"
            results["num_plots_generated"] = 0
            return results
    else:
        association_results = Path(association_results)
        try:
            # Load association results from file
            if association_results.suffix.lower() == ".csv":
                import pandas as pd

                df = pd.read_csv(association_results)
            elif association_results.suffix.lower() == ".tsv":
                import pandas as pd

                df = pd.read_csv(association_results, sep="\t")
            else:
                logger.error(f"Unsupported file format: {association_results.suffix}")
                results["status"] = "failed"
                results["num_plots_generated"] = 0
                return results
        except Exception as e:
            logger.error(f"Failed to load association results: {e}")
            results["status"] = "failed"
            results["num_plots_generated"] = 0
            return results

    try:
        # Resolve flexible column names
        p_col = _get_p_value_column(df)
        chrom_col = _get_chrom_column(df)
        pos_col = _get_pos_column(df)

        # Convert DataFrame rows to list-of-dicts format for general.py plot functions
        records = []
        for _, row in df.iterrows():
            record: Dict[str, Any] = {
                "chrom": str(row[chrom_col]),
                "pos": int(row[pos_col]),
                "p_value": float(row[p_col]),
            }
            # Include optional columns if present
            if "BETA" in df.columns:
                record["beta"] = float(row["BETA"])
            if "SE" in df.columns:
                record["se"] = float(row["SE"])
            records.append(record)

        # Generate Manhattan plot
        manhattan_file = output_dir / "manhattan_plot.png"
        try:
            from .general import manhattan_plot

            plot = manhattan_plot(records, output_path=manhattan_file, significance_threshold=significance_threshold)
            if plot:
                results["plots"]["manhattan"] = str(manhattan_file)
                results["files_created"].append(str(manhattan_file))
        except Exception as e:
            logger.warning(f"Failed to create Manhattan plot: {e}")

        # Generate Q-Q plot
        qq_file = output_dir / "qq_plot.png"
        try:
            from .general import qq_plot

            p_values = df[p_col].dropna().values
            plot = qq_plot(p_values.tolist(), output_path=qq_file)
            if plot:
                results["plots"]["qq"] = str(qq_file)
                results["files_created"].append(str(qq_file))
        except Exception as e:
            logger.warning(f"Failed to create Q-Q plot: {e}")

        # Generate regional plots for top hits
        try:
            from .general import regional_plot

            top_hits = df.nsmallest(5, p_col)

            for i, (_, hit) in enumerate(top_hits.iterrows()):
                chrom = str(hit[chrom_col])
                pos = int(hit[pos_col])
                regional_file = output_dir / f"regional_plot_chr{chrom}_{pos}.png"

                plot = regional_plot(records, chrom, pos - 500000, pos + 500000, output_path=regional_file)
                if plot:
                    results["plots"][f"regional_{i+1}"] = str(regional_file)
                    results["files_created"].append(str(regional_file))
        except Exception as e:
            logger.warning(f"Failed to create regional plots: {e}")

        # Generate effect size plot if BETA column present
        if "BETA" in df.columns:
            effect_file = output_dir / "effect_size_plot.png"
            try:
                from .general import effect_size_plot

                plot = effect_size_plot(records, output_path=effect_file)
                if plot:
                    results["plots"]["effect_size"] = str(effect_file)
                    results["files_created"].append(str(effect_file))
            except Exception as e:
                logger.warning(f"Failed to create effect size plot: {e}")

        # Generate PCA plot if PCA data available
        if pca_file and Path(pca_file).exists():
            try:
                pca_file_data = np.load(pca_file)
                pca_plot_file = output_dir / "pca_plot.png"

                from .general import pca_plot

                plot = pca_plot(pca_file_data, output_path=pca_plot_file)
                if plot:
                    results["plots"]["pca"] = str(pca_plot_file)
                    results["files_created"].append(str(pca_plot_file))
            except Exception as e:
                logger.warning(f"Failed to create PCA plot: {e}")

        # Generate kinship heatmap if kinship data available
        if kinship_file and Path(kinship_file).exists():
            try:
                kinship_data = np.load(kinship_file)
                kinship_plot_file = output_dir / "kinship_heatmap.png"

                from .general import kinship_heatmap

                plot = kinship_heatmap(kinship_data, output_path=kinship_plot_file)
                if plot:
                    results["plots"]["kinship"] = str(kinship_plot_file)
                    results["files_created"].append(str(kinship_plot_file))
            except Exception as e:
                logger.warning(f"Failed to create kinship heatmap: {e}")

        # Generate enhanced visualizations from new modules
        # Phenotype distribution (if phenotype data embedded in results)
        try:
            from .visualization_phenotype import phenotype_distribution

            if "phenotype" in df.columns:
                pheno_file = output_dir / "phenotype_distribution.png"
                pheno_values = df["phenotype"].dropna().tolist()
                if pheno_values:
                    result = phenotype_distribution(pheno_values, output_file=pheno_file)
                    if result.get("status") == "success":
                        results["plots"]["phenotype_distribution"] = str(pheno_file)
                        results["files_created"].append(str(pheno_file))
        except Exception as e:
            logger.warning(f"Failed to create phenotype distribution: {e}")

        # Composite summary panel
        try:
            from .visualization_composite import gwas_summary_panel

            summary_file = output_dir / "gwas_summary_panel.png"
            result = gwas_summary_panel(
                records,
                output_file=summary_file,
                significance_threshold=significance_threshold,
            )
            if result.get("status") == "success":
                results["plots"]["summary_panel"] = str(summary_file)
                results["files_created"].append(str(summary_file))
        except Exception as e:
            logger.warning(f"Failed to create summary panel: {e}")

        # Interactive Manhattan (HTML)
        try:
            from .visualization_interactive import interactive_manhattan

            interactive_file = output_dir / "manhattan_interactive.html"
            result = interactive_manhattan(records, output_file=interactive_file)
            if result.get("status") == "success":
                results["plots"]["manhattan_interactive"] = str(interactive_file)
                results["files_created"].append(str(interactive_file))
        except Exception as e:
            logger.warning(f"Failed to create interactive Manhattan: {e}")

        # ---------------------------------------------------------------
        # LD Decay and LD Heatmap (requires genotype_matrix + positions)
        # ---------------------------------------------------------------
        if genotype_matrix is not None and len(genotype_matrix) >= 2:
            positions_list = [int(row[pos_col]) for _, row in df.iterrows()]

            # LD Decay curve
            try:
                from .visualization_ld import compute_ld_decay, ld_decay_plot

                logger.info("Computing LD decay curve...")
                ld_decay_data = compute_ld_decay(
                    genotypes_by_variant=genotype_matrix,
                    positions=positions_list,
                )
                if ld_decay_data.get("status") == "success":
                    ld_decay_file = output_dir / "ld_decay_plot.png"
                    ld_result = ld_decay_plot(ld_decay_data, output_file=ld_decay_file)
                    if ld_result.get("status") == "success":
                        results["plots"]["ld_decay"] = str(ld_decay_file)
                        results["files_created"].append(str(ld_decay_file))
            except Exception as e:
                logger.warning(f"Failed to create LD decay plot: {e}")

            # LD Heatmap for top hit region
            try:
                from .visualization_ld import ld_heatmap_region

                top_hit_row = df.nsmallest(1, p_col).iloc[0]
                top_pos = int(top_hit_row[pos_col])
                region_start = max(0, top_pos - 500000)
                region_end = top_pos + 500000

                ld_heatmap_file = output_dir / "ld_heatmap_top_region.png"
                ld_hm_result = ld_heatmap_region(
                    genotypes_by_variant=genotype_matrix,
                    positions=positions_list,
                    output_file=ld_heatmap_file,
                    region_start=region_start,
                    region_end=region_end,
                    title=f"LD Heatmap (top hit region: {region_start}-{region_end})",
                )
                if ld_hm_result.get("status") == "success":
                    results["plots"]["ld_heatmap"] = str(ld_heatmap_file)
                    results["files_created"].append(str(ld_heatmap_file))
            except Exception as e:
                logger.warning(f"Failed to create LD heatmap: {e}")

        # ---------------------------------------------------------------
        # Geographic plots (requires metadata with lat/lon)
        # ---------------------------------------------------------------
        if metadata is not None and len(metadata) > 0:
            # Sample map
            try:
                from .visualization_geography import sample_map

                sample_map_file = output_dir / "sample_map.png"
                sm_result = sample_map(metadata=metadata, output_file=sample_map_file)
                if sm_result.get("status") == "success":
                    results["plots"]["sample_map"] = str(sample_map_file)
                    results["files_created"].append(str(sample_map_file))
            except Exception as e:
                logger.warning(f"Failed to create sample map: {e}")

            # Population count map
            try:
                from .visualization_geography import population_count_map

                pop_count_file = output_dir / "population_count_map.png"
                pcm_result = population_count_map(metadata=metadata, output_file=pop_count_file)
                if pcm_result.get("status") == "success":
                    results["plots"]["population_count_map"] = str(pop_count_file)
                    results["files_created"].append(str(pop_count_file))
            except Exception as e:
                logger.warning(f"Failed to create population count map: {e}")

        # ---------------------------------------------------------------
        # Kinship dendrogram and clustermap (requires kinship_matrix)
        # ---------------------------------------------------------------
        if kinship_matrix is not None:
            # Kinship dendrogram
            try:
                from .visualization_population import kinship_dendrogram

                dendro_file = output_dir / "kinship_dendrogram.png"
                dendro_result = kinship_dendrogram(
                    kinship_matrix=kinship_matrix,
                    output_file=dendro_file,
                )
                if dendro_result.get("status") == "success":
                    results["plots"]["kinship_dendrogram"] = str(dendro_file)
                    results["files_created"].append(str(dendro_file))
            except Exception as e:
                logger.warning(f"Failed to create kinship dendrogram: {e}")

            # Kinship clustermap
            try:
                from .visualization_population import kinship_clustermap

                clustermap_file = output_dir / "kinship_clustermap.png"
                cm_result = kinship_clustermap(
                    kinship_matrix=kinship_matrix,
                    output_file=clustermap_file,
                )
                if cm_result.get("status") == "success":
                    results["plots"]["kinship_clustermap"] = str(clustermap_file)
                    results["files_created"].append(str(clustermap_file))
            except Exception as e:
                logger.warning(f"Failed to create kinship clustermap: {e}")

        # ---------------------------------------------------------------
        # Population structure panel (requires pca_data + kinship_matrix)
        # ---------------------------------------------------------------
        if pca_data is not None and kinship_matrix is not None:
            try:
                from .visualization_composite import population_structure_panel

                pop_panel_file = output_dir / "population_structure_panel.png"
                psp_result = population_structure_panel(
                    pca_data=pca_data,
                    kinship_matrix=kinship_matrix,
                    metadata=metadata,
                    output_file=pop_panel_file,
                )
                if psp_result.get("status") == "success":
                    results["plots"]["population_structure_panel"] = str(pop_panel_file)
                    results["files_created"].append(str(pop_panel_file))
            except Exception as e:
                logger.warning(f"Failed to create population structure panel: {e}")

        # ---------------------------------------------------------------
        # Fine-mapping: credible set plot (requires assoc_results)
        # ---------------------------------------------------------------
        try:
            from .visualization_finemapping import credible_set_plot

            cs_input = finemapping_results if finemapping_results is not None else records
            if cs_input and len(cs_input) >= 2:
                cs_file = output_dir / "credible_set_plot.png"
                cs_result = credible_set_plot(
                    assoc_results=cs_input,
                    output_file=cs_file,
                )
                if cs_result.get("status") == "success":
                    results["plots"]["credible_set"] = str(cs_file)
                    results["files_created"].append(str(cs_file))
        except Exception as e:
            logger.warning(f"Failed to create credible set plot: {e}")

        # ---------------------------------------------------------------
        # Phenotype correlation matrix (requires phenotypes_dict with 2+ traits)
        # ---------------------------------------------------------------
        if phenotypes_dict is not None and len(phenotypes_dict) >= 2:
            try:
                from .visualization_phenotype import phenotype_correlation_matrix

                corr_file = output_dir / "phenotype_correlation_matrix.png"
                corr_result = phenotype_correlation_matrix(
                    phenotypes_dict=phenotypes_dict,
                    output_file=corr_file,
                )
                if corr_result.get("status") == "success":
                    results["plots"]["phenotype_correlation_matrix"] = str(corr_file)
                    results["files_created"].append(str(corr_file))
            except Exception as e:
                logger.warning(f"Failed to create phenotype correlation matrix: {e}")

        # ---------------------------------------------------------------
        # Genotype-phenotype boxplot for top hit (requires genotype_matrix + phenotype_values)
        # ---------------------------------------------------------------
        if genotype_matrix is not None and phenotype_values is not None and len(genotype_matrix) > 0:
            try:
                from .visualization_phenotype import genotype_phenotype_boxplot

                # Use the top hit by p-value
                top_idx = int(df[p_col].idxmin())
                if top_idx < len(genotype_matrix):
                    top_genotypes = genotype_matrix[top_idx]
                    top_variant_id = f"chr{df.iloc[top_idx][chrom_col]}:{df.iloc[top_idx][pos_col]}"
                    gp_file = output_dir / "genotype_phenotype_boxplot.png"
                    gp_result = genotype_phenotype_boxplot(
                        genotypes=top_genotypes,
                        phenotypes=phenotype_values,
                        variant_id=top_variant_id,
                        output_file=gp_file,
                    )
                    if gp_result.get("status") == "success":
                        results["plots"]["genotype_phenotype_boxplot"] = str(gp_file)
                        results["files_created"].append(str(gp_file))
            except Exception as e:
                logger.warning(f"Failed to create genotype-phenotype boxplot: {e}")

        # ---------------------------------------------------------------
        # Top hits genotype-phenotype (requires assoc_results + genotype_matrix + phenotype_values)
        # ---------------------------------------------------------------
        if genotype_matrix is not None and phenotype_values is not None and len(records) > 0:
            try:
                from .visualization_phenotype import top_hits_genotype_phenotype

                top_hits_dir = output_dir / "top_hits_genotype_phenotype"
                th_result = top_hits_genotype_phenotype(
                    assoc_results=records,
                    genotype_matrix=genotype_matrix,
                    phenotypes=phenotype_values,
                    output_dir=top_hits_dir,
                    n_top=min(5, len(records)),
                )
                if th_result.get("status") == "success":
                    results["plots"]["top_hits_genotype_phenotype"] = str(top_hits_dir)
                    for fpath in th_result.get("generated_files", []):
                        results["files_created"].append(fpath)
            except Exception as e:
                logger.warning(f"Failed to create top hits genotype-phenotype plots: {e}")

        # ---------------------------------------------------------------
        # Phenotype-PCA correlation (requires pca_data + phenotype_values)
        # ---------------------------------------------------------------
        if pca_data is not None and phenotype_values is not None:
            try:
                from .visualization_phenotype import phenotype_pca_correlation

                pcs_raw = pca_data.get("pcs")
                if pcs_raw is not None:
                    pcs_list = pcs_raw if isinstance(pcs_raw, list) else pcs_raw.tolist()
                    if len(pcs_list) == len(phenotype_values):
                        pca_corr_file = output_dir / "phenotype_pca_correlation.png"
                        pca_corr_result = phenotype_pca_correlation(
                            pcs=pcs_list,
                            phenotype_values=phenotype_values,
                            output_file=pca_corr_file,
                        )
                        if pca_corr_result.get("status") == "success":
                            results["plots"]["phenotype_pca_correlation"] = str(pca_corr_file)
                            results["files_created"].append(str(pca_corr_file))
            except Exception as e:
                logger.warning(f"Failed to create phenotype-PCA correlation plot: {e}")

        # ---------------------------------------------------------------
        # Interactive PCA (requires pca_data with >= 3 components)
        # ---------------------------------------------------------------
        if pca_data is not None:
            try:
                from .visualization_interactive import interactive_pca

                pcs_raw = pca_data.get("pcs")
                if pcs_raw is not None:
                    pcs_arr = np.asarray(pcs_raw)
                    if pcs_arr.ndim == 2 and pcs_arr.shape[1] >= 3:
                        interactive_pca_file = output_dir / "pca_interactive.html"
                        ipca_result = interactive_pca(
                            pca_data=pca_data,
                            metadata=metadata,
                            output_file=interactive_pca_file,
                        )
                        if ipca_result.get("status") == "success":
                            results["plots"]["pca_interactive"] = str(interactive_pca_file)
                            results["files_created"].append(str(interactive_pca_file))
            except Exception as e:
                logger.warning(f"Failed to create interactive PCA: {e}")

        # ---------------------------------------------------------------
        # Interactive Volcano (requires records with beta)
        # ---------------------------------------------------------------
        if any(r.get("beta") is not None for r in records):
            try:
                from .visualization_interactive import interactive_volcano

                volcano_file = output_dir / "volcano_interactive.html"
                iv_result = interactive_volcano(
                    assoc_results=records,
                    output_file=volcano_file,
                    significance_threshold=significance_threshold,
                )
                if iv_result.get("status") == "success":
                    results["plots"]["volcano_interactive"] = str(volcano_file)
                    results["files_created"].append(str(volcano_file))
            except Exception as e:
                logger.warning(f"Failed to create interactive volcano plot: {e}")

        # ---------------------------------------------------------------
        # Heritability bar chart (requires h2_data)
        # ---------------------------------------------------------------
        if h2_data is not None:
            try:
                from ..analysis.heritability import heritability_bar_chart

                h2_file = output_dir / "heritability_bar_chart.png"
                h2_result = heritability_bar_chart(
                    h2_data=h2_data,
                    output_file=h2_file,
                )
                if h2_result.get("status") == "success":
                    results["plots"]["heritability_bar_chart"] = str(h2_file)
                    results["files_created"].append(str(h2_file))
            except Exception as e:
                logger.warning(f"Failed to create heritability bar chart: {e}")

        # Calculate basic statistics
        try:
            p_values = df[p_col].dropna().astype(float).values
            results["statistics"]["total_variants"] = len(df)
            results["statistics"]["significant_variants"] = int((p_values < significance_threshold).sum())
            if len(p_values) > 0:
                results["statistics"]["lambda_gc"] = float(np.median(p_values) / 0.456)
                results["statistics"]["min_p"] = float(np.min(p_values))
        except Exception as e:
            logger.warning(f"Failed to calculate statistics: {e}")

        num_plots = len(results["files_created"])
        logger.info(f"Generated {num_plots} plots in {output_dir}")
        results["status"] = "completed"
        results["num_plots_generated"] = num_plots

    except Exception as e:
        logger.error(f"Failed to generate plots: {e}")
        results["status"] = "failed"
        results["num_plots_generated"] = len(results["files_created"])

    return results
