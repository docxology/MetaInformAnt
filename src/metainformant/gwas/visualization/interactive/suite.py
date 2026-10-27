"""Comprehensive GWAS visualization suite.

This module provides functions to generate complete sets of GWAS plots.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, TextIO, Union

import numpy as np

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

#: Field delimiter of each text table format (``None`` splits on any run of whitespace).
_TEXT_TABLE_DELIMITERS: Dict[str, Optional[str]] = {".tsv": "\t", ".csv": ",", ".txt": None}
_SUPPORTED_TABLE_FORMATS = ".npy, .npz, .tsv, .csv, .txt or .json"


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


def _is_number(cell: str) -> bool:
    try:
        float(cell)
    except ValueError:
        return False
    return True


def _read_numeric_table(path: Path) -> np.ndarray:
    """Read a ``.tsv``/``.csv``/``.txt`` table of numbers into a 2-D float array.

    A header row (a first row with any non-numeric cell) and a column of sample labels (a first column that
    is non-numeric in every data row) are dropped, so both a bare matrix and a pandas-style table with
    sample names load the same way. Blank lines and ``#`` comment lines are ignored.

    Raises:
        ValueError: If the file has no data rows, is not rectangular, or has non-numeric cells.
    """
    delimiter = _TEXT_TABLE_DELIMITERS[path.suffix.lower()]
    rows: List[List[str]] = []
    with open(path, newline="") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            rows.append([cell.strip() for cell in (text.split(delimiter) if delimiter else text.split())])

    if rows and not all(_is_number(cell) for cell in rows[0]):
        rows = rows[1:]  # header
    if not rows:
        raise ValueError(f"{path.name} contains no data rows")
    if all(not _is_number(row[0]) for row in rows):
        rows = [row[1:] for row in rows]  # sample labels
    width = len(rows[0])
    if width == 0 or any(len(row) != width for row in rows):
        raise ValueError(f"{path.name} is not a rectangular table of numbers")
    try:
        return np.array([[float(cell) for cell in row] for row in rows], dtype=float)
    except ValueError as exc:
        raise ValueError(f"{path.name} contains non-numeric cells") from exc


def _first_nonempty(payload: Dict[str, Any], keys: List[str]) -> Optional[List[float]]:
    """First non-empty sequence among ``payload[key]`` for ``keys``, as a list of floats."""
    for key in keys:
        values = payload.get(key)
        if values is not None and len(values) > 0:
            return [float(v) for v in np.ravel(values)]
    return None


def _load_pca_file(path: Path) -> Dict[str, Any]:
    """Load principal components from a ``.npy``, ``.npz``, ``.tsv``, ``.csv``, ``.txt`` or ``.json`` file.

    Orientation is always ``(n_samples, n_components)``: arrays and tables are read as they are, an ``.npz``
    prefers its ``pcs`` array, and the ``components`` list the GWAS workflow writes to JSON (one row per
    component) is transposed.

    Returns:
        ``{"pcs": ndarray}``, plus ``"explained_variance_ratio"`` (list of floats) when the file carries it.

    Raises:
        ValueError: If the format is unsupported or the content is not a 2-D numeric table.
    """
    suffix = path.suffix.lower()
    explained: Optional[List[float]] = None
    if suffix == ".npy":
        pcs = np.load(path, allow_pickle=False)
    elif suffix == ".npz":
        with np.load(path, allow_pickle=False) as archive:
            names = list(archive.files)
            if not names:
                raise ValueError(f"{path.name} contains no arrays")
            if "pcs" in names:
                pcs = archive["pcs"]
            elif "components" in names:
                pcs = archive["components"].T
            else:
                pcs = archive[names[0]]
            variance_keys = ["explained_variance_ratio", "explained_variance"]
            explained = _first_nonempty({k: archive[k] for k in variance_keys if k in names}, variance_keys)
    elif suffix == ".json":
        payload = json.loads(path.read_text())
        if isinstance(payload, dict):
            if "pcs" in payload:
                pcs = np.asarray(payload["pcs"], dtype=float)
            elif "components" in payload:
                pcs = np.asarray(payload["components"], dtype=float).T
            else:
                raise ValueError(f"{path.name} has neither a 'pcs' nor a 'components' entry")
            explained = _first_nonempty(payload, ["explained_variance_ratio", "explained_variance", "variance"])
        else:
            pcs = np.asarray(payload, dtype=float)
    elif suffix in _TEXT_TABLE_DELIMITERS:
        pcs = _read_numeric_table(path)
    else:
        raise ValueError(f"Unsupported file format {suffix!r} for {path.name}; use {_SUPPORTED_TABLE_FORMATS}")

    pcs = np.asarray(pcs, dtype=float)
    if pcs.ndim != 2 or pcs.size == 0:
        raise ValueError(
            f"PCA data in {path.name} must be a 2-D (n_samples, n_components) table, got shape {pcs.shape}"
        )
    loaded: Dict[str, Any] = {"pcs": pcs}
    if explained is not None:
        loaded["explained_variance_ratio"] = explained
    return loaded


def _load_kinship_file(path: Path) -> np.ndarray:
    """Load a square kinship matrix from a ``.npy``, ``.npz``, ``.tsv``, ``.csv``, ``.txt`` or ``.json`` file.

    An ``.npz`` is read from its ``kinship_matrix``/``kinship``/``matrix`` array (else its first array); a
    JSON file may be the matrix itself or a result dict with ``kinship_matrix`` or ``matrix``.

    Raises:
        ValueError: If the format is unsupported or the content is not a square numeric matrix.
    """
    suffix = path.suffix.lower()
    if suffix == ".npy":
        matrix = np.load(path, allow_pickle=False)
    elif suffix == ".npz":
        with np.load(path, allow_pickle=False) as archive:
            names = list(archive.files)
            if not names:
                raise ValueError(f"{path.name} contains no arrays")
            key = next((k for k in ("kinship_matrix", "kinship", "matrix") if k in names), names[0])
            matrix = archive[key]
    elif suffix == ".json":
        payload = json.loads(path.read_text())
        if isinstance(payload, dict):
            payload = payload.get("kinship_matrix", payload.get("matrix"))
            if payload is None:
                raise ValueError(f"{path.name} has neither a 'kinship_matrix' nor a 'matrix' entry")
        matrix = np.asarray(payload, dtype=float)
    elif suffix in _TEXT_TABLE_DELIMITERS:
        matrix = _read_numeric_table(path)
    else:
        raise ValueError(f"Unsupported file format {suffix!r} for {path.name}; use {_SUPPORTED_TABLE_FORMATS}")

    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"Kinship matrix in {path.name} must be square, got shape {matrix.shape}")
    return matrix


def _read_vcf_sites(path: Path) -> List[Dict[str, Any]]:
    """Read CHROM, POS, REF and ALT of every record of a plain or gzip-compressed VCF.

    ``ALT`` is kept verbatim (comma-separated for a multi-allelic site), so the transition/transversion
    count, which only classifies single-base REF/ALT pairs, never guesses which allele to use.
    """
    handle: TextIO = gzip.open(path, "rt") if path.suffix.lower() == ".gz" else open(path, "rt")
    sites: List[Dict[str, Any]] = []
    with handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t", 5)
            if len(fields) < 5:
                continue
            chrom, pos, _variant_id, ref, alt = fields[:5]
            try:
                position = int(pos)
            except ValueError:
                continue
            sites.append({"CHROM": chrom, "POS": position, "REF": ref, "ALT": alt})
    return sites


def _close_figure(figure: Any) -> None:
    """Close a figure returned by a plotting helper once the helper has saved it.

    The helpers save their figure and hand it back still open; closing it here keeps repeated suite runs from
    accumulating open pyplot figures (matplotlib warns at 20, and this repository's test configuration turns
    that warning into an error that silently drops the plots drawn afterwards). Non-figures are ignored.
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.figure import Figure
    except ImportError:
        return
    if isinstance(figure, Figure):
        plt.close(figure)


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
        pca_file: Optional principal-component file (``.npy``, ``.npz``, ``.tsv``, ``.csv``, ``.txt`` or
            ``.json``; one row per sample, one column per component). It gets its own PCA plot and, unless
            ``pca_data`` is given, feeds every PCA-dependent plot including the summary panel.
        kinship_file: Optional square kinship matrix file in the same formats. It gets its own heatmap and,
            unless ``kinship_matrix`` is given, feeds every kinship-dependent plot including the summary panel.
        vcf_file: Optional VCF file (plain or ``.gz``). Produces a variant-density plot and a
            transition/transversion plot, and records ``vcf_variants`` and ``ts_tv_ratio`` in the statistics.
        significance_threshold: P-value threshold for significance
        metadata: Optional sample metadata dict {sample_id: {latitude, longitude,
            population, ...}} for geographic and population-aware plots.
        kinship_matrix: Optional kinship matrix (ndarray or list of lists) for
            dendrogram, clustermap and summary-panel visualizations. Wins over ``kinship_file``.
        genotype_matrix: Optional genotype matrix (variants x samples) for
            LD decay, genotype-phenotype, and top-hits plots.
        pca_data: Optional PCA data dict with 'pcs', 'explained_variance_ratio',
            and optionally 'sample_ids' for interactive, composite and summary-panel PCA plots.
            Wins over ``pca_file``.
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
        Every figure the call opens is saved and then closed, so it leaves no pyplot figure open.

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
            from ..general import manhattan_plot

            plot = manhattan_plot(
                records,
                output_path=manhattan_file,
                significance_threshold=significance_threshold,
            )
            if plot:
                results["plots"]["manhattan"] = str(manhattan_file)
                results["files_created"].append(str(manhattan_file))
            _close_figure(plot)
        except Exception as e:
            logger.warning(f"Failed to create Manhattan plot: {e}")

        # Generate Q-Q plot
        qq_file = output_dir / "qq_plot.png"
        try:
            from ..general import qq_plot

            p_values = df[p_col].dropna().values
            plot = qq_plot(p_values.tolist(), output_path=qq_file)
            if plot:
                results["plots"]["qq"] = str(qq_file)
                results["files_created"].append(str(qq_file))
            _close_figure(plot)
        except Exception as e:
            logger.warning(f"Failed to create Q-Q plot: {e}")

        # Generate regional plots for top hits
        try:
            from ..general import regional_plot

            top_hits = df.nsmallest(5, p_col)

            for i, (_, hit) in enumerate(top_hits.iterrows()):
                chrom = str(hit[chrom_col])
                pos = int(hit[pos_col])
                regional_file = output_dir / f"regional_plot_chr{chrom}_{pos}.png"

                plot = regional_plot(
                    records,
                    chrom,
                    pos - 500000,
                    pos + 500000,
                    output_path=regional_file,
                )
                if plot:
                    results["plots"][f"regional_{i + 1}"] = str(regional_file)
                    results["files_created"].append(str(regional_file))
                _close_figure(plot)
        except Exception as e:
            logger.warning(f"Failed to create regional plots: {e}")

        # Generate effect size plot if BETA column present
        if "BETA" in df.columns:
            effect_file = output_dir / "effect_size_plot.png"
            try:
                from ..general import effect_size_plot

                plot = effect_size_plot(records, output_path=effect_file)
                if plot:
                    results["plots"]["effect_size"] = str(effect_file)
                    results["files_created"].append(str(effect_file))
                _close_figure(plot)
            except Exception as e:
                logger.warning(f"Failed to create effect size plot: {e}")

        # Generate PCA plot if PCA data available. The file is read in whichever supported format it is in
        # (see _load_pca_file); a file that cannot be read is reported with the reason, not silently skipped.
        pca_file_data: Optional[Dict[str, Any]] = None
        if pca_file:
            if Path(pca_file).exists():
                try:
                    pca_file_data = _load_pca_file(Path(pca_file))
                    pca_plot_file = output_dir / "pca_plot.png"

                    from ..general import pca_plot

                    # pca_plot takes (components, variance, loadings) with one row per component
                    plot = pca_plot(
                        (pca_file_data["pcs"].T, [], []),
                        output_path=pca_plot_file,
                        explained_var=pca_file_data.get("explained_variance_ratio"),
                    )
                    if plot:
                        results["plots"]["pca"] = str(pca_plot_file)
                        results["files_created"].append(str(pca_plot_file))
                    _close_figure(plot)
                except Exception as e:
                    logger.warning(f"Failed to create PCA plot from {pca_file}: {e}")
            else:
                logger.warning(f"PCA file not found: {pca_file}")

        # Generate kinship heatmap if kinship data available
        kinship_file_data: Optional[np.ndarray] = None
        if kinship_file:
            if Path(kinship_file).exists():
                try:
                    kinship_file_data = _load_kinship_file(Path(kinship_file))
                    kinship_plot_file = output_dir / "kinship_heatmap.png"

                    from ..general import kinship_heatmap

                    plot = kinship_heatmap(kinship_file_data, output_path=kinship_plot_file)
                    if plot:
                        results["plots"]["kinship"] = str(kinship_plot_file)
                        results["files_created"].append(str(kinship_plot_file))
                    _close_figure(plot)
                except Exception as e:
                    logger.warning(f"Failed to create kinship heatmap from {kinship_file}: {e}")
            else:
                logger.warning(f"Kinship file not found: {kinship_file}")

        # Data passed in memory wins; the files fill in whatever was not passed, so the summary panel and
        # every other PCA/kinship-dependent plot below see the same inputs the standalone plots used
        if pca_data is None:
            pca_data = pca_file_data
        if kinship_matrix is None:
            kinship_matrix = kinship_file_data

        # Variant plots from the VCF: variant density along the genome and the transition/transversion spectrum
        if vcf_file:
            vcf_path = Path(vcf_file)
            if not vcf_path.exists():
                logger.warning(f"VCF file not found, no variant plots made: {vcf_path}")
            else:
                try:
                    from ..genomic.variants import transition_transversion_plot, variant_density_plot

                    vcf_sites = _read_vcf_sites(vcf_path)
                    if not vcf_sites:
                        logger.warning(f"No variant records in {vcf_path}, no variant plots made")
                    else:
                        results["statistics"]["vcf_variants"] = len(vcf_sites)

                        density_file = output_dir / "variant_density.png"
                        if variant_density_plot(vcf_sites, output_file=density_file).get("status") == "success":
                            results["plots"]["variant_density"] = str(density_file)
                            results["files_created"].append(str(density_file))

                        tstv_file = output_dir / "transition_transversion.png"
                        tstv_result = transition_transversion_plot(vcf_sites, output_file=tstv_file)
                        if tstv_result.get("status") == "success":
                            results["plots"]["transition_transversion"] = str(tstv_file)
                            results["files_created"].append(str(tstv_file))
                            results["statistics"]["ts_tv_ratio"] = tstv_result["ts_tv_ratio"]
                except Exception as e:
                    logger.warning(f"Failed to create VCF variant plots from {vcf_path}: {e}")

        # Generate enhanced visualizations from new modules
        # Phenotype distribution (if phenotype data embedded in results)
        try:
            from .phenotype import phenotype_distribution

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
            from .composite import gwas_summary_panel

            summary_file = output_dir / "gwas_summary_panel.png"
            result = gwas_summary_panel(
                records,
                pca_data=pca_data,
                kinship_matrix=kinship_matrix,
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
            from .interactive import interactive_manhattan

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
                from ..genomic.ld import compute_ld_decay, ld_decay_plot

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
                from ..genomic.ld import ld_heatmap_region

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
                from ..population.geography import sample_map

                sample_map_file = output_dir / "sample_map.png"
                sm_result = sample_map(metadata=metadata, output_file=sample_map_file)
                if sm_result.get("status") == "success":
                    results["plots"]["sample_map"] = str(sample_map_file)
                    results["files_created"].append(str(sample_map_file))
            except Exception as e:
                logger.warning(f"Failed to create sample map: {e}")

            # Population count map
            try:
                from ..population.geography import population_count_map

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
                from ..population.population import kinship_dendrogram

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
                from ..population.population import kinship_clustermap

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
                from .composite import population_structure_panel

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
            from .finemapping import credible_set_plot

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
                from .phenotype import phenotype_correlation_matrix

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
                from .phenotype import genotype_phenotype_boxplot

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
                from .phenotype import top_hits_genotype_phenotype

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
                from .phenotype import phenotype_pca_correlation

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
                from .interactive import interactive_pca

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
                from .interactive import interactive_volcano

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
                from ...analysis.heritability import heritability_bar_chart

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
                from ...analysis.correction import lambda_gc_from_p_values

                lambda_gc = lambda_gc_from_p_values(p_values.tolist())
                if lambda_gc is not None:
                    results["statistics"]["lambda_gc"] = float(lambda_gc)
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
