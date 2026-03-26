"""Strain-specific variant analysis for population-structured GWAS.

Provides per-variant Fst computation, strain-private allele detection,
and within-strain allele frequency tables for Apis mellifera population
groups (C, I, M, R lineages).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

try:
    import numpy as np  # noqa: F401

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

# ── Constants ───────────────────────────────────────────────────────────
STRAIN_PALETTE = {"C": "#2196F3", "I": "#FF9800", "M": "#4CAF50", "R": "#F44336"}
STRAIN_ORDER = ["C", "I", "M", "R"]


def _extract_strain(sample_id: str) -> str:
    """Extract strain letter from a sample ID like 'C15ITQ' -> 'C'."""
    if sample_id and sample_id[0].isalpha():
        return sample_id[0].upper()
    return "?"


def build_strain_map(sample_ids: List[str]) -> Dict[str, List[int]]:
    """Map strain letter to list of sample indices.

    Args:
        sample_ids: List of sample IDs (e.g. ['C15ITQ', 'I12G', ...])

    Returns:
        Dict mapping strain -> [indices]
    """
    strain_map: Dict[str, List[int]] = {}
    for i, sid in enumerate(sample_ids):
        strain = _extract_strain(sid)
        strain_map.setdefault(strain, []).append(i)
    return strain_map


def compute_allele_frequencies_by_strain(
    genotype_matrix: List[List[int]],
    sample_ids: List[str],
) -> Dict[str, List[float]]:
    """Compute per-variant allele frequencies within each strain.

    Args:
        genotype_matrix: Variant-major (n_variants x n_samples), values 0/1/2/-1
        sample_ids: Sample IDs in column order

    Returns:
        Dict mapping strain -> list of allele frequencies per variant
    """
    strain_map = build_strain_map(sample_ids)
    result: Dict[str, List[float]] = {}

    for strain, indices in sorted(strain_map.items()):
        freqs = []
        for geno_row in genotype_matrix:
            valid = [geno_row[i] for i in indices if 0 <= geno_row[i] <= 2]
            if valid:
                af = sum(valid) / (2.0 * len(valid))
                freqs.append(af)
            else:
                freqs.append(float("nan"))
        result[strain] = freqs

    return result


def compute_fst_per_variant(
    genotype_matrix: List[List[int]],
    sample_ids: List[str],
    strain_pairs: Optional[List[tuple[str, str]]] = None,
    method: str = "weir_cockerham",
) -> Dict[str, List[float]]:
    """Compute per-variant Fst between strain pairs.

    Supports two estimators:
      - **weir_cockerham** (default): Weir & Cockerham (1984) Fst with finite
        sample size correction. Preferred for population-structured GWAS with
        unequal strain sizes. Reference: Weir & Cockerham (1984) Evolution 38:1358.
      - **hudson**: Simplified Hudson estimator: Fst = 1 - (Hw / Hb).
        Reference: Bhatia et al. (2013) Genome Research 23:1514.

    Args:
        genotype_matrix: Variant-major (n_variants x n_samples)
        sample_ids: Sample IDs
        strain_pairs: Pairs to compute (default: all pairwise)
        method: 'weir_cockerham' or 'hudson'

    Returns:
        Dict mapping "A_vs_B" -> list of Fst values per variant
    """
    strain_map = build_strain_map(sample_ids)
    strains_present = sorted(strain_map.keys())

    if strain_pairs is None:
        strain_pairs = [
            (a, b)
            for i, a in enumerate(strains_present)
            for b in strains_present[i + 1 :]
        ]

    af_by_strain = compute_allele_frequencies_by_strain(genotype_matrix, sample_ids)
    n_variants = len(genotype_matrix)
    result: Dict[str, List[float]] = {}

    for s1, s2 in strain_pairs:
        key = f"{s1}_vs_{s2}"
        fst_values = []
        af1 = af_by_strain.get(s1, [])
        af2 = af_by_strain.get(s2, [])
        n1 = len(strain_map.get(s1, []))
        n2 = len(strain_map.get(s2, []))

        for v in range(n_variants):
            p1 = af1[v] if v < len(af1) else float("nan")
            p2 = af2[v] if v < len(af2) else float("nan")

            if math.isnan(p1) or math.isnan(p2):
                fst_values.append(float("nan"))
                continue

            if method == "weir_cockerham":
                fst = _weir_cockerham_fst(p1, p2, n1, n2)
            else:
                fst = _hudson_fst(p1, p2, n1, n2)

            fst_values.append(fst)

        result[key] = fst_values

    return result


def _hudson_fst(p1: float, p2: float, n1: int, n2: int) -> float:
    """Hudson Fst estimator: Fst = 1 - (Hw / Hb)."""
    p_mean = (n1 * p1 + n2 * p2) / (n1 + n2)
    h_total = 2 * p_mean * (1 - p_mean)
    h_within = (n1 * 2 * p1 * (1 - p1) + n2 * 2 * p2 * (1 - p2)) / (n1 + n2)
    if h_total > 0:
        return max(0.0, (h_total - h_within) / h_total)
    return 0.0


def _weir_cockerham_fst(p1: float, p2: float, n1: int, n2: int) -> float:
    """Weir & Cockerham (1984) Fst with finite sample correction.

    Uses the ANOVA-based estimator theta-hat that properly accounts
    for unequal sample sizes.

    Reference: Weir & Cockerham (1984) Evolution 38:1358-1370.
    """
    r = 2  # number of populations
    n_total = n1 + n2
    if n_total < 2 or n1 < 1 or n2 < 1:
        return 0.0

    # Sample sizes (diploids: allele count = 2n)
    n_bar = n_total / r
    # Coefficient of variation of sample sizes
    n_c = (n_total - (n1**2 + n2**2) / n_total) / (r - 1)

    # Sample allele frequencies (already computed as allele freqs)
    p_bar = (n1 * p1 + n2 * p2) / n_total

    # Mean square components
    # s² = sample variance of allele frequencies across populations
    s_sq = (n1 * (p1 - p_bar) ** 2 + n2 * (p2 - p_bar) ** 2) / ((r - 1) * n_bar)

    # h_bar = average within-population heterozygosity
    h1 = 2 * p1 * (1 - p1)
    h2 = 2 * p2 * (1 - p2)
    h_bar = (n1 * h1 + n2 * h2) / n_total

    # Variance components (Weir & Cockerham eq. 2, 3, 4)
    a = (n_bar / n_c) * (
        s_sq
        - (1 / (n_bar - 1)) * (p_bar * (1 - p_bar) - ((r - 1) / r) * s_sq - h_bar / 4)
    )
    b = (n_bar / (n_bar - 1)) * (
        p_bar * (1 - p_bar)
        - ((r - 1) / r) * s_sq
        - ((2 * n_bar - 1) / (4 * n_bar)) * h_bar
    )
    c = h_bar / 2

    denom = a + b + c
    if denom <= 0:
        return 0.0

    return max(0.0, a / denom)


def strain_specific_variants(
    genotype_matrix: List[List[int]],
    sample_ids: List[str],
    maf_threshold: float = 0.05,
) -> Dict[str, List[int]]:
    """Identify variants that are private to or highly differentiated in each strain.

    A variant is considered strain-private if its MAF > threshold in that strain
    and MAF < threshold in all other strains.

    Returns:
        Dict mapping strain -> list of variant indices
    """
    af_by_strain = compute_allele_frequencies_by_strain(genotype_matrix, sample_ids)
    strains = sorted(af_by_strain.keys())
    n_variants = len(genotype_matrix)
    result: Dict[str, List[int]] = {s: [] for s in strains}

    for v in range(n_variants):
        for target_strain in strains:
            target_af = af_by_strain[target_strain][v]
            target_maf = (
                min(target_af, 1 - target_af) if not math.isnan(target_af) else 0
            )

            if target_maf < maf_threshold:
                continue

            # Check that all other strains have low MAF
            is_private = True
            for other_strain in strains:
                if other_strain == target_strain:
                    continue
                other_af = af_by_strain[other_strain][v]
                other_maf = (
                    min(other_af, 1 - other_af) if not math.isnan(other_af) else 0
                )
                if other_maf >= maf_threshold:
                    is_private = False
                    break

            if is_private:
                result[target_strain].append(v)

    for strain, variants in result.items():
        logger.info(
            f"Strain {strain}: {len(variants)} private variants (MAF>{maf_threshold})"
        )

    return result


def compute_global_fst(
    genotype_matrix: List[List[int]],
    sample_ids: List[str],
) -> Dict[str, Any]:
    """Compute genome-wide mean Fst across all strain pairs.

    Returns:
        Dict with per-pair mean Fst and overall mean
    """
    pairwise_fst = compute_fst_per_variant(genotype_matrix, sample_ids)
    summary: Dict[str, Any] = {"pairwise": {}}

    all_means = []
    for pair_key, fst_values in pairwise_fst.items():
        valid = [f for f in fst_values if not math.isnan(f)]
        if valid:
            mean_fst = sum(valid) / len(valid)
            summary["pairwise"][pair_key] = {
                "mean_fst": round(mean_fst, 6),
                "n_variants": len(valid),
                "max_fst": round(max(valid), 6),
            }
            all_means.append(mean_fst)

    summary["global_mean_fst"] = (
        round(sum(all_means) / len(all_means), 6) if all_means else 0.0
    )
    return summary


def run_strain_analysis(
    gwas_results: Dict[str, Any],
    *,
    vcf_path: str,
    output_dir: str,
    results_dir: str,
    sample_ids: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Full strain-aware analysis orchestration.

    Performs:
      1. Strain-aware PCA grid
      2. Dendrogram + clustered kinship heatmap
      3. Pairwise Fst (per-variant + global mean)
      4. Allele frequency heatmap (top N differentiating loci)
      5. Strain-specific (private) variant detection

    Imported from ``run_partial_gwas.py`` and ``assess_gwas_partial.py``
    to eliminate duplication.

    Args:
        gwas_results: Dict returned by ``execute_gwas_workflow`` (must have
            ``outputs`` key with pca/kinship/association_results).
        vcf_path: Path to the VCF used for the GWAS run.
        output_dir: Directory containing pca_results.json, kinship_results.json.
        results_dir: Base results directory for strain_analysis/ subdirectory.
        sample_ids: Sample IDs (extracted from VCF if ``None``).

    Returns:
        Dict with paths to all generated outputs and summary statistics.
    """
    import json as _json
    from pathlib import Path as _Path

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as _plt
    except ImportError:
        logger.warning("matplotlib not available — skipping strain plots")
        return {"error": "matplotlib not available"}

    output_dir_p = _Path(output_dir)
    strain_dir = _Path(results_dir) / "strain_analysis"
    strain_dir.mkdir(parents=True, exist_ok=True)

    gwas_results.get("outputs", {})
    report: Dict[str, Any] = {"strain_dir": str(strain_dir)}

    # Resolve sample IDs from VCF header if not provided
    if sample_ids is None:
        from metainformant.gwas.data.vcf_utils import extract_sample_ids

        sample_ids = extract_sample_ids(vcf_path)
    logger.info(f"Strain analysis: {len(sample_ids)} samples from VCF")

    # ── 1. Strain-aware PCA ──
    pca_file = output_dir_p / "pca_results.json"
    if pca_file.exists() and sample_ids:
        try:
            from metainformant.gwas.visualization.strain_plots import strain_pca_plot

            with open(pca_file) as f:
                pca_data = _json.load(f)

            components = pca_data.get("components", [])
            variance = pca_data.get("variance", pca_data.get("explained_variance", []))

            if components and isinstance(components[0], list):
                n_pcs = len(components)
                n_samples = len(components[0])
                pcs_sample_major = [
                    [components[pc][s] for pc in range(n_pcs)] for s in range(n_samples)
                ]
            else:
                pcs_sample_major = pca_data.get("pcs", [])
                n_samples = len(pcs_sample_major)

            effective_ids = sample_ids[:n_samples]
            pca_dict = {"pcs": pcs_sample_major, "explained_variance_ratio": variance}
            fig = strain_pca_plot(
                pca_dict, effective_ids, output_path=strain_dir / "strain_pca_grid.png"
            )
            if fig:
                _plt.close(fig)
                report["strain_pca"] = str(strain_dir / "strain_pca_grid.png")
                logger.info("Saved strain PCA grid")
        except Exception as e:
            logger.warning(f"Strain PCA failed: {e}", exc_info=True)

    # ── 2. Dendrogram + Clustered Kinship ──
    kinship_file = output_dir_p / "kinship_results.json"
    if kinship_file.exists() and sample_ids:
        try:
            from metainformant.gwas.visualization.strain_plots import (
                dendrogram_plot,
                kinship_heatmap_clustered,
            )

            with open(kinship_file) as f:
                kinship_matrix = _json.load(f)

            n_kin = len(kinship_matrix)
            effective_ids = sample_ids[:n_kin]

            fig = dendrogram_plot(
                kinship_matrix, effective_ids, output_path=strain_dir / "dendrogram.png"
            )
            if fig:
                _plt.close(fig)
                report["dendrogram"] = str(strain_dir / "dendrogram.png")

            fig = kinship_heatmap_clustered(
                kinship_matrix,
                effective_ids,
                output_path=strain_dir / "kinship_clustered.png",
            )
            if fig:
                _plt.close(fig)
                report["kinship_clustered"] = str(strain_dir / "kinship_clustered.png")
        except Exception as e:
            logger.warning(f"Dendrogram/kinship failed: {e}", exc_info=True)

    # ── 3. Fst + AF Heatmap + Private Variants ──
    if sample_ids:
        try:
            from metainformant.gwas.visualization.strain_plots import (
                fst_manhattan_plot,
                allele_frequency_heatmap,
            )
            from metainformant.gwas.analysis.quality import parse_vcf_full

            logger.info(f"Re-parsing VCF for Fst: {vcf_path}")
            vcf_data = parse_vcf_full(vcf_path)
            genotypes = vcf_data.get("genotypes", [])
            variants_info = vcf_data.get("variants", [])
            vcf_samples = vcf_data.get("samples", [])

            if genotypes:
                # Transpose if needed: (samples × variants) → (variants × samples)
                if genotypes and len(genotypes[0]) != len(vcf_samples):
                    geno_vm = genotypes
                else:
                    n_s = len(genotypes)
                    n_v = len(genotypes[0]) if genotypes else 0
                    geno_vm = [
                        [genotypes[s][v] for s in range(n_s)] for v in range(n_v)
                    ]

                effective_ids = vcf_samples if vcf_samples else sample_ids

                # Global Fst
                fst_summary = compute_global_fst(geno_vm, effective_ids)
                with open(strain_dir / "fst_summary.json", "w") as f:
                    _json.dump(fst_summary, f, indent=2)
                report["fst_summary"] = fst_summary
                logger.info(
                    f"Global mean Fst: {fst_summary.get('global_mean_fst', 'N/A')}"
                )

                # Fst Manhattan
                fst_data = compute_fst_per_variant(geno_vm, effective_ids)
                fig = fst_manhattan_plot(
                    fst_data,
                    variants_info,
                    output_path=strain_dir / "fst_manhattan.png",
                )
                if fig:
                    _plt.close(fig)
                    report["fst_manhattan"] = str(strain_dir / "fst_manhattan.png")

                # AF heatmap
                af_data = compute_allele_frequencies_by_strain(geno_vm, effective_ids)
                fig = allele_frequency_heatmap(
                    af_data,
                    variants_info,
                    output_path=strain_dir / "af_heatmap.png",
                    top_n=50,
                )
                if fig:
                    _plt.close(fig)
                    report["af_heatmap"] = str(strain_dir / "af_heatmap.png")

                # Private variants
                private = strain_specific_variants(geno_vm, effective_ids)
                with open(strain_dir / "strain_private_variants.json", "w") as f:
                    _json.dump(private, f, indent=2)
                report["private_variants"] = {k: len(v) for k, v in private.items()}

        except Exception as e:
            logger.warning(f"Fst analysis failed: {e}", exc_info=True)

    logger.info(f"Strain analysis outputs in {strain_dir}")
    return report
