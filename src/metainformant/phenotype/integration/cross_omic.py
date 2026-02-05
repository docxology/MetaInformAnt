"""Cross-omic integration functions connecting phenotype to DNA, RNA, GWAS.

This module provides functions for integrating phenotypic data with:
- Genotypic variants (GWAS-style associations)
- Gene expression data (trait-expression correlations)
- Environmental factors (GxE interactions)
- Multiple phenotype domains (multi-trait analysis)
"""

from __future__ import annotations

import math
import statistics
from typing import Dict, List, Any, Optional, Tuple, Union
from collections import defaultdict

from metainformant.core.utils.errors import ValidationError


def phenotype_genotype_association(
    phenotypes: Dict[str, float],
    genotypes: Dict[str, List[int]],
    method: str = "linear_regression",
    covariates: Optional[Dict[str, List[float]]] = None,
) -> Dict[str, Any]:
    """Test association between phenotypes and genetic variants.

    Implements basic association testing similar to GWAS.
    For production use, integrate with metainformant.gwas module.

    Args:
        phenotypes: Dict mapping sample_id to phenotype value.
        genotypes: Dict mapping variant_id to list of genotype calls (0/1/2) per sample.
        method: 'linear_regression' or 'correlation'.
        covariates: Optional covariates to include.

    Returns:
        Dict with per-variant association statistics.
    """
    if not phenotypes or not genotypes:
        return {"error": "Empty input data", "associations": {}}

    sample_ids = list(phenotypes.keys())
    n_samples = len(sample_ids)
    pheno_values = [phenotypes[s] for s in sample_ids]

    results = {"n_samples": n_samples, "n_variants": len(genotypes), "method": method, "associations": {}}

    for variant_id, geno_calls in genotypes.items():
        if len(geno_calls) != n_samples:
            continue

        # Filter to non-missing
        valid = [(pheno_values[i], geno_calls[i]) for i in range(n_samples) if geno_calls[i] >= 0]
        if len(valid) < 3:
            continue

        y_vals = [v[0] for v in valid]
        x_vals = [float(v[1]) for v in valid]

        if method == "correlation":
            r, p = _pearson_correlation(x_vals, y_vals)
            results["associations"][variant_id] = {
                "correlation": r,
                "p_value": p,
                "n_valid": len(valid),
            }
        else:
            # Linear regression
            slope, intercept, r_sq, p = _linear_regression_stats(x_vals, y_vals)
            results["associations"][variant_id] = {
                "beta": slope,
                "intercept": intercept,
                "r_squared": r_sq,
                "p_value": p,
                "n_valid": len(valid),
            }

    return results


def trait_expression_correlation(
    trait_values: Dict[str, float],
    expression_matrix: Dict[str, Dict[str, float]],
    gene_list: Optional[List[str]] = None,
    method: str = "pearson",
) -> Dict[str, Any]:
    """Correlate phenotypic traits with gene expression levels.

    Args:
        trait_values: Dict mapping sample_id to trait value.
        expression_matrix: Dict mapping gene_id to {sample_id: expression}.
        gene_list: Optional subset of genes to test.
        method: 'pearson' or 'spearman'.

    Returns:
        Dict with per-gene correlation statistics.
    """
    if not trait_values or not expression_matrix:
        return {"error": "Empty input data", "correlations": {}}

    sample_ids = list(trait_values.keys())
    genes_to_test = gene_list if gene_list else list(expression_matrix.keys())

    results = {
        "n_samples": len(sample_ids),
        "n_genes_tested": len(genes_to_test),
        "method": method,
        "correlations": {},
    }

    for gene_id in genes_to_test:
        if gene_id not in expression_matrix:
            continue

        gene_expr = expression_matrix[gene_id]

        # Align samples
        paired = [(trait_values[s], gene_expr[s]) for s in sample_ids if s in gene_expr]
        if len(paired) < 3:
            continue

        traits = [p[0] for p in paired]
        exprs = [p[1] for p in paired]

        if method == "spearman":
            r, p = _spearman_correlation(traits, exprs)
        else:
            r, p = _pearson_correlation(traits, exprs)

        results["correlations"][gene_id] = {
            "correlation": r,
            "p_value": p,
            "n_samples": len(paired),
        }

    # Sort by absolute correlation
    sorted_genes = sorted(results["correlations"].items(), key=lambda x: abs(x[1]["correlation"]), reverse=True)
    results["top_genes"] = [g[0] for g in sorted_genes[:10]]

    return results


def multi_phenotype_integration(
    phenotype_matrices: Dict[str, Dict[str, float]],
    method: str = "correlation",
) -> Dict[str, Any]:
    """Integrate multiple phenotype domains to find cross-phenotype relationships.

    Args:
        phenotype_matrices: Dict mapping phenotype_name to {sample_id: value}.
        method: 'correlation' or 'pca'.

    Returns:
        Dict with cross-phenotype correlation matrix and clusters.
    """
    if len(phenotype_matrices) < 2:
        return {"error": "Need at least 2 phenotypes", "correlations": {}}

    phenotype_names = list(phenotype_matrices.keys())
    n_phenotypes = len(phenotype_names)

    # Find common samples
    all_samples = [set(pm.keys()) for pm in phenotype_matrices.values()]
    common_samples = set.intersection(*all_samples)

    if len(common_samples) < 3:
        return {"error": "Insufficient common samples", "n_common": len(common_samples)}

    common_samples = sorted(common_samples)

    # Build correlation matrix
    correlation_matrix: List[List[float]] = [[0.0] * n_phenotypes for _ in range(n_phenotypes)]

    for i in range(n_phenotypes):
        for j in range(i, n_phenotypes):
            if i == j:
                correlation_matrix[i][j] = 1.0
            else:
                vals_i = [phenotype_matrices[phenotype_names[i]][s] for s in common_samples]
                vals_j = [phenotype_matrices[phenotype_names[j]][s] for s in common_samples]
                r, _ = _pearson_correlation(vals_i, vals_j)
                correlation_matrix[i][j] = r
                correlation_matrix[j][i] = r

    # Find strongly correlated pairs
    strong_pairs = []
    for i in range(n_phenotypes):
        for j in range(i + 1, n_phenotypes):
            if abs(correlation_matrix[i][j]) > 0.5:
                strong_pairs.append(
                    {
                        "phenotype_a": phenotype_names[i],
                        "phenotype_b": phenotype_names[j],
                        "correlation": correlation_matrix[i][j],
                    }
                )

    return {
        "n_phenotypes": n_phenotypes,
        "n_common_samples": len(common_samples),
        "phenotype_names": phenotype_names,
        "correlation_matrix": correlation_matrix,
        "strong_pairs": strong_pairs,
    }


def phenotype_environment_interaction(
    phenotypes: Dict[str, float],
    genotypes: Dict[str, List[int]],
    environment: Dict[str, float],
    interaction_model: str = "multiplicative",
) -> Dict[str, Any]:
    """Test for genotype-by-environment (GxE) interactions on phenotype.

    Args:
        phenotypes: Dict mapping sample_id to phenotype value.
        genotypes: Dict mapping variant_id to list of genotype calls per sample.
        environment: Dict mapping sample_id to environmental variable.
        interaction_model: 'multiplicative' or 'additive'.

    Returns:
        Dict with GxE interaction statistics.
    """
    if not phenotypes or not genotypes or not environment:
        return {"error": "Empty input data", "interactions": {}}

    sample_ids = [s for s in phenotypes.keys() if s in environment]
    n_samples = len(sample_ids)

    if n_samples < 5:
        return {"error": "Insufficient samples with both phenotype and environment"}

    results = {
        "n_samples": n_samples,
        "n_variants": len(genotypes),
        "model": interaction_model,
        "interactions": {},
    }

    pheno_vals = [phenotypes[s] for s in sample_ids]
    env_vals = [environment[s] for s in sample_ids]

    for variant_id, geno_calls in genotypes.items():
        if len(geno_calls) != len(phenotypes):
            continue

        # Align genotypes with sample_ids
        geno_map = dict(zip(phenotypes.keys(), geno_calls))
        geno_aligned = [geno_map.get(s, -1) for s in sample_ids]

        # Filter valid
        valid_idx = [i for i in range(n_samples) if geno_aligned[i] >= 0]
        if len(valid_idx) < 5:
            continue

        y = [pheno_vals[i] for i in valid_idx]
        g = [float(geno_aligned[i]) for i in valid_idx]
        e = [env_vals[i] for i in valid_idx]

        # Create interaction term
        if interaction_model == "multiplicative":
            gxe = [g[i] * e[i] for i in range(len(g))]
        else:
            gxe = [g[i] + e[i] for i in range(len(g))]

        # Test interaction term
        r, p = _pearson_correlation(gxe, y)

        results["interactions"][variant_id] = {
            "gxe_correlation": r,
            "p_value": p,
            "n_valid": len(valid_idx),
        }

    # Find significant interactions
    sig_interactions = [{"variant": v, **data} for v, data in results["interactions"].items() if data["p_value"] < 0.05]
    results["significant_interactions"] = sorted(sig_interactions, key=lambda x: x["p_value"])

    return results


# ============== Helper Statistics Functions ==============


def _pearson_correlation(x: List[float], y: List[float]) -> Tuple[float, float]:
    """Pearson correlation coefficient with p-value approximation."""
    n = len(x)
    if n < 3 or len(y) != n:
        return 0.0, 1.0

    mean_x = sum(x) / n
    mean_y = sum(y) / n

    cov_xy = sum((x[i] - mean_x) * (y[i] - mean_y) for i in range(n))
    var_x = sum((xi - mean_x) ** 2 for xi in x)
    var_y = sum((yi - mean_y) ** 2 for yi in y)

    if var_x == 0 or var_y == 0:
        return 0.0, 1.0

    r = cov_xy / math.sqrt(var_x * var_y)

    # P-value approximation using t-distribution
    if abs(r) >= 1.0:
        return r, 0.0 if r != 0 else 1.0

    t_stat = r * math.sqrt((n - 2) / (1 - r**2))
    # Rough p-value approximation
    p_value = 2 * (1 - _t_cdf(abs(t_stat), n - 2))

    return r, max(0.0, min(1.0, p_value))


def _spearman_correlation(x: List[float], y: List[float]) -> Tuple[float, float]:
    """Spearman rank correlation."""
    n = len(x)
    if n < 3 or len(y) != n:
        return 0.0, 1.0

    # Convert to ranks
    rank_x = _to_ranks(x)
    rank_y = _to_ranks(y)

    return _pearson_correlation(rank_x, rank_y)


def _to_ranks(values: List[float]) -> List[float]:
    """Convert values to ranks."""
    indexed = [(v, i) for i, v in enumerate(values)]
    sorted_vals = sorted(indexed, key=lambda x: x[0])
    ranks = [0.0] * len(values)
    for rank, (_, orig_idx) in enumerate(sorted_vals):
        ranks[orig_idx] = float(rank + 1)
    return ranks


def _linear_regression_stats(x: List[float], y: List[float]) -> Tuple[float, float, float, float]:
    """Simple linear regression with stats."""
    n = len(x)
    if n < 3:
        return 0.0, 0.0, 0.0, 1.0

    mean_x = sum(x) / n
    mean_y = sum(y) / n

    ss_xy = sum((x[i] - mean_x) * (y[i] - mean_y) for i in range(n))
    ss_xx = sum((xi - mean_x) ** 2 for xi in x)
    ss_yy = sum((yi - mean_y) ** 2 for yi in y)

    if ss_xx == 0:
        return 0.0, mean_y, 0.0, 1.0

    slope = ss_xy / ss_xx
    intercept = mean_y - slope * mean_x

    r_squared = (ss_xy**2) / (ss_xx * ss_yy) if ss_yy > 0 else 0.0

    # P-value from F-statistic
    if r_squared < 1.0 and r_squared > 0:
        f_stat = (r_squared / 1) / ((1 - r_squared) / (n - 2))
        p_value = 1 - _f_cdf(f_stat, 1, n - 2)
    else:
        p_value = 0.0 if r_squared > 0 else 1.0

    return slope, intercept, r_squared, max(0.0, min(1.0, p_value))


def _t_cdf(t: float, df: int) -> float:
    """Approximate t-distribution CDF."""
    # Simple approximation using normal for large df
    if df > 30:
        return _normal_cdf(t)
    # Very rough approximation for small df
    x = df / (df + t**2)
    return 0.5 + 0.5 * (1 - x ** (df / 2)) * (1 if t >= 0 else -1)


def _f_cdf(f: float, df1: int, df2: int) -> float:
    """Approximate F-distribution CDF."""
    if f <= 0:
        return 0.0
    # Simple approximation
    x = df2 / (df2 + df1 * f)
    return 1 - x ** (df2 / 2)


def _normal_cdf(z: float) -> float:
    """Approximate standard normal CDF using error function approximation."""
    return 0.5 * (1 + math.erf(z / math.sqrt(2)))
