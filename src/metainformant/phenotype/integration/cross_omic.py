"""Cross-omic integration functions connecting phenotype to DNA, RNA, GWAS.

This module provides functions for integrating phenotypic data with:
- Genotypic variants (GWAS-style associations)
- Gene expression data (trait-expression correlations)
- Environmental factors (GxE interactions)
- Multiple phenotype domains (multi-trait analysis)
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy import stats as scipy_stats


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
        genotypes: Dict mapping variant_id to list of genotype calls (0/1/2) per sample, in the order of
            ``phenotypes``. A negative call is missing.
        method: 'linear_regression' or 'correlation'.
        covariates: Optional covariates (``method='linear_regression'`` only): dict mapping covariate name to its
            list of values per sample, in the order of ``phenotypes``. They enter the model
            (``phenotype ~ 1 + covariates + genotype``) with the genotype, whose p-value then uses the degrees
            of freedom they leave; ``r_squared`` is that of the full model. A sample with a missing covariate
            (NaN) is left out.

    Returns:
        Dict with per-variant association statistics under ``"associations"`` and, under ``"skipped_variants"``,
        the reason for each variant that could not be tested (a genotype list of the wrong length, fewer than
        3 complete samples, or no variation in the genotype among them, which leaves no slope to estimate).

    Raises:
        ValueError: If ``method`` is not one of the two above, a covariate does not have one value per sample, or
            ``covariates`` is combined with ``method='correlation'``.
    """
    if method not in ("linear_regression", "correlation"):
        raise ValueError(f"Unknown method {method!r}. Use 'linear_regression' or 'correlation'.")
    if covariates and method == "correlation":
        raise ValueError("covariates are only supported with method='linear_regression'")

    if not phenotypes or not genotypes:
        return {"error": "Empty input data", "associations": {}}

    sample_ids = list(phenotypes.keys())
    n_samples = len(sample_ids)
    pheno_values = [phenotypes[s] for s in sample_ids]
    covariate_rows: Optional[List[List[float]]] = None
    if covariates:
        for name, values in covariates.items():
            if len(values) != n_samples:
                raise ValueError(f"Covariate {name!r} has {len(values)} values for {n_samples} samples")
        covariate_rows = [[float(covariates[name][i]) for name in covariates] for i in range(n_samples)]

    results: Dict[str, Any] = {
        "n_samples": n_samples,
        "n_variants": len(genotypes),
        "method": method,
        "associations": {},
        "skipped_variants": {},
    }

    for variant_id, geno_calls in genotypes.items():
        if len(geno_calls) != n_samples:
            results["skipped_variants"][variant_id] = f"{len(geno_calls)} genotype calls for {n_samples} samples"
            continue

        # Filter to non-missing
        keep = [
            i
            for i in range(n_samples)
            if geno_calls[i] >= 0 and (covariate_rows is None or not any(math.isnan(c) for c in covariate_rows[i]))
        ]
        if len(keep) < 3:
            results["skipped_variants"][variant_id] = f"{len(keep)} complete samples"
            continue

        y_vals = [pheno_values[i] for i in keep]
        x_vals = [float(geno_calls[i]) for i in keep]
        if len(set(x_vals)) < 2:
            results["skipped_variants"][variant_id] = "no variation in the genotype among the complete samples"
            continue

        if method == "correlation":
            r, p = _pearson_correlation(x_vals, y_vals)
            results["associations"][variant_id] = {
                "correlation": r,
                "p_value": p,
                "n_valid": len(keep),
            }
        else:
            # Linear regression
            if covariate_rows is None:
                slope, intercept, r_sq, p = _linear_regression_stats(x_vals, y_vals)
            else:
                fit = _covariate_adjusted_fit(y_vals, x_vals, [covariate_rows[i] for i in keep])
                if fit is None:
                    results["skipped_variants"][
                        variant_id
                    ] = "the genotype adds nothing to the covariates, or too few samples for them"
                    continue
                slope, intercept, r_sq, p, _ = fit
            results["associations"][variant_id] = {
                "beta": slope,
                "intercept": intercept,
                "r_squared": r_sq,
                "p_value": p,
                "n_valid": len(keep),
            }

    return results


def _covariate_adjusted_fit(
    y: List[float], x: List[float], covariates: List[List[float]]
) -> Optional[Tuple[float, float, float, float, float]]:
    """OLS of y on an intercept, the covariates and x.

    Returns ``(slope of x, intercept, model R-squared, p-value of x, standard error of x)``.

    The p-value is the two-sided Student's t test of x's coefficient on ``n - rank`` degrees of freedom. Returns
    ``None`` when x adds no information beyond the intercept and covariates or there are no spare degrees of freedom.
    """
    n = len(y)
    target = np.asarray(y, dtype=float)
    design = np.column_stack([np.ones(n), np.asarray(covariates, dtype=float), np.asarray(x, dtype=float)])
    rank = int(np.linalg.matrix_rank(design))
    reduced_rank = int(np.linalg.matrix_rank(design[:, :-1]))
    df = n - rank
    if df < 1 or rank == reduced_rank:
        return None
    coefficients, *_ = np.linalg.lstsq(design, target, rcond=None)
    residuals = target - design @ coefficients
    ss_res = float(residuals @ residuals)
    ss_tot = float(((target - target.mean()) ** 2).sum())
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    # Factor the design directly: forming X.T @ X squares its condition number.
    design_inverse = np.linalg.pinv(design)
    cov_beta = (ss_res / df) * (design_inverse @ design_inverse.T)
    se = math.sqrt(max(float(cov_beta[-1, -1]), 0.0))
    slope = float(coefficients[-1])
    if se == 0.0:
        p_value = 0.0 if slope != 0.0 else 1.0
    else:
        p_value = float(2.0 * scipy_stats.t.sf(abs(slope / se), df=df))
    return slope, float(coefficients[0]), r_squared, max(0.0, min(1.0, p_value)), se


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
        Dict with per-gene correlation statistics under ``"correlations"``, the ten genes with the largest
        absolute correlation under ``"top_genes"``, and under ``"skipped_genes"`` the reason for each gene that
        could not be correlated: fewer than 3 samples shared with the trait, or no variation in the trait or in
        the gene's expression among them (a correlation with a constant is undefined, not 0 with p = 1).

    Raises:
        ValueError: If ``method`` is neither 'pearson' nor 'spearman'.
    """
    if method not in ("pearson", "spearman"):
        raise ValueError(f"Unknown method {method!r}. Use 'pearson' or 'spearman'.")
    if not trait_values or not expression_matrix:
        return {"error": "Empty input data", "correlations": {}}

    sample_ids = list(trait_values.keys())
    genes_to_test = gene_list if gene_list else list(expression_matrix.keys())

    results: Dict[str, Any] = {
        "n_samples": len(sample_ids),
        "n_genes_tested": len(genes_to_test),
        "method": method,
        "correlations": {},
        "skipped_genes": {},
    }

    for gene_id in genes_to_test:
        if gene_id not in expression_matrix:
            results["skipped_genes"][gene_id] = "not in the expression matrix"
            continue

        gene_expr = expression_matrix[gene_id]

        # Align samples
        paired = [(trait_values[s], gene_expr[s]) for s in sample_ids if s in gene_expr]
        if len(paired) < 3:
            results["skipped_genes"][gene_id] = f"{len(paired)} samples shared with the trait"
            continue

        traits = [p[0] for p in paired]
        exprs = [p[1] for p in paired]
        if len(set(traits)) < 2 or len(set(exprs)) < 2:
            which = "trait" if len(set(traits)) < 2 else "expression"
            results["skipped_genes"][gene_id] = f"no variation in the {which} among the shared samples"
            continue

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
        method: 'correlation' or 'pca'. Both return the pairwise Pearson correlations and the strongly
            correlated pairs (|r| > 0.5); 'pca' adds a principal component analysis of the standardised
            phenotypes under ``"pca"``: ``eigenvalues`` and ``explained_variance_ratio`` of the components,
            ``loadings`` (correlation of each phenotype with each component) and ``scores`` per sample.

    Returns:
        Dict with the cross-phenotype correlation matrix (``nan`` for a pair with a constant phenotype, whose
        correlation is undefined) and the strongly correlated pairs.

    Raises:
        ValueError: If ``method`` is neither 'correlation' nor 'pca', or (for 'pca') a phenotype is constant
            over the common samples, so it cannot be standardised.
    """
    if method not in ("correlation", "pca"):
        raise ValueError(f"Unknown method {method!r}. Use 'correlation' or 'pca'.")
    if len(phenotype_matrices) < 2:
        return {"error": "Need at least 2 phenotypes", "correlations": {}}

    phenotype_names = list(phenotype_matrices.keys())
    n_phenotypes = len(phenotype_names)

    # Find common samples
    all_samples = [set(pm.keys()) for pm in phenotype_matrices.values()]
    common_samples = set.intersection(*all_samples)

    if len(common_samples) < 3:
        return {"error": "Insufficient common samples", "n_common": len(common_samples)}

    common_samples_list = sorted(common_samples)

    # Build correlation matrix
    correlation_matrix: List[List[float]] = [[0.0] * n_phenotypes for _ in range(n_phenotypes)]

    for i in range(n_phenotypes):
        for j in range(i, n_phenotypes):
            if i == j:
                correlation_matrix[i][j] = 1.0
            else:
                vals_i = [phenotype_matrices[phenotype_names[i]][s] for s in common_samples_list]
                vals_j = [phenotype_matrices[phenotype_names[j]][s] for s in common_samples_list]
                if len(set(vals_i)) < 2 or len(set(vals_j)) < 2:
                    r = math.nan  # no correlation with a constant
                else:
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

    result: Dict[str, Any] = {
        "n_phenotypes": n_phenotypes,
        "n_common_samples": len(common_samples),
        "phenotype_names": phenotype_names,
        "correlation_matrix": correlation_matrix,
        "strong_pairs": strong_pairs,
    }

    if method == "pca":
        data = np.array(
            [[phenotype_matrices[name][sample] for name in phenotype_names] for sample in common_samples_list],
            dtype=float,
        )
        spread = data.std(axis=0, ddof=1)
        if np.any(spread == 0.0):
            constant = [name for name, sd in zip(phenotype_names, spread) if sd == 0.0]
            raise ValueError(f"Cannot run PCA: phenotype(s) {constant} are constant over the common samples")
        standardised = (data - data.mean(axis=0)) / spread
        _, singular_values, components = np.linalg.svd(standardised, full_matrices=False)
        eigenvalues = singular_values**2 / (len(common_samples_list) - 1)
        loadings = components.T * np.sqrt(eigenvalues)
        scores = standardised @ components.T
        result["pca"] = {
            "eigenvalues": eigenvalues.tolist(),
            "explained_variance_ratio": (eigenvalues / eigenvalues.sum()).tolist(),
            "loadings": {name: loadings[k].tolist() for k, name in enumerate(phenotype_names)},
            "scores": {sample: scores[k].tolist() for k, sample in enumerate(common_samples_list)},
        }

    return result


def phenotype_environment_interaction(
    phenotypes: Dict[str, float],
    genotypes: Dict[str, List[int]],
    environment: Dict[str, float],
    interaction_model: str = "multiplicative",
) -> Dict[str, Any]:
    """Test for genotype-by-environment (GxE) interactions on phenotype.

    For each variant the phenotype is regressed on an intercept, the genotype, the environment and their
    product, ``y ~ 1 + g + e + g*e``, and the interaction is the coefficient of the product, tested with
    Student's t on ``n - 4`` degrees of freedom. Main effects of the genotype or the environment alone
    therefore do not register as interaction. (Correlating the product with the phenotype, which this function
    used to report as the interaction test, is significant whenever the genotype or the environment has a main
    effect, because the product is correlated with both.)

    Args:
        phenotypes: Dict mapping sample_id to phenotype value.
        genotypes: Dict mapping variant_id to list of genotype calls per sample, in the order of ``phenotypes``
            (a negative call is missing).
        environment: Dict mapping sample_id to environmental variable.
        interaction_model: ``'multiplicative'`` (the product term above). An ``'additive'`` combination ``g + e``
            contains no interaction term (it is the sum of the two main effects) and is rejected.

    Returns:
        Dict with ``"interactions"``: per variant ``gxe_beta`` and ``gxe_se`` (the interaction coefficient and its
        standard error), ``gxe_correlation`` (the descriptive Pearson correlation of the product with the
        phenotype), ``p_value`` (the interaction test) and ``n_valid``; ``"significant_interactions"`` (p < 0.05,
        sorted); and ``"skipped_variants"`` with the reason for each variant that could not be tested (a genotype
        list of the wrong length, fewer than 5 complete samples, or a design in which the product adds nothing).

    Raises:
        ValueError: If ``interaction_model`` is not ``'multiplicative'``.
    """
    if interaction_model != "multiplicative":
        raise ValueError(
            f"Unsupported interaction_model {interaction_model!r}: only 'multiplicative' (the product g*e) defines "
            "an interaction; 'additive' (g + e) would just sum the main effects"
        )
    if not phenotypes or not genotypes or not environment:
        return {"error": "Empty input data", "interactions": {}}

    sample_ids = [s for s in phenotypes.keys() if s in environment]
    n_samples = len(sample_ids)

    if n_samples < 5:
        return {"error": "Insufficient samples with both phenotype and environment"}

    results: Dict[str, Any] = {
        "n_samples": n_samples,
        "n_variants": len(genotypes),
        "model": interaction_model,
        "interactions": {},
        "skipped_variants": {},
    }

    pheno_vals = [phenotypes[s] for s in sample_ids]
    env_vals = [environment[s] for s in sample_ids]

    for variant_id, geno_calls in genotypes.items():
        if len(geno_calls) != len(phenotypes):
            results["skipped_variants"][variant_id] = f"{len(geno_calls)} genotype calls for {len(phenotypes)} samples"
            continue

        # Align genotypes with sample_ids
        geno_map = dict(zip(phenotypes.keys(), geno_calls))
        geno_aligned = [geno_map.get(s, -1) for s in sample_ids]

        # Filter valid
        valid_idx = [i for i in range(n_samples) if geno_aligned[i] >= 0]
        if len(valid_idx) < 5:
            results["skipped_variants"][variant_id] = f"{len(valid_idx)} complete samples"
            continue

        y = [pheno_vals[i] for i in valid_idx]
        g = [float(geno_aligned[i]) for i in valid_idx]
        e = [env_vals[i] for i in valid_idx]
        gxe = [g[i] * e[i] for i in range(len(g))]

        # Test the interaction term in y ~ 1 + g + e + g*e
        fit = _covariate_adjusted_fit(y, gxe, [[g[i], e[i]] for i in range(len(g))])
        if fit is None:
            results["skipped_variants"][
                variant_id
            ] = "the product of genotype and environment adds nothing to the two main effects"
            continue
        beta, _, _, p, se = fit
        r = _pearson_correlation(gxe, y)[0] if len(set(gxe)) > 1 else math.nan

        results["interactions"][variant_id] = {
            "gxe_beta": beta,
            "gxe_se": se,
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
    """Pearson correlation coefficient with its two-sided p-value.

    The p-value is that of the t statistic ``r * sqrt((n - 2) / (1 - r**2))`` under Student's t distribution
    with ``n - 2`` degrees of freedom (the test ``scipy.stats.pearsonr`` performs), whatever the sample size.
    """
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

    # Two-sided p-value from Student's t with n - 2 degrees of freedom
    if abs(r) >= 1.0:
        return r, 0.0 if r != 0 else 1.0

    t_stat = r * math.sqrt((n - 2) / (1 - r**2))
    p_value = float(2.0 * scipy_stats.t.sf(abs(t_stat), df=n - 2))

    return r, max(0.0, min(1.0, p_value))


def _spearman_correlation(x: List[float], y: List[float]) -> Tuple[float, float]:
    """Spearman rank correlation and its two-sided p-value (Student's t on the ranks, as ``scipy.stats.spearmanr``).

    Tied values get average ranks, so the result equals ``scipy.stats.spearmanr`` on data with ties.
    """
    n = len(x)
    if n < 3 or len(y) != n:
        return 0.0, 1.0

    # Convert to ranks
    rank_x = _to_ranks(x)
    rank_y = _to_ranks(y)

    return _pearson_correlation(rank_x, rank_y)


def _to_ranks(values: List[float]) -> List[float]:
    """Convert values to ranks, 1 for the smallest; tied values share the mean of the ranks they span.

    Averaging makes a rank depend on the value alone: ranking tied values by their position in the list (1, 2, 3
    for three equal values) would make Spearman's coefficient depend on the order of the samples, and on data
    with many ties, such as expression with many zeros, it would not be Spearman's coefficient at all.
    """
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    first = 0
    while first < len(order):
        last = first
        while last + 1 < len(order) and values[order[last + 1]] == values[order[first]]:
            last += 1
        average_rank = (first + last) / 2.0 + 1.0  # the mean of the 1-based ranks first + 1 .. last + 1
        for position in range(first, last + 1):
            ranks[order[position]] = average_rank
        first = last + 1
    return ranks


def _linear_regression_stats(x: List[float], y: List[float]) -> Tuple[float, float, float, float]:
    """Simple linear regression: slope, intercept, R-squared and the p-value of the slope.

    The p-value is that of the F statistic ``R**2 / ((1 - R**2) / (n - 2))`` with 1 and ``n - 2`` degrees of
    freedom, which equals the two-sided t-test of the slope (the p-value of ``scipy.stats.linregress``).
    """
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

    # P-value from the F distribution with (1, n - 2) degrees of freedom
    if r_squared < 1.0 and r_squared > 0:
        f_stat = (r_squared / 1) / ((1 - r_squared) / (n - 2))
        p_value = float(scipy_stats.f.sf(f_stat, 1, n - 2))
    else:
        p_value = 0.0 if r_squared > 0 else 1.0

    return slope, intercept, r_squared, max(0.0, min(1.0, p_value))
