"""Phenome-wide association study (PheWAS) and GWAS-phenotype integration.

Provides functions for phenome-wide association scans, phenotype correlation
analysis, polygenic/genetic risk score computation, heritability screening,
and phenotype categorisation. All statistical computations use pure Python
implementations with optional NumPy acceleration.

PheWAS tests a single genetic variant against many phenotypes simultaneously,
making it the inverse of a traditional GWAS which tests many variants against
a single phenotype.
"""

from __future__ import annotations

import math
import random
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _mean(values: list[float]) -> float:
    """Compute arithmetic mean."""
    if not values:
        return 0.0
    return sum(values) / len(values)


def _variance(values: list[float], ddof: int = 1) -> float:
    """Compute sample variance."""
    n = len(values)
    if n < 2:
        return 0.0
    mu = _mean(values)
    return sum((x - mu) ** 2 for x in values) / (n - ddof)


def _std(values: list[float], ddof: int = 1) -> float:
    """Compute sample standard deviation."""
    return math.sqrt(max(0.0, _variance(values, ddof)))


def _pearson_r(x: list[float], y: list[float]) -> tuple[float, float]:
    """Compute Pearson correlation and approximate p-value.

    Args:
        x: First variable values.
        y: Second variable values.

    Returns:
        Tuple of (correlation, p_value).
    """
    n = len(x)
    if n < 3:
        return 0.0, 1.0
    mx = _mean(x)
    my = _mean(y)
    sx = _std(x)
    sy = _std(y)
    if sx == 0.0 or sy == 0.0:
        return 0.0, 1.0
    r = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y)) / ((n - 1) * sx * sy)
    r = max(-1.0, min(1.0, r))
    # t-test approximation
    if abs(r) >= 1.0:
        return r, 0.0
    t_stat = r * math.sqrt((n - 2) / (1.0 - r * r))
    # Approximate two-tailed p from t using normal approximation for large n
    p_value = 2.0 * _normal_sf(abs(t_stat), df=n - 2)
    return r, p_value


def _normal_sf(t: float, df: int = 100) -> float:
    """Approximate survival function for t-distribution.

    Uses the normal approximation for df > 30; otherwise a rough
    Beta-function based approximation.
    """
    if df > 30:
        # Normal approximation
        return 0.5 * math.erfc(t / math.sqrt(2.0))
    # Rough approximation for smaller df
    x = df / (df + t * t)
    return 0.5 * _regularized_incomplete_beta(df / 2.0, 0.5, x)


def _regularized_incomplete_beta(a: float, b: float, x: float) -> float:
    """Rough continued-fraction approximation of I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    # Use a simple series expansion (first 60 terms)
    result = 0.0
    term = 1.0
    for n in range(60):
        if n == 0:
            term = (x**a) * ((1 - x) ** b) / a
            coeff = 1.0
        else:
            coeff *= (n - 1 - b) / n
            term *= x
        contrib = coeff * term / (a + n)
        result += contrib
        if abs(contrib) < 1e-12:
            break
    try:
        beta_val = math.gamma(a) * math.gamma(b) / math.gamma(a + b)
    except (OverflowError, ValueError):
        beta_val = 1.0
    if beta_val == 0:
        return 0.5
    return min(1.0, max(0.0, result / beta_val))


def _simple_linear_regression(
    y: list[float],
    x: list[float],
    covariates: list[list[float]] | None = None,
) -> tuple[float, float, float]:
    """Simple linear regression of y on x.

    Returns (beta, se, p_value). Covariates are residualised out first
    if provided.
    """
    n = len(y)
    if n < 3:
        return 0.0, 1.0, 1.0

    # If covariates, residualise y and x on covariates via OLS
    dep = list(y)
    pred = list(x)

    if covariates and len(covariates) > 0 and len(covariates[0]) > 0:
        dep = _residualise(dep, covariates)
        pred = _residualise(pred, covariates)

    mx = _mean(pred)
    my = _mean(dep)
    ss_xx = sum((xi - mx) ** 2 for xi in pred)
    if ss_xx == 0:
        return 0.0, 1.0, 1.0
    ss_xy = sum((xi - mx) * (yi - my) for xi, yi in zip(pred, dep))
    beta = ss_xy / ss_xx
    y_hat = [my + beta * (xi - mx) for xi in pred]
    residuals = [yi - yhi for yi, yhi in zip(dep, y_hat)]
    ss_res = sum(r * r for r in residuals)
    mse = ss_res / max(n - 2, 1)
    se = math.sqrt(max(mse / ss_xx, 1e-20))
    if se == 0:
        return beta, 1.0, 1.0
    t_stat = beta / se
    p_value = 2.0 * _normal_sf(abs(t_stat), df=n - 2)
    return beta, se, p_value


def _residualise(y: list[float], covariates: list[list[float]]) -> list[float]:
    """Residualise y on covariates using OLS (pure Python)."""
    n = len(y)
    k = len(covariates[0]) if covariates else 0
    if k == 0 or n < k + 1:
        return y

    # Build design matrix X (with intercept)
    # Use simple sequential residualisation for pure Python
    residuals = list(y)
    for col_idx in range(k):
        col = [covariates[i][col_idx] for i in range(n)]
        mc = _mean(col)
        mr = _mean(residuals)
        ss = sum((c - mc) ** 2 for c in col)
        if ss == 0:
            continue
        sp = sum((c - mc) * (r - mr) for c, r in zip(col, residuals))
        b = sp / ss
        residuals = [r - b * (c - mc) for r, c in zip(residuals, col)]
    return residuals


# ---------------------------------------------------------------------------
# PheWAS
# ---------------------------------------------------------------------------


def run_phewas(
    genotype: list[int],
    phenotypes: dict,
    covariates: list[list[float]] | None = None,
) -> list[dict]:
    """Run a phenome-wide association scan.

    Tests one genetic variant (encoded as 0/1/2 additive dosage) against
    multiple phenotypes using linear regression.

    Args:
        genotype: List of genotype values (0, 1, 2) for each sample.
        phenotypes: Dictionary mapping phenotype name to list of values
            (continuous or binary). Each list must have the same length as
            ``genotype``.
        covariates: Optional list of covariate vectors. Outer list is
            per-sample, inner list is per-covariate.

    Returns:
        List of result dicts sorted by p-value, each with keys:
            - phenotype: Name of the phenotype.
            - beta: Regression coefficient.
            - se: Standard error of beta.
            - p_value: Two-sided p-value.
            - n_samples: Number of non-missing samples.
            - phenotype_category: Category assigned by :func:`categorize_phenotypes`
              (``"uncategorized"`` if not categorised).

    Raises:
        ValueError: If genotype and phenotype lengths do not match.
    """
    n = len(genotype)
    results: list[dict] = []

    for pheno_name, pheno_values in phenotypes.items():
        if len(pheno_values) != n:
            logger.warning(
                "Phenotype '%s' length %d != genotype length %d; skipping",
                pheno_name,
                len(pheno_values),
                n,
            )
            continue

        # Filter missing values (None or NaN)
        valid_indices = []
        for i in range(n):
            g = genotype[i]
            p = pheno_values[i]
            if g is not None and p is not None:
                try:
                    if not math.isnan(float(p)) and not math.isnan(float(g)):
                        valid_indices.append(i)
                except (TypeError, ValueError):
                    continue

        if len(valid_indices) < 3:
            logger.debug("Phenotype '%s': too few valid samples (%d)", pheno_name, len(valid_indices))
            continue

        geno_valid = [float(genotype[i]) for i in valid_indices]
        pheno_valid = [float(pheno_values[i]) for i in valid_indices]
        cov_valid = None
        if covariates is not None:
            cov_valid = [[covariates[i][j] for j in range(len(covariates[0]))] for i in valid_indices]

        beta, se, p_value = _simple_linear_regression(pheno_valid, geno_valid, cov_valid)

        results.append(
            {
                "phenotype": pheno_name,
                "beta": beta,
                "se": se,
                "p_value": p_value,
                "n_samples": len(valid_indices),
                "phenotype_category": "uncategorized",
            }
        )

    # Sort by p-value
    results.sort(key=lambda r: r["p_value"])
    logger.info("PheWAS complete: tested %d phenotypes", len(results))
    return results


# ---------------------------------------------------------------------------
# Phenotype correlation
# ---------------------------------------------------------------------------


def phenotype_correlation_matrix(phenotypes: dict) -> dict:
    """Compute pairwise Pearson correlation matrix between all phenotypes.

    Args:
        phenotypes: Dictionary mapping phenotype name to list of numeric values.
            All lists must have the same length.

    Returns:
        Dictionary with keys:
            - correlation_matrix: 2D list of correlation coefficients.
            - p_value_matrix: 2D list of p-values.
            - phenotype_names: Ordered list of phenotype names.
            - n_samples: Number of samples used.

    Raises:
        ValueError: If fewer than 2 phenotypes or inconsistent lengths.
    """
    names = list(phenotypes.keys())
    if len(names) < 2:
        raise ValueError("Need at least 2 phenotypes for correlation matrix")

    values = [phenotypes[name] for name in names]
    n = len(values[0])
    for i, v in enumerate(values):
        if len(v) != n:
            raise ValueError(f"Phenotype '{names[i]}' has {len(v)} samples, expected {n}")

    k = len(names)
    corr_matrix = [[0.0] * k for _ in range(k)]
    p_matrix = [[0.0] * k for _ in range(k)]

    for i in range(k):
        corr_matrix[i][i] = 1.0
        p_matrix[i][i] = 0.0
        for j in range(i + 1, k):
            # Filter out None/NaN pairs
            xi = []
            xj = []
            for idx in range(n):
                vi = values[i][idx]
                vj = values[j][idx]
                if vi is not None and vj is not None:
                    try:
                        fi = float(vi)
                        fj = float(vj)
                        if not math.isnan(fi) and not math.isnan(fj):
                            xi.append(fi)
                            xj.append(fj)
                    except (TypeError, ValueError):
                        continue

            if len(xi) < 3:
                r, p = 0.0, 1.0
            else:
                r, p = _pearson_r(xi, xj)

            corr_matrix[i][j] = r
            corr_matrix[j][i] = r
            p_matrix[i][j] = p
            p_matrix[j][i] = p

    logger.info("Computed %dx%d phenotype correlation matrix", k, k)

    return {
        "correlation_matrix": corr_matrix,
        "p_value_matrix": p_matrix,
        "phenotype_names": names,
        "n_samples": n,
    }


# ---------------------------------------------------------------------------
# Genetic risk score
# ---------------------------------------------------------------------------


def genetic_risk_score(
    genotypes: Any,
    effect_sizes: list[float],
    weights: list[float] | None = None,
) -> dict:
    """Compute polygenic / genetic risk score (GRS) across individuals.

    The GRS for individual *i* is the weighted sum of genotype dosages
    times effect sizes:  ``GRS_i = sum(w_j * beta_j * g_ij)``

    Args:
        genotypes: 2D structure (list of lists or numpy array), shape
            ``(n_samples, n_variants)``. Each element is 0, 1, or 2.
        effect_sizes: List of effect sizes (beta values) per variant.
        weights: Optional per-variant weights. If ``None``, all weights
            are 1.0.

    Returns:
        Dictionary with keys:
            - risk_scores: List of GRS values per individual.
            - mean: Mean GRS.
            - std: Standard deviation of GRS.
            - percentiles: Dict of {5, 25, 50, 75, 95} percentile values.
            - n_variants: Number of variants used.
            - n_samples: Number of individuals.

    Raises:
        ValueError: If dimensions are inconsistent.
    """
    # Convert to list-of-lists if needed
    if HAS_NUMPY and isinstance(genotypes, np.ndarray):
        geno_list = genotypes.tolist()
    else:
        geno_list = [list(row) for row in genotypes]

    n_samples = len(geno_list)
    if n_samples == 0:
        raise ValueError("genotypes must contain at least one sample")

    n_variants = len(geno_list[0])
    if len(effect_sizes) != n_variants:
        raise ValueError(f"effect_sizes length {len(effect_sizes)} != n_variants {n_variants}")

    if weights is None:
        weights = [1.0] * n_variants
    elif len(weights) != n_variants:
        raise ValueError(f"weights length {len(weights)} != n_variants {n_variants}")

    scores: list[float] = []
    for i in range(n_samples):
        row = geno_list[i]
        grs = sum(weights[j] * effect_sizes[j] * float(row[j]) for j in range(n_variants))
        scores.append(grs)

    mu = _mean(scores)
    sd = _std(scores) if len(scores) > 1 else 0.0

    # Compute percentiles
    sorted_scores = sorted(scores)

    def _percentile(data: list[float], p: float) -> float:
        idx = (p / 100.0) * (len(data) - 1)
        low = int(math.floor(idx))
        high = int(math.ceil(idx))
        if low == high:
            return data[low]
        frac = idx - low
        return data[low] * (1 - frac) + data[high] * frac

    percentiles = {pct: _percentile(sorted_scores, pct) for pct in [5, 25, 50, 75, 95]}

    logger.info(
        "Computed GRS for %d samples x %d variants: mean=%.4f, std=%.4f",
        n_samples,
        n_variants,
        mu,
        sd,
    )

    return {
        "risk_scores": scores,
        "mean": mu,
        "std": sd,
        "percentiles": percentiles,
        "n_variants": n_variants,
        "n_samples": n_samples,
    }


# ---------------------------------------------------------------------------
# Heritability screening
# ---------------------------------------------------------------------------


def phenotype_heritability_screen(
    phenotypes: dict,
    kinship_matrix: Any,
) -> list[dict]:
    """Screen multiple phenotypes for heritability using Haseman-Elston regression.

    The Haseman-Elston method regresses squared phenotypic differences
    on kinship coefficients.  h^2 is estimated from the regression slope.

    Args:
        phenotypes: Dictionary mapping phenotype name to list of values.
        kinship_matrix: 2D kinship matrix (list of lists or numpy array),
            shape ``(n_samples, n_samples)``.

    Returns:
        List of dicts sorted by p-value, each with keys:
            - phenotype: Phenotype name.
            - h2: Estimated heritability (clamped to [0, 1]).
            - se: Standard error of h2.
            - p_value: P-value for h2 > 0.

    Raises:
        ValueError: If dimensions do not match.
    """
    # Convert kinship to list-of-lists
    if HAS_NUMPY and isinstance(kinship_matrix, np.ndarray):
        kin = kinship_matrix.tolist()
    else:
        kin = [list(row) for row in kinship_matrix]

    n = len(kin)
    results: list[dict] = []

    for pheno_name, pheno_values in phenotypes.items():
        if len(pheno_values) != n:
            logger.warning(
                "Phenotype '%s' length %d != kinship dimension %d; skipping",
                pheno_name,
                len(pheno_values),
                n,
            )
            continue

        # Filter valid values
        valid = [
            i
            for i in range(n)
            if pheno_values[i] is not None and not (isinstance(pheno_values[i], float) and math.isnan(pheno_values[i]))
        ]

        if len(valid) < 5:
            continue

        vals = [float(pheno_values[i]) for i in valid]
        pheno_var = _variance(vals)
        if pheno_var == 0:
            results.append(
                {
                    "phenotype": pheno_name,
                    "h2": 0.0,
                    "se": 0.0,
                    "p_value": 1.0,
                }
            )
            continue

        # Haseman-Elston: regress (Yi - Yj)^2 on kinship_ij
        diffs_sq: list[float] = []
        kin_vals: list[float] = []

        for ii in range(len(valid)):
            for jj in range(ii + 1, len(valid)):
                i_idx = valid[ii]
                j_idx = valid[jj]
                diff = vals[ii] - vals[jj]
                diffs_sq.append(diff * diff)
                kin_vals.append(float(kin[i_idx][j_idx]))

        if len(diffs_sq) < 3:
            continue

        beta, se, p_value = _simple_linear_regression(diffs_sq, kin_vals)

        # h2 = -beta / (2 * Var(Y)) for squared-difference HE regression
        h2 = -beta / (2.0 * pheno_var) if pheno_var > 0 else 0.0
        h2 = max(0.0, min(1.0, h2))
        h2_se = abs(se / (2.0 * pheno_var)) if pheno_var > 0 else 0.0

        results.append(
            {
                "phenotype": pheno_name,
                "h2": h2,
                "se": h2_se,
                "p_value": p_value,
            }
        )

    results.sort(key=lambda r: r["p_value"])
    logger.info("Heritability screening: %d phenotypes assessed", len(results))
    return results


# ---------------------------------------------------------------------------
# Phenotype categorisation
# ---------------------------------------------------------------------------

_PHECODE_CATEGORIES: dict[str, list[str]] = {
    "infectious": ["infection", "viral", "bacterial", "fungal", "parasit", "sepsis", "pneumonia"],
    "neoplasms": ["cancer", "tumor", "neoplasm", "carcinoma", "lymphoma", "leukemia", "melanoma"],
    "endocrine": ["diabetes", "thyroid", "adrenal", "pituitary", "hormone", "metabolic", "obesity"],
    "hematologic": ["anemia", "bleeding", "coagul", "platelet", "leukocyte", "hemoglobin"],
    "mental_health": ["depression", "anxiety", "schizophren", "bipolar", "psycho", "adhd", "autism"],
    "neurological": ["epilepsy", "seizure", "migraine", "headache", "neuropathy", "parkinson", "alzheimer"],
    "circulatory": ["hypertension", "heart", "cardiac", "coronary", "stroke", "atheroscler", "arrhythmia"],
    "respiratory": ["asthma", "copd", "bronch", "pulmonary", "lung", "respiratory"],
    "digestive": ["gastric", "hepat", "liver", "colitis", "crohn", "celiac", "pancreat", "bowel"],
    "genitourinary": ["renal", "kidney", "bladder", "prostate", "ovarian", "uterine"],
    "musculoskeletal": ["arthritis", "osteopor", "fracture", "bone", "joint", "muscle", "back_pain"],
    "dermatologic": ["dermat", "eczema", "psoriasis", "skin", "rash", "acne"],
    "congenital": ["congenital", "birth_defect", "chromosom", "genetic_disorder"],
    "laboratory": ["cholesterol", "glucose", "triglyceride", "creatinine", "alt", "ast", "bmi", "weight", "height"],
}


def categorize_phenotypes(
    phenotypes: dict,
    method: str = "phecode",
) -> dict:
    """Group phenotype names into categories using keyword matching.

    Args:
        phenotypes: Dictionary mapping phenotype name to any value (values
            are ignored; only keys are categorised).
        method: Categorisation method. Currently only ``"phecode"`` (keyword-
            based) is implemented.

    Returns:
        Dictionary with keys:
            - categories: Dict mapping category name to list of phenotype
              names in that category.
            - uncategorized: List of phenotype names that could not be
              categorised.
            - method: Method used.
            - n_categorized: Number of phenotypes successfully categorised.
    """
    categories: dict[str, list[str]] = {cat: [] for cat in _PHECODE_CATEGORIES}
    uncategorized: list[str] = []

    for pheno_name in phenotypes:
        name_lower = pheno_name.lower().replace(" ", "_")
        found = False
        for cat, keywords in _PHECODE_CATEGORIES.items():
            for kw in keywords:
                if kw in name_lower:
                    categories[cat].append(pheno_name)
                    found = True
                    break
            if found:
                break
        if not found:
            uncategorized.append(pheno_name)

    # Remove empty categories
    categories = {k: v for k, v in categories.items() if v}

    n_categorized = sum(len(v) for v in categories.values())
    logger.info(
        "Categorised %d/%d phenotypes into %d categories",
        n_categorized,
        len(phenotypes),
        len(categories),
    )

    return {
        "categories": categories,
        "uncategorized": uncategorized,
        "method": method,
        "n_categorized": n_categorized,
    }
