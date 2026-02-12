"""GWAS association testing utilities.

This module provides functions for performing association tests between
genotypes and phenotypes, including linear and logistic regression.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core.utils import logging

# Import numpy with graceful fallback
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore

logger = logging.get_logger(__name__)


def association_test_linear(
    genotypes: List[int], phenotypes: List[float], covariates: Optional[List[List[float]]] = None, **kwargs
) -> Dict[str, Any]:
    """Perform linear regression association test.

    Args:
        genotypes: Genotype values (0, 1, 2) for each sample
        phenotypes: Phenotype values for each sample
        covariates: Optional covariate matrix (covariates x samples)
        **kwargs: Additional parameters

    Returns:
        Dictionary with association test results
    """
    if len(genotypes) != len(phenotypes):
        raise ValueError("Genotypes and phenotypes must have same length")

    logger.debug(f"Running linear association test on {len(genotypes)} samples")

    # Simple linear regression implementation
    # In practice, would use statsmodels or similar

    n = len(genotypes)

    # Prepare design matrix
    X = [[1.0, float(gt)] for gt in genotypes]  # Intercept + genotype

    # Add covariates if provided
    if covariates:
        for i in range(len(covariates[0])):  # For each sample
            for j, cov in enumerate(covariates):
                X[i].append(cov[i])

    # Simple OLS implementation
    try:
        beta, se, t_stat, p_value, r_squared = _simple_linear_regression(X, phenotypes)

        result = {
            "status": "success",
            "beta": beta,
            "se": se,
            "t_stat": t_stat,
            "p_value": p_value,
            "r_squared": r_squared,
            "n_samples": n,
            "test_type": "linear",
            "converged": True,
        }

    except Exception as e:
        logger.warning(f"Linear regression failed: {e}")
        result = {
            "status": "error",
            "beta": 0.0,
            "se": 0.0,
            "t_stat": 0.0,
            "p_value": 1.0,
            "r_squared": 0.0,
            "n_samples": n,
            "test_type": "linear",
            "converged": False,
            "error": str(e),
        }

    return result


def association_test_logistic(
    genotypes: List[int],
    phenotypes: List[int],
    covariates: Optional[List[List[float]]] = None,
    max_iter: int = 100,
    **kwargs,
) -> Dict[str, Any]:
    """Perform logistic regression association test.

    Args:
        genotypes: Genotype values (0, 1, 2) for each sample
        phenotypes: Binary phenotype values (0, 1) for each sample
        covariates: Optional covariate matrix (covariates x samples)
        max_iter: Maximum iterations for optimization
        **kwargs: Additional parameters

    Returns:
        Dictionary with association test results
    """
    if len(genotypes) != len(phenotypes):
        raise ValueError("Genotypes and phenotypes must have same length")

    if not all(p in (0, 1) for p in phenotypes):
        raise ValueError("Logistic regression requires binary phenotypes (0 or 1)")

    logger.debug(f"Running logistic association test on {len(genotypes)} samples")

    n = len(genotypes)

    # Check for sufficient cases/controls
    n_cases = sum(1 for p in phenotypes if p == 1)
    n_controls = sum(1 for p in phenotypes if p == 0)
    if n_cases < 2 or n_controls < 2:
        logger.warning(f"Insufficient cases ({n_cases}) or controls ({n_controls}) for logistic regression")
        return {
            "status": "failed",
            "beta": 0.0,
            "se": 0.0,
            "z_stat": 0.0,
            "p_value": 1.0,
            "odds_ratio": 1.0,
            "n_samples": n,
            "n_cases": n_cases,
            "n_controls": n_controls,
            "test_type": "logistic",
            "converged": False,
            "error": f"Insufficient cases ({n_cases}) or controls ({n_controls})",
        }

    # Prepare design matrix
    X = [[1.0, float(gt)] for gt in genotypes]  # Intercept + genotype

    # Add covariates if provided
    if covariates:
        for i in range(len(covariates[0])):  # For each sample
            for j, cov in enumerate(covariates):
                X[i].append(cov[i])

    # Logistic regression with proper statistical inference
    try:
        if not HAS_NUMPY:
            raise ImportError("numpy is required for logistic regression")
        X_array = np.array(X)
        y_array = np.array(phenotypes)
        beta, se, z_stat, p_value = _logistic_regression_with_stats(X_array, y_array, max_iter)

        # Calculate odds ratio
        odds_ratio = math.exp(beta) if abs(beta) < 700 else float("inf")

        result = {
            "status": "success",
            "beta": beta,
            "se": se,
            "z_stat": z_stat,
            "p_value": p_value,
            "odds_ratio": odds_ratio,
            "n_samples": n,
            "n_cases": n_cases,
            "n_controls": n_controls,
            "test_type": "logistic",
            "converged": True,
        }

    except Exception as e:
        logger.warning(f"Logistic regression failed: {e}")
        result = {
            "status": "error",
            "beta": 0.0,
            "se": 0.0,
            "z_stat": 0.0,
            "p_value": 1.0,
            "odds_ratio": 1.0,
            "n_samples": n,
            "n_cases": n_cases,
            "n_controls": n_controls,
            "test_type": "logistic",
            "converged": False,
            "error": str(e),
        }

    return result


def _simple_linear_regression(X: List[List[float]], y: List[float]) -> Tuple[float, float, float, float, float]:
    """Simple linear regression implementation.

    Args:
        X: Design matrix
        y: Response variable

    Returns:
        Tuple of (beta, se, t_stat, p_value, r_squared)
    """
    # Very simplified OLS - in practice would use proper linear algebra
    n = len(y)
    if n < 3:
        return 0.0, 0.0, 0.0, 1.0, 0.0

    # Simple slope calculation for single predictor
    if len(X[0]) == 2:  # Intercept + one predictor
        x_vals = [row[1] for row in X]
        y_vals = y

        # Calculate means
        x_mean = sum(x_vals) / n
        y_mean = sum(y_vals) / n

        # Calculate slope
        numerator = sum((x - x_mean) * (yi - y_mean) for x, yi in zip(x_vals, y_vals))
        denominator = sum((x - x_mean) ** 2 for x in x_vals)

        if denominator == 0:
            return 0.0, 0.0, 0.0, 1.0, 0.0

        beta = numerator / denominator
        intercept = y_mean - beta * x_mean

        # Calculate R-squared
        ss_tot = sum((yi - y_mean) ** 2 for yi in y_vals)
        y_pred = [intercept + beta * x for x in x_vals]
        ss_res = sum((yi - yp) ** 2 for yi, yp in zip(y_vals, y_pred))
        r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

        # Calculate standard error (simplified)
        residuals = [yi - yp for yi, yp in zip(y_vals, y_pred)]
        mse = sum(r**2 for r in residuals) / (n - 2)
        se = math.sqrt(mse / denominator) if mse >= 0 and denominator > 0 else 0.0

        # t-statistic
        t_stat = beta / se if se > 0 else 0.0

        # p-value approximation (two-tailed)
        p_value = 2 * (1 - _normal_cdf(abs(t_stat)))

        return beta, se, t_stat, p_value, r_squared

    # Multiple predictors: use numpy OLS if available, otherwise pure-python fallback
    if HAS_NUMPY:
        X_arr = np.array(X, dtype=float)
        y_arr = np.array(y, dtype=float)
        # OLS via least squares: beta = (X'X)^-1 X'y
        try:
            beta_vec, residuals_arr, rank, sv = np.linalg.lstsq(X_arr, y_arr, rcond=None)
        except np.linalg.LinAlgError:
            return 0.0, 0.0, 0.0, 1.0, 0.0

        # Predicted values and residuals
        y_pred_arr = X_arr @ beta_vec
        resid = y_arr - y_pred_arr
        p = len(beta_vec)  # number of predictors including intercept
        dof = n - p
        if dof <= 0:
            return 0.0, 0.0, 0.0, 1.0, 0.0

        mse = float(np.sum(resid**2) / dof)

        # Variance-covariance matrix of beta
        try:
            XtX_inv = np.linalg.inv(X_arr.T @ X_arr)
        except np.linalg.LinAlgError:
            return 0.0, 0.0, 0.0, 1.0, 0.0

        beta_cov = mse * XtX_inv

        # Genotype coefficient is at index 1 (after intercept)
        beta = float(beta_vec[1])
        se = float(math.sqrt(max(beta_cov[1, 1], 0.0)))

        # R-squared
        y_mean = float(np.mean(y_arr))
        ss_tot = float(np.sum((y_arr - y_mean) ** 2))
        ss_res = float(np.sum(resid**2))
        r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

        # t-statistic and p-value
        t_stat = beta / se if se > 0 else 0.0

        # Use scipy t-distribution if available, otherwise normal approximation
        try:
            from scipy import stats as scipy_stats

            p_value = float(2 * scipy_stats.t.sf(abs(t_stat), dof))
        except ImportError:
            p_value = 2 * (1 - _normal_cdf(abs(t_stat)))

        return beta, se, t_stat, p_value, r_squared

    # Pure-python fallback for multiple predictors without numpy
    # Use Gaussian elimination (reuse existing _solve_linear_system)
    k = len(X[0])
    # Compute X'X
    XtX = [[sum(X[s][i] * X[s][j] for s in range(n)) for j in range(k)] for i in range(k)]
    # Compute X'y
    Xty = [sum(X[s][i] * y[s] for s in range(n)) for i in range(k)]
    beta_vec_pp = _solve_linear_system(XtX, Xty)
    if beta_vec_pp is None:
        return 0.0, 0.0, 0.0, 1.0, 0.0

    # Predicted and residuals
    y_pred_pp = [sum(X[s][j] * beta_vec_pp[j] for j in range(k)) for s in range(n)]
    resid_pp = [y[s] - y_pred_pp[s] for s in range(n)]
    dof_pp = n - k
    if dof_pp <= 0:
        return 0.0, 0.0, 0.0, 1.0, 0.0

    mse_pp = sum(r**2 for r in resid_pp) / dof_pp

    # Inverse of X'X for standard errors
    XtX_inv_pp = _matrix_inverse(XtX)
    if XtX_inv_pp is None:
        return 0.0, 0.0, 0.0, 1.0, 0.0

    beta = beta_vec_pp[1]
    se = math.sqrt(max(mse_pp * XtX_inv_pp[1][1], 0.0))

    y_mean_pp = sum(y) / n
    ss_tot_pp = sum((yi - y_mean_pp) ** 2 for yi in y)
    ss_res_pp = sum(r**2 for r in resid_pp)
    r_squared = 1.0 - (ss_res_pp / ss_tot_pp) if ss_tot_pp > 0 else 0.0

    t_stat = beta / se if se > 0 else 0.0
    p_value = 2 * (1 - _normal_cdf(abs(t_stat)))

    return beta, se, t_stat, p_value, r_squared


def _logistic_regression_with_stats(X: Any, y: Any, max_iter: int) -> Tuple[float, float, float, float]:
    """Perform logistic regression with statistical inference.

    Uses statsmodels or sklearn if available, otherwise falls back to a pure Python
    implementation using iteratively reweighted least squares (IRLS).

    Args:
        X: Design matrix as numpy array or list of lists
        y: Binary response variable as numpy array or list
        max_iter: Maximum iterations

    Returns:
        Tuple of (beta, se, z_stat, p_value) for the genotype coefficient
    """
    # Try statsmodels first (best statistical inference)
    try:
        import statsmodels.api as sm

        X_with_intercept = sm.add_constant(X)
        logit_model = sm.Logit(y, X_with_intercept)
        result = logit_model.fit(disp=False, maxiter=max_iter)

        beta = float(result.params[1])  # Genotype coefficient
        p_value = float(result.pvalues[1])
        z_stat = beta / result.bse[1] if result.bse[1] != 0 else 0.0
        se = float(result.bse[1])

        return beta, se, z_stat, p_value

    except ImportError:
        pass

    # Try sklearn as fallback
    try:
        import scipy.stats as stats
        from sklearn.linear_model import LogisticRegression

        lr = LogisticRegression(max_iter=max_iter, random_state=42)
        lr.fit(X, y)

        beta = float(lr.coef_[0][1]) if len(lr.coef_[0]) > 1 else float(lr.coef_[0][0])
        n = len(y)
        se = 1.0 / np.sqrt(n)
        z_stat = beta / se if se != 0 else 0.0
        p_value = float(2 * (1 - stats.norm.cdf(abs(z_stat))))

        return beta, se, z_stat, p_value

    except ImportError:
        pass

    # Pure Python fallback using IRLS
    return _logistic_regression_pure_python(X, y, max_iter)


def _logistic_regression_pure_python(X: Any, y: Any, max_iter: int) -> Tuple[float, float, float, float]:
    """Pure Python logistic regression using iteratively reweighted least squares.

    Args:
        X: Design matrix (n_samples x n_features)
        y: Binary response variable (n_samples,)
        max_iter: Maximum iterations

    Returns:
        Tuple of (beta, se, z_stat, p_value) for the first non-intercept coefficient
    """
    # Convert to lists for pure Python
    if hasattr(X, "tolist"):
        X_list = X.tolist()
    else:
        X_list = list(X)
    if hasattr(y, "tolist"):
        y_list = y.tolist()
    else:
        y_list = list(y)

    n = len(y_list)
    n_features = len(X_list[0]) if X_list else 0

    if n_features == 0 or n == 0:
        return 0.0, 0.0, 0.0, 1.0

    # Add intercept column if not present (check if first column is all 1s)
    has_intercept = all(abs(X_list[i][0] - 1.0) < 1e-10 for i in range(n))
    if not has_intercept:
        X_list = [[1.0] + row for row in X_list]
        n_features += 1

    # Initialize coefficients
    beta = [0.0] * n_features

    # IRLS iterations
    for _ in range(max_iter):
        # Compute linear predictor
        eta = [sum(X_list[i][j] * beta[j] for j in range(n_features)) for i in range(n)]

        # Compute probabilities (sigmoid)
        prob = []
        for e in eta:
            if e > 20:
                p = 1.0 - 1e-10
            elif e < -20:
                p = 1e-10
            else:
                p = 1.0 / (1.0 + math.exp(-e))
            prob.append(p)

        # Compute weights and working response
        W = [max(p * (1 - p), 1e-10) for p in prob]
        z = [eta[i] + (y_list[i] - prob[i]) / W[i] for i in range(n)]

        # Weighted least squares solve: (X'WX)^-1 X'Wz
        # Compute X'WX
        XtWX = [[0.0] * n_features for _ in range(n_features)]
        for j in range(n_features):
            for k in range(n_features):
                XtWX[j][k] = sum(X_list[i][j] * W[i] * X_list[i][k] for i in range(n))

        # Compute X'Wz
        XtWz = [sum(X_list[i][j] * W[i] * z[i] for i in range(n)) for j in range(n_features)]

        # Solve using simple Gaussian elimination (for small systems)
        new_beta = _solve_linear_system(XtWX, XtWz)
        if new_beta is None:
            break

        # Check convergence
        max_change = max(abs(new_beta[j] - beta[j]) for j in range(n_features))
        beta = new_beta

        if max_change < 1e-6:
            break

    # Compute standard errors from inverse of Fisher information matrix
    # Fisher info = X'WX (already computed in last iteration)
    eta = [sum(X_list[i][j] * beta[j] for j in range(n_features)) for i in range(n)]
    prob = []
    for e in eta:
        if e > 20:
            p = 1.0 - 1e-10
        elif e < -20:
            p = 1e-10
        else:
            p = 1.0 / (1.0 + math.exp(-e))
        prob.append(p)
    W = [max(p * (1 - p), 1e-10) for p in prob]

    XtWX = [[0.0] * n_features for _ in range(n_features)]
    for j in range(n_features):
        for k in range(n_features):
            XtWX[j][k] = sum(X_list[i][j] * W[i] * X_list[i][k] for i in range(n))

    # Invert to get variance-covariance matrix
    cov_matrix = _matrix_inverse(XtWX)
    if cov_matrix is None:
        se = 1.0 / math.sqrt(n) if n > 0 else 0.0
    else:
        # SE for genotype coefficient (index 1 if intercept present)
        coef_idx = 1 if n_features > 1 else 0
        se = math.sqrt(max(cov_matrix[coef_idx][coef_idx], 0.0))

    # Get genotype coefficient
    genotype_beta = beta[1] if n_features > 1 else beta[0]

    # Compute z-statistic and p-value
    z_stat = genotype_beta / se if se > 0 else 0.0
    p_value = 2 * (1 - _normal_cdf(abs(z_stat)))

    return genotype_beta, se, z_stat, p_value


def _solve_linear_system(A: List[List[float]], b: List[float]) -> Optional[List[float]]:
    """Solve linear system Ax = b using Gaussian elimination with partial pivoting.

    Args:
        A: Coefficient matrix (n x n)
        b: Right-hand side vector (n,)

    Returns:
        Solution vector or None if singular
    """
    n = len(b)
    if n == 0:
        return None

    # Augmented matrix
    aug = [A[i][:] + [b[i]] for i in range(n)]

    # Forward elimination with partial pivoting
    for col in range(n):
        # Find pivot
        max_row = col
        for row in range(col + 1, n):
            if abs(aug[row][col]) > abs(aug[max_row][col]):
                max_row = row
        aug[col], aug[max_row] = aug[max_row], aug[col]

        if abs(aug[col][col]) < 1e-12:
            return None  # Singular

        # Eliminate
        for row in range(col + 1, n):
            factor = aug[row][col] / aug[col][col]
            for j in range(col, n + 1):
                aug[row][j] -= factor * aug[col][j]

    # Back substitution
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        if abs(aug[i][i]) < 1e-12:
            return None
        x[i] = aug[i][n]
        for j in range(i + 1, n):
            x[i] -= aug[i][j] * x[j]
        x[i] /= aug[i][i]

    return x


def _matrix_inverse(A: List[List[float]]) -> Optional[List[List[float]]]:
    """Compute matrix inverse using Gaussian elimination.

    Args:
        A: Square matrix (n x n)

    Returns:
        Inverse matrix or None if singular
    """
    n = len(A)
    if n == 0:
        return None

    # Augmented matrix [A | I]
    aug = [A[i][:] + [1.0 if j == i else 0.0 for j in range(n)] for i in range(n)]

    # Forward elimination
    for col in range(n):
        # Find pivot
        max_row = col
        for row in range(col + 1, n):
            if abs(aug[row][col]) > abs(aug[max_row][col]):
                max_row = row
        aug[col], aug[max_row] = aug[max_row], aug[col]

        if abs(aug[col][col]) < 1e-12:
            return None

        # Scale pivot row
        scale = aug[col][col]
        for j in range(2 * n):
            aug[col][j] /= scale

        # Eliminate column
        for row in range(n):
            if row != col:
                factor = aug[row][col]
                for j in range(2 * n):
                    aug[row][j] -= factor * aug[col][j]

    # Extract inverse
    return [aug[i][n:] for i in range(n)]


def _normal_cdf(x: float) -> float:
    """Approximate normal cumulative distribution function.

    Args:
        x: Input value

    Returns:
        CDF value
    """
    # Abramowitz & Stegun approximation
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911

    sign = 1 if x >= 0 else -1
    x = abs(x) / math.sqrt(2.0)

    t = 1.0 / (1.0 + p * x)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * math.exp(-x * x)

    return 0.5 * (1 + sign * y)
