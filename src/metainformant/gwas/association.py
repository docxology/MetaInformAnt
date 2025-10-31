"""Association testing for GWAS."""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any

import numpy as np

from ..core.io import ensure_directory, read_delimited, write_delimited
from .quality import parse_vcf_full

logger = logging.getLogger(__name__)

# Try importing scipy for statistical tests
try:
    from scipy.stats import norm

    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    norm = None


def association_test_linear(
    genotypes: list[int],
    phenotypes: list[float],
    covariates: list[list[float]] | None = None,
) -> dict[str, Any]:
    """Test association using linear regression.

    Model: phenotype ~ genotype + covariates

    Args:
        genotypes: Genotype vector for one variant (0/1/2/-1 encoding)
        phenotypes: Phenotype values
        covariates: Optional covariates matrix (list of lists, one per covariate)

    Returns:
        Dictionary with:
        - beta: Effect size (regression coefficient)
        - se: Standard error
        - p_value: P-value
        - r_squared: R-squared for the model
    """
    if len(genotypes) != len(phenotypes):
        return {
            "status": "failed",
            "error": "Genotype and phenotype lengths do not match",
        }

    # Remove missing data
    valid_indices = [i for i in range(len(genotypes)) if genotypes[i] != -1 and not (isinstance(phenotypes[i], float) and (np.isnan(phenotypes[i]) or np.isinf(phenotypes[i])))]
    
    if len(valid_indices) < 3:
        return {
            "status": "failed",
            "error": "Insufficient non-missing data",
        }

    genos = np.array([genotypes[i] for i in valid_indices], dtype=float)
    phenos = np.array([phenotypes[i] for i in valid_indices], dtype=float)

    n = len(genos)

    # Build design matrix
    X = np.ones((n, 1))  # Intercept
    X = np.column_stack([X, genos])  # Add genotype

    # Add covariates if provided
    if covariates and len(covariates) > 0:
        num_covariates = len(covariates[0]) if isinstance(covariates[0], (list, tuple)) else 0
        for cov_idx in range(num_covariates):
            cov_vector = [covariates[i][cov_idx] for i in valid_indices if i < len(covariates)]
            if len(cov_vector) == n:
                X = np.column_stack([X, cov_vector])

    try:
        # Fit linear regression using normal equation
        # beta = (X'X)^(-1) X'y
        XtX = np.dot(X.T, X)
        Xty = np.dot(X.T, phenos)

        # Add small regularization for numerical stability
        reg = 1e-8 * np.eye(XtX.shape[0])
        beta = np.linalg.solve(XtX + reg, Xty)

        # Predictions
        y_pred = np.dot(X, beta)
        residuals = phenos - y_pred

        # Residual sum of squares
        rss = np.sum(residuals ** 2)
        df_residual = n - X.shape[1]

        # Standard errors
        if df_residual > 0:
            mse = rss / df_residual
            var_beta = mse * np.linalg.inv(XtX + reg)
            se = np.sqrt(np.diag(var_beta))
        else:
            se = np.zeros(len(beta))

        # T-test for genotype coefficient (first non-intercept coefficient)
        beta_geno = beta[1]  # Genotype coefficient
        se_geno = se[1] if len(se) > 1 else 0.0

        if se_geno > 0:
            t_stat = beta_geno / se_geno
            # Two-tailed p-value
            if SCIPY_AVAILABLE and norm:
                p_value = 2.0 * (1.0 - norm.cdf(abs(t_stat)))
            else:
                # Simple approximation: p ≈ 2*(1 - erf(|t|/√2))
                # For large |t|, use exponential approximation
                if abs(t_stat) > 6:
                    p_value = 2.0 * math.exp(-t_stat * t_stat / 2.0) / (abs(t_stat) * math.sqrt(2 * math.pi))
                else:
                    # Rough approximation
                    p_value = 2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(t_stat) / math.sqrt(2.0))))
        else:
            p_value = 1.0

        # R-squared
        ss_res = rss
        ss_tot = np.sum((phenos - np.mean(phenos)) ** 2)
        r_squared = 1.0 - (ss_res / (ss_tot + 1e-10))

        return {
            "status": "success",
            "beta": float(beta_geno),
            "se": float(se_geno),
            "p_value": float(p_value),
            "r_squared": float(r_squared),
            "n": n,
        }

    except Exception as exc:
        logger.debug(f"test_association_linear: Regression failed: {exc}")
        return {
            "status": "failed",
            "error": str(exc),
        }


def association_test_logistic(
    genotypes: list[int],
    phenotypes: list[int],
    covariates: list[list[float]] | None = None,
    max_iter: int = 100,
) -> dict[str, Any]:
    """Test association using logistic regression for binary traits.

    Model: logit(P(Y=1)) = beta0 + beta_geno * genotype + covariates

    Args:
        genotypes: Genotype vector for one variant (0/1/2/-1 encoding)
        phenotypes: Binary phenotype values (0/1)
        covariates: Optional covariates matrix
        max_iter: Maximum iterations for logistic regression

    Returns:
        Dictionary with:
        - beta: Log odds ratio (coefficient)
        - se: Standard error
        - p_value: P-value
        - odds_ratio: Odds ratio (exp(beta))
    """
    if len(genotypes) != len(phenotypes):
        return {
            "status": "failed",
            "error": "Genotype and phenotype lengths do not match",
        }

    # Remove missing data
    valid_indices = [i for i in range(len(genotypes)) if genotypes[i] != -1 and phenotypes[i] in [0, 1]]
    
    if len(valid_indices) < 3:
        return {
            "status": "failed",
            "error": "Insufficient non-missing data",
        }

    genos = np.array([genotypes[i] for i in valid_indices], dtype=float)
    phenos = np.array([phenotypes[i] for i in valid_indices], dtype=int)

    n = len(genos)

    # Check for sufficient cases and controls
    if np.sum(phenos) < 2 or np.sum(1 - phenos) < 2:
        return {
            "status": "failed",
            "error": "Insufficient cases or controls",
        }

    # Build design matrix
    X = np.ones((n, 1))  # Intercept
    X = np.column_stack([X, genos])  # Add genotype

    # Add covariates if provided
    if covariates and len(covariates) > 0:
        num_covariates = len(covariates[0]) if isinstance(covariates[0], (list, tuple)) else 0
        for cov_idx in range(num_covariates):
            cov_vector = [covariates[i][cov_idx] for i in valid_indices if i < len(covariates)]
            if len(cov_vector) == n:
                X = np.column_stack([X, cov_vector])

    try:
        # Simple logistic regression using Newton-Raphson
        beta = np.zeros(X.shape[1])

        for iteration in range(max_iter):
            # Compute predictions
            logits = np.dot(X, beta)
            # Clip logits to prevent overflow
            logits = np.clip(logits, -500, 500)
            probs = 1.0 / (1.0 + np.exp(-logits))

            # Gradient
            gradient = np.dot(X.T, (phenos - probs))

            # Hessian (negative for Newton-Raphson in minimization context)
            W = probs * (1.0 - probs)
            hessian = -np.dot(X.T, W[:, np.newaxis] * X)

            # Update
            try:
                beta_update = np.linalg.solve(-hessian + 1e-8 * np.eye(hessian.shape[0]), gradient)
                # Limit step size to prevent divergence
                max_step = 10.0
                beta_update = np.clip(beta_update, -max_step, max_step)
                beta = beta - beta_update
            except np.linalg.LinAlgError:
                break

            # Check convergence
            if np.max(np.abs(beta_update)) < 1e-6:
                break

        # Standard errors from Hessian
        logits = np.dot(X, beta)
        logits = np.clip(logits, -500, 500)
        probs = 1.0 / (1.0 + np.exp(-logits))
        W = probs * (1.0 - probs)
        hessian = -np.dot(X.T, W[:, np.newaxis] * X)

        try:
            var_beta = np.linalg.inv(-hessian + 1e-8 * np.eye(hessian.shape[0]))
            se = np.sqrt(np.diag(var_beta))
        except np.linalg.LinAlgError:
            se = np.ones(len(beta)) * np.inf

        # Genotype coefficient
        beta_geno = beta[1] if len(beta) > 1 else 0.0
        se_geno = se[1] if len(se) > 1 else np.inf

        # Wald test
        if se_geno > 0 and not np.isinf(se_geno):
            z_stat = beta_geno / se_geno
            if SCIPY_AVAILABLE and norm:
                p_value = 2.0 * (1.0 - norm.cdf(abs(z_stat)))
            else:
                # Approximation
                if abs(z_stat) > 6:
                    p_value = 2.0 * math.exp(-z_stat * z_stat / 2.0) / (abs(z_stat) * math.sqrt(2 * math.pi))
                else:
                    p_value = 2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(z_stat) / math.sqrt(2.0))))
        else:
            p_value = 1.0

        # Calculate odds ratio, handling extreme values
        if not np.isnan(beta_geno) and not np.isinf(beta_geno):
            # Clip beta to prevent overflow in exp
            beta_clipped = np.clip(beta_geno, -50, 50)
            odds_ratio = math.exp(beta_clipped)
        else:
            odds_ratio = 1.0

        return {
            "status": "success",
            "beta": float(beta_geno),
            "se": float(se_geno),
            "p_value": float(p_value),
            "odds_ratio": float(odds_ratio),
            "n": n,
        }

    except Exception as exc:
        logger.debug(f"test_association_logistic: Regression failed: {exc}")
        return {
            "status": "failed",
            "error": str(exc),
        }


def run_gwas(
    vcf_path: str | Path,
    phenotype_path: str | Path,
    config: dict[str, Any],
    output_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Run association tests for all variants in VCF.

    Args:
        vcf_path: Path to VCF file
        phenotype_path: Path to phenotype file (TSV with sample_id, trait columns)
        config: Association configuration
        output_dir: Optional directory for results

    Returns:
        Dictionary with association test results
    """
    logger.info(f"run_gwas: Running association tests for {vcf_path}")

    model = config.get("model", "linear")
    trait = config.get("trait")
    covariate_names = config.get("covariates", [])
    min_sample_size = config.get("min_sample_size", 50)

    if not trait:
        return {
            "status": "failed",
            "error": "Trait name not specified in config",
        }

    # Parse VCF
    vcf_data = parse_vcf_full(vcf_path)
    samples = vcf_data["samples"]
    variants = vcf_data["variants"]
    genotypes = vcf_data["genotypes"]

    if not variants:
        return {
            "status": "failed",
            "error": "No variants found in VCF",
        }

    # Load phenotypes
    phenotype_path_obj = Path(phenotype_path)
    if not phenotype_path_obj.exists():
        return {
            "status": "failed",
            "error": f"Phenotype file not found: {phenotype_path}",
        }
    
    try:
        phenotype_data_list = list(read_delimited(phenotype_path, delimiter="\t"))
    except (FileNotFoundError, IOError) as e:
        return {
            "status": "failed",
            "error": f"Cannot read phenotype file: {e}",
        }
    if not phenotype_data_list or len(phenotype_data_list) < 1:
        return {
            "status": "failed",
            "error": "Invalid phenotype file",
        }

    # Parse phenotype header (first row has keys)
    if not phenotype_data_list:
        return {
            "status": "failed",
            "error": "Empty phenotype file",
        }
    
    header = list(phenotype_data_list[0].keys())
    
    if trait not in header:
        return {
            "status": "failed",
            "error": f"Trait '{trait}' not found in phenotype file",
        }

    # Extract phenotype values
    phenotype_dict: dict[str, float] = {}
    for row in phenotype_data_list:
        sample_id = row.get("sample_id", "")
        if sample_id:
            try:
                trait_value = float(row.get(trait, float("nan")))
                if not (isinstance(trait_value, float) and (np.isnan(trait_value) or np.isinf(trait_value))):
                    phenotype_dict[sample_id] = trait_value
            except (ValueError, TypeError):
                continue

    # Match samples to phenotypes
    phenotypes = [phenotype_dict.get(sample_id, float("nan")) for sample_id in samples]

    # Load covariates if specified
    covariates_list: list[list[float]] | None = None
    if covariate_names:
        covariate_dicts: dict[str, dict[str, float]] = {}
        for cov_name in covariate_names:
            if cov_name in header:
                cov_dict = {}
                for row in phenotype_data_list:
                    sample_id = row.get("sample_id", "")
                    if sample_id:
                        try:
                            cov_value = float(row.get(cov_name, 0.0))
                            cov_dict[sample_id] = cov_value
                        except (ValueError, TypeError):
                            continue
                covariate_dicts[cov_name] = cov_dict

        # Build covariates matrix
        if covariate_dicts:
            covariates_list = []
            for sample_id in samples:
                cov_row = [covariate_dicts[cov_name].get(sample_id, 0.0) for cov_name in covariate_names]
                covariates_list.append(cov_row)

    # Run association tests
    results: list[dict[str, Any]] = []
    num_variants = len(variants)

    logger.info(f"run_gwas: Testing {num_variants} variants")

    for var_idx in range(num_variants):
        variant = variants[var_idx]
        var_genotypes = [genotypes[sample_idx][var_idx] for sample_idx in range(len(samples))]

        # Skip if insufficient sample size
        valid_samples = sum(1 for i, g in enumerate(var_genotypes) if g != -1 and not (isinstance(phenotypes[i], float) and (np.isnan(phenotypes[i]) or np.isinf(phenotypes[i]))))
        if valid_samples < min_sample_size:
            continue

        # Test association
        if model == "logistic":
            # Convert phenotypes to binary (0/1)
            binary_phenos = [1 if p > 0.5 else 0 for p in phenotypes]
            test_result = association_test_logistic(var_genotypes, binary_phenos, covariates_list)
        else:
            test_result = association_test_linear(var_genotypes, phenotypes, covariates_list)

        if test_result.get("status") == "success":
            result_row = {
                "CHROM": variant["CHROM"],
                "POS": variant["POS"],
                "ID": variant["ID"],
                "REF": variant["REF"],
                "ALT": variant["ALT"],
                "beta": test_result["beta"],
                "se": test_result["se"],
                "p_value": test_result["p_value"],
                "n": test_result.get("n", valid_samples),
            }
            if "r_squared" in test_result:
                result_row["r_squared"] = test_result["r_squared"]
            if "odds_ratio" in test_result:
                result_row["odds_ratio"] = test_result["odds_ratio"]

            results.append(result_row)

        if (var_idx + 1) % 1000 == 0:
            logger.info(f"run_gwas: Processed {var_idx + 1}/{num_variants} variants")

    logger.info(f"run_gwas: Completed {len(results)} association tests")

    # Prepare output
    output_data: dict[str, Any] = {
        "status": "success",
        "num_variants_tested": len(results),
        "model": model,
        "trait": trait,
        "results": results,
    }

    # Write results if output directory provided
    if output_dir:
        out_dir = ensure_directory(output_dir)

        # Write results table
        if results:
            results_table = [["CHROM", "POS", "ID", "REF", "ALT", "beta", "se", "p_value", "n"]]
            for result in results:
                row = [
                    str(result.get("CHROM", "")),
                    str(result.get("POS", "")),
                    str(result.get("ID", "")),
                    str(result.get("REF", "")),
                    str(result.get("ALT", "")),
                    f"{result.get('beta', 0.0):.6e}",
                    f"{result.get('se', 0.0):.6e}",
                    f"{result.get('p_value', 1.0):.6e}",
                    str(result.get("n", 0)),
                ]
                if "r_squared" in result:
                    results_table[0].append("r_squared")
                    row.append(f"{result['r_squared']:.6f}")
                if "odds_ratio" in result:
                    results_table[0].append("odds_ratio")
                    row.append(f"{result['odds_ratio']:.6f}")

                results_table.append(row)

            results_path = out_dir / "association_results.tsv"
            from ..core.io import write_tsv
            write_tsv(results_table, results_path)
            logger.info(f"run_gwas: Wrote results to {results_path}")

    return output_data

