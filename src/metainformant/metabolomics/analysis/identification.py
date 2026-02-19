"""Metabolite identification and quantification.

Provides methods for matching observed m/z values to known metabolite databases,
normalizing intensity data, and performing differential abundance analysis.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class MetaboliteMatch:
    """Result of a metabolite identification match.

    Attributes:
        query_mz: Queried mass-to-charge ratio.
        matched_name: Name of matched metabolite.
        matched_mz: Database m/z value.
        delta_ppm: Mass error in parts per million.
        score: Match confidence score (0-1).
    """

    query_mz: float
    matched_name: str
    matched_mz: float
    delta_ppm: float
    score: float


def identify_metabolites(
    observed_mz: np.ndarray,
    database: dict[str, float],
    ppm_tolerance: float = 10.0,
) -> list[list[MetaboliteMatch]]:
    """Identify metabolites by matching observed m/z to a reference database.

    For each observed m/z, finds all database entries within the specified
    mass tolerance (in ppm).

    Args:
        observed_mz: 1D array of observed m/z values.
        database: Dict mapping metabolite names to their exact m/z values.
        ppm_tolerance: Maximum mass error in ppm for a match.

    Returns:
        List of lists; each inner list contains MetaboliteMatch entries
        for the corresponding observed m/z.
    """
    results: list[list[MetaboliteMatch]] = []
    db_names = list(database.keys())
    db_mz = np.array(list(database.values()))

    for mz in observed_mz:
        ppm_errors = np.abs(db_mz - mz) / mz * 1e6
        matches = []
        for i in np.where(ppm_errors <= ppm_tolerance)[0]:
            score = max(0.0, 1.0 - ppm_errors[i] / ppm_tolerance)
            matches.append(
                MetaboliteMatch(
                    query_mz=float(mz),
                    matched_name=db_names[i],
                    matched_mz=float(db_mz[i]),
                    delta_ppm=float(ppm_errors[i]),
                    score=float(score),
                )
            )
        matches.sort(key=lambda m: m.score, reverse=True)
        results.append(matches)

    return results


def normalize_intensities(
    intensities: np.ndarray,
    method: str = "total_ion_count",
) -> np.ndarray:
    """Normalize metabolite intensity matrix.

    Args:
        intensities: 2D array (metabolites × samples).
        method: Normalization method. One of:
            - 'total_ion_count': Scale by total intensity per sample.
            - 'median': Scale by median intensity per sample.
            - 'log2': Log2(x + 1) transformation.
            - 'pareto': Pareto scaling (divide by sqrt of std dev).

    Returns:
        Normalized intensity matrix.
    """
    data = intensities.copy().astype(float)

    if method == "total_ion_count":
        totals = data.sum(axis=0, keepdims=True)
        totals = np.where(totals > 0, totals, 1.0)
        return data / totals * np.median(totals)

    elif method == "median":
        medians = np.median(data, axis=0, keepdims=True)
        medians = np.where(medians > 0, medians, 1.0)
        return data / medians * np.median(medians)

    elif method == "log2":
        return np.log2(data + 1)

    elif method == "pareto":
        mean = data.mean(axis=1, keepdims=True)
        std = data.std(axis=1, keepdims=True)
        std = np.where(std > 0, std, 1.0)
        return (data - mean) / np.sqrt(std)

    else:
        raise ValueError(f"Unknown normalization method: {method}")


def fold_change(
    intensities: np.ndarray,
    group_a: list[int],
    group_b: list[int],
) -> np.ndarray:
    """Compute log2 fold change between two groups.

    Args:
        intensities: 2D array (metabolites × samples).
        group_a: Sample indices for group A.
        group_b: Sample indices for group B.

    Returns:
        1D array of log2 fold changes per metabolite.
    """
    mean_a = intensities[:, group_a].mean(axis=1)
    mean_b = intensities[:, group_b].mean(axis=1)

    safe_a = np.where(mean_a > 0, mean_a, 1e-10)
    safe_b = np.where(mean_b > 0, mean_b, 1e-10)

    return np.log2(safe_a / safe_b)


def differential_abundance(
    intensities: np.ndarray,
    group_a: list[int],
    group_b: list[int],
) -> tuple[np.ndarray, np.ndarray]:
    """Compute differential abundance statistics between two groups.

    Performs a two-sample t-test for each metabolite comparing group A vs B.

    Args:
        intensities: 2D array (metabolites × samples).
        group_a: Sample indices for group A.
        group_b: Sample indices for group B.

    Returns:
        Tuple of (t_statistics, p_values), each 1D array per metabolite.
    """
    data_a = intensities[:, group_a]
    data_b = intensities[:, group_b]

    n_a = data_a.shape[1]
    n_b = data_b.shape[1]

    mean_a = data_a.mean(axis=1)
    mean_b = data_b.mean(axis=1)

    var_a = data_a.var(axis=1, ddof=1) if n_a > 1 else np.zeros(intensities.shape[0])
    var_b = data_b.var(axis=1, ddof=1) if n_b > 1 else np.zeros(intensities.shape[0])

    pooled_se = np.sqrt(var_a / max(n_a, 1) + var_b / max(n_b, 1))
    pooled_se = np.where(pooled_se > 0, pooled_se, 1e-10)

    t_stats = (mean_a - mean_b) / pooled_se

    # Welch's degrees of freedom approximation
    num = (var_a / max(n_a, 1) + var_b / max(n_b, 1)) ** 2
    denom_a = (var_a / max(n_a, 1)) ** 2 / max(n_a - 1, 1) if n_a > 1 else np.zeros_like(var_a)
    denom_b = (var_b / max(n_b, 1)) ** 2 / max(n_b - 1, 1) if n_b > 1 else np.zeros_like(var_b)
    denom = denom_a + denom_b
    with np.errstate(divide="ignore", invalid="ignore"):
        df = np.where(denom > 0, num / denom, 1.0)

    # Approximate p-value using normal distribution for large df
    p_values = 2.0 * _normal_sf(np.abs(t_stats))

    return t_stats, p_values


def _normal_sf(x: np.ndarray) -> np.ndarray:
    """Survival function (1 - CDF) of standard normal, approximation."""
    # Abramowitz and Stegun approximation 7.1.26
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911

    sign = np.sign(x)
    x = np.abs(x) / np.sqrt(2.0)

    t = 1.0 / (1.0 + p * x)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * np.exp(-(x**2))

    return 0.5 * (1.0 - sign * y)


def cosine_spectral_similarity(
    spectrum_a: np.ndarray,
    spectrum_b: np.ndarray,
) -> float:
    """Compute cosine similarity between two mass spectra.

    Used for spectral library matching where each spectrum is represented
    as a vector of intensities at discrete m/z bins.

    Args:
        spectrum_a: 1D intensity vector for spectrum A.
        spectrum_b: 1D intensity vector for spectrum B.

    Returns:
        Cosine similarity score between 0 and 1.
    """
    dot = np.dot(spectrum_a, spectrum_b)
    norm_a = np.linalg.norm(spectrum_a)
    norm_b = np.linalg.norm(spectrum_b)

    if norm_a == 0 or norm_b == 0:
        return 0.0

    return float(dot / (norm_a * norm_b))


@dataclass
class AdductMatch:
    """Result of adduct-aware identification.

    Attributes:
        query_mz: Observed m/z.
        neutral_mass: Inferred neutral mass after adduct removal.
        adduct_type: Matched adduct type (e.g., '[M+H]+').
        matched_name: Matched metabolite name.
        delta_ppm: Mass error in ppm.
    """

    query_mz: float
    neutral_mass: float
    adduct_type: str
    matched_name: str
    delta_ppm: float


# Common ESI adducts: name -> mass offset from neutral mass
COMMON_ADDUCTS: dict[str, float] = {
    "[M+H]+": 1.007276,
    "[M+Na]+": 22.989218,
    "[M+K]+": 38.963158,
    "[M+NH4]+": 18.034164,
    "[M-H]-": -1.007276,
    "[M+Cl]-": 34.969402,
    "[M+FA-H]-": 44.998201,
}


def identify_with_adducts(
    observed_mz: np.ndarray,
    database: dict[str, float],
    adducts: dict[str, float] | None = None,
    ppm_tolerance: float = 10.0,
    ion_mode: str = "positive",
) -> list[list[AdductMatch]]:
    """Identify metabolites considering common adduct ions.

    For each observed m/z and each adduct type, computes the putative
    neutral mass and matches against the database.

    Args:
        observed_mz: 1D array of observed m/z values.
        database: Dict mapping metabolite names to neutral monoisotopic masses.
        adducts: Dict of adduct name -> mass delta. If None, uses common adducts
            filtered by ion_mode.
        ppm_tolerance: Maximum mass error in ppm.
        ion_mode: 'positive' or 'negative' to filter default adducts.

    Returns:
        List of lists of AdductMatch per observed m/z.
    """
    if adducts is None:
        if ion_mode == "positive":
            adducts = {k: v for k, v in COMMON_ADDUCTS.items() if "+" in k}
        else:
            adducts = {k: v for k, v in COMMON_ADDUCTS.items() if "-" in k}

    db_names = list(database.keys())
    db_masses = np.array(list(database.values()))

    results: list[list[AdductMatch]] = []
    for mz in observed_mz:
        matches: list[AdductMatch] = []
        for adduct_name, adduct_mass in adducts.items():
            neutral = mz - adduct_mass
            if neutral <= 0:
                continue
            ppm_errors = np.abs(db_masses - neutral) / neutral * 1e6
            for i in np.where(ppm_errors <= ppm_tolerance)[0]:
                matches.append(
                    AdductMatch(
                        query_mz=float(mz),
                        neutral_mass=float(neutral),
                        adduct_type=adduct_name,
                        matched_name=db_names[i],
                        delta_ppm=float(ppm_errors[i]),
                    )
                )
        matches.sort(key=lambda m: m.delta_ppm)
        results.append(matches)

    return results


def missing_value_imputation(
    intensities: np.ndarray,
    method: str = "min_half",
) -> np.ndarray:
    """Impute missing (zero or NaN) values in metabolomics intensity data.

    Args:
        intensities: 2D array (metabolites × samples). Zeros and NaNs treated as missing.
        method: Imputation method:
            - 'min_half': Replace with half the minimum non-zero value per metabolite.
            - 'knn': Simple k-nearest neighbors using row (metabolite) correlation.
            - 'median': Replace with row median.

    Returns:
        Imputed intensity matrix.
    """
    data = intensities.copy().astype(float)
    data[data == 0] = np.nan

    if method == "min_half":
        for i in range(data.shape[0]):
            row = data[i]
            non_missing = row[~np.isnan(row)]
            if len(non_missing) > 0:
                fill_val = non_missing.min() / 2.0
            else:
                fill_val = 0.0
            data[i] = np.where(np.isnan(row), fill_val, row)

    elif method == "knn":
        # Simple approach: for each missing value, use mean of k=5 most correlated rows
        k = min(5, data.shape[0] - 1)
        filled = data.copy()
        for i in range(data.shape[0]):
            missing_mask = np.isnan(data[i])
            if not np.any(missing_mask):
                continue
            # Correlation with other rows (using available pairwise values)
            correlations = np.zeros(data.shape[0])
            for j in range(data.shape[0]):
                if i == j:
                    correlations[j] = -np.inf
                    continue
                valid = ~(np.isnan(data[i]) | np.isnan(data[j]))
                if valid.sum() < 2:
                    correlations[j] = -np.inf
                    continue
                ci = data[i, valid] - data[i, valid].mean()
                cj = data[j, valid] - data[j, valid].mean()
                denom = np.sqrt(np.sum(ci**2) * np.sum(cj**2))
                correlations[j] = np.sum(ci * cj) / denom if denom > 0 else 0

            top_k = np.argsort(correlations)[-k:]
            for col in np.where(missing_mask)[0]:
                neighbor_vals = data[top_k, col]
                valid_neighbors = neighbor_vals[~np.isnan(neighbor_vals)]
                filled[i, col] = valid_neighbors.mean() if len(valid_neighbors) > 0 else 0.0
            data[i] = filled[i]

    elif method == "median":
        for i in range(data.shape[0]):
            row = data[i]
            non_missing = row[~np.isnan(row)]
            fill_val = np.median(non_missing) if len(non_missing) > 0 else 0.0
            data[i] = np.where(np.isnan(row), fill_val, row)

    else:
        raise ValueError(f"Unknown imputation method: {method}")

    return np.nan_to_num(data, nan=0.0)
