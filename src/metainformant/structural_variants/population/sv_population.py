"""Population-scale structural variant analysis.

Provides population genotyping, allele frequency computation, association
testing, population structure via PCA, linkage disequilibrium analysis between
SVs and SNPs, and multi-sample callset merging.

Designed for cohort-level SV analysis where the same structural variants are
genotyped across multiple samples and tested for association with phenotypes
or population structure.

Optional dependencies:
    - numpy: Numerical computation (required for most functions)
    - scipy: Statistical tests, distance computations
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional dependencies
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]

try:
    from scipy import stats as scipy_stats

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    scipy_stats = None  # type: ignore[assignment]


def genotype_sv_population(
    sv_calls: list[dict[str, Any]],
    samples: list[str],
    method: str = "depth",
) -> dict[str, Any]:
    """Population-scale SV genotyping across multiple samples.

    For each SV, determines the genotype (0/0, 0/1, or 1/1) in each sample
    using the specified method. The depth method uses read depth signals,
    while the split method uses split-read and discordant-pair evidence.

    Args:
        sv_calls: List of SV call dicts, each containing:
            - ``chrom`` (str): Chromosome.
            - ``start`` (int): Start position.
            - ``end`` (int): End position.
            - ``sv_type`` (str): SV type (DEL, DUP, INV, etc.).
            - ``samples`` (dict): Mapping sample name to evidence dict with
              ``depth_ratio`` (float), ``split_reads`` (int),
              ``discordant_pairs`` (int).
        samples: List of sample names to genotype.
        method: Genotyping method, ``"depth"`` or ``"split"``.

    Returns:
        Dictionary with keys:
            - ``genotype_matrix``: 2D list (n_svs x n_samples) of genotype
              integers (0, 1, or 2 for number of alt alleles).
            - ``quality_scores``: 2D list (n_svs x n_samples) of genotype
              quality scores (0-100).
            - ``allele_frequencies``: List of allele frequencies per SV.

    Raises:
        ValueError: If method is unrecognized.
    """
    if method not in ("depth", "split"):
        raise ValueError(f"Unknown genotyping method: {method}. Use 'depth' or 'split'.")

    n_svs = len(sv_calls)
    n_samples = len(samples)

    genotype_matrix: list[list[int]] = []
    quality_matrix: list[list[float]] = []

    for sv in sv_calls:
        sv_genotypes: list[int] = []
        sv_qualities: list[float] = []
        sample_evidence = sv.get("samples", {})

        for sample in samples:
            evidence = sample_evidence.get(sample, {})

            if method == "depth":
                gt, qual = _genotype_by_depth(
                    evidence.get("depth_ratio", 1.0),
                    sv.get("sv_type", "DEL"),
                )
            else:
                gt, qual = _genotype_by_split(
                    evidence.get("split_reads", 0),
                    evidence.get("discordant_pairs", 0),
                    evidence.get("total_reads", 30),
                )

            sv_genotypes.append(gt)
            sv_qualities.append(qual)

        genotype_matrix.append(sv_genotypes)
        quality_matrix.append(sv_qualities)

    # Compute allele frequencies
    allele_freqs: list[float] = []
    for sv_gt in genotype_matrix:
        total_alleles = len(sv_gt) * 2
        alt_count = sum(sv_gt)
        af = alt_count / total_alleles if total_alleles > 0 else 0.0
        allele_freqs.append(af)

    logger.info(
        "Genotyped %d SVs across %d samples (method=%s)",
        n_svs,
        n_samples,
        method,
    )

    return {
        "genotype_matrix": genotype_matrix,
        "quality_scores": quality_matrix,
        "allele_frequencies": allele_freqs,
    }


def sv_allele_frequency(
    genotype_matrix: Any,
    sample_labels: list[str] | None = None,
) -> dict[str, Any]:
    """Compute SV allele frequencies, optionally stratified by population.

    Calculates minor allele frequency (MAF), site frequency spectrum (SFS),
    and optionally per-population allele frequencies.

    Args:
        genotype_matrix: Genotype matrix (n_svs x n_samples) with values
            0, 1, or 2 (number of alt alleles).
        sample_labels: Optional population labels for each sample (length
            n_samples). If provided, frequencies are computed per population.

    Returns:
        Dictionary with keys:
            - ``frequencies``: List of allele frequencies per SV.
            - ``maf``: List of minor allele frequencies per SV.
            - ``n_polymorphic``: Number of polymorphic SVs (MAF > 0).
            - ``sfs``: Site frequency spectrum as a list of counts.
            - ``population_frequencies``: Dict mapping population label to
              list of per-SV allele frequencies (only if sample_labels given).

    Raises:
        ImportError: If numpy is not available.
    """
    if not HAS_NUMPY:
        raise ImportError("numpy is required: uv pip install numpy")

    gt = np.asarray(genotype_matrix, dtype=np.int32)
    n_svs, n_samples = gt.shape

    # Overall allele frequencies
    alt_counts = gt.sum(axis=1)
    total_alleles = n_samples * 2
    frequencies = (alt_counts / total_alleles).tolist()

    # MAF
    maf = [min(f, 1.0 - f) for f in frequencies]
    n_polymorphic = sum(1 for m in maf if m > 0)

    # Site frequency spectrum (binned alt allele count)
    sfs_bins = total_alleles + 1
    sfs = [0] * sfs_bins
    for ac in alt_counts:
        sfs[int(ac)] += 1

    result: dict[str, Any] = {
        "frequencies": frequencies,
        "maf": maf,
        "n_polymorphic": n_polymorphic,
        "sfs": sfs,
    }

    # Per-population frequencies
    if sample_labels is not None:
        pop_freqs: dict[str, list[float]] = {}
        unique_pops = sorted(set(sample_labels))

        for pop in unique_pops:
            pop_indices = [i for i, l in enumerate(sample_labels) if l == pop]
            if not pop_indices:
                pop_freqs[pop] = [0.0] * n_svs
                continue

            pop_gt = gt[:, pop_indices]
            pop_alt = pop_gt.sum(axis=1)
            pop_total = len(pop_indices) * 2
            pop_freqs[pop] = (pop_alt / pop_total).tolist()

        result["population_frequencies"] = pop_freqs

    logger.info(
        "Computed allele frequencies for %d SVs: %d polymorphic",
        n_svs,
        n_polymorphic,
    )

    return result


def sv_association_test(
    genotypes: list[int],
    phenotypes: list[float],
    covariates: list[list[float]] | None = None,
) -> dict[str, Any]:
    """Test association between SV genotype and phenotype.

    Performs linear regression (continuous phenotype) or logistic regression
    (binary phenotype) of phenotype on genotype, optionally adjusting for
    covariates.

    Args:
        genotypes: Genotype values per sample (0, 1, or 2).
        phenotypes: Phenotype values per sample (continuous or binary 0/1).
        covariates: Optional covariate matrix (n_samples x n_covariates).

    Returns:
        Dictionary with keys:
            - ``beta``: Effect size (regression coefficient for genotype).
            - ``se``: Standard error of beta.
            - ``p_value``: Two-sided p-value.
            - ``odds_ratio``: Odds ratio (only for binary phenotype, else None).
            - ``n_samples``: Number of samples analyzed.

    Raises:
        ImportError: If numpy is not available.
        ValueError: If genotypes and phenotypes have different lengths.
    """
    if not HAS_NUMPY:
        raise ImportError("numpy is required: uv pip install numpy")

    gt = np.asarray(genotypes, dtype=np.float64)
    pheno = np.asarray(phenotypes, dtype=np.float64)

    if len(gt) != len(pheno):
        raise ValueError(f"Length mismatch: {len(gt)} genotypes vs {len(pheno)} phenotypes")

    n = len(gt)

    # Determine if binary phenotype
    unique_pheno = set(pheno.tolist())
    is_binary = unique_pheno.issubset({0.0, 1.0})

    # Build design matrix
    if covariates is not None:
        cov = np.asarray(covariates, dtype=np.float64)
        design = np.column_stack([np.ones(n), gt, cov])
    else:
        design = np.column_stack([np.ones(n), gt])

    if is_binary:
        beta, se, p_value = _logistic_regression(design, pheno)
    else:
        beta, se, p_value = _linear_regression(design, pheno)

    # Genotype coefficient is at index 1
    gt_beta = float(beta[1]) if len(beta) > 1 else 0.0
    gt_se = float(se[1]) if len(se) > 1 else 0.0
    gt_pvalue = float(p_value[1]) if len(p_value) > 1 else 1.0

    odds_ratio = math.exp(gt_beta) if is_binary else None

    logger.info(
        "Association test: beta=%.4f, se=%.4f, p=%.2e, n=%d",
        gt_beta,
        gt_se,
        gt_pvalue,
        n,
    )

    return {
        "beta": gt_beta,
        "se": gt_se,
        "p_value": gt_pvalue,
        "odds_ratio": odds_ratio,
        "n_samples": n,
    }


def sv_population_structure(
    genotype_matrix: Any,
    n_components: int = 10,
) -> dict[str, Any]:
    """PCA on SV genotype matrix for population structure analysis.

    Centers and scales the genotype matrix, then computes the top principal
    components to visualize population structure driven by structural
    variant genotypes.

    Args:
        genotype_matrix: Genotype matrix (n_svs x n_samples) with values
            0, 1, or 2.
        n_components: Number of principal components to compute.

    Returns:
        Dictionary with keys:
            - ``pcs``: 2D array (n_samples x n_components) of PC coordinates.
            - ``eigenvalues``: List of eigenvalues for each component.
            - ``variance_explained``: List of fraction of variance explained
              per component.
            - ``loadings``: 2D array (n_svs x n_components) of SV loadings.

    Raises:
        ImportError: If numpy is not available.
    """
    if not HAS_NUMPY:
        raise ImportError("numpy is required: uv pip install numpy")

    gt = np.asarray(genotype_matrix, dtype=np.float64)

    # Transpose: PCA on samples, so we need samples as rows
    # gt is (n_svs x n_samples), transpose to (n_samples x n_svs)
    data = gt.T
    n_samples, n_svs = data.shape

    n_components = min(n_components, min(n_samples, n_svs))

    # Center and scale
    mean = data.mean(axis=0)
    std = data.std(axis=0)
    std[std == 0] = 1.0
    data_scaled = (data - mean) / std

    # Covariance matrix (sample x sample)
    if n_samples < n_svs:
        # Compute n_samples x n_samples covariance
        cov = data_scaled @ data_scaled.T / max(n_svs - 1, 1)
        eigenvalues, eigenvectors = np.linalg.eigh(cov)

        # Sort by descending eigenvalue
        idx = np.argsort(eigenvalues)[::-1][:n_components]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]

        pcs = eigenvectors  # (n_samples x n_components)

        # Loadings: project back to SV space
        loadings = data_scaled.T @ pcs  # (n_svs x n_components)
        for j in range(n_components):
            norm = np.linalg.norm(loadings[:, j])
            if norm > 0:
                loadings[:, j] /= norm
    else:
        # Standard PCA on SV x SV covariance
        cov = data_scaled.T @ data_scaled / max(n_samples - 1, 1)
        eigenvalues, eigenvectors = np.linalg.eigh(cov)

        idx = np.argsort(eigenvalues)[::-1][:n_components]
        eigenvalues = eigenvalues[idx]
        loadings = eigenvectors[:, idx]

        pcs = data_scaled @ loadings

    total_var = float(np.sum(np.abs(eigenvalues)))
    variance_explained = [float(abs(ev)) / total_var if total_var > 0 else 0.0 for ev in eigenvalues]

    logger.info(
        "PCA on %d SVs x %d samples: top %d PCs explain %.1f%% variance",
        n_svs,
        n_samples,
        n_components,
        sum(variance_explained) * 100,
    )

    return {
        "pcs": pcs,
        "eigenvalues": eigenvalues.tolist(),
        "variance_explained": variance_explained,
        "loadings": loadings,
    }


def sv_ld_analysis(
    genotype_matrix_sv: Any,
    genotype_matrix_snp: Any,
    sv_positions: list[int],
    snp_positions: list[int],
) -> dict[str, Any]:
    """Compute LD between SVs and nearby SNPs.

    For each SV, computes the squared Pearson correlation (r^2) with all
    SNPs within a window, identifying tag SNPs (SNPs in highest LD with
    each SV).

    Args:
        genotype_matrix_sv: SV genotype matrix (n_svs x n_samples), values
            0/1/2.
        genotype_matrix_snp: SNP genotype matrix (n_snps x n_samples),
            values 0/1/2.
        sv_positions: Genomic position for each SV.
        snp_positions: Genomic position for each SNP.

    Returns:
        Dictionary with keys:
            - ``ld_matrix``: 2D list (n_svs x n_snps) of r^2 values.
            - ``tag_snps``: List of dicts per SV with ``snp_index``,
              ``snp_position``, ``r2``.
            - ``max_r2_per_sv``: List of maximum r^2 value for each SV.

    Raises:
        ImportError: If numpy is not available.
    """
    if not HAS_NUMPY:
        raise ImportError("numpy is required: uv pip install numpy")

    sv_gt = np.asarray(genotype_matrix_sv, dtype=np.float64)
    snp_gt = np.asarray(genotype_matrix_snp, dtype=np.float64)

    n_svs = sv_gt.shape[0]
    n_snps = snp_gt.shape[0]

    ld_matrix: list[list[float]] = []
    tag_snps: list[dict[str, Any]] = []
    max_r2_per_sv: list[float] = []

    for i in range(n_svs):
        sv_row = sv_gt[i, :]
        sv_var = float(np.var(sv_row))

        row_r2: list[float] = []
        best_r2 = 0.0
        best_snp_idx = 0

        for j in range(n_snps):
            snp_row = snp_gt[j, :]
            snp_var = float(np.var(snp_row))

            if sv_var == 0 or snp_var == 0:
                r2 = 0.0
            else:
                cov = float(np.mean((sv_row - np.mean(sv_row)) * (snp_row - np.mean(snp_row))))
                r2 = (cov**2) / (sv_var * snp_var)

            row_r2.append(r2)

            if r2 > best_r2:
                best_r2 = r2
                best_snp_idx = j

        ld_matrix.append(row_r2)
        max_r2_per_sv.append(best_r2)
        tag_snps.append(
            {
                "snp_index": best_snp_idx,
                "snp_position": snp_positions[best_snp_idx] if best_snp_idx < len(snp_positions) else 0,
                "r2": best_r2,
            }
        )

    logger.info(
        "LD analysis: %d SVs x %d SNPs, mean max_r2=%.3f",
        n_svs,
        n_snps,
        sum(max_r2_per_sv) / len(max_r2_per_sv) if max_r2_per_sv else 0.0,
    )

    return {
        "ld_matrix": ld_matrix,
        "tag_snps": tag_snps,
        "max_r2_per_sv": max_r2_per_sv,
    }


def merge_sv_callsets(
    callsets: list[list[dict[str, Any]]],
    min_overlap: float = 0.5,
    max_breakpoint_distance: int = 500,
) -> list[dict[str, Any]]:
    """Merge SV calls across samples using reciprocal overlap and breakpoint proximity.

    Compares SV calls from multiple samples/callers and merges variants that
    represent the same event based on reciprocal overlap fraction and
    breakpoint distance criteria.

    Args:
        callsets: List of callsets, each a list of SV call dicts containing:
            - ``chrom`` (str): Chromosome.
            - ``start`` (int): Start position.
            - ``end`` (int): End position.
            - ``sv_type`` (str): SV type.
            - ``sample`` (str, optional): Source sample name.
        min_overlap: Minimum reciprocal overlap fraction (0-1) for merging.
        max_breakpoint_distance: Maximum distance between breakpoints for
            merging.

    Returns:
        List of merged SV dicts, each containing:
            - ``chrom``, ``start``, ``end``, ``sv_type``: Consensus call.
            - ``n_samples``: Number of samples supporting this SV.
            - ``samples``: List of sample names.
            - ``support_count``: Total number of supporting calls.
            - ``merged_calls``: List of original calls merged into this entry.
    """
    # Flatten all calls with source tracking
    all_calls: list[dict[str, Any]] = []
    for cs_idx, callset in enumerate(callsets):
        for call in callset:
            entry = dict(call)
            if "sample" not in entry:
                entry["sample"] = f"sample_{cs_idx}"
            entry["_source_idx"] = cs_idx
            all_calls.append(entry)

    if not all_calls:
        return []

    # Sort by chromosome and start position
    all_calls.sort(key=lambda c: (c.get("chrom", ""), c.get("start", 0)))

    # Greedy merging
    merged: list[dict[str, Any]] = []
    used: set[int] = set()

    for i in range(len(all_calls)):
        if i in used:
            continue

        cluster = [all_calls[i]]
        used.add(i)

        for j in range(i + 1, len(all_calls)):
            if j in used:
                continue

            if _should_merge(all_calls[i], all_calls[j], min_overlap, max_breakpoint_distance):
                cluster.append(all_calls[j])
                used.add(j)

        # Build consensus call
        consensus = _build_consensus(cluster)
        merged.append(consensus)

    logger.info(
        "Merged %d calls from %d callsets into %d consensus SVs",
        len(all_calls),
        len(callsets),
        len(merged),
    )

    return merged


# ---------------------------------------------------------------------------
# Internal helper functions
# ---------------------------------------------------------------------------


def _genotype_by_depth(
    depth_ratio: float,
    sv_type: str,
) -> tuple[int, float]:
    """Genotype an SV based on read depth ratio.

    For deletions, depth_ratio < 0.75 suggests het, < 0.25 suggests hom.
    For duplications, depth_ratio > 1.25 suggests het, > 1.75 suggests hom.

    Args:
        depth_ratio: Observed/expected read depth ratio.
        sv_type: SV type string.

    Returns:
        Tuple of (genotype, quality) where genotype is 0/1/2.
    """
    if sv_type in ("DEL", "deletion"):
        if depth_ratio < 0.25:
            quality = min(99.0, (0.25 - depth_ratio) * 400)
            return 2, quality
        elif depth_ratio < 0.75:
            quality = min(99.0, abs(0.5 - depth_ratio) * 200)
            return 1, quality
        else:
            quality = min(99.0, (depth_ratio - 0.75) * 200)
            return 0, quality
    elif sv_type in ("DUP", "duplication"):
        if depth_ratio > 1.75:
            quality = min(99.0, (depth_ratio - 1.75) * 200)
            return 2, quality
        elif depth_ratio > 1.25:
            quality = min(99.0, abs(depth_ratio - 1.5) * 200)
            return 1, quality
        else:
            quality = min(99.0, (1.25 - depth_ratio) * 200)
            return 0, quality
    else:
        # For INV/TRA, depth is less informative
        if depth_ratio < 0.5:
            return 1, 30.0
        return 0, 30.0


def _genotype_by_split(
    split_reads: int,
    discordant_pairs: int,
    total_reads: int,
) -> tuple[int, float]:
    """Genotype an SV based on split-read and discordant-pair evidence.

    Uses the fraction of supporting reads to estimate genotype.

    Args:
        split_reads: Number of split reads supporting the SV.
        discordant_pairs: Number of discordant pairs.
        total_reads: Total reads in the region.

    Returns:
        Tuple of (genotype, quality).
    """
    support = split_reads + discordant_pairs
    total = max(total_reads, 1)
    fraction = support / total

    if fraction >= 0.85:
        return 2, min(99.0, support * 5.0)
    elif fraction >= 0.15:
        return 1, min(99.0, support * 5.0)
    else:
        return 0, min(99.0, (total - support) * 2.0)


def _linear_regression(
    design: Any,
    y: Any,
) -> tuple[Any, Any, Any]:
    """Ordinary least squares linear regression.

    Args:
        design: Design matrix (n x p) including intercept.
        y: Response vector (n,).

    Returns:
        Tuple of (coefficients, standard_errors, p_values).
    """
    n, p = design.shape

    # OLS: beta = (X^T X)^{-1} X^T y
    xtx = design.T @ design
    xtx += np.eye(p) * 1e-8  # Regularization
    xtx_inv = np.linalg.inv(xtx)
    beta = xtx_inv @ design.T @ y

    # Residuals
    residuals = y - design @ beta
    rss = float(np.sum(residuals**2))
    sigma2 = rss / max(n - p, 1)

    # Standard errors
    se = np.sqrt(np.diag(xtx_inv) * sigma2)
    se[se == 0] = 1e-10

    # T-statistics and p-values
    t_stats = beta / se
    p_values = np.array([2.0 * (1.0 - _normal_cdf_np(abs(t))) for t in t_stats])

    return beta, se, p_values


def _logistic_regression(
    design: Any,
    y: Any,
    max_iter: int = 50,
) -> tuple[Any, Any, Any]:
    """Iteratively reweighted least squares logistic regression.

    Args:
        design: Design matrix (n x p) including intercept.
        y: Binary response vector (n,).
        max_iter: Maximum IRLS iterations.

    Returns:
        Tuple of (coefficients, standard_errors, p_values).
    """
    n, p = design.shape
    beta = np.zeros(p, dtype=np.float64)

    for _iteration in range(max_iter):
        eta = design @ beta
        # Clip to avoid overflow
        eta = np.clip(eta, -20, 20)
        mu = 1.0 / (1.0 + np.exp(-eta))

        # Weight matrix diagonal
        w = mu * (1.0 - mu)
        w = np.maximum(w, 1e-10)

        # Weighted least squares update
        z = eta + (y - mu) / w
        w_diag = np.diag(w)

        xtwx = design.T @ w_diag @ design
        xtwx += np.eye(p) * 1e-8
        xtwx_inv = np.linalg.inv(xtwx)
        beta_new = xtwx_inv @ design.T @ w_diag @ z

        # Check convergence
        if np.max(np.abs(beta_new - beta)) < 1e-6:
            beta = beta_new
            break
        beta = beta_new

    # Standard errors from final information matrix
    eta = design @ beta
    eta = np.clip(eta, -20, 20)
    mu = 1.0 / (1.0 + np.exp(-eta))
    w = mu * (1.0 - mu)
    w = np.maximum(w, 1e-10)

    info_matrix = design.T @ np.diag(w) @ design
    info_matrix += np.eye(p) * 1e-8

    try:
        cov_matrix = np.linalg.inv(info_matrix)
        se = np.sqrt(np.maximum(np.diag(cov_matrix), 0.0))
    except np.linalg.LinAlgError:
        se = np.ones(p) * 1e10

    se[se == 0] = 1e-10

    # Wald test p-values
    z_stats = beta / se
    p_values = np.array([2.0 * (1.0 - _normal_cdf_np(abs(z))) for z in z_stats])

    return beta, se, p_values


def _normal_cdf_np(x: float) -> float:
    """Standard normal CDF approximation."""
    return 0.5 * (1.0 + math.erf(float(x) / math.sqrt(2.0)))


def _should_merge(
    call_a: dict[str, Any],
    call_b: dict[str, Any],
    min_overlap: float,
    max_bp_dist: int,
) -> bool:
    """Determine if two SV calls should be merged.

    Checks chromosome, SV type, reciprocal overlap, and breakpoint distance.

    Args:
        call_a: First SV call dict.
        call_b: Second SV call dict.
        min_overlap: Minimum reciprocal overlap fraction.
        max_bp_dist: Maximum breakpoint distance.

    Returns:
        True if calls should be merged.
    """
    # Must be same chromosome and type
    if call_a.get("chrom", "") != call_b.get("chrom", ""):
        return False
    if call_a.get("sv_type", "") != call_b.get("sv_type", ""):
        return False

    start_a = call_a.get("start", 0)
    end_a = call_a.get("end", 0)
    start_b = call_b.get("start", 0)
    end_b = call_b.get("end", 0)

    # Breakpoint distance check
    if abs(start_a - start_b) > max_bp_dist:
        return False
    if abs(end_a - end_b) > max_bp_dist:
        return False

    # Reciprocal overlap
    overlap_start = max(start_a, start_b)
    overlap_end = min(end_a, end_b)
    overlap = max(0, overlap_end - overlap_start)

    size_a = max(end_a - start_a, 1)
    size_b = max(end_b - start_b, 1)

    recip_a = overlap / size_a
    recip_b = overlap / size_b

    return min(recip_a, recip_b) >= min_overlap


def _build_consensus(
    cluster: list[dict[str, Any]],
) -> dict[str, Any]:
    """Build a consensus SV call from a cluster of merged calls.

    Uses median breakpoints and majority SV type.

    Args:
        cluster: List of SV call dicts to merge.

    Returns:
        Consensus SV call dict.
    """
    starts = [c.get("start", 0) for c in cluster]
    ends = [c.get("end", 0) for c in cluster]

    starts.sort()
    ends.sort()
    n = len(cluster)

    consensus_start = starts[n // 2]
    consensus_end = ends[n // 2]

    # Majority SV type
    type_counts: dict[str, int] = defaultdict(int)
    for c in cluster:
        type_counts[c.get("sv_type", "UNKNOWN")] += 1
    consensus_type = max(type_counts, key=type_counts.get)  # type: ignore[arg-type]

    # Collect samples
    samples: list[str] = []
    for c in cluster:
        sample = c.get("sample", "")
        if sample and sample not in samples:
            samples.append(sample)

    return {
        "chrom": cluster[0].get("chrom", ""),
        "start": consensus_start,
        "end": consensus_end,
        "sv_type": consensus_type,
        "n_samples": len(samples),
        "samples": samples,
        "support_count": len(cluster),
        "merged_calls": cluster,
    }
