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
import statistics
from collections import defaultdict
from typing import Any

from metainformant.core.utils.logging import get_logger
from metainformant.structural_variants._breakends import (
    canonical_junction,
    consensus_junction,
    is_translocation,
    junctions_agree,
)

logger = get_logger(__name__)

# Optional dependencies
np: Any  # None when numpy is unavailable (optional dependency)
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None

scipy_stats: Any  # None when scipy is unavailable (optional dependency)
try:
    from scipy import stats as scipy_stats

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    scipy_stats = None


def genotype_sv_population(
    sv_calls: list[dict[str, Any]],
    samples: list[str],
    method: str = "depth",
) -> dict[str, Any]:
    """Population-scale SV genotyping across multiple samples.

    For each SV, determines the genotype (0/0, 0/1, or 1/1) in each sample
    using the specified method. The depth method uses read depth signals,
    while the split method uses split-read and discordant-pair evidence.

    A sample with no usable evidence for an SV (no entry in ``samples``, or no
    ``depth_ratio`` for the depth method / no ``total_reads`` for the split
    method) is called 0 (reference) with quality 0.0: nothing was observed, so it
    counts as reference in the allele frequency (which therefore understates the
    frequency when many samples lack data) but carries no confidence. The number
    of such calls is logged as a warning.

    Args:
        sv_calls: List of SV call dicts, each containing:
            - ``chrom`` (str): Chromosome.
            - ``start`` (int): Start position.
            - ``end`` (int): End position.
            - ``sv_type`` (str): SV type (DEL, DUP, INV, etc.); required by the
              depth method, whose thresholds depend on it.
            - ``samples`` (dict): Mapping sample name to evidence dict with
              ``depth_ratio`` (float), ``split_reads`` (int),
              ``discordant_pairs`` (int), ``total_reads`` (int). The split
              method treats missing ``split_reads``/``discordant_pairs`` as 0.
        samples: List of sample names to genotype.
        method: Genotyping method, ``"depth"`` or ``"split"``.

    Returns:
        Dictionary with keys:
            - ``genotype_matrix``: 2D list (n_svs x n_samples) of genotype
              integers (0, 1, or 2 for number of alt alleles).
            - ``quality_scores``: 2D list (n_svs x n_samples) of genotype
              quality scores (0-100; 0.0 where there was no evidence).
            - ``allele_frequencies``: List of allele frequencies per SV.

    Raises:
        ValueError: If method is unrecognized, or an SV has no ``sv_type`` under
            the depth method.
    """
    if method not in ("depth", "split"):
        raise ValueError(f"Unknown genotyping method: {method}. Use 'depth' or 'split'.")

    n_svs = len(sv_calls)
    n_samples = len(samples)

    genotype_matrix: list[list[int]] = []
    quality_matrix: list[list[float]] = []
    n_without_evidence = 0

    for index, sv in enumerate(sv_calls):
        sv_genotypes: list[int] = []
        sv_qualities: list[float] = []
        sample_evidence = sv.get("samples", {})
        if method == "depth" and sv.get("sv_type") is None:
            raise ValueError(f"sv_calls[{index}] has no 'sv_type': the depth method needs it to pick thresholds")

        for sample in samples:
            evidence = sample_evidence.get(sample, {})

            if method == "depth":
                depth_ratio = evidence.get("depth_ratio")
                observed = depth_ratio is not None
                gt, qual = _genotype_by_depth(depth_ratio, sv["sv_type"]) if observed else (0, 0.0)
            else:
                total_reads = evidence.get("total_reads")
                observed = total_reads is not None
                gt, qual = (
                    _genotype_by_split(evidence.get("split_reads", 0), evidence.get("discordant_pairs", 0), total_reads)
                    if observed
                    else (0, 0.0)
                )

            n_without_evidence += not observed
            sv_genotypes.append(gt)
            sv_qualities.append(qual)

        genotype_matrix.append(sv_genotypes)
        quality_matrix.append(sv_qualities)

    if n_without_evidence:
        logger.warning(
            "%d of %d sample genotypes had no %s evidence: called 0 (reference) with quality 0.0",
            n_without_evidence,
            n_svs * n_samples,
            "depth_ratio" if method == "depth" else "total_reads",
        )

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
        ValueError: If a genotype is not 0, 1 or 2 (a missing-genotype code such
            as -1 would otherwise be summed as an allele count), or
            ``sample_labels`` does not have one label per sample.
    """
    if not HAS_NUMPY:
        raise ImportError("numpy is required: uv pip install numpy")

    raw = np.asarray(genotype_matrix, dtype=np.float64)
    if not np.all(np.isin(raw, (0.0, 1.0, 2.0))):
        raise ValueError("genotype_matrix must contain only 0, 1 and 2 (remove or impute missing genotypes first)")
    gt = raw.astype(np.int32)
    n_svs, n_samples = gt.shape
    if sample_labels is not None and len(sample_labels) != n_samples:
        raise ValueError(f"sample_labels has {len(sample_labels)} entries for {n_samples} samples")

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
            pop_indices = [i for i, label in enumerate(sample_labels) if label == pop]
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

    A test that cannot be carried out is reported as missing, not as "no
    association": ``beta``, ``se`` and ``p_value`` are NaN (and ``odds_ratio``
    NaN for a binary phenotype) and ``error`` gives the reason. That happens
    when the design matrix is singular (a monomorphic SV, or covariates
    collinear with the genotype; checked by its rank, never by whether the
    solver raises), when a linear fit has no residual degrees of freedom or no
    residual variance (a constant phenotype or an exact fit), and when the
    logistic fit does not converge (complete or quasi-complete separation, or
    a constant phenotype).

    Args:
        genotypes: Genotype values per sample (0, 1, or 2).
        phenotypes: Phenotype values per sample (continuous or binary 0/1).
        covariates: Optional covariate matrix (n_samples x n_covariates).

    Returns:
        Dictionary with keys:
            - ``beta``: Effect size (regression coefficient for genotype).
            - ``se``: Standard error of beta.
            - ``p_value``: Two-sided p-value: the t distribution with the
              residual degrees of freedom for a continuous phenotype, the
              normal distribution (Wald test) for a binary one. Without SciPy
              the normal distribution stands in for the t distribution, which
              understates the p-value for small samples; ``p_value_method`` says
              which was used.
            - ``p_value_method``: ``"t"``, ``"normal_approximation_of_t"`` or
              ``"wald_normal"``.
            - ``odds_ratio``: Odds ratio (only for binary phenotype, else None).
            - ``n_samples``: Number of samples analyzed.
            - ``error``: None, or why the test could not be carried out.

    Raises:
        ImportError: If numpy is not available.
        ValueError: If genotypes and phenotypes have different lengths, or any
            value is not finite.
    """
    if not HAS_NUMPY:
        raise ImportError("numpy is required: uv pip install numpy")

    gt = np.asarray(genotypes, dtype=np.float64)
    pheno = np.asarray(phenotypes, dtype=np.float64)

    if len(gt) != len(pheno):
        raise ValueError(f"Length mismatch: {len(gt)} genotypes vs {len(pheno)} phenotypes")

    n = len(gt)

    # Build design matrix
    if covariates is not None:
        cov = np.asarray(covariates, dtype=np.float64)
        design = np.column_stack([np.ones(n), gt, cov])
    else:
        design = np.column_stack([np.ones(n), gt])

    if not (np.all(np.isfinite(design)) and np.all(np.isfinite(pheno))):
        raise ValueError("genotypes, phenotypes and covariates must be finite (remove or impute missing values)")

    # Determine if binary phenotype
    unique_pheno = set(pheno.tolist())
    is_binary = unique_pheno.issubset({0.0, 1.0})

    n_params = design.shape[1]
    # The rank test catches a monomorphic SV or a collinear covariate on every platform, including those whose
    # solver would still invert the nearly singular matrix without raising.
    if np.linalg.matrix_rank(design) < n_params:
        return _failed_association(n, is_binary, "design matrix is singular (monomorphic SV or collinear covariates)")

    if is_binary:
        fit = _logistic_regression(design, pheno)
        method = "wald_normal"
        df = None
    else:
        df = n - n_params
        if df <= 0:
            return _failed_association(
                n, is_binary, f"no residual degrees of freedom ({n} samples, {n_params} parameters)"
            )
        fit = _linear_regression(design, pheno)
        method = "t" if HAS_SCIPY else "normal_approximation_of_t"
    if isinstance(fit, str):
        return _failed_association(n, is_binary, fit)
    beta, se = fit

    # Genotype coefficient is at index 1
    gt_beta = float(beta[1])
    gt_se = float(se[1])
    gt_pvalue = _two_sided_p_value(gt_beta / gt_se, df)

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
        "p_value_method": method,
        "odds_ratio": odds_ratio,
        "n_samples": n,
        "error": None,
    }


def sv_population_structure(
    genotype_matrix: Any,
    n_components: int = 10,
) -> dict[str, Any]:
    """PCA on SV genotype matrix for population structure analysis.

    Centers each SV to mean 0 and scales it to unit (population) variance, giving
    the standardised matrix ``X`` (n_samples x n_svs). With the singular value
    decomposition ``X = U S Vt``, every returned quantity is defined the same way
    whatever the shape of the data:

    * ``pcs``: the sample scores ``X V = U S`` (principal-component coordinates);
    * ``eigenvalues``: the eigenvalues of the sample covariance matrix
      ``Xt X / (n_samples - 1)``, i.e. ``S**2 / (n_samples - 1)``;
    * ``loadings``: ``V``, the unit-length eigenvectors of that covariance matrix.

    When there are fewer samples than SVs the decomposition is computed from the
    smaller ``n_samples x n_samples`` matrix ``X Xt`` (same eigenvalues, after
    dividing by ``n_samples - 1``), and the scores and loadings are recovered from
    its eigenvectors; the result is the same as for the covariance matrix.

    The sign of each component is arbitrary (it depends on the LAPACK build), so
    compare components up to sign. A component beyond the rank of the data (the
    centered data has rank at most ``n_samples - 1``) has eigenvalue 0, scores 0
    and loadings 0: its direction is not defined.

    Args:
        genotype_matrix: Genotype matrix (n_svs x n_samples) with values
            0, 1, or 2.
        n_components: Number of principal components to compute (at least 1; at
            most ``min(n_samples, n_svs)`` are returned).

    Returns:
        Dictionary with keys:
            - ``pcs``: 2D array (n_samples x n_components) of PC coordinates
              (sample scores).
            - ``eigenvalues``: List of covariance-matrix eigenvalues for each
              component, in descending order.
            - ``variance_explained``: List of the fraction of the total variance
              (the sum of all eigenvalues of the covariance matrix, not just the
              returned components) that each component explains. It sums to less
              than 1 when ``n_components`` is smaller than the rank of the data.
            - ``loadings``: 2D array (n_svs x n_components) of SV loadings.

    Raises:
        ImportError: If numpy is not available.
        ValueError: If the matrix is not 2D or has non-finite values, or
            ``n_components`` is below 1.
    """
    if not HAS_NUMPY:
        raise ImportError("numpy is required: uv pip install numpy")

    if n_components < 1:
        raise ValueError(f"n_components must be at least 1, got {n_components}")

    gt = np.asarray(genotype_matrix, dtype=np.float64)
    if gt.ndim != 2:
        raise ValueError(f"genotype_matrix must be 2D (n_svs x n_samples), got shape {gt.shape}")
    if not np.all(np.isfinite(gt)):
        raise ValueError("genotype_matrix has non-finite values: remove or impute missing genotypes first")

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

    dof = max(n_samples - 1, 1)

    if n_samples < n_svs:
        # Gram matrix X Xt (n_samples x n_samples): eigenvalues S**2, eigenvectors U.
        squared_singular, u = np.linalg.eigh(data_scaled @ data_scaled.T)
        order = np.argsort(squared_singular)[::-1]
        squared_singular = np.clip(squared_singular[order], 0.0, None)
        u = u[:, order]
        all_eigenvalues = squared_singular / dof  # eigenvalues of the covariance matrix Xt X / dof
        keep = all_eigenvalues[:n_components]
        singular = np.sqrt(squared_singular[:n_components])
        pcs = u[:, :n_components] * singular  # U S
        defined = _defined_components(keep, all_eigenvalues, n_samples, n_svs)
        # V = Xt U / S, where the component has a direction
        loadings = data_scaled.T @ u[:, :n_components]
        loadings[:, defined] /= singular[defined]
        loadings[:, ~defined] = 0.0
        pcs[:, ~defined] = 0.0
    else:
        # Covariance matrix Xt X / dof (n_svs x n_svs).
        values, vectors = np.linalg.eigh(data_scaled.T @ data_scaled / dof)
        order = np.argsort(values)[::-1]
        all_eigenvalues = np.clip(values[order], 0.0, None)
        vectors = vectors[:, order]
        keep = all_eigenvalues[:n_components]
        defined = _defined_components(keep, all_eigenvalues, n_samples, n_svs)
        loadings = vectors[:, :n_components].copy()
        loadings[:, ~defined] = 0.0
        pcs = data_scaled @ loadings  # X V

    total_var = float(np.sum(all_eigenvalues))
    eigenvalues = np.where(defined, keep, 0.0)

    # Each component's share of the total variance, which comes from all eigenvalues
    # (taken above, before the top components were selected).
    variance_explained = [float(ev) / total_var if total_var > 0 else 0.0 for ev in eigenvalues]

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
    """Compute LD between SVs and SNPs.

    For each SV, computes the squared Pearson correlation (r^2) with every
    SNP in ``genotype_matrix_snp`` (there is no window: the positions are not
    used to select SNPs, so pass only the SNPs near the SVs), identifying the
    tag SNP (the SNP in highest LD with each SV).

    Args:
        genotype_matrix_sv: SV genotype matrix (n_svs x n_samples), values
            0/1/2.
        genotype_matrix_snp: SNP genotype matrix (n_snps x n_samples),
            values 0/1/2.
        sv_positions: Genomic position for each SV (one per SV; checked, not
            used in the computation).
        snp_positions: Genomic position for each SNP (one per SNP); the
            position of the tag SNP is reported from it.

    Returns:
        Dictionary with keys:
            - ``ld_matrix``: 2D list (n_svs x n_snps) of r^2 values.
            - ``tag_snps``: List of dicts per SV with ``snp_index``,
              ``snp_position``, ``r2``. When no SNP has r^2 above 0 (for example
              a monomorphic SV) there is no tag SNP: the first SNP is listed
              with ``r2`` 0.0 (``snp_index`` and ``snp_position`` are None when
              there are no SNPs at all).
            - ``max_r2_per_sv``: List of maximum r^2 value for each SV.

    Raises:
        ImportError: If numpy is not available.
        ValueError: If a position list does not have one entry per SV or SNP.
    """
    if not HAS_NUMPY:
        raise ImportError("numpy is required: uv pip install numpy")

    sv_gt = np.asarray(genotype_matrix_sv, dtype=np.float64)
    snp_gt = np.asarray(genotype_matrix_snp, dtype=np.float64)

    n_svs = sv_gt.shape[0]
    n_snps = snp_gt.shape[0]

    if len(sv_positions) != n_svs:
        raise ValueError(f"sv_positions has {len(sv_positions)} entries for {n_svs} SVs")
    if len(snp_positions) != n_snps:
        raise ValueError(f"snp_positions has {len(snp_positions)} entries for {n_snps} SNPs")

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
                "snp_index": best_snp_idx if n_snps else None,
                "snp_position": snp_positions[best_snp_idx] if n_snps else None,
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

    Calls are clustered greedily in chromosome and start order: the first call
    not yet merged anchors a cluster, and a call joins it only when it matches
    the anchor itself (a chain of calls that each match the next is not
    followed).

    A translocation (``TRA``/``BND``) has no interval, so reciprocal overlap
    cannot match it. It is matched on both breakpoints instead: ``chrom`` and
    ``start`` of one call against the other's, and ``chrom2`` and ``pos2``
    likewise, each pair of positions within ``max_breakpoint_distance``,
    whichever breakpoint each call wrote first. ``min_overlap`` does not apply
    to it, and it never merges with an interval variant. The consensus of
    translocations takes the median of each breakpoint, in (chromosome,
    position) order, with ``end == start`` and the second breakpoint in
    ``chrom2``/``pos2``.

    Args:
        callsets: List of callsets, each a list of SV call dicts containing:
            - ``chrom`` (str): Chromosome.
            - ``start`` (int): Start position.
            - ``end`` (int): End position.
            - ``sv_type`` (str): SV type.
            - ``chrom2`` (str), ``pos2`` (int): Second breakpoint (translocations).
            - ``sample`` (str, optional): Source sample name; a call without one
              is named ``sample_<index of its callset>``.
        min_overlap: Minimum reciprocal overlap fraction (0-1) for merging.
        max_breakpoint_distance: Maximum distance between breakpoints for
            merging.

    Returns:
        List of merged SV dicts, each containing:
            - ``chrom``, ``start``, ``end``, ``sv_type``: Consensus call (and
              ``chrom2``, ``pos2`` for a translocation).
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


def _defined_components(kept: Any, all_eigenvalues: Any, n_samples: int, n_svs: int) -> Any:
    """Which of the kept components have a direction: eigenvalue above the numerical noise of the decomposition.

    The eigenvalues of a covariance matrix cannot be negative; the tiny ones that rounding produces (and the
    ``n_samples - 1`` rank limit of centered data) are noise, in the sense of :func:`numpy.linalg.matrix_rank`.

    Args:
        kept: Eigenvalues of the returned components.
        all_eigenvalues: Every eigenvalue of the covariance matrix, descending (not only the leading ones).
        n_samples: Number of samples.
        n_svs: Number of SVs.

    Returns:
        Boolean array, True for each kept component that is above the tolerance.
    """
    if all_eigenvalues.size == 0 or float(all_eigenvalues[0]) <= 0.0:
        return np.zeros(len(kept), dtype=bool)
    tolerance = max(n_samples, n_svs) * np.finfo(np.float64).eps * float(all_eigenvalues[0])
    return kept > tolerance


# Largest linear predictor the logistic fit works with: exp(-20) is below 1e-8, so the fitted probability is 0 or 1.
_MAX_LINEAR_PREDICTOR = 20.0


def _failed_association(n: int, is_binary: bool, error: str) -> dict[str, Any]:
    """Result of an association test that could not be carried out: NaN statistics and the reason."""
    nan = float("nan")
    return {
        "beta": nan,
        "se": nan,
        "p_value": nan,
        "p_value_method": None,
        "odds_ratio": nan if is_binary else None,
        "n_samples": n,
        "error": error,
    }


def _two_sided_p_value(statistic: float, df: int | None) -> float:
    """Two-sided p-value of a Wald or t statistic.

    Args:
        statistic: The statistic.
        df: Residual degrees of freedom for a t statistic; None for a Wald (normal) statistic.

    Returns:
        The t distribution's tail when ``df`` is given and SciPy is available, else the normal tail. The normal
        tail is ``erfc(|z| / sqrt(2))``, exact far into the tail (``2 * (1 - cdf)`` underflows to 0 beyond about 8).
    """
    if df is not None and HAS_SCIPY:
        return float(2.0 * scipy_stats.t.sf(abs(statistic), df))
    return math.erfc(abs(statistic) / math.sqrt(2.0))


def _linear_regression(
    design: Any,
    y: Any,
) -> tuple[Any, Any] | str:
    """Ordinary least squares linear regression.

    Args:
        design: Design matrix (n x p) of full column rank, including intercept.
        y: Response vector (n,).

    Returns:
        ``(coefficients, standard_errors)``, or the reason why the standard errors are undefined (a constant
        response, or an exact fit).
    """
    n, p = design.shape

    # Solve through the QR decomposition rather than the normal equations: squaring the condition number of a
    # nearly collinear design would make the answer depend on the BLAS build. The caller has already rejected
    # a design that is rank deficient.
    q, r = np.linalg.qr(design)
    try:
        beta = np.linalg.solve(r, q.T @ y)
        r_inv = np.linalg.inv(r)
    except np.linalg.LinAlgError:
        return "design matrix is numerically singular"
    xtx_inv = r_inv @ r_inv.T  # (Xt X)^-1 = R^-1 R^-T

    # Residuals
    residuals = y - design @ beta
    rss = float(np.sum(residuals**2))
    tss = float(np.sum((y - np.mean(y)) ** 2))
    if tss == 0.0:
        return "phenotype is constant"
    if rss <= 1e-12 * tss:
        return "zero residual variance (exact fit): the test statistic is undefined"
    sigma2 = rss / (n - p)

    # Standard errors
    se = np.sqrt(np.diag(xtx_inv) * sigma2)
    return beta, se


def _logistic_regression(
    design: Any,
    y: Any,
    max_iter: int = 50,
) -> tuple[Any, Any] | str:
    """Iteratively reweighted least squares logistic regression.

    Args:
        design: Design matrix (n x p) of full column rank, including intercept.
        y: Binary response vector (n,).
        max_iter: Maximum IRLS iterations.

    Returns:
        ``(coefficients, standard_errors)``, or the reason why the fit is not usable (no convergence, which is
        what complete or quasi-complete separation and a constant response look like).
    """
    n, p = design.shape
    beta = np.zeros(p, dtype=np.float64)
    converged = False

    for _iteration in range(max_iter):
        eta = design @ beta
        # Clip to avoid overflow
        eta = np.clip(eta, -_MAX_LINEAR_PREDICTOR, _MAX_LINEAR_PREDICTOR)
        mu = 1.0 / (1.0 + np.exp(-eta))

        # Weight matrix diagonal
        w = mu * (1.0 - mu)
        w = np.maximum(w, 1e-10)

        # Weighted least squares update
        z = eta + (y - mu) / w
        xtw = design.T * w
        try:
            beta_new = np.linalg.solve(xtw @ design, xtw @ z)
        except np.linalg.LinAlgError:
            return "information matrix is singular (separation or a constant phenotype)"

        # Check convergence
        converged = bool(np.all(np.isfinite(beta_new)) and np.max(np.abs(beta_new - beta)) < 1e-6)
        beta = beta_new
        if converged:
            break

    if not converged:
        return (
            f"logistic regression did not converge in {max_iter} iterations "
            "(complete or quasi-complete separation, or a constant phenotype)"
        )

    # Under separation the linear predictor runs off to infinity: the clip above then stalls the iteration at a
    # finite value that looks converged, but the fitted probabilities are numerically 0 or 1 and the Wald test is
    # meaningless (the standard error explodes with the coefficient).
    if np.any(np.abs(design @ beta) >= _MAX_LINEAR_PREDICTOR):
        return "fitted probabilities are numerically 0 or 1 (complete or quasi-complete separation)"

    # Standard errors from final information matrix
    eta = design @ beta
    mu = 1.0 / (1.0 + np.exp(-eta))
    w = np.maximum(mu * (1.0 - mu), 1e-10)

    try:
        cov_matrix = np.linalg.inv((design.T * w) @ design)
    except np.linalg.LinAlgError:
        return "information matrix is singular (separation or a constant phenotype)"
    se = np.sqrt(np.maximum(np.diag(cov_matrix), 0.0))
    if not np.all(se > 0.0):
        return "standard error of a coefficient is zero or undefined"
    return beta, se


def _should_merge(
    call_a: dict[str, Any],
    call_b: dict[str, Any],
    min_overlap: float,
    max_bp_dist: int,
) -> bool:
    """Determine if two SV calls should be merged.

    Checks chromosome, SV type, reciprocal overlap, and breakpoint distance; two
    translocations are compared on both breakpoints instead (see
    :func:`merge_sv_callsets`).

    Args:
        call_a: First SV call dict.
        call_b: Second SV call dict.
        min_overlap: Minimum reciprocal overlap fraction.
        max_bp_dist: Maximum breakpoint distance.

    Returns:
        True if calls should be merged.
    """
    # A translocation is compared on both breakpoints, never with an interval variant
    translocation_a, translocation_b = is_translocation(call_a), is_translocation(call_b)
    if translocation_a or translocation_b:
        return (
            translocation_a
            and translocation_b
            and junctions_agree(canonical_junction(call_a), canonical_junction(call_b), max_bp_dist)
        )

    # Must be same chromosome and type
    if call_a.get("chrom", "") != call_b.get("chrom", ""):
        return False
    if call_a.get("sv_type", "") != call_b.get("sv_type", ""):
        return False

    start_a: int = call_a.get("start", 0)
    end_a: int = call_a.get("end", 0)
    start_b: int = call_b.get("start", 0)
    end_b: int = call_b.get("end", 0)

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
        Consensus SV call dict. A cluster of translocations has ``end == start`` and
        carries its second breakpoint as ``chrom2``/``pos2`` (``''`` and None for a
        single breakend).
    """
    # True median breakpoints. Positions are integer coordinates, so the half-integer
    # median of an even-sized cluster is truncated, as in filtering.merge._median_position.
    translocation = is_translocation(cluster[0])
    chrom = cluster[0].get("chrom", "")
    second = None
    if translocation:
        (chrom, consensus_start), second = consensus_junction(
            [canonical_junction(c) for c in cluster], lambda positions: int(statistics.median(positions))
        )
        consensus_end = consensus_start
    else:
        consensus_start = int(statistics.median([c.get("start", 0) for c in cluster]))
        consensus_end = int(statistics.median([c.get("end", 0) for c in cluster]))

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

    consensus: dict[str, Any] = {
        "chrom": chrom,
        "start": consensus_start,
        "end": consensus_end,
        "sv_type": consensus_type,
        "n_samples": len(samples),
        "samples": samples,
        "support_count": len(cluster),
        "merged_calls": cluster,
    }
    if translocation:
        consensus["chrom2"], consensus["pos2"] = ("", None) if second is None else second
    return consensus
