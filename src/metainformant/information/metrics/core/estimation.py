"""Bias-corrected entropy and information estimation methods.

This module implements various estimation methods for information-theoretic
quantities with bias correction to improve accuracy on finite samples.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Any, Dict, List, Optional, Union

import numpy as np

from metainformant.core.data import validation
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def entropy_estimator(
    counts: "Union[Dict[Any, int], List[int], np.ndarray]",
    method: str = "plugin",
    bias_correction: bool = True,
) -> float:
    """Estimate Shannon entropy with various methods and bias correction.

    Args:
        counts: Either dict mapping items to counts or list of counts
        method: Estimation method ('plugin', 'miller_madow', 'chao_shen', 'jackknife')
        bias_correction: Whether to apply bias correction

    Returns:
        Entropy estimate in bits

    Raises:
        ValueError: If invalid method or insufficient data
    """
    if isinstance(counts, dict):
        count_values = list(counts.values())
    else:
        count_values = list(counts)

    validation.validate_type(count_values, list, "counts")

    count_array = np.array(count_values, dtype=int)
    if np.any(count_array < 0):
        raise ValueError("Counts cannot be negative")

    total = np.sum(count_array)
    if total == 0:
        return 0.0

    if method == "plugin":
        return _plugin_entropy_estimator(count_array, total, bias_correction)
    elif method == "miller_madow":
        return _miller_madow_entropy_estimator(count_array, total)
    elif method == "chao_shen":
        return _chao_shen_entropy_estimator(count_array, total)
    elif method == "jackknife":
        return _jackknife_entropy_estimator(count_array, total)
    else:
        raise ValueError(f"Unknown entropy estimation method: {method}")


def _plugin_entropy_estimator(
    counts: np.ndarray, total: int, bias_correction: bool
) -> float:
    """Plugin (maximum-likelihood) entropy estimator, in bits."""
    # Convert to probabilities
    probs = counts / total
    probs = probs[probs > 0]  # Remove zeros

    if len(probs) == 0:
        return 0.0

    # Plugin entropy: -sum(p * log2(p))
    entropy = -np.sum(probs * np.log2(probs))

    # Bias correction (Miller, 1955): the plugin estimator is biased LOW by
    # approximately (k-1)/(2n) nats, i.e. (k-1)/(2n * ln 2) bits; add it.
    if bias_correction and total > 1:
        k = len(probs)  # Number of non-zero categories
        correction = (k - 1) / (2 * total * math.log(2))
        entropy += correction

    return max(0.0, float(entropy))


def _miller_madow_entropy_estimator(counts: np.ndarray, total: int) -> float:
    """Miller-Madow entropy estimator (bias-corrected plugin)."""
    # Convert to probabilities
    probs = counts / total
    probs = probs[probs > 0]  # Remove zeros

    if len(probs) == 0:
        return 0.0

    # Miller-Madow correction: the plugin estimator is biased LOW by
    # approximately (k-1)/(2n) nats where k is the number of non-zero
    # categories; in bits the additive correction is (k-1)/(2n * ln 2).
    k = len(probs)
    correction = (k - 1) / (2 * total * math.log(2)) if total > 1 else 0

    entropy = -np.sum(probs * np.log2(probs)) + correction

    return max(0.0, float(entropy))


def _chao_shen_entropy_estimator(counts: np.ndarray, total: int) -> float:
    """Chao-Shen (2003) coverage-adjusted Horvitz-Thompson entropy estimator (bits).

    Implements the published estimator (Chao & Shen 2003; identical to
    ``entropy.ChaoShen`` in the R `entropy` package):

      1. sample coverage  ``C = 1 - f1 / n`` (f1 = number of singleton counts)
      2. adjusted probs   ``pa_i = C * p_i`` (p_i = observed frequencies)
      3. inclusion probs  ``la_i = 1 - (1 - pa_i) ** n`` (Horvitz-Thompson)
      4. ``H = -sum_i (pa_i / la_i) * log2(pa_i)``

    When every observed category is a singleton the coverage estimate
    ``C <= 0`` and the estimator is undefined; this implementation then
    returns 0.0.

    Args:
        counts: Count per category (zeros ignored)
        total: Total number of observations n

    Returns:
        Chao-Shen entropy estimate in bits
    """
    positive = counts[counts > 0]
    if len(positive) == 0 or total <= 0:
        return 0.0

    f1 = int(np.sum(positive == 1))
    coverage = 1.0 - f1 / total
    if coverage <= 0.0:
        return 0.0

    probs = positive / total
    adjusted = coverage * probs
    inclusion = 1.0 - (1.0 - adjusted) ** total

    h = -np.sum(adjusted * np.log2(adjusted) / inclusion)
    return max(0.0, float(h))


def _jackknife_entropy_estimator(counts: np.ndarray, total: int) -> float:
    """Jackknife (Zahl 1977 category jackknife) entropy estimator, in bits."""
    counts = counts[counts > 0]  # Remove zeros
    k = len(counts)

    if k <= 1:
        return 0.0

    # Plugin entropy
    h_plugin = entropy_estimator(counts, method="plugin", bias_correction=False)

    # Jackknife: remove one category at a time
    h_jackknife_terms = []

    for i in range(k):
        # Remove i-th category
        reduced_counts = np.delete(counts, i)
        reduced_total = total - counts[i]

        if reduced_total > 0 and len(reduced_counts) > 0:
            reduced_probs = reduced_counts / reduced_total
            h_reduced = -np.sum(reduced_probs * np.log2(reduced_probs))
            h_jackknife_terms.append(h_reduced)

    if not h_jackknife_terms:
        return h_plugin

    # Jackknife estimate: k * h_plugin - ((k-1)/k) * sum(h_{-i})
    h_avg_reduced = np.mean(h_jackknife_terms)
    h_jackknife = k * h_plugin - (k - 1) * h_avg_reduced

    return float(max(0.0, h_jackknife))


def mutual_information_estimator(
    x: List[Any], y: List[Any], method: str = "plugin", bias_correction: bool = True
) -> float:
    """Estimate mutual information I(X;Y) with first-order bias correction, in bits.

    The plugin estimate ``I = H(X) + H(Y) - H(X,Y)`` is biased LOW by
    (Miller 1955; Panzeri & Treves 1996, first order)::

        E[I_plugin] - I = -(|X|*|Y| - |X| - |Y| + 1) / (2n)  nats

    where |X| and |Y| are the numbers of POSSIBLE marginal states, so the
    possible joint alphabet has ``|X|*|Y|`` states. The correction is
    therefore computed from ``k_x * k_y`` (the full possible joint alphabet),
    NOT from the number of jointly occupied states, and added in bits. At
    small n this first-order correction can over-shoot the true MI for
    strongly dependent variables; the result is clipped at 0 because
    MI >= 0 by definition (the clip does not alter the correction itself).

    Args:
        x: Samples from first variable
        y: Samples from second variable
        method: Estimation method ('plugin', 'miller_madow'); both apply the
            identical joint-alphabet first-order MI correction. Other methods
            (e.g. 'chao_shen', 'jackknife') apply their own internal
            per-entropy corrections.
        bias_correction: Whether to apply the first-order MI bias correction
            (plugin method only; 'miller_madow' always corrects)

    Returns:
        Mutual information estimate in bits (>= 0)

    Raises:
        ValueError: If sequences have different lengths
    """
    validation.validate_type(x, list, "x")
    validation.validate_type(y, list, "y")

    if len(x) != len(y):
        raise ValueError("Sequences must have the same length")

    n = len(x)
    if n == 0:
        return 0.0

    # Joint counts
    joint_counts = Counter(zip(x, y))

    # Marginal counts
    x_counts = Counter(x)
    y_counts = Counter(y)

    if method in ("plugin", "miller_madow"):
        # Raw (uncorrected) plugin entropies in bits...
        h_x = _plugin_entropy_estimator(
            np.array(list(x_counts.values()), dtype=int), n, False
        )
        h_y = _plugin_entropy_estimator(
            np.array(list(y_counts.values()), dtype=int), n, False
        )
        h_xy = _plugin_entropy_estimator(
            np.array(list(joint_counts.values()), dtype=int), n, False
        )

        mi = h_x + h_y - h_xy

        # ...plus the exact first-order MI bias over the POSSIBLE joint
        # alphabet |X|*|Y|: (k_x - 1)(k_y - 1) / (2n) nats, converted to bits.
        if method == "miller_madow" or bias_correction:
            k_x = len(x_counts)
            k_y = len(y_counts)
            mi += (k_x * k_y - k_x - k_y + 1) / (2 * n * math.log(2))
    else:
        # Other estimators (e.g. chao_shen, jackknife) correct each entropy
        # term internally.
        h_x = entropy_estimator(
            x_counts, method=method, bias_correction=bias_correction
        )
        h_y = entropy_estimator(
            y_counts, method=method, bias_correction=bias_correction
        )
        h_xy = entropy_estimator(
            joint_counts, method=method, bias_correction=bias_correction
        )
        mi = h_x + h_y - h_xy

    return max(0.0, mi)


def kl_divergence_estimator(
    p: List[Any], q: List[Any], method: str = "plugin", bias_correction: bool = True
) -> float:
    """Estimate KL divergence D_KL(P||Q) with bias correction, in bits.

    Args:
        p: Samples from distribution P
        q: Samples from distribution Q
        method: Estimation method
        bias_correction: Whether to apply bias correction

    Returns:
        KL divergence estimate in bits (inf if q assigns zero probability
        to an outcome observed under p)

    Raises:
        ValueError: If sample lists have different lengths
    """
    if len(p) != len(q):
        raise ValueError("Sample lists must have the same length")

    # Convert to probability distributions using counts
    p_counts = Counter(p)
    q_counts = Counter(q)

    total = len(p)

    # Convert to probabilities
    p_probs = {item: count / total for item, count in p_counts.items()}
    q_probs = {item: count / total for item, count in q_counts.items()}

    # All possible items
    all_items = set(p_probs.keys()) | set(q_probs.keys())

    # KL divergence: sum(p * log(p/q))
    kl_div = 0.0
    for item in all_items:
        p_prob = p_probs.get(item, 0.0)
        q_prob = q_probs.get(item, 0.0)

        if p_prob > 0:
            if q_prob > 0:
                kl_div += p_prob * math.log2(p_prob / q_prob)
            else:
                return float("inf")  # Infinite divergence

    return max(0.0, kl_div)


def bias_correction(entropy: float, sample_size: int, alphabet_size: int) -> float:
    """Add the first-order (Miller-Madow) entropy bias correction, in bits.

    The plugin Shannon entropy estimator is biased LOW by approximately
    ``(alphabet_size - 1) / (2 * sample_size)`` nats (Miller 1955). This
    helper adds the bit-valued equivalent
    ``(alphabet_size - 1) / (2 * sample_size * ln 2)`` to the supplied
    estimate. :func:`effective_sample_size_correction` is a
    backward-compatible alias of this function.

    The correction is unreliable for sample_size <= 1 (no distribution can
    be estimated from a single draw); the entropy is returned unchanged
    there, consistent with the plugin estimator which skips the correction
    when n <= 1.

    Args:
        entropy: Raw entropy estimate (bits)
        sample_size: Number of samples (n)
        alphabet_size: Size of alphabet (d)

    Returns:
        Bias-corrected entropy estimate (bits)

    Raises:
        ValueError: If parameters are invalid
    """
    if sample_size <= 0:
        raise ValueError("Sample size must be positive")
    if alphabet_size <= 0:
        raise ValueError("Alphabet size must be positive")
    if sample_size <= 1:
        return entropy

    correction = (alphabet_size - 1) / (2 * sample_size * math.log(2))

    return max(0.0, entropy + correction)


def entropy_bootstrap_confidence(
    counts: Union[Dict[Any, int], List[int]],
    method: str = "plugin",
    n_bootstraps: int = 1000,
    confidence_level: float = 0.95,
    random_state: Optional[int] = None,
) -> Dict[str, float]:
    """Calculate bootstrap confidence interval for entropy estimate.

    Args:
        counts: Count data
        method: Estimation method
        n_bootstraps: Number of bootstrap samples
        confidence_level: Confidence level (0-1)
        random_state: Random state for reproducibility

    Returns:
        Dictionary with entropy estimate and confidence interval

    Raises:
        ValueError: If invalid parameters
    """
    if not (0 < confidence_level < 1):
        raise ValueError("Confidence level must be between 0 and 1")

    if isinstance(counts, dict):
        items = []
        for item, count in counts.items():
            items.extend([item] * count)
    else:
        items = []
        for i, count in enumerate(counts):
            items.extend([i] * count)

    if not items:
        return {
            "entropy": 0.0,
            "ci_lower": 0.0,
            "ci_upper": 0.0,
            "confidence_level": confidence_level,
        }

    np.random.seed(random_state)
    bootstrap_entropies = []

    # Generate bootstrap samples
    for _ in range(n_bootstraps):
        # Bootstrap resampling
        bootstrap_sample = np.random.choice(items, size=len(items), replace=True)

        # Convert back to counts
        bootstrap_counts = Counter(bootstrap_sample)

        # Estimate entropy
        entropy_est = entropy_estimator(
            bootstrap_counts, method=method, bias_correction=True
        )
        bootstrap_entropies.append(entropy_est)

    bootstrap_entropies_arr = np.array(bootstrap_entropies)

    # Calculate confidence interval
    alpha = 1 - confidence_level
    ci_lower = np.percentile(bootstrap_entropies_arr, alpha / 2 * 100)
    ci_upper = np.percentile(bootstrap_entropies_arr, (1 - alpha / 2) * 100)

    # Main estimate (using original data)
    if isinstance(counts, dict):
        main_counts = counts
    else:
        main_counts = Counter(items)

    main_entropy = entropy_estimator(main_counts, method=method, bias_correction=True)

    return {
        "entropy": main_entropy,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
        "confidence_level": confidence_level,
        "n_bootstraps": n_bootstraps,
    }


def effective_sample_size_correction(
    entropy: float, sample_size: int, alphabet_size: int
) -> float:
    """Backward-compatible alias of :func:`bias_correction`.

    Applies the identical additive Miller-Madow first-order correction
    ``(alphabet_size - 1) / (2 * sample_size * ln 2)`` bits (sample_size <= 1
    returns the entropy unchanged). The two functions used to carry
    duplicated implementations; the logic now lives only in
    :func:`bias_correction`. New code should call that function.

    Args:
        entropy: Raw entropy estimate (bits)
        sample_size: Number of samples (n)
        alphabet_size: Alphabet size (d)

    Returns:
        Bias-corrected entropy estimate (bits)
    """
    return bias_correction(entropy, sample_size, alphabet_size)


def _panzeri_treves_support(probabilities: np.ndarray, sample_size: int) -> float:
    """Panzeri-Treves (1996) Bayesian estimate of the response support R.

    ``probabilities`` is the probability vector over the FULL alphabet
    (unobserved responses carry probability 0). Starting from the number of
    occupied responses, the estimate grows the support one extra ("quasi")
    response at a time: occupied responses receive quasi-Bayes add-one
    probabilities over ``n + R`` pseudo-counts, the extra responses share the
    remaining weight ``gamma = x * (1 - (n / (n + R)) ** (1/n))``, and the
    support is grown as long as the expected number of responses seen at
    least once in ``n`` draws moves closer to the observed occupied count.

    Faithful to ``pt_bayescount`` in pyentropy (Ince et al. 2009), the
    reference implementation of the Panzeri-Treves correction.
    """
    eps = np.finfo(float).eps
    non_zero = probabilities[probabilities > eps]
    r_naive = non_zero.size

    if r_naive >= probabilities.size:
        return float(r_naive)

    r_expected = r_naive - float(((1.0 - non_zero) ** sample_size).sum())
    delta_prev = float(probabilities.size)
    delta = abs(r_naive - r_expected)
    extra = 0.0
    while (delta < delta_prev) and (r_naive + extra) < probabilities.size:
        extra += 1.0
        # Occupied responses: quasi-Bayes (add-one) probabilities.
        gamma = extra * (
            1.0 - (sample_size / (sample_size + r_naive)) ** (1.0 / sample_size)
        )
        p_bayes = ((1.0 - gamma) / (sample_size + r_naive)) * (
            non_zero * sample_size + 1.0
        )
        r_expected = float((1.0 - (1.0 - p_bayes) ** sample_size).sum())
        # The extra, so-far-unobserved quasi-responses.
        p_bayes = gamma / extra
        r_expected += extra * (1.0 - (1.0 - p_bayes) ** sample_size)
        delta_prev = delta
        delta = abs(r_naive - r_expected)

    r_naive = r_naive + extra - 1.0
    if delta < delta_prev:
        r_naive += 1.0
    return float(r_naive)


def panzeri_treves_bias_correction(
    entropy: float,
    sample_size: int,
    alphabet_size: int,
    response_frequencies: Optional[np.ndarray] = None,
) -> float:
    """Apply the Panzeri-Treves (1996) sampling-bias correction, in bits.

    The plugin entropy is biased LOW because responses that exist in the
    alphabet but go unobserved (and singletons seen only once) contribute
    entropy that a finite sample misses. PT (1996) estimate this
    analytically: the effective support ``R`` of the response space
    (occupied responses plus a Bayesian estimate of the number of unobserved
    "quasi-sample" responses, see :func:`_panzeri_treves_support`) replaces
    the occupied count in the first-order bias term, giving::

        H_PT = H_plugin + (R - 1) / (2 * n * ln 2)   bits

    The correction is ADDED (the plugin estimate underestimates the true
    entropy). When every response of the alphabet is observed
    (``R = alphabet_size``) this reduces exactly to the Miller-Madow
    correction of :func:`bias_correction`.

    Args:
        entropy: Raw (plugin) entropy estimate (bits)
        sample_size: Sample size n (number of draws)
        alphabet_size: Size of the full response alphabet (may exceed the
            number of observed responses)
        response_frequencies: Observed response counts (non-negative,
            summing to sample_size, at most alphabet_size entries). If None,
            a uniform response distribution over the full alphabet is assumed.

    Returns:
        Panzeri-Treves bias-corrected entropy estimate (bits)

    Raises:
        ValueError: If parameters are invalid

    References:
        Panzeri & Treves (1996). Analytical estimates of limited sampling
        biases in different information measures. Network 7, 87-107.
    """
    if sample_size <= 1:
        return entropy
    if sample_size <= 0:
        raise ValueError("Sample size must be positive")
    if alphabet_size <= 0:
        raise ValueError("Alphabet size must be positive")

    if response_frequencies is None:
        # Uniform response fallback: identical counts over the full alphabet.
        freq_array = np.full(alphabet_size, sample_size / alphabet_size, dtype=float)
    else:
        freq_array = np.asarray(response_frequencies, dtype=float)
        if freq_array.ndim != 1:
            raise ValueError("response_frequencies must be a 1D array of counts")
        if freq_array.size > alphabet_size:
            raise ValueError("response_frequencies cannot exceed the alphabet size")
        if np.any(freq_array < 0):
            raise ValueError("Response counts cannot be negative")
        if not math.isclose(
            float(freq_array.sum()), float(sample_size), rel_tol=1e-9, abs_tol=1e-9
        ):
            raise ValueError("Response counts must sum to the sample size")

    probabilities = np.zeros(alphabet_size, dtype=float)
    probabilities[: freq_array.size] = freq_array / sample_size

    support = _panzeri_treves_support(probabilities, sample_size)
    correction = (support - 1.0) / (2.0 * sample_size * math.log(2))

    return max(0.0, entropy + correction)


def entropy_rate_estimator(
    sequence: List[Any], order: int = 1, method: str = "plugin"
) -> float:
    """Estimate entropy rate of a sequence.

    The entropy rate is the limit of n-block entropy divided by n as n→∞.
    For Markov chains, this equals the conditional entropy H(X_{n+1}|X_n).
    Args:
        sequence: Input sequence
        order: Markov order (1 for first-order Markov)
        method: Entropy estimation method

    Returns:
        Entropy rate estimate in bits (>= 0)

    Raises:
        ValueError: If sequence is too short or order is invalid
    """
    validation.validate_type(sequence, list, "sequence")

    n = len(sequence)
    if n <= order + 1:
        raise ValueError(f"Sequence too short for order {order}")

    if order < 1:
        raise ValueError("Order must be >= 1")

    # For entropy rate, we estimate H(X_{n+1} | X_1^n)
    # Using the chain rule: H(X_1^{n+1}) - H(X_1^n)

    # Build n+1 blocks
    blocks_n1 = [tuple(sequence[i : i + order + 1]) for i in range(n - order)]
    # Build n blocks
    blocks_n = [tuple(sequence[i : i + order]) for i in range(n - order + 1)]

    # Estimate entropies
    h_n1 = entropy_estimator(Counter(blocks_n1), method=method)
    h_n = entropy_estimator(Counter(blocks_n), method=method)

    # Entropy rate: H(X_{n+1}|X_1^n) = H(X_1^{n+1}) - H(X_1^n)
    entropy_rate = h_n1 - h_n

    return max(0.0, entropy_rate)
