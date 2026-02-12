"""Partial Information Decomposition (PID) measures for biological data.

This module implements information decomposition measures including co-information,
partial information decomposition (PID), unique/redundant/synergistic information,
dual total correlation, and O-information. These quantities are essential for
understanding how multiple biological variables (e.g., gene expression, epigenetic
marks, environmental factors) jointly encode information about a target phenotype.

References:
    - Williams, P. L., & Beer, R. D. (2010). Nonnegative decomposition of
      multivariate information.
    - Barrett, A. B. (2015). Exploration of synergistic and redundant information
      sharing in static and dynamical Gaussian systems. Physical Review E.
    - Rosas, F. E., et al. (2019). Quantifying high-order interdependencies via
      multivariate extensions of the mutual information. Physical Review A.
"""

from __future__ import annotations

import math
from collections import Counter
from itertools import combinations
from typing import Any, Dict, List, Sequence

import numpy as np

from metainformant.core.data import validation
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _entropy_from_sequence(data: Sequence[Any], base: float = 2.0) -> float:
    """Compute Shannon entropy of a single discrete variable from raw observations.

    Args:
        data: Sequence of observed values.
        base: Logarithm base (2.0 for bits, math.e for nats).

    Returns:
        Shannon entropy H(X) in the requested base.
    """
    counts = Counter(data)
    total = len(data)
    probs = [c / total for c in counts.values()]
    return shannon_entropy(probs, base=base)


def _joint_entropy_n(variables: List[Sequence[Any]], base: float = 2.0) -> float:
    """Compute joint entropy H(X1, X2, ..., Xn) for N discrete variables.

    Args:
        variables: List of sequences (all must have the same length).
        base: Logarithm base.

    Returns:
        Joint entropy of all variables combined.
    """
    joint_counts = Counter(zip(*variables))
    total = len(variables[0])
    probs = [c / total for c in joint_counts.values()]
    return shannon_entropy(probs, base=base)


def _mutual_information_joint_target(sources: List[Sequence[Any]], target: Sequence[Any], base: float = 2.0) -> float:
    """Compute I(S1,S2,...,Sn ; T) -- mutual information between joint sources and target.

    Args:
        sources: List of source variable sequences.
        target: Target variable sequence.
        base: Logarithm base.

    Returns:
        Mutual information between joint sources and target.
    """
    # Create joint source variable by zipping all sources into tuples
    joint_source = list(zip(*sources))
    # I(S_joint; T) = H(S_joint) + H(T) - H(S_joint, T)
    h_s = _entropy_from_sequence(joint_source, base=base)
    h_t = _entropy_from_sequence(target, base=base)
    # Joint of (S_joint, T) -- combine the joint-source tuple with target
    combined = [(s, t) for s, t in zip(joint_source, target)]
    h_st = _entropy_from_sequence(combined, base=base)
    return max(0.0, h_s + h_t - h_st)


def _validate_variables(variables: List[Sequence[Any]], min_count: int = 2) -> None:
    """Validate a list of variable sequences for multivariate information measures.

    Args:
        variables: List of sequences to validate.
        min_count: Minimum number of variables required.

    Raises:
        TypeError: If variables is not a list.
        ValueError: If fewer than min_count variables or lengths differ.
    """
    validation.validate_type(variables, list, "variables")

    if len(variables) < min_count:
        raise ValueError(f"At least {min_count} variables required, got {len(variables)}")

    lengths = [len(v) for v in variables]
    if len(set(lengths)) != 1:
        raise ValueError(f"All variable sequences must have the same length, got lengths {lengths}")

    if lengths[0] == 0:
        raise ValueError("Variable sequences must not be empty")


def _validate_base(base: float) -> None:
    """Validate logarithm base parameter.

    Args:
        base: Logarithm base to validate.

    Raises:
        ValueError: If base is not positive or equals 1.
    """
    if base <= 0:
        raise ValueError(f"Logarithm base must be positive, got {base}")
    if base == 1.0:
        raise ValueError("Logarithm base must not be 1.0")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def co_information(variables: List[Sequence[Any]], base: float = 2.0) -> float:
    """Compute co-information (interaction information) for N variables.

    Co-information generalises mutual information to an arbitrary number of
    variables via the inclusion-exclusion principle on entropies:

        CI(X1; X2; ...; Xn) = -sum_{S subset {1,...,n}, S != empty}
                                (-1)^{|S|} H(X_S)

    For two variables this equals the mutual information I(X;Y).
    For three variables:
        CI(X;Y;Z) = H(X) + H(Y) + H(Z) - H(X,Y) - H(X,Z) - H(Y,Z) + H(X,Y,Z)
                   = I(X;Y) - I(X;Y|Z)

    Positive co-information indicates redundancy among the variables;
    negative co-information indicates synergy.

    Args:
        variables: List of at least 2 sequences (all the same length) representing
            discrete random variables observed point-wise.
        base: Logarithm base (2.0 for bits, math.e for nats).

    Returns:
        Co-information value. Positive means redundancy-dominated,
        negative means synergy-dominated.

    Raises:
        TypeError: If variables is not a list.
        ValueError: If fewer than 2 variables, sequences differ in length,
            or sequences are empty.
    """
    _validate_variables(variables, min_count=2)
    _validate_base(base)

    n = len(variables)

    # Inclusion-exclusion: CI = sum over non-empty subsets S of (-1)^(|S|+1) H(X_S)
    # Equivalently: CI = -sum over non-empty subsets S of (-1)^|S| H(X_S)
    ci = 0.0
    for size in range(1, n + 1):
        sign = (-1.0) ** (size + 1)
        for subset_indices in combinations(range(n), size):
            subset_vars = [variables[i] for i in subset_indices]
            if len(subset_vars) == 1:
                h = _entropy_from_sequence(subset_vars[0], base=base)
            else:
                h = _joint_entropy_n(subset_vars, base=base)
            ci += sign * h

    return float(ci)


def partial_information_decomposition(
    x: Sequence[Any],
    y: Sequence[Any],
    target: Sequence[Any],
    base: float = 2.0,
) -> Dict[str, float]:
    """Decompose I(X,Y;T) into redundant, unique, and synergistic components.

    Uses the minimum mutual information (MMI) approach (Barrett 2015) to
    decompose the total information that two source variables X and Y carry
    about a target variable T:

        I(X,Y;T) = Redundancy + Unique_X + Unique_Y + Synergy

    where:
        Redundancy = min(I(X;T), I(Y;T))
        Unique_X   = I(X;T) - Redundancy
        Unique_Y   = I(Y;T) - Redundancy
        Synergy    = I(X,Y;T) - I(X;T) - I(Y;T) + Redundancy

    Args:
        x: First source variable sequence.
        y: Second source variable sequence.
        target: Target variable sequence (same length as x and y).
        base: Logarithm base (2.0 for bits, math.e for nats).

    Returns:
        Dictionary with keys:
            ``redundancy`` -- shared information both sources carry about target.
            ``unique_x``   -- information only X carries about target.
            ``unique_y``   -- information only Y carries about target.
            ``synergy``    -- information available only from joint observation.
            ``total``      -- total mutual information I(X,Y;T).

    Raises:
        TypeError: If any input is not a list or tuple.
        ValueError: If sequences differ in length or are empty.
    """
    validation.validate_type(x, (list, tuple), "x")
    validation.validate_type(y, (list, tuple), "y")
    validation.validate_type(target, (list, tuple), "target")
    _validate_base(base)

    if not (len(x) == len(y) == len(target)):
        raise ValueError(f"All sequences must have the same length, got {len(x)}, {len(y)}, {len(target)}")
    if len(x) == 0:
        raise ValueError("Sequences must not be empty")

    # Individual mutual informations
    i_x_t = mutual_information(list(x), list(target), base=base)
    i_y_t = mutual_information(list(y), list(target), base=base)

    # Joint mutual information I(X,Y;T)
    i_xy_t = _mutual_information_joint_target([list(x), list(y)], list(target), base=base)

    # MMI decomposition
    redundancy = min(i_x_t, i_y_t)
    unique_x = max(0.0, i_x_t - redundancy)
    unique_y = max(0.0, i_y_t - redundancy)
    synergy = max(0.0, i_xy_t - i_x_t - i_y_t + redundancy)

    # The total should match I(X,Y;T)
    total = i_xy_t

    logger.debug(
        "PID decomposition: redundancy=%.4f, unique_x=%.4f, unique_y=%.4f, " "synergy=%.4f, total=%.4f",
        redundancy,
        unique_x,
        unique_y,
        synergy,
        total,
    )

    return {
        "redundancy": float(redundancy),
        "unique_x": float(unique_x),
        "unique_y": float(unique_y),
        "synergy": float(synergy),
        "total": float(total),
    }


def unique_information(
    source: Sequence[Any],
    target: Sequence[Any],
    other_source: Sequence[Any],
    base: float = 2.0,
) -> float:
    """Compute unique information that *source* provides about *target* beyond *other_source*.

    UI(source -> target \\ other_source) = I(source; target) - Redundancy

    where Redundancy is the MMI redundancy between source and other_source
    with respect to target.

    Args:
        source: The source whose unique contribution is measured.
        target: The target variable.
        other_source: The competing source variable.
        base: Logarithm base (2.0 for bits, math.e for nats).

    Returns:
        Unique information (non-negative) that source provides about target
        beyond what other_source also provides.

    Raises:
        TypeError: If any input is not a list or tuple.
        ValueError: If sequences differ in length or are empty.
    """
    validation.validate_type(source, (list, tuple), "source")
    validation.validate_type(target, (list, tuple), "target")
    validation.validate_type(other_source, (list, tuple), "other_source")
    _validate_base(base)

    if not (len(source) == len(target) == len(other_source)):
        raise ValueError(
            f"All sequences must have the same length, got {len(source)}, " f"{len(target)}, {len(other_source)}"
        )
    if len(source) == 0:
        raise ValueError("Sequences must not be empty")

    i_source_target = mutual_information(list(source), list(target), base=base)
    i_other_target = mutual_information(list(other_source), list(target), base=base)

    # MMI redundancy
    red = min(i_source_target, i_other_target)

    return max(0.0, float(i_source_target - red))


def redundant_information(
    sources: List[Sequence[Any]],
    target: Sequence[Any],
    base: float = 2.0,
) -> float:
    """Compute redundant information that all sources share about target.

    Uses the MMI (minimum mutual information) measure:

        Redundancy = min_i I(S_i ; T)

    This captures the common information floor -- the minimum amount of
    information that every single source individually provides about the target.

    Args:
        sources: List of at least 2 source variable sequences.
        target: Target variable sequence (same length as each source).
        base: Logarithm base (2.0 for bits, math.e for nats).

    Returns:
        Redundant information (non-negative).

    Raises:
        TypeError: If sources is not a list or target is not a list/tuple.
        ValueError: If fewer than 2 sources, lengths differ, or sequences are empty.
    """
    validation.validate_type(sources, list, "sources")
    validation.validate_type(target, (list, tuple), "target")
    _validate_base(base)

    if len(sources) < 2:
        raise ValueError(f"At least 2 sources required, got {len(sources)}")

    target_len = len(target)
    for i, src in enumerate(sources):
        if len(src) != target_len:
            raise ValueError(f"Source {i} length ({len(src)}) differs from target length ({target_len})")
    if target_len == 0:
        raise ValueError("Sequences must not be empty")

    # I(Si; T) for each source
    mis = [mutual_information(list(src), list(target), base=base) for src in sources]

    return float(min(mis))


def synergistic_information(
    sources: List[Sequence[Any]],
    target: Sequence[Any],
    base: float = 2.0,
) -> float:
    """Compute synergistic information about target available only from joint sources.

    Synergy is the information about target that can only be obtained by
    observing all sources jointly -- it is not present in any individual source.

        Synergy = I(S1, S2, ..., Sn ; T) - sum_i Unique_i - Redundancy

    For the 2-source case with MMI:
        Synergy = I(S1, S2; T) - I(S1; T) - I(S2; T) + Redundancy

    For the N-source case the generalised formula is:
        Synergy = I(S_joint; T) - sum_i I(Si; T) + (n-1) * Redundancy

    This arises because each I(Si; T) = Unique_i + Redundancy, so
    sum_i I(Si; T) = sum_i Unique_i + n * Redundancy.  Rearranging the
    PID identity I_joint = Synergy + sum_i Unique_i + Redundancy gives
    the formula above.

    Args:
        sources: List of at least 2 source variable sequences.
        target: Target variable sequence (same length as each source).
        base: Logarithm base (2.0 for bits, math.e for nats).

    Returns:
        Synergistic information (non-negative under the MMI decomposition).

    Raises:
        TypeError: If sources is not a list or target is not a list/tuple.
        ValueError: If fewer than 2 sources, lengths differ, or sequences are empty.
    """
    validation.validate_type(sources, list, "sources")
    validation.validate_type(target, (list, tuple), "target")
    _validate_base(base)

    if len(sources) < 2:
        raise ValueError(f"At least 2 sources required, got {len(sources)}")

    target_len = len(target)
    for i, src in enumerate(sources):
        if len(src) != target_len:
            raise ValueError(f"Source {i} length ({len(src)}) differs from target length ({target_len})")
    if target_len == 0:
        raise ValueError("Sequences must not be empty")

    n = len(sources)
    sources_list = [list(s) for s in sources]

    # Joint MI: I(S_joint; T)
    i_joint_t = _mutual_information_joint_target(sources_list, list(target), base=base)

    # Individual MIs
    mis = [mutual_information(src, list(target), base=base) for src in sources_list]

    # MMI Redundancy
    red = min(mis)

    # Synergy = I_joint - sum(I_i) + (n-1)*Redundancy
    synergy = i_joint_t - sum(mis) + (n - 1) * red

    return max(0.0, float(synergy))


def dual_total_correlation(variables: List[Sequence[Any]], base: float = 2.0) -> float:
    """Compute the dual total correlation (binding information) for N variables.

    The dual total correlation (Abdallah & Plumbley 2012, Han 1978) measures
    the shared or redundant multivariate dependence:

        DTC(X1, ..., Xn) = H(X1, ..., Xn) - sum_i H(Xi | X_{\\i})

    where X_{\\i} denotes all variables except Xi.

    Equivalently:
        DTC = H(X_joint) - sum_i [H(X_joint) - H(X_{\\i})]
            = (1 - n) * H(X_joint) + sum_i H(X_{\\i})

    DTC is non-negative and equals zero if and only if at least one variable
    is independent of the rest.

    Args:
        variables: List of at least 2 sequences (all the same length).
        base: Logarithm base (2.0 for bits, math.e for nats).

    Returns:
        Dual total correlation (non-negative).

    Raises:
        TypeError: If variables is not a list.
        ValueError: If fewer than 2 variables, lengths differ, or sequences empty.
    """
    _validate_variables(variables, min_count=2)
    _validate_base(base)

    n = len(variables)

    # H(X1, ..., Xn)
    h_joint = _joint_entropy_n(variables, base=base)

    # For each i, compute H(Xi | X_{\i}) = H(X_joint) - H(X_{\i})
    # So DTC = H(X_joint) - sum_i [H(X_joint) - H(X_{\i})]
    #        = H(X_joint) - n*H(X_joint) + sum_i H(X_{\i})
    #        = (1-n)*H(X_joint) + sum_i H(X_{\i})
    sum_h_others = 0.0
    for i in range(n):
        others = [variables[j] for j in range(n) if j != i]
        if len(others) == 1:
            h_others = _entropy_from_sequence(others[0], base=base)
        else:
            h_others = _joint_entropy_n(others, base=base)
        sum_h_others += h_others

    dtc = (1 - n) * h_joint + sum_h_others

    return max(0.0, float(dtc))


def o_information(variables: List[Sequence[Any]], base: float = 2.0) -> float:
    """Compute the O-information for N variables (Rosas et al. 2019).

    O-information quantifies the balance between redundancy and synergy
    in a multivariate system:

        O(X1, ..., Xn) = TC(X1, ..., Xn) - DTC(X1, ..., Xn)

    where:
        TC  = sum_i H(Xi)  - H(X1, ..., Xn)          (total correlation)
        DTC = H(X1, ..., Xn) - sum_i H(Xi | X_{\\i})  (dual total correlation)

    Expanding:
        O = sum_i H(Xi) - H(X_joint) - H(X_joint) + sum_i H(Xi | X_{\\i})

    But H(Xi | X_{\\i}) = H(X_joint) - H(X_{\\i}), so:
        O = sum_i H(Xi) - 2*H(X_joint) + n*H(X_joint) - sum_i H(X_{\\i})

    Simplifying:
        O = (n - 2)*H(X_joint) + sum_i H(Xi) - sum_i H(X_{\\i})

    Interpretation:
        O > 0 : redundancy-dominated system
        O < 0 : synergy-dominated system
        O = 0 : balanced

    Args:
        variables: List of at least 2 sequences (all the same length).
        base: Logarithm base (2.0 for bits, math.e for nats).

    Returns:
        O-information value (positive = redundancy, negative = synergy).

    Raises:
        TypeError: If variables is not a list.
        ValueError: If fewer than 2 variables, lengths differ, or sequences empty.
    """
    _validate_variables(variables, min_count=2)
    _validate_base(base)

    n = len(variables)

    # H(X_joint)
    h_joint = _joint_entropy_n(variables, base=base)

    # sum_i H(Xi)
    sum_h_individual = sum(_entropy_from_sequence(variables[i], base=base) for i in range(n))

    # sum_i H(X_{\i})
    sum_h_others = 0.0
    for i in range(n):
        others = [variables[j] for j in range(n) if j != i]
        if len(others) == 1:
            h_others = _entropy_from_sequence(others[0], base=base)
        else:
            h_others = _joint_entropy_n(others, base=base)
        sum_h_others += h_others

    # O = (n - 2) * H(X_joint) + sum H(Xi) - sum H(X_{\i})
    o_val = (n - 2) * h_joint + sum_h_individual - sum_h_others

    return float(o_val)
