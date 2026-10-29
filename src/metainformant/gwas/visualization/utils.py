"""Shared utility functions for GWAS visualization.

Provides common helpers used across multiple visualization modules,
extracted to eliminate code duplication. Includes p-value key detection,
safe log10 transforms, chromosome handling (key aliases, labels, natural
ordering and cumulative genome-wide offsets) and the standard spellings of
the chromosome, position and p-value columns of a results table.
"""

from __future__ import annotations

import math
import numbers
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

# Floor value for p-value clamping to avoid log(0)
_P_VALUE_FLOOR: float = 1e-300

#: Result keys that may hold the chromosome, in priority order. This is the same tuple the interactive
#: Manhattan plot uses (``interactive.interactive._CHROMOSOME_KEYS``; a test pins the two together) so
#: every plot reads the same records the same way.
CHROMOSOME_KEYS: Tuple[str, ...] = ("chromosome", "chrom", "chr", "CHROM", "CHR")

#: Result keys that may hold the base-pair position, in priority order: the interactive Manhattan's keys
#: (``interactive.interactive._POSITION_KEYS``) plus PLINK's ``BP``.
POSITION_KEYS: Tuple[str, ...] = ("position", "pos", "POS", "BP")

#: Result keys that may hold the p-value, in priority order: the interactive Manhattan's keys
#: (``interactive.interactive._P_VALUE_KEYS``) plus ``pval``. A bare lowercase ``p`` is left out on purpose:
#: it is just as often an allele frequency.
P_VALUE_KEYS: Tuple[str, ...] = ("p_value", "pvalue", "pval", "P")

#: Result keys that may hold the effect size, in priority order.
BETA_KEYS: Tuple[str, ...] = ("beta", "BETA")

#: Each standard column of a results table (see :func:`standardize_result_columns`), with the name of what it
#: holds and the spellings that may stand for it.
_STANDARD_COLUMNS: Dict[str, Tuple[str, Tuple[str, ...]]] = {
    "CHR": ("chromosome", CHROMOSOME_KEYS),
    "BP": ("position", POSITION_KEYS),
    "P": ("p-value", P_VALUE_KEYS),
}


def detect_p_value_key(result_dict: Dict[str, Any]) -> Optional[str]:
    """Detect which key holds the p-value in a GWAS result dictionary.

    Checks for common p-value key names in priority order:
    p_value > pvalue > pval > P > p

    Args:
        result_dict: A single GWAS result dictionary.

    Returns:
        The key name holding the p-value, or None if no recognized key is found.
    """
    # Priority order matches the most common conventions across the codebase:
    # metainformant internal format, then PLINK/BOLT-LMM format, then minimal
    for key in ("p_value", "pvalue", "pval", "P", "p"):
        if key in result_dict:
            return key
    return None


def safe_log10_p(p_value: float) -> float:
    """Compute -log10(p) with clamping for zero and tiny values.

    Clamps p_value to [1e-300, +inf) before applying -log10, preventing
    -inf or math domain errors. Negative p-values are also clamped to
    the floor.

    Args:
        p_value: The p-value to transform. Should be in (0, 1] but
            zero, negative, and very small values are handled safely.

    Returns:
        -log10(clamped_p) as a float. Always finite and >= 0 for
        valid p-values.
    """
    clamped = max(p_value, _P_VALUE_FLOOR)
    return float(-math.log10(clamped))


def _first_present(record: Mapping[str, Any], keys: Sequence[str], default: Any) -> Any:
    """Value of the first of ``keys`` that ``record`` has with a value other than ``None``."""
    for key in keys:
        value = record.get(key)
        if value is not None:
            return value
    return default


def get_chromosome(record: Mapping[str, Any], default: Any = None) -> Any:
    """Return the chromosome of a result record under any of the :data:`CHROMOSOME_KEYS` spellings.

    The first key whose value is not ``None`` wins; ``default`` is returned when none is present.

    Args:
        record: A single GWAS result dictionary.
        default: Value to return when the record has no chromosome.
    """
    return _first_present(record, CHROMOSOME_KEYS, default)


def get_position(record: Mapping[str, Any], default: Any = None) -> Any:
    """Return the base-pair position of a result record under any of the :data:`POSITION_KEYS` spellings.

    The first key whose value is not ``None`` wins; ``default`` is returned when none is present.
    """
    return _first_present(record, POSITION_KEYS, default)


def finite_float(value: Any) -> float:
    """``value`` as a float, or NaN when it is missing (None), not a number (such as ``"NA"``) or not finite.

    A boolean is not a number here.
    """
    if value is None or isinstance(value, bool):
        return float("nan")
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if math.isfinite(number) else float("nan")


def get_p_value(record: Mapping[str, Any]) -> float:
    """Return the p-value of a result record under any of the :data:`P_VALUE_KEYS` spellings.

    The first spelling that holds a value decides (as for :func:`get_chromosome` and :func:`get_position`), so a
    failed test whose ``p_value`` is NaN stays missing instead of falling through to another spelling. The
    result is NaN when the record has no p-value, or one that is not a finite number (``None``, ``"NA"``,
    infinity). Callers decide what a missing p-value means; it is never turned into 1.
    """
    for key in P_VALUE_KEYS:
        value = record.get(key)
        if value is not None:
            return finite_float(value)
    return float("nan")


def get_beta(record: Mapping[str, Any]) -> float:
    """Return the effect size of a result record under any of the :data:`BETA_KEYS` spellings (NaN if none).

    The first spelling that holds a value decides; a value that is not a finite number is NaN, never 0.
    """
    for key in BETA_KEYS:
        value = record.get(key)
        if value is not None:
            return finite_float(value)
    return float("nan")


def _is_missing(value: Any) -> bool:
    """Whether a chromosome or position value is absent: None, NaN or blank text."""
    if value is None:
        return True
    if isinstance(value, str):
        return not value.strip()
    if isinstance(value, bool):
        return False
    return isinstance(value, numbers.Real) and math.isnan(float(value))


def _accepted_keys(column: str) -> str:
    """``chromosome (accepted keys: 'chromosome', 'chrom', ...)`` for a standard column."""
    name, keys = _STANDARD_COLUMNS[column]
    return f"{name} (accepted keys: {', '.join(repr(key) for key in keys)})"


def require_result_columns(df: Any, columns: Sequence[str], *, plot: str) -> None:
    """Raise ``ValueError`` naming the accepted keys of each standard column that ``df`` lacks.

    Call it on the frame :func:`standardize_result_columns` returns, with the subset of ``"CHR"``, ``"BP"``
    and ``"P"`` that ``plot`` needs. Without it a missing column surfaces as a bare ``KeyError`` from deep
    inside the plot.

    Args:
        df: Results DataFrame whose columns have been standardized.
        columns: Standard columns the plot cannot do without.
        plot: Name of the plot, for the message.
    """
    missing = [column for column in columns if column not in df.columns]
    if missing:
        needed = "; ".join(_accepted_keys(column) for column in missing)
        present = df.attrs.get("input_columns", [str(column) for column in df.columns])
        raise ValueError(f"{plot}: the results have no {needed}; the keys present are {sorted(present)}")


def record_placement(record: Mapping[str, Any], index: int, *, plot: str) -> Tuple[str, float]:
    """Return the ``(chromosome label, position)`` of a result record, or raise ``ValueError``.

    The chromosome is read under any of :data:`CHROMOSOME_KEYS` and the position under any of
    :data:`POSITION_KEYS`. A record that has neither cannot be placed on a genome-wide axis, so it is an error
    that names the accepted keys: it is not moved to chromosome 1 or to position 0.

    Args:
        record: A single GWAS result dictionary.
        index: Position of the record in its list, for the message.
        plot: Name of the plot, for the message.
    """
    chromosome = get_chromosome(record)
    position = get_position(record)
    absent = [column for column, value in (("CHR", chromosome), ("BP", position)) if _is_missing(value)]
    if absent:
        needed = " and no ".join(_accepted_keys(column) for column in absent)
        raise ValueError(f"{plot}: result {index} has no {needed}")
    number = finite_float(position)
    if math.isnan(number):
        raise ValueError(f"{plot}: result {index} has a position that is not a finite number: {position!r}")
    return chromosome_label(chromosome), number


def standardize_result_columns(df: Any) -> Any:
    """Return a copy of a results DataFrame whose chromosome, position and p-value columns are CHR, BP and P.

    In each family (:data:`CHROMOSOME_KEYS`, :data:`POSITION_KEYS`, :data:`P_VALUE_KEYS`) the first alias
    present is renamed. A canonical column that already exists is never overwritten and its aliases are left
    alone, so the pipeline's own records (``chrom``/``pos``/``p_value``) and PLINK-style frames
    (``CHR``/``BP``/``P``) are read alike. The input is not modified. The copy remembers the column names it was
    given in ``attrs["input_columns"]``, which :func:`require_result_columns` quotes in its error message.
    """
    renames: Dict[str, str] = {}
    for canonical, aliases in (("CHR", CHROMOSOME_KEYS), ("BP", POSITION_KEYS), ("P", P_VALUE_KEYS)):
        if canonical in df.columns:
            continue
        for alias in aliases:
            if alias in df.columns:
                renames[alias] = canonical
                break
    standardized = df.rename(columns=renames) if renames else df.copy()
    standardized.attrs["input_columns"] = [str(column) for column in df.columns]
    return standardized


def chromosome_label(value: Any) -> str:
    """Normalise a chromosome value to text so ``1``, ``1.0`` and ``"1"`` name the same chromosome.

    Integral numbers (including NumPy scalars, as read from a numeric column with gaps) become their
    integer text; anything else is stripped text. A leading ``chr`` is kept, so ``"chr1"`` and ``"1"``
    stay different labels.
    """
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, numbers.Real):
        number = float(value)
        if math.isfinite(number) and number.is_integer():
            return str(int(number))
    return str(value).strip()


def chromosome_sort_key(label: str) -> Tuple[int, str]:
    """Sort key placing numeric chromosomes numerically, then X, Y, MT, ... by name.

    A ``chr`` prefix is ignored, so ``"chr2"`` sorts before ``"chr10"``.
    """
    name = label.replace("chr", "")
    try:
        return (0, str(int(name)).zfill(5))
    except ValueError:
        return (1, name)


def cumulative_chromosome_offsets(extents: Mapping[str, float], order: Sequence[str]) -> Dict[str, float]:
    """Offset of each chromosome when chromosomes are laid end to end in ``order``.

    Each chromosome starts one base after the largest position of the previous one, so the genome-wide
    x position of a variant is ``offsets[chromosome] + position``.

    Args:
        extents: Largest position observed on each chromosome.
        order: Chromosome labels in plotting order.
    """
    offsets: Dict[str, float] = {}
    running: float = 0
    for chrom in order:
        offsets[chrom] = running
        running += extents[chrom] + 1
    return offsets
