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


def standardize_result_columns(df: Any) -> Any:
    """Return a copy of a results DataFrame whose chromosome, position and p-value columns are CHR, BP and P.

    In each family (:data:`CHROMOSOME_KEYS`, :data:`POSITION_KEYS`, :data:`P_VALUE_KEYS`) the first alias
    present is renamed. A canonical column that already exists is never overwritten and its aliases are left
    alone, so the pipeline's own records (``chrom``/``pos``/``p_value``) and PLINK-style frames
    (``CHR``/``BP``/``P``) are read alike. The input is not modified.
    """
    renames: Dict[str, str] = {}
    for canonical, aliases in (("CHR", CHROMOSOME_KEYS), ("BP", POSITION_KEYS), ("P", P_VALUE_KEYS)):
        if canonical in df.columns:
            continue
        for alias in aliases:
            if alias in df.columns:
                renames[alias] = canonical
                break
    return df.rename(columns=renames) if renames else df.copy()


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
