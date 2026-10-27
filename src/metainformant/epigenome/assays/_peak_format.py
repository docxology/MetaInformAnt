"""Shared reader and writer for the standard peak-file layouts.

ChIP-seq and ATAC-seq peaks are stored in the same three tab-separated layouts. The
column order is defined once, here, and both assay modules read and write through these
helpers, so a file written with one layout can always be read back with that layout.

Layouts (ENCODE / UCSC):

* ``narrowpeak`` -- ten columns: ``chrom``, ``chromStart``, ``chromEnd``, ``name``, ``score``,
  ``strand``, ``signalValue``, ``pValue``, ``qValue``, ``peak``. ``peak`` is the summit as an
  offset from ``chromStart``.
* ``broadpeak`` -- the first nine narrowPeak columns (there is no ``peak`` column).
* ``bed`` -- the first six columns: ``chrom``, ``chromStart``, ``chromEnd``, ``name``, ``score``,
  ``strand``. Only three columns are required when reading.

A value that is not available is written as ``-1`` in the ``pValue``, ``qValue`` and ``peak``
columns, as the specification requires, and is read back from ``-1`` or ``.`` as ``None``.

The specification (and MACS2) store ``pValue`` and ``qValue`` as ``-log10`` values, while peak
objects carry probabilities. The conversion happens here, at the file boundary: reading turns
a column value ``v`` into ``10 ** -v``, and writing turns a probability ``p`` into ``-log10(p)``.
A probability of exactly 0 is written as the ``-log10`` of the smallest positive float.
"""

from __future__ import annotations

import math
from typing import NamedTuple, Optional, Sequence, Tuple

PEAK_FORMATS: Tuple[str, ...] = ("narrowpeak", "broadpeak", "bed")

# Columns a line needs for each layout before it can be parsed.
MIN_COLUMNS = {"narrowpeak": 10, "broadpeak": 9, "bed": 3}


class PeakFields(NamedTuple):
    """The fields parsed from one peak line; absent optional columns are ``None``."""

    chromosome: str
    start: int
    end: int
    score: float
    strand: str
    signal_value: Optional[float]
    p_value: Optional[float]
    q_value: Optional[float]
    summit: Optional[int]


def normalise_peak_format(format: str) -> str:
    """Return the lower-case layout name, rejecting anything that is not a supported layout.

    Args:
        format: ``"narrowPeak"``, ``"broadPeak"`` or ``"bed"`` in any letter case.

    Returns:
        ``"narrowpeak"``, ``"broadpeak"`` or ``"bed"``.

    Raises:
        ValueError: If the layout is not supported.
    """
    key = str(format).strip().lower()
    if key not in PEAK_FORMATS:
        raise ValueError(f"Unsupported peak format {format!r}; expected one of {', '.join(PEAK_FORMATS)}")
    return key


def is_skipped_line(line: str) -> bool:
    """Return True for blank, comment, ``track`` and ``browser`` lines (already stripped)."""
    return not line or line.startswith(("#", "track", "browser"))


# -log10 of the smallest positive float: how a probability of exactly 0 is written.
_MAX_MINUS_LOG10 = -math.log10(5e-324)


def _probability_from_column(token: str, column: str) -> Optional[float]:
    """Read a ``-log10`` significance column as a probability; ``.`` and ``-1`` mean "not available"."""
    if token == ".":
        return None
    value = float(token)
    if value == -1.0:
        return None
    if not value >= 0.0:  # also rejects NaN
        raise ValueError(f"{column} column must hold -log10 of a probability (>= 0) or -1, found {token!r}")
    return float(10.0**-value)  # float() because typeshed types float ** float as Any (the result may be complex)


def _column_from_probability(value: Optional[float], column: str) -> str:
    """Write a probability as its ``-log10`` column value, or ``-1`` when it is not available."""
    if value is None:
        return "-1"
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{column} must be a probability in [0, 1] (or None), got {value!r}")
    return _format_number(_MAX_MINUS_LOG10 if value == 0.0 else -math.log10(value))


def _optional_offset(token: str) -> Optional[int]:
    """Parse the ``peak`` column: an integer offset, with ``.`` and ``-1`` meaning "no summit"."""
    if token == ".":
        return None
    value = int(token)
    return None if value == -1 else value


def parse_peak_fields(parts: Sequence[str], format: str) -> PeakFields:
    """Parse the tab-split columns of one peak line.

    Args:
        parts: The columns of the line.
        format: Layout of the line (``"narrowpeak"``, ``"broadpeak"`` or ``"bed"``).

    Returns:
        The parsed fields. ``p_value`` and ``q_value`` are probabilities converted from the
        ``-log10`` columns. ``signal_value``, ``p_value``, ``q_value`` and ``summit`` are ``None``
        when the layout has no such column or the column is not available.

    Raises:
        ValueError: If the line has too few columns for the layout or a column does not parse.
    """
    fmt = normalise_peak_format(format)
    needed = MIN_COLUMNS[fmt]
    if len(parts) < needed:
        raise ValueError(f"{fmt} lines need at least {needed} columns, found {len(parts)}")

    start = int(parts[1])
    end = int(parts[2])
    score = 0.0 if len(parts) < 5 or parts[4] == "." else float(parts[4])
    strand = parts[5] if len(parts) > 5 else "."

    signal_value: Optional[float] = None
    p_value: Optional[float] = None
    q_value: Optional[float] = None
    summit: Optional[int] = None

    if fmt in ("narrowpeak", "broadpeak"):
        signal_value = float(parts[6])
        p_value = _probability_from_column(parts[7], "pValue")
        q_value = _probability_from_column(parts[8], "qValue")
    if fmt == "narrowpeak":
        offset = _optional_offset(parts[9])
        summit = None if offset is None else start + offset

    return PeakFields(parts[0], start, end, score, strand, signal_value, p_value, q_value, summit)


def _format_number(value: float) -> str:
    """Write a number without a spurious fraction and without losing precision."""
    if math.isfinite(value) and value == int(value) and abs(value) < 1e15:
        return str(int(value))
    return repr(float(value))


def format_peak_line(
    *,
    name: str,
    chromosome: str,
    start: int,
    end: int,
    score: float,
    strand: str,
    signal_value: float,
    p_value: Optional[float],
    q_value: Optional[float],
    summit: Optional[int],
    format: str,
) -> str:
    """Format one peak as a line (without newline) of the requested layout.

    Args:
        name: Peak name (column 4).
        chromosome: Chromosome name.
        start: Peak start.
        end: Peak end.
        score: Peak score.
        strand: ``"+"``, ``"-"`` or ``"."``.
        signal_value: Signal value (``narrowpeak`` and ``broadpeak`` only).
        p_value: p-value as a probability, written as ``-log10``; ``None`` for ``-1``.
        q_value: q-value as a probability, written as ``-log10``; ``None`` for ``-1``.
        summit: Absolute summit position, written as an offset from ``start`` (``narrowpeak`` only), or
            ``None`` for ``-1``.
        format: Layout to write.

    Returns:
        The tab-separated line.

    Raises:
        ValueError: If ``p_value`` or ``q_value`` is not a probability in [0, 1].
    """
    fmt = normalise_peak_format(format)
    columns = [chromosome, str(start), str(end), name, _format_number(score), strand or "."]
    if fmt != "bed":
        columns.append(_format_number(signal_value))
        columns.append(_column_from_probability(p_value, "p_value"))
        columns.append(_column_from_probability(q_value, "q_value"))
    if fmt == "narrowpeak":
        columns.append("-1" if summit is None else str(summit - start))
    return "\t".join(columns)
