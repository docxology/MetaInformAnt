"""Breakpoint pairs of translocations (internal helpers).

An intra-chromosomal variant (DEL, DUP, INV, INS) occupies the interval ``chrom:start-end``. A translocation
(``TRA``, or a ``BND`` breakend pair) joins two loci instead, so it has no interval: ``chrom``/``start`` is its
first breakpoint, ``chrom2``/``pos2`` its second, and ``end`` repeats ``start`` (see
:class:`~metainformant.structural_variants.detection.sv_calling.StructuralVariant`).

Reciprocal overlap of two such records is always 0, and the same junction can be written from either of its two
sides (``chr1:100 -> chr2:500`` or ``chr2:500 -> chr1:100``). Translocations are therefore compared breakpoint by
breakpoint: two of them are the same event when their breakpoints agree on both chromosomes and each pair of
positions lies within a distance, whichever side each caller wrote first.

Comparison works on the *canonical* form of a junction, its two breakpoints ordered by ``(chromosome, position)``.
For two breakpoints on different chromosomes that order is the same for every caller. For two on one chromosome,
a match in the swapped orientation implies a match in the direct one once both pairs are sorted by position, so
comparing canonical forms finds every match.

The second breakpoint is ``pos2`` on ``chrom2`` (absent, ``None`` or ``""`` means ``chrom``). A record that gives
``chrom2`` on another chromosome without ``pos2`` is read the way VCF-style callers write it: ``end`` is the
position on ``chrom2``. A record with neither (``end == start`` on one chromosome) has no known partner: it is a
single breakend, which agrees only with another single breakend at the same locus.
"""

from __future__ import annotations

import numbers
from collections.abc import Callable, Iterator, Mapping, Sequence
from typing import Any

#: SV types that describe a junction between two loci rather than an interval.
TRANSLOCATION_TYPES = frozenset({"TRA", "BND"})

#: A breakpoint: chromosome and 0-based position.
Breakend = tuple[str, int]

#: Canonical junction: the two breakpoints in ``(chromosome, position)`` order (second None: single breakend).
Junction = tuple[Breakend, "Breakend | None"]


def sv_type_name(variant: Mapping[str, Any]) -> str:
    """Return the variant's SV type as a plain string (enum members are unwrapped; "" when absent)."""
    sv_type = variant.get("sv_type", "")
    sv_type = getattr(sv_type, "value", sv_type)
    return sv_type if isinstance(sv_type, str) else ""


def is_translocation(variant: Mapping[str, Any]) -> bool:
    """Return True for a ``TRA`` or ``BND`` record."""
    return sv_type_name(variant) in TRANSLOCATION_TYPES


def _position(value: Any, field: str) -> int:
    """Return ``value`` as an int coordinate.

    Raises:
        ValueError: If it is not an integral number.
    """
    if not isinstance(value, bool):
        if isinstance(value, numbers.Integral):
            return int(value)
        if isinstance(value, float) and value.is_integer():
            return int(value)
    raise ValueError(f"{field} must be an integer position, got {value!r}")


def breakends(variant: Mapping[str, Any]) -> tuple[Breakend, Breakend | None]:
    """Return the first and second breakpoint of a translocation record, as the record gives them.

    Args:
        variant: Dictionary with ``chrom``, ``start`` and the second breakpoint as ``chrom2``/``pos2`` (or
            ``chrom2``/``end`` in the VCF style; see the module docstring).

    Returns:
        ``((chrom, start), (chrom2, pos2))``; the second breakpoint is None when the record has no known partner.

    Raises:
        ValueError: If a coordinate is not an integral number.
    """
    chrom = variant.get("chrom", "")
    start = _position(variant.get("start", 0), "start")
    chrom2 = variant.get("chrom2") or chrom
    field, raw = "pos2", variant.get("pos2")
    if raw is None:
        # VCF style: END is a position on CHR2 (or, on one chromosome, a second breakpoint after the first).
        end = variant.get("end")
        if end is not None and (chrom2 != chrom or _position(end, "end") != start):
            field, raw = "end", end
    second = None if raw is None else (chrom2, _position(raw, field))
    return (chrom, start), second


def canonical_junction(variant: Mapping[str, Any]) -> Junction:
    """Return a translocation's breakpoints ordered by ``(chromosome, position)``.

    Two records of the same junction written from opposite sides have the same canonical junction.
    """
    first, second = breakends(variant)
    if second is not None and second < first:
        return second, first
    return first, second


def junctions_agree(a: Junction, b: Junction, max_distance: int) -> bool:
    """Return True when two canonical junctions are the same event.

    Both breakpoints must lie on the same chromosomes and within ``max_distance`` bp of each other. A single
    breakend (no known partner) agrees only with another single breakend.
    """
    (a_first, a_second), (b_first, b_second) = a, b
    if a_first[0] != b_first[0] or abs(a_first[1] - b_first[1]) > max_distance:
        return False
    if a_second is None or b_second is None:
        return a_second is None and b_second is None
    return a_second[0] == b_second[0] and abs(a_second[1] - b_second[1]) <= max_distance


def junction_links(junctions: Sequence[Junction], max_distance: int) -> Iterator[tuple[int, int]]:
    """Yield the index pairs ``(i, j)``, ``i < j``, of junctions that are the same event.

    Junctions are grouped by chromosome signature and swept in position order, so the cost is that of sorting
    plus the number of junctions that lie within ``max_distance`` of each other.

    Args:
        junctions: Canonical junctions (see :func:`canonical_junction`).
        max_distance: Largest distance between corresponding breakpoints.
    """
    groups: dict[tuple[str, str | None], list[int]] = {}
    for index, (first, second) in enumerate(junctions):
        groups.setdefault((first[0], None if second is None else second[0]), []).append(index)

    for members in groups.values():
        members.sort(key=lambda k: junctions[k][0][1])
        for rank, i in enumerate(members):
            for j in members[rank + 1 :]:
                if junctions[j][0][1] - junctions[i][0][1] > max_distance:
                    break
                if junctions_agree(junctions[i], junctions[j], max_distance):
                    yield (i, j) if i < j else (j, i)


def consensus_junction(junctions: Sequence[Junction], median: Callable[[list[int]], int]) -> Junction:
    """Return the consensus of junctions that agree: per breakpoint, the median position.

    Args:
        junctions: Canonical junctions of one group (same chromosomes, as :func:`junction_links` guarantees).
        median: Function returning the median of a list of positions as an int.
    """
    first_chrom = junctions[0][0][0]
    first_pos = median([first[1] for first, _ in junctions])
    seconds = [second for _, second in junctions if second is not None]
    if not seconds:
        return (first_chrom, first_pos), None
    return (first_chrom, first_pos), (seconds[0][0], median([second[1] for second in seconds]))


def variant_regions(variant: Mapping[str, Any]) -> list[tuple[str, int, int]]:
    """Return the genomic regions a variant touches.

    An interval variant touches ``(chrom, start, end)``. A translocation touches one zero-length region per
    breakpoint: the junction itself on each chromosome (a single breakend has just the one).
    """
    if is_translocation(variant):
        first, second = breakends(variant)
        regions = [(first[0], first[1], first[1])]
        if second is not None:
            regions.append((second[0], second[1], second[1]))
        return regions
    return [(variant.get("chrom", ""), variant.get("start", 0), variant.get("end", 0))]
