"""Long-read BAM file processing with methylation and supplementary alignment support.

Provides real parsing of BAM files from long-read sequencers with support for
modified base tags (MM/ML), supplementary alignment extraction for structural
variant detection, and comprehensive alignment statistics.

Modified base tags (MM/ML) follow the SAM optional-tags specification (hts-specs ``SAMtags``,
"Base modifications"), checked against the htslib reference implementation:

* An MM entry ``<base><strand><codes>[<flag>],<skip>,<skip>,...;`` counts the letter ``<base>`` in the
  read **as it was sequenced**. For an alignment with FLAG 0x10 the SEQ field is the reverse
  complement of that read, so SEQ is reverse-complemented before counting and the positions are mapped
  back to the stored (SEQ) orientation afterwards. Each skip counts the letters passed over before the
  next modified one.
* The strand says on which strand the modification sits: ``+`` the sequenced strand, ``-`` the opposite
  one. The counted letter is the same either way, so ``G-m`` counts the Gs of the read and marks the
  5-methylcytosines opposite them (``modified_base`` of such an entry is the complement of the letter).
* An entry with several codes (``C+mh``) carries one ML value per code at every listed position,
  interleaved in code order; ML is consumed entry after entry in MM order.
* The flag after the codes says what the skipped bases are: ``.`` unmodified, ``?`` unknown. Without
  a flag the specification's default, ``.``, applies.
* MM without ML leaves the probabilities unknown (not zero).
* When present, MN must equal the stored SEQ length; a mismatch makes MM/ML stale and decoding raises.

Optional dependencies:
    - pysam: For reading BAM/CRAM files
"""

from __future__ import annotations

import operator
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, NamedTuple, Sequence

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

try:
    import pysam
except ImportError:
    pysam = None


@dataclass
class LongReadAlignment:
    """Represents a single alignment from a long-read BAM file.

    Attributes:
        read_name: Query read name.
        query_sequence: The aligned read sequence.
        query_qualities: Per-base quality scores (list of ints).
        reference_name: Reference contig name.
        reference_start: 0-based start position on reference.
        reference_end: 0-based end position on reference.
        mapping_quality: MAPQ score.
        cigar_string: CIGAR alignment string.
        cigar_tuples: CIGAR as list of (operation, length) tuples.
        is_supplementary: Whether this is a supplementary alignment.
        is_secondary: Whether this is a secondary alignment.
        is_reverse: Whether the read is reverse complemented.
        is_unmapped: Whether the read is unmapped.
        tags: Dictionary of BAM tags.
        methylation_tags: Parsed methylation data from MM/ML tags.
        query_length: Total length of the original query.
        aligned_pairs: List of (query_pos, ref_pos) alignment pairs.
    """

    read_name: str
    query_sequence: str = ""
    query_qualities: list[int] = field(default_factory=list)
    reference_name: str = ""
    reference_start: int = 0
    reference_end: int = 0
    mapping_quality: int = 0
    cigar_string: str = ""
    cigar_tuples: list[tuple[int, int]] = field(default_factory=list)
    is_supplementary: bool = False
    is_secondary: bool = False
    is_reverse: bool = False
    is_unmapped: bool = False
    tags: dict[str, Any] = field(default_factory=dict)
    methylation_tags: dict[str, Any] = field(default_factory=dict)
    query_length: int = 0
    aligned_pairs: list[tuple[int | None, int | None]] = field(default_factory=list)


@dataclass
class AlignmentStats:
    """Summary statistics for a set of long-read alignments.

    Attributes:
        total_reads: Total number of reads.
        mapped_reads: Number of mapped reads.
        unmapped_reads: Number of unmapped reads.
        supplementary_count: Number of supplementary alignments.
        secondary_count: Number of secondary alignments.
        mean_mapping_quality: Average MAPQ across mapped reads.
        mean_read_length: Average query length.
        mean_aligned_length: Average aligned length.
        mean_identity: Average alignment identity (matches / aligned_length).
        total_bases: Total bases in all reads.
        total_aligned_bases: Total aligned bases.
        coverage_depth: Estimated coverage depth (if reference length provided).
    """

    total_reads: int = 0
    mapped_reads: int = 0
    unmapped_reads: int = 0
    supplementary_count: int = 0
    secondary_count: int = 0
    mean_mapping_quality: float = 0.0
    mean_read_length: float = 0.0
    mean_aligned_length: float = 0.0
    mean_identity: float = 0.0
    total_bases: int = 0
    total_aligned_bases: int = 0
    coverage_depth: float = 0.0


def read_long_read_bam(
    filepath: str | Path,
    region: str | None = None,
    include_supplementary: bool = True,
    include_secondary: bool = False,
    min_mapping_quality: int = 0,
) -> list[LongReadAlignment]:
    """Read a BAM file with long-read specific tag extraction.

    Parses BAM files from ONT or PacBio sequencers, extracting methylation
    tags (MM/ML), haplotype tags (HP), and other long-read specific information.

    Args:
        filepath: Path to the BAM file (must be indexed with .bai).
        region: Optional genomic region in format "chr:start-end".
        include_supplementary: Whether to include supplementary alignments.
        include_secondary: Whether to include secondary alignments.
        min_mapping_quality: Minimum MAPQ filter.

    Returns:
        List of LongReadAlignment objects.

    Raises:
        FileNotFoundError: If the BAM file does not exist.
        ImportError: If pysam is not installed.
    """
    if pysam is None:
        raise ImportError("pysam is required for reading BAM files. " "Install it with: uv pip install pysam")

    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"BAM file not found: {filepath}")

    alignments: list[LongReadAlignment] = []

    with pysam.AlignmentFile(str(filepath), "rb") as bam:
        iterator = bam.fetch(region=region) if region else bam.fetch(until_eof=True)

        for read in iterator:
            if read.is_secondary and not include_secondary:
                continue
            if read.is_supplementary and not include_supplementary:
                continue
            if not read.is_unmapped and read.mapping_quality < min_mapping_quality:
                continue

            # Extract tags as dict
            tags_dict: dict[str, Any] = {}
            try:
                tags_dict = dict(read.get_tags())
            except Exception:
                pass

            # Parse methylation tags
            meth_tags = _parse_methylation_from_tags(tags_dict)

            # Get CIGAR
            cigar_tuples = list(read.cigartuples) if read.cigartuples else []
            cigar_string = read.cigarstring or ""

            # Get aligned pairs (limited to avoid memory issues on very long reads)
            aligned_pairs: list[tuple[int | None, int | None]] = []
            try:
                pairs = read.get_aligned_pairs()
                aligned_pairs = [(p[0], p[1]) for p in pairs]
            except Exception:
                pass

            aln = LongReadAlignment(
                read_name=read.query_name or "",
                query_sequence=read.query_sequence or "",
                query_qualities=list(read.query_qualities) if read.query_qualities is not None else [],
                reference_name=read.reference_name or "",
                reference_start=read.reference_start,
                reference_end=read.reference_end,
                mapping_quality=read.mapping_quality,
                cigar_string=cigar_string,
                cigar_tuples=cigar_tuples,
                is_supplementary=read.is_supplementary,
                is_secondary=read.is_secondary,
                is_reverse=read.is_reverse,
                is_unmapped=read.is_unmapped,
                tags=tags_dict,
                methylation_tags=meth_tags,
                query_length=read.query_length or 0,
                aligned_pairs=aligned_pairs,
            )
            alignments.append(aln)

    logger.info(
        "Read %d alignments from %s%s",
        len(alignments),
        filepath.name,
        f" region={region}" if region else "",
    )
    return alignments


#: Scale of the ML values: 0 (certainly unmodified) to 255; the specification quantises a probability in
#: 1/256ths, and this module reports value / 255 on the 0-1 scale.
_ML_SCALE = 255

_MM_ENTRY = re.compile(r"([ACGTUN])([+-])([A-Za-z]+|[0-9]+)([.?]?)((?:\s*,\s*[0-9]+)*)")
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


class _MmEntry(NamedTuple):
    """One parsed entry of an MM tag: ``C+mh?,5,12,0;`` is base C, strand +, codes (m, h), flag ?, skips (5, 12, 0)."""

    base: str
    strand: str
    codes: tuple[str, ...]
    flag: str
    skips: tuple[int, ...]


def _complement_base(base: str) -> str:
    """Complement of one A/C/G/T/U/N letter (U and T pair with A; N stays N)."""
    return "A" if base == "U" else base.translate(_COMPLEMENT)


def _parse_mm_tag(mm_tag: Any) -> list[_MmEntry]:
    """Parse the text of an MM tag into entries.

    Grammar (SAM tags specification): ``([ACGTUN][-+]([a-z]+|[0-9]+)[.?]?(,[0-9]+)*;)*``. A run of letters
    is several codes (``mh`` is ``m`` and ``h``), a run of digits is one ChEBI code, and an entry may list
    no positions (``C+m;``). A missing final ``;`` is accepted.

    Raises:
        ValueError: If an entry does not follow the grammar. The skips of a malformed entry cannot be told
            apart from the ML values that belong to it, so nothing after it can be decoded reliably.
    """
    text = mm_tag.decode("ascii", errors="replace") if isinstance(mm_tag, (bytes, bytearray)) else str(mm_tag)
    entries: list[_MmEntry] = []
    for raw in text.split(";"):
        raw = raw.strip()
        if not raw:
            continue
        match = _MM_ENTRY.fullmatch(raw)
        if match is None:
            raise ValueError(f"malformed MM entry {raw!r}")
        base, strand, codes, flag, skips = match.groups()
        entries.append(
            _MmEntry(
                base=base,
                strand=strand,
                codes=(codes,) if codes.isdigit() else tuple(codes),
                flag=flag,
                skips=tuple(int(token) for token in skips.split(",") if token.strip()),
            )
        )
    return entries


def _parse_methylation_from_tags(tags: dict[str, Any]) -> dict[str, Any]:
    """Parse MM/ML base-modification tags from BAM alignment tags.

    The MM tag lists, per base and strand, how many bases to skip before each modified one::

        MM:Z:C+m,0,1,0;   (5mC at the 1st, 3rd and 4th C of the read as sequenced)
        MM:Z:C+mh?,2;A+a,0;

    and the ML tag holds one probability (0-255) per listed position and code, entry after entry in MM
    order, interleaved in code order within an entry (``C+mh,3,1`` with ML ``10,20,30,40`` is m=10, h=20 at
    the first position and m=30, h=40 at the second). See the module docstring for the full rules.

    The optional MN value is retained for validation against the read sequence.
    Positions are not decoded here: that needs the read sequence and its orientation, see
    :func:`extract_methylation_tags`.

    Args:
        tags: Dictionary of BAM tags. ``ML`` may be a list, tuple, bytes, numpy
            array or the ``array('B')`` that pysam returns for B-type tags. The
            draft names ``Mm`` and ``Ml`` are read when ``MM`` and ``ML`` are absent.

    Returns:
        ``{}`` without an MM tag, or when the tag is malformed (a warning says why). Otherwise
        ``{"modifications": [...]}`` with one dictionary per entry and code, in MM order, holding ``base``
        (the letter counted), ``strand``, ``modification`` (one code), ``modified_base`` (the base that
        carries the modification: ``base`` for ``+``, its complement for ``-``), ``implicit`` (True when the
        skipped bases are unmodified, the ``.`` flag or no flag; False for ``?``, where they are unknown),
        ``delta_positions``, ``probabilities`` (this code's raw 0-255 ML values, one per position; empty
        when there is no ML tag or its length does not match the MM tag, which a warning reports) and
        ``probability_scale``.
    """
    result: dict[str, Any] = {}

    # Explicit None checks: ``or`` would call bool() on an array
    mm_tag = tags.get("MM")
    if mm_tag is None:
        mm_tag = tags.get("Mm")
    ml_tag = tags.get("ML")
    if ml_tag is None:
        ml_tag = tags.get("Ml")

    if mm_tag is None:
        return result

    try:
        entries = _parse_mm_tag(mm_tag)
    except ValueError as exc:
        logger.warning("Base modifications not decoded: %s", exc)
        return result

    ml_probs = _ml_values(ml_tag)
    expected = sum(len(entry.skips) * len(entry.codes) for entry in entries)
    if ml_probs and len(ml_probs) != expected:
        logger.warning(
            "Base modification probabilities not used: the ML tag has %d values but the MM tag lists %d "
            "(positions x codes)",
            len(ml_probs),
            expected,
        )
        ml_probs = []

    modifications: list[dict[str, Any]] = []
    offset = 0
    for entry in entries:
        n_codes = len(entry.codes)
        n_positions = len(entry.skips)
        for code_index, code in enumerate(entry.codes):
            # ML interleaves the codes of an entry: position p, code c is at offset + p * n_codes + c
            probs = ml_probs[offset + code_index : offset + n_positions * n_codes : n_codes] if ml_probs else []
            modifications.append(
                {
                    "base": entry.base,
                    "strand": entry.strand,
                    "modification": code,
                    "modified_base": entry.base if entry.strand == "+" else _complement_base(entry.base),
                    "implicit": entry.flag != "?",
                    "delta_positions": list(entry.skips),
                    "probabilities": probs,
                    "probability_scale": _ML_SCALE,
                }
            )
        offset += n_positions * n_codes

    if modifications:
        result["modifications"] = modifications
        if "MN" in tags:
            result["MN"] = tags["MN"]

    return result


def _ml_values(ml_tag: Any) -> list[int]:
    """Return the values of an ML tag as a list of ints.

    Accepts what BAM readers hand back for the B-type ML tag: a list or tuple,
    bytes, a numpy array, pysam's ``array('B')``, or a single value.
    """
    if ml_tag is None:
        return []
    if isinstance(ml_tag, (bytes, bytearray)):
        return list(ml_tag)
    if isinstance(ml_tag, str):
        return [int(ml_tag)]
    try:
        return [int(value) for value in ml_tag]
    except TypeError:  # a scalar
        return [int(ml_tag)]


def _original_orientation(sequence: str, is_reverse: bool) -> str:
    """The read as it was sequenced: SEQ, reverse-complemented for an alignment on the reverse strand (FLAG 0x10).

    Upper-cased, with RNA ``U`` read as ``T`` (SEQ stores T). Letters other than A/C/G/T stay as they are
    (their complement is not needed: only the four bases and ``N`` are ever counted).
    """
    seq = sequence.upper().replace("U", "T")
    if not is_reverse:
        return seq
    return seq[::-1].translate(_COMPLEMENT)


def _locate_sites(
    sequence: str, base: str, skips: Sequence[int], is_reverse: bool
) -> tuple[list[tuple[int, int]], list[int]]:
    """Place the listed positions of one MM entry on the stored read.

    Counts the letter ``base`` (any base for ``N``) in the read as it was sequenced, passing over ``skip``
    letters before each modified one, and maps the positions back to indices of ``sequence`` as stored
    (the reverse complement's index ``i`` is stored index ``len - 1 - i``).

    Returns:
        ``(placed, skipped)``: ``placed`` holds ``(stored position, index into skips)`` in ascending
        stored order for the positions found; a skip that runs past the last letter of the read ends the
        list, so the positions after it are missing. ``skipped`` holds the stored positions of every
        counted letter that is not a listed one, ascending.
    """
    length = len(sequence)
    original = _original_orientation(sequence, is_reverse)
    letter = "T" if base == "U" else base
    counted = [i for i, c in enumerate(original) if letter == "N" or c == letter]

    found: list[tuple[int, int]] = []  # (index into counted letters of the sequenced read, index into skips)
    cursor = 0
    for index, skip in enumerate(skips):
        cursor += skip
        if cursor >= len(counted):
            break
        found.append((counted[cursor], index))
        cursor += 1

    def stored(i: int) -> int:
        return length - 1 - i if is_reverse else i

    listed = {position for position, _ in found}
    placed = sorted((stored(position), index) for position, index in found)
    skipped = sorted(stored(i) for i in counted if i not in listed)
    return placed, skipped


def _aligned_pairs_from_cigar(
    cigar_tuples: Sequence[tuple[int, int]], reference_start: int
) -> list[tuple[int | None, int | None]]:
    """``(query position, reference position)`` pairs of an alignment, like pysam's ``get_aligned_pairs()``.

    Query positions index SEQ as stored. Matches (``M``, ``=``, ``X``) pair both coordinates, insertions and
    soft clips (``I``, ``S``) have no reference position, deletions and skipped regions (``D``, ``N``) have no
    query position, and hard clips and padding (``H``, ``P``) are not in SEQ and give no pair.
    """
    pairs: list[tuple[int | None, int | None]] = []
    query = 0
    reference = reference_start
    for operation, length in cigar_tuples:
        if operation in (0, 7, 8):
            pairs.extend((query + i, reference + i) for i in range(length))
            query += length
            reference += length
        elif operation in (1, 4):
            pairs.extend((query + i, None) for i in range(length))
            query += length
        elif operation in (2, 3):
            pairs.extend((None, reference + i) for i in range(length))
            reference += length
    return pairs


_CIGAR_OPERATIONS = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5, "P": 6, "=": 7, "X": 8}


def _alignment_from_read_dict(read: dict[str, Any]) -> LongReadAlignment:
    """A LongReadAlignment for a read dictionary as the workflow loaders produce it.

    Reads the sequence from ``query_sequence`` or ``sequence``, and the tags, ``is_reverse`` and
    ``reference_*`` fields under their usual names. ``aligned_pairs`` are taken from the dictionary or,
    when it has none, built from ``cigar_tuples`` / ``cigar_string`` and ``reference_start``; a read with
    neither has no reference coordinates.
    """
    cigar_tuples = list(read.get("cigar_tuples") or [])
    cigar_string = read.get("cigar_string") or ""
    if not cigar_tuples and cigar_string:
        cigar_tuples = [(_CIGAR_OPERATIONS[op], int(n)) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar_string)]

    aligned_pairs = read.get("aligned_pairs")
    if aligned_pairs is None or len(aligned_pairs) == 0:
        aligned_pairs = _aligned_pairs_from_cigar(cigar_tuples, int(read.get("reference_start") or 0))

    return LongReadAlignment(
        read_name=str(read.get("read_name") or read.get("read_id") or ""),
        query_sequence=read.get("query_sequence") or read.get("sequence") or "",
        reference_name=read.get("reference_name") or "",
        reference_start=int(read.get("reference_start") or 0),
        cigar_string=cigar_string,
        cigar_tuples=cigar_tuples,
        is_reverse=bool(read.get("is_reverse", False)),
        is_unmapped=bool(read.get("is_unmapped", False)),
        tags=dict(read.get("tags") or {}),
        aligned_pairs=list(aligned_pairs),
    )


def _validate_mn_length(tags: dict[str, Any], sequence_length: int, read_name: str) -> None:
    """Reject modification coordinates produced for a different SEQ length."""
    if "MN" not in tags:
        return
    value = tags["MN"]
    try:
        mn_length = operator.index(value)
    except TypeError:
        raise ValueError(f"Read {read_name}: MN tag must be a nonnegative integer") from None
    if isinstance(value, bool) or mn_length < 0:
        raise ValueError(f"Read {read_name}: MN tag must be a nonnegative integer")
    if mn_length != sequence_length:
        raise ValueError(
            f"Read {read_name}: MN tag length {mn_length} does not match SEQ length {sequence_length}; "
            "MM/ML values are out-of-date"
        )


def extract_methylation_tags(alignment: LongReadAlignment) -> dict[str, Any]:
    """Extract and place the base modifications of a LongReadAlignment.

    Decodes the MM/ML tags into positions on the read and on the reference, following the SAM tags
    specification (see the module docstring): the skip counts refer to the read as sequenced, so
    for a reverse-aligned read (``is_reverse``) the counting runs over the reverse complement of
    ``query_sequence`` and the positions are mapped back to the stored sequence. Covers CpG (5mC), 6mA and
    every other code the tags carry.

    Args:
        alignment: A LongReadAlignment with its ``tags`` (or ``methylation_tags``), ``query_sequence``,
            ``is_reverse`` and, for reference coordinates, ``aligned_pairs``.

    Returns:
        Dictionary keyed by base, strand and modification code as written in the MM tag
        (``"C+m"`` and ``"C-m"`` are separate entries; ``C+mh`` gives ``"C+m"`` and ``"C+h"``), each holding:

        * ``base``, ``strand``, ``modification``, ``modified_base`` and ``implicit`` as in
          :func:`_parse_methylation_from_tags`.
        * ``read_positions``: 0-based indices into ``query_sequence`` as stored, ascending.
        * ``genomic_positions``: ``(reference position or None, probability or None)`` per listed site, in
          the order of ``read_positions``; the reference position is None where the base is not aligned
          (insertion, soft clip, unaligned read), and the probability is None without an ML tag.
        * ``probabilities``: the ML values on a 0-1 scale (value / 255), in the order of ``read_positions``;
          empty without usable ML values.
        * ``skipped_read_positions`` and ``skipped_genomic_positions``: the counted bases the entry does not
          list. They are unmodified when ``implicit`` is True and of unknown status when it is False.
        * ``n_unplaced``: listed positions that could not be placed because the read ran out of counted
          bases (a warning reports them) or has no sequence.

    Raises:
        ValueError: If an MN tag is invalid or differs from the stored SEQ length,
            including an MN value retained in prematerialized methylation tags.
    """
    result: dict[str, Any] = {}
    seq = alignment.query_sequence or ""
    _validate_mn_length(alignment.tags, len(seq), alignment.read_name)
    _validate_mn_length(alignment.methylation_tags, len(seq), alignment.read_name)

    if not alignment.methylation_tags or "modifications" not in alignment.methylation_tags:
        # Try parsing from raw tags
        parsed = _parse_methylation_from_tags(alignment.tags)
        if not parsed or "modifications" not in parsed:
            return result
        mods = parsed["modifications"]
    else:
        mods = alignment.methylation_tags["modifications"]

    to_reference = {
        query: reference for query, reference in alignment.aligned_pairs if query is not None and reference is not None
    }

    for mod_info in mods:
        base_char = mod_info["base"]
        strand = mod_info.get("strand", "+")
        mod_code = mod_info["modification"]
        deltas = mod_info["delta_positions"]
        probs = mod_info.get("probabilities", [])
        if probs and len(probs) != len(deltas):
            logger.warning(
                "Read %s: %s%s%s lists %d positions but has %d probabilities; the probabilities are not used",
                alignment.read_name,
                base_char,
                strand,
                mod_code,
                len(deltas),
                len(probs),
            )
            probs = []

        placed, skipped = _locate_sites(seq, base_char, deltas, alignment.is_reverse) if seq else ([], [])
        n_unplaced = len(deltas) - len(placed)
        if n_unplaced and not seq:
            logger.warning(
                "Read %s: %s%s%s lists %d positions but the alignment has no sequence to place them on",
                alignment.read_name,
                base_char,
                strand,
                mod_code,
                n_unplaced,
            )
        elif n_unplaced:
            logger.warning(
                "Read %s: %s%s%s lists %d positions but only %d could be placed on the %d bp read; %d lie beyond "
                "its last %s",
                alignment.read_name,
                base_char,
                strand,
                mod_code,
                len(deltas),
                len(placed),
                len(seq),
                n_unplaced,
                base_char,
            )

        read_positions = [position for position, _ in placed]
        site_probs = [probs[index] / _ML_SCALE for _, index in placed] if probs else []
        genomic_positions: list[tuple[int | None, float | None]] = [
            (to_reference.get(position), site_probs[i] if site_probs else None)
            for i, position in enumerate(read_positions)
        ]

        mod_key = f"{base_char}{strand}{mod_code}"
        result[mod_key] = {
            "base": base_char,
            "strand": strand,
            "modification": mod_code,
            "modified_base": mod_info.get("modified_base")
            or (base_char if strand == "+" else _complement_base(base_char)),
            "implicit": mod_info.get("implicit", True),
            "read_positions": read_positions,
            "genomic_positions": genomic_positions,
            "probabilities": site_probs,
            "skipped_read_positions": skipped,
            "skipped_genomic_positions": [to_reference.get(position) for position in skipped],
            "n_unplaced": n_unplaced,
        }

    return result


def get_supplementary_alignments(alignment: LongReadAlignment) -> list[dict[str, Any]]:
    """Extract supplementary alignment information for split-read SV detection.

    Parses the SA tag which encodes supplementary alignments in the format:
        SA:Z:chr,pos,strand,CIGAR,mapQ,NM;

    Args:
        alignment: A LongReadAlignment object.

    Returns:
        List of dictionaries, each containing supplementary alignment details
        including reference name, position, strand, CIGAR, and mapping quality.
    """
    sa_tag = alignment.tags.get("SA", "")
    if not sa_tag:
        return []

    supplementary: list[dict[str, Any]] = []

    entries = [e.strip() for e in str(sa_tag).split(";") if e.strip()]
    for entry in entries:
        parts = entry.split(",")
        if len(parts) < 6:
            continue

        ref_name = parts[0]
        try:
            pos = int(parts[1]) - 1  # Convert 1-based SAM to 0-based
        except ValueError:
            continue
        strand = parts[2]
        cigar = parts[3]
        try:
            mapq = int(parts[4])
        except ValueError:
            mapq = 0
        try:
            nm = int(parts[5])
        except ValueError:
            nm = 0

        # Parse CIGAR to get aligned length
        aligned_length = _cigar_aligned_length(cigar)

        supplementary.append(
            {
                "reference_name": ref_name,
                "reference_start": pos,
                "strand": strand,
                "cigar": cigar,
                "mapping_quality": mapq,
                "edit_distance": nm,
                "aligned_length": aligned_length,
                "is_reverse": strand == "-",
            }
        )

    return supplementary


def _cigar_aligned_length(cigar_string: str) -> int:
    """Calculate the aligned reference length from a CIGAR string.

    CIGAR operations that consume reference:
        M (alignment match), D (deletion), N (skipped region), = (match), X (mismatch)

    Args:
        cigar_string: CIGAR alignment string.

    Returns:
        Total reference-consuming length.
    """
    ref_consuming = {"M", "D", "N", "=", "X"}
    total = 0
    for match in re.finditer(r"(\d+)([MIDNSHP=X])", cigar_string):
        length = int(match.group(1))
        op = match.group(2)
        if op in ref_consuming:
            total += length
    return total


def _cigar_query_length(cigar_string: str) -> int:
    """Calculate the query-consuming length from a CIGAR string.

    CIGAR operations that consume query:
        M (alignment match), I (insertion), S (soft clip), = (match), X (mismatch)

    Args:
        cigar_string: CIGAR alignment string.

    Returns:
        Total query-consuming length.
    """
    query_consuming = {"M", "I", "S", "=", "X"}
    total = 0
    for match in re.finditer(r"(\d+)([MIDNSHP=X])", cigar_string):
        length = int(match.group(1))
        op = match.group(2)
        if op in query_consuming:
            total += length
    return total


def calculate_alignment_stats(
    alignments: Sequence[LongReadAlignment],
    reference_length: int | None = None,
) -> AlignmentStats:
    """Calculate comprehensive alignment statistics for a set of long-read alignments.

    Computes mapping quality, identity, coverage, and read length distributions
    from a collection of alignments.

    Args:
        alignments: Sequence of LongReadAlignment objects.
        reference_length: Total reference genome length (for coverage calculation).

    Returns:
        AlignmentStats with aggregated metrics.
    """
    if not alignments:
        return AlignmentStats()

    total = len(alignments)
    mapped = 0
    unmapped = 0
    supplementary = 0
    secondary = 0
    mapq_sum = 0.0
    read_lengths: list[int] = []
    aligned_lengths: list[int] = []
    identities: list[float] = []
    total_bases = 0
    total_aligned = 0

    for aln in alignments:
        if aln.is_unmapped:
            unmapped += 1
            continue

        if aln.is_supplementary:
            supplementary += 1
        if aln.is_secondary:
            secondary += 1

        mapped += 1
        mapq_sum += aln.mapping_quality

        ql = aln.query_length or len(aln.query_sequence)
        read_lengths.append(ql)
        total_bases += ql

        # Calculate aligned length from CIGAR
        if aln.cigar_string:
            al = _cigar_aligned_length(aln.cigar_string)
        elif aln.reference_end > aln.reference_start:
            al = aln.reference_end - aln.reference_start
        else:
            al = ql

        aligned_lengths.append(al)
        total_aligned += al

        # Calculate identity from NM (edit distance) tag or CIGAR
        nm = aln.tags.get("NM", aln.tags.get("nM", None))
        if nm is not None and al > 0:
            identity = 1.0 - (int(nm) / al)
            identities.append(max(0.0, identity))
        elif aln.cigar_tuples:
            # Estimate identity from CIGAR: matches / (matches + mismatches + indels)
            matches = sum(length for op, length in aln.cigar_tuples if op in (0, 7))  # M or =
            total_ops = sum(length for op, length in aln.cigar_tuples if op in (0, 1, 2, 7, 8))
            if total_ops > 0:
                identities.append(matches / total_ops)

    mean_mapq = mapq_sum / mapped if mapped > 0 else 0.0
    mean_read_len = sum(read_lengths) / len(read_lengths) if read_lengths else 0.0
    mean_aligned_len = sum(aligned_lengths) / len(aligned_lengths) if aligned_lengths else 0.0
    mean_identity = sum(identities) / len(identities) if identities else 0.0

    coverage = 0.0
    if reference_length and reference_length > 0:
        coverage = total_aligned / reference_length

    return AlignmentStats(
        total_reads=total,
        mapped_reads=mapped,
        unmapped_reads=unmapped,
        supplementary_count=supplementary,
        secondary_count=secondary,
        mean_mapping_quality=mean_mapq,
        mean_read_length=mean_read_len,
        mean_aligned_length=mean_aligned_len,
        mean_identity=mean_identity,
        total_bases=total_bases,
        total_aligned_bases=total_aligned,
        coverage_depth=coverage,
    )
