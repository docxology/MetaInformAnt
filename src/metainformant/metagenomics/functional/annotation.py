"""Functional annotation of metagenomic sequences.

Implements HMM-based functional annotation, open reading frame (ORF) prediction,
and gene family classification. Supports COG, KEGG, and Pfam-style databases.

Annotation approach:
1. Predict ORFs from assembled contigs or reads using start/stop codon scanning.
2. Translate predicted ORFs to protein sequences.
3. Compare protein sequences against functional databases using profile HMM scoring.
4. Assign gene family classifications based on best-hit or domain architecture.
"""

from __future__ import annotations

import math
import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from functools import lru_cache

from Bio.Data import CodonTable

from metainformant.core.sequence import reverse_complement
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

_ENV_PREFIX = "META_"

# Genetic code (Standard, NCBI Table 1)
_CODON_TABLE: dict[str, str] = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

# Start codons used for ORF calling: the NCBI table's starts restricted to the canonical
# ATG/GTG/TTG. The rare alternative starts of table 11 (CTG, ATT, ATC, ATA) are excluded
# because they would add large numbers of spurious nested ORFs.
_CANONICAL_STARTS = frozenset({"ATG", "GTG", "TTG"})


@lru_cache(maxsize=None)
def _genetic_code(table_id: int) -> tuple[dict[str, str], frozenset[str], frozenset[str]]:
    """Return ``(codon -> amino acid, start codons, stop codons)`` for an NCBI genetic code.

    Args:
        table_id: NCBI translation table id (1 = standard, 11 = bacterial, 4 = Mycoplasma, ...).

    Returns:
        Codon table (stop codons map to ``"*"``), usable start codons and stop codons.

    Raises:
        ValueError: If the table id is not an NCBI genetic code.
    """
    try:
        table = CodonTable.unambiguous_dna_by_id[table_id]
    except (KeyError, TypeError):
        supported = sorted(CodonTable.unambiguous_dna_by_id)
        raise ValueError(f"Unsupported translation_table {table_id!r}; NCBI table ids: {supported}") from None
    codons: dict[str, str] = dict(table.forward_table)
    stops = frozenset(table.stop_codons)
    for stop in stops:
        codons[stop] = "*"
    all_starts = frozenset(table.start_codons)
    return codons, (all_starts & _CANONICAL_STARTS) or all_starts, stops


@dataclass
class ORF:
    """Predicted open reading frame."""

    sequence_id: str
    start: int
    end: int
    strand: str  # '+' or '-'
    nucleotide_seq: str
    protein_seq: str
    frame: int  # 0, 1, or 2
    partial: bool = False  # True if ORF runs to edge of contig

    @property
    def length_nt(self) -> int:
        """Nucleotide length."""
        return len(self.nucleotide_seq)

    @property
    def length_aa(self) -> int:
        """Amino acid length."""
        return len(self.protein_seq)


@dataclass
class HMMHit:
    """Result from HMM-based database search.

    ``e_value``, ``bias`` and ``domain_e_value`` are NaN: the position-specific scoring
    matrices used here have no calibrated null model, so no E-value (or composition
    bias) is computed. ``query_start``/``query_end`` (0-based, end exclusive) bound the
    part of the query protein covered by the best ungapped placement of the profile, and
    ``target_start``/``target_end`` the profile positions scored in it.
    """

    query_id: str
    target_id: str
    target_name: str
    score: float
    e_value: float
    bias: float
    domain_score: float
    domain_e_value: float
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    description: str = ""


@dataclass
class FunctionalAnnotation:
    """Functional annotation for a gene/ORF.

    ``profile_coverage`` and ``consensus_identity`` describe the best hit (``best_hit``) and
    are None while there is no hit:

    - ``profile_coverage``: fraction of the profile's positions that the protein covers in the
      best ungapped placement (1.0 when the whole profile is aligned, less for a protein that
      is truncated at an end of the profile).
    - ``consensus_identity``: fraction of the aligned positions where the protein carries the
      profile's best-scoring residue (the profile consensus).

    Neither is a probability or a significance; there is no calibrated null model here.
    """

    orf_id: str
    gene_families: list[str] = field(default_factory=list)
    kegg_orthologs: list[str] = field(default_factory=list)
    cog_categories: list[str] = field(default_factory=list)
    pfam_domains: list[str] = field(default_factory=list)
    ec_numbers: list[str] = field(default_factory=list)
    best_hit: HMMHit | None = None
    description: str = ""
    profile_coverage: float | None = None
    consensus_identity: float | None = None


def _translate(
    nucleotide_seq: str,
    stop_at_stop: bool = False,
    codon_table: Mapping[str, str] | None = None,
) -> str:
    """Translate nucleotide sequence to protein.

    Args:
        nucleotide_seq: DNA sequence (must be uppercase).
        stop_at_stop: If True, stop translation at first stop codon.
        codon_table: Codon -> amino acid mapping (stops as ``"*"``); the standard
            code when None.

    Returns:
        Amino acid sequence.
    """
    table = _CODON_TABLE if codon_table is None else codon_table
    protein = []
    seq = nucleotide_seq.upper()
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        if len(codon) < 3:
            break
        aa = table.get(codon, "X")
        if aa == "*" and stop_at_stop:
            break
        if aa != "*":
            protein.append(aa)
    return "".join(protein)


def predict_orfs(
    sequence: str,
    sequence_id: str = "seq",
    min_length: int = 100,
    allow_partial: bool = True,
    translation_table: int = 11,
) -> list[ORF]:
    """Predict open reading frames in a nucleotide sequence.

    Scans all six reading frames (3 forward, 3 reverse) for ORFs defined
    by start and stop codons. Supports both complete ORFs and partial
    ORFs that extend to contig edges (common in metagenomic assemblies).

    Args:
        sequence: Input nucleotide sequence.
        sequence_id: Identifier for the source sequence.
        min_length: Minimum ORF length in nucleotides (default 100).
        allow_partial: Allow ORFs that extend to sequence edges.
        translation_table: NCBI genetic code table id (default 11 = Bacterial). It sets
            the amino-acid assignments and the stop codons (for example TGA is
            tryptophan, not a stop, in table 4) and, through the table's start codons,
            which of ATG/GTG/TTG can open an ORF.

    Returns:
        List of predicted ORFs sorted by position. Each (strand, start, end) appears
        once, even when the sequence begins with a start codon.

    Raises:
        ValueError: If ``translation_table`` is not an NCBI genetic code id.
    """
    codon_table, start_codons, stop_codons = _genetic_code(translation_table)
    seq_upper = sequence.upper().replace(" ", "").replace("\n", "")
    seq_len = len(seq_upper)
    orfs: list[ORF] = []

    if seq_len < min_length:
        logger.debug(f"Sequence too short ({seq_len} bp) for ORF prediction with min_length={min_length}")
        return orfs

    for strand, working_seq in [("+", seq_upper), ("-", reverse_complement(seq_upper))]:
        for frame in range(3):
            current_starts: list[int] = []
            # If allowing partial ORFs, treat position 0 as a potential start
            if allow_partial and frame == 0:
                current_starts.append(0)

            pos = frame
            while pos + 2 < seq_len:
                codon = working_seq[pos : pos + 3]

                if codon in start_codons:
                    # Position 0 may already be there as the implicit partial start
                    if pos not in current_starts:
                        current_starts.append(pos)

                elif codon in stop_codons:
                    for start_pos in current_starts:
                        orf_seq = working_seq[start_pos : pos + 3]
                        if len(orf_seq) >= min_length:
                            protein = _translate(orf_seq, stop_at_stop=True, codon_table=codon_table)
                            if strand == "+":
                                genome_start = start_pos
                                genome_end = pos + 3
                            else:
                                genome_start = seq_len - (pos + 3)
                                genome_end = seq_len - start_pos

                            orfs.append(
                                ORF(
                                    sequence_id=sequence_id,
                                    start=genome_start,
                                    end=genome_end,
                                    strand=strand,
                                    nucleotide_seq=orf_seq,
                                    protein_seq=protein,
                                    frame=frame,
                                    partial=start_pos == 0 and allow_partial,
                                )
                            )
                    current_starts = []

                pos += 3

            # Handle partial ORFs at sequence end
            if allow_partial and current_starts:
                for start_pos in current_starts:
                    orf_seq = working_seq[start_pos:pos]
                    if len(orf_seq) >= min_length:
                        protein = _translate(orf_seq, stop_at_stop=False, codon_table=codon_table)
                        if strand == "+":
                            genome_start = start_pos
                            genome_end = pos
                        else:
                            genome_start = seq_len - pos
                            genome_end = seq_len - start_pos

                        orfs.append(
                            ORF(
                                sequence_id=sequence_id,
                                start=genome_start,
                                end=genome_end,
                                strand=strand,
                                nucleotide_seq=orf_seq,
                                protein_seq=protein,
                                frame=frame,
                                partial=True,
                            )
                        )

    orfs.sort(key=lambda o: (o.sequence_id, o.start))
    logger.info(f"Predicted {len(orfs)} ORFs from sequence '{sequence_id}' ({seq_len} bp)")
    return orfs


@dataclass(frozen=True)
class _ProfileAlignment:
    """Best ungapped placement of a position-specific scoring matrix on a protein.

    Attributes:
        score: Mean score per scored position.
        offset: Protein index aligned with profile position 0.
        scored_positions: Number of profile positions scored in the window.
        consensus_matches: Scored positions where the residue is the column's best-scoring one.
        profile_positions: Number of positions the profile defines.
        first_position: First scored profile position.
        last_position: Last scored profile position.
    """

    score: float
    offset: int
    scored_positions: int
    consensus_matches: int
    profile_positions: int
    first_position: int
    last_position: int

    @property
    def coverage(self) -> float:
        """Fraction of the profile's positions that were scored."""
        return self.scored_positions / self.profile_positions

    @property
    def identity(self) -> float:
        """Fraction of the scored positions that carry the profile consensus residue."""
        return self.consensus_matches / self.scored_positions


def _align_profile(
    query_profile: dict[int, dict[str, float]],
    target_seq: str,
) -> _ProfileAlignment | None:
    """Place a simple position-specific scoring matrix on a protein without gaps.

    Uses a simplified HMM-like scoring: the mean log-odds score over the profile positions
    that fall on the protein, maximised over the offsets of the profile along the protein.
    Real HMM implementations (HMMER) use full profile HMMs with match/insert/delete states.
    A residue a profile column does not list scores -1.0.

    Args:
        query_profile: Position-specific scoring matrix {pos: {aa: score}}.
        target_seq: Target amino acid sequence.

    Returns:
        The best placement (the earliest offset on ties), or None when the profile or the
        sequence is empty and so nothing can be aligned.
    """
    if not query_profile or not target_seq:
        return None

    profile_len = max(query_profile.keys()) + 1
    seq_len = len(target_seq)

    best: _ProfileAlignment | None = None
    for offset in range(max(1, seq_len - profile_len + 1)):
        score = 0.0
        matched = 0
        consensus = 0
        first_position = -1
        last_position = -1
        for pos in range(min(profile_len, seq_len - offset)):
            column = query_profile.get(pos)
            if column is None:
                continue
            aa = target_seq[offset + pos]
            column_score = column.get(aa, -1.0)
            score += column_score
            matched += 1
            if aa in column and column_score >= max(column.values()):
                consensus += 1
            if first_position < 0:
                first_position = pos
            last_position = pos
        if matched > 0 and (best is None or score / matched > best.score):
            best = _ProfileAlignment(
                score=score / matched,
                offset=offset,
                scored_positions=matched,
                consensus_matches=consensus,
                profile_positions=len(query_profile),
                first_position=first_position,
                last_position=last_position,
            )

    return best


def annotate_genes(
    sequences: dict[str, str],
    hmm_db: dict[str, dict[int, dict[str, float]]] | None = None,
    e_value_threshold: float | None = None,
    min_score: float = 25.0,
) -> list[FunctionalAnnotation]:
    """Annotate gene sequences using HMM-based database search.

    Compares protein sequences against a database of profile HMMs
    (position-specific scoring matrices) and returns functional annotations
    based on significant hits.

    Args:
        sequences: Dict mapping gene/ORF IDs to protein sequences.
        hmm_db: HMM database as {profile_id: {pos: {aa: score}}}.
            If None, only basic sequence statistics are computed.
        e_value_threshold: Must be None. The position-specific scoring used here has no
            calibrated null model, so E-values cannot be computed and hits are filtered
            by ``min_score`` only; a numeric threshold would be silently ignored, so it
            is rejected instead.
        min_score: Minimum alignment score threshold (mean score per matched position).

    Returns:
        List of FunctionalAnnotation objects. ``best_hit.e_value`` is NaN (not computed);
        ``profile_coverage`` and ``consensus_identity`` describe the best hit and are None
        for a gene without one.

    Raises:
        ValueError: If ``e_value_threshold`` is given.
    """
    if e_value_threshold is not None:
        raise ValueError(
            "e_value_threshold is not supported: E-values cannot be computed for position-specific "
            "scoring matrices without a calibrated null model. Filter hits with min_score instead."
        )

    annotations: list[FunctionalAnnotation] = []

    for seq_id, protein_seq in sequences.items():
        annotation = FunctionalAnnotation(orf_id=seq_id)

        if hmm_db is not None:
            hits: list[tuple[str, float, _ProfileAlignment]] = []
            for profile_id, profile in hmm_db.items():
                alignment = _align_profile(profile, protein_seq)
                if alignment is not None and alignment.score >= min_score:
                    hits.append((profile_id, alignment.score, alignment))

            hits.sort(key=lambda h: h[1], reverse=True)

            for profile_id, _, _ in hits:
                # Parse family/database from profile ID convention
                if profile_id.startswith("KO:"):
                    annotation.kegg_orthologs.append(profile_id.removeprefix("KO:"))
                elif profile_id.startswith("COG:"):
                    annotation.cog_categories.append(profile_id.removeprefix("COG:"))
                elif profile_id.startswith("PF"):
                    annotation.pfam_domains.append(profile_id)
                else:
                    annotation.gene_families.append(profile_id)

            if hits:
                best_id, best_score, best_alignment = hits[0]
                annotation.best_hit = HMMHit(
                    query_id=seq_id,
                    target_id=best_id,
                    target_name=best_id,
                    score=best_score,
                    e_value=math.nan,
                    bias=math.nan,
                    domain_score=best_score,
                    domain_e_value=math.nan,
                    query_start=best_alignment.offset + best_alignment.first_position,
                    query_end=best_alignment.offset + best_alignment.last_position + 1,
                    target_start=best_alignment.first_position,
                    target_end=best_alignment.last_position + 1,
                )
                annotation.profile_coverage = best_alignment.coverage
                annotation.consensus_identity = best_alignment.identity

        # Extract EC numbers from description patterns (if embedded in IDs)
        ec_pattern = re.compile(r"EC:(\d+\.\d+\.\d+\.\d+)")
        for family in annotation.gene_families:
            ec_match = ec_pattern.search(family)
            if ec_match:
                annotation.ec_numbers.append(ec_match.group(1))

        annotations.append(annotation)

    logger.info(
        f"Annotated {len(annotations)} sequences; "
        f"{sum(1 for a in annotations if a.best_hit is not None)} with database hits"
    )
    return annotations


def classify_gene_families(
    genes: dict[str, str],
    database: str = "COG",
    reference_profiles: dict[str, dict[int, dict[str, float]]] | None = None,
    min_score: float = 20.0,
) -> dict[str, list[str]]:
    """Classify genes into functional families.

    Assigns each gene to one or more functional families based on
    profile HMM similarity. Supports COG, KEGG, and Pfam databases.

    Args:
        genes: Dict mapping gene IDs to protein sequences.
        database: Target database ('COG', 'KEGG', 'Pfam').
        reference_profiles: Profile database used to assign families. A family can only
            be assigned by comparison with reference profiles, so when this is None
            no family is invented: every gene maps to an empty list (unannotated).
        min_score: Minimum score for family assignment.

    Returns:
        Dict mapping gene IDs to lists of family assignments (empty = unannotated).
    """
    prefix_map = {"COG": "COG:", "KEGG": "KO:", "Pfam": "PF"}
    prefix = prefix_map.get(database, "")

    classifications: dict[str, list[str]] = {}

    if reference_profiles is not None:
        # Filter profiles by database prefix
        db_profiles = {pid: profile for pid, profile in reference_profiles.items() if pid.startswith(prefix)}

        for gene_id, protein_seq in genes.items():
            families: list[str] = []
            for profile_id, profile in db_profiles.items():
                alignment = _align_profile(profile, protein_seq)
                if alignment is not None and alignment.score >= min_score:
                    families.append(profile_id)
            classifications[gene_id] = families
    else:
        # Amino-acid composition cannot identify a COG/KEGG/Pfam family, so without reference
        # profiles every gene is reported as unannotated rather than given an invented label.
        logger.warning(
            f"No reference_profiles supplied: {len(genes)} genes left unannotated for {database} "
            "(composition alone cannot assign a family)"
        )
        for gene_id in genes:
            classifications[gene_id] = []

    assigned = sum(1 for fams in classifications.values() if fams)
    logger.info(f"Classified {assigned}/{len(genes)} genes into {database} families")
    return classifications


__all__ = [
    "ORF",
    "HMMHit",
    "FunctionalAnnotation",
    "predict_orfs",
    "annotate_genes",
    "classify_gene_families",
]
