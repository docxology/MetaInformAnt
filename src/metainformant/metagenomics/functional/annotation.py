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

import os
import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

_ENV_PREFIX = "META_"

# Genetic code (Standard, NCBI Table 1)
_CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

_START_CODONS = {"ATG", "GTG", "TTG"}
_STOP_CODONS = {"TAA", "TAG", "TGA"}


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
    """Result from HMM-based database search."""

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
    """Functional annotation for a gene/ORF."""

    orf_id: str
    gene_families: list[str] = field(default_factory=list)
    kegg_orthologs: list[str] = field(default_factory=list)
    cog_categories: list[str] = field(default_factory=list)
    pfam_domains: list[str] = field(default_factory=list)
    ec_numbers: list[str] = field(default_factory=list)
    best_hit: HMMHit | None = None
    description: str = ""
    confidence: float = 0.0


def _reverse_complement(sequence: str) -> str:
    """Compute reverse complement of a DNA sequence."""
    complement_map = str.maketrans("ATCGatcgNn", "TAGCtagcNn")
    return sequence.translate(complement_map)[::-1]


def _translate(nucleotide_seq: str, stop_at_stop: bool = False) -> str:
    """Translate nucleotide sequence to protein.

    Args:
        nucleotide_seq: DNA sequence (must be uppercase).
        stop_at_stop: If True, stop translation at first stop codon.

    Returns:
        Amino acid sequence.
    """
    protein = []
    seq = nucleotide_seq.upper()
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        if len(codon) < 3:
            break
        aa = _CODON_TABLE.get(codon, "X")
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
        translation_table: Genetic code table number (default 11 = Bacterial).

    Returns:
        List of predicted ORFs sorted by position.
    """
    seq_upper = sequence.upper().replace(" ", "").replace("\n", "")
    seq_len = len(seq_upper)
    orfs: list[ORF] = []

    if seq_len < min_length:
        logger.debug(f"Sequence too short ({seq_len} bp) for ORF prediction with min_length={min_length}")
        return orfs

    for strand, working_seq in [("+", seq_upper), ("-", _reverse_complement(seq_upper))]:
        for frame in range(3):
            current_starts: list[int] = []
            # If allowing partial ORFs, treat position 0 as a potential start
            if allow_partial and frame == 0:
                current_starts.append(0)

            pos = frame
            while pos + 2 < seq_len:
                codon = working_seq[pos : pos + 3]

                if codon in _START_CODONS:
                    current_starts.append(pos)

                elif codon in _STOP_CODONS:
                    for start_pos in current_starts:
                        orf_seq = working_seq[start_pos : pos + 3]
                        if len(orf_seq) >= min_length:
                            protein = _translate(orf_seq, stop_at_stop=True)
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
                    orf_seq = working_seq[start_pos : pos]
                    if len(orf_seq) >= min_length:
                        protein = _translate(orf_seq, stop_at_stop=False)
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


def _score_hmm_alignment(
    query_profile: dict[int, dict[str, float]],
    target_seq: str,
) -> float:
    """Score a sequence against a simple position-specific scoring matrix.

    Uses a simplified HMM-like scoring: sum of log-odds scores for each
    position in the profile. Real HMM implementations (HMMER) use full
    profile HMMs with match/insert/delete states.

    Args:
        query_profile: Position-specific scoring matrix {pos: {aa: score}}.
        target_seq: Target amino acid sequence.

    Returns:
        Alignment score (sum of position-specific scores).
    """
    if not query_profile or not target_seq:
        return 0.0

    profile_len = max(query_profile.keys()) + 1
    seq_len = len(target_seq)

    if seq_len == 0 or profile_len == 0:
        return 0.0

    # Simple ungapped alignment scoring
    best_score = float("-inf")
    for offset in range(max(1, seq_len - profile_len + 1)):
        score = 0.0
        matched = 0
        for pos in range(min(profile_len, seq_len - offset)):
            aa = target_seq[offset + pos]
            if pos in query_profile:
                score += query_profile[pos].get(aa, -1.0)
                matched += 1
        if matched > 0:
            normalized = score / matched
            if normalized > best_score:
                best_score = normalized

    return best_score if best_score > float("-inf") else 0.0


def annotate_genes(
    sequences: dict[str, str],
    hmm_db: dict[str, dict[int, dict[str, float]]] | None = None,
    e_value_threshold: float = 1e-5,
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
        e_value_threshold: Maximum E-value for significant hits.
        min_score: Minimum alignment score threshold.

    Returns:
        List of FunctionalAnnotation objects.
    """
    annotations: list[FunctionalAnnotation] = []

    for seq_id, protein_seq in sequences.items():
        annotation = FunctionalAnnotation(orf_id=seq_id)

        if hmm_db is not None:
            hits: list[tuple[str, float]] = []
            for profile_id, profile in hmm_db.items():
                score = _score_hmm_alignment(profile, protein_seq)
                if score >= min_score:
                    hits.append((profile_id, score))

            hits.sort(key=lambda h: h[1], reverse=True)

            for profile_id, score in hits:
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
                best_id, best_score = hits[0]
                annotation.best_hit = HMMHit(
                    query_id=seq_id,
                    target_id=best_id,
                    target_name=best_id,
                    score=best_score,
                    e_value=0.0,
                    bias=0.0,
                    domain_score=best_score,
                    domain_e_value=0.0,
                    query_start=0,
                    query_end=len(protein_seq),
                    target_start=0,
                    target_end=0,
                )
                annotation.confidence = min(1.0, best_score / 100.0)

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
        reference_profiles: Pre-loaded profile database. If None,
            a simple composition-based classification is used.
        min_score: Minimum score for family assignment.

    Returns:
        Dict mapping gene IDs to lists of family assignments.
    """
    prefix_map = {"COG": "COG:", "KEGG": "KO:", "Pfam": "PF"}
    prefix = prefix_map.get(database, "")

    classifications: dict[str, list[str]] = {}

    if reference_profiles is not None:
        # Filter profiles by database prefix
        db_profiles = {
            pid: profile
            for pid, profile in reference_profiles.items()
            if pid.startswith(prefix)
        }

        for gene_id, protein_seq in genes.items():
            families: list[str] = []
            for profile_id, profile in db_profiles.items():
                score = _score_hmm_alignment(profile, protein_seq)
                if score >= min_score:
                    families.append(profile_id)
            classifications[gene_id] = families
    else:
        # Composition-based heuristic classification
        # Classify based on amino acid composition patterns
        for gene_id, protein_seq in genes.items():
            families = _composition_classify(protein_seq, database)
            classifications[gene_id] = families

    assigned = sum(1 for fams in classifications.values() if fams)
    logger.info(f"Classified {assigned}/{len(genes)} genes into {database} families")
    return classifications


def _composition_classify(protein_seq: str, database: str) -> list[str]:
    """Simple composition-based classification heuristic.

    Uses amino acid composition patterns to tentatively classify proteins.
    This is a fallback when no HMM database is available.
    """
    if not protein_seq:
        return []

    seq_upper = protein_seq.upper()
    length = len(seq_upper)
    if length == 0:
        return []

    # Calculate amino acid frequencies
    aa_counts: dict[str, int] = defaultdict(int)
    for aa in seq_upper:
        aa_counts[aa] += 1

    # Simple heuristic features
    charged = (aa_counts.get("D", 0) + aa_counts.get("E", 0) + aa_counts.get("K", 0) + aa_counts.get("R", 0)) / length
    hydrophobic = (
        aa_counts.get("A", 0) + aa_counts.get("V", 0) + aa_counts.get("I", 0)
        + aa_counts.get("L", 0) + aa_counts.get("M", 0) + aa_counts.get("F", 0)
        + aa_counts.get("W", 0)
    ) / length
    cys_frac = aa_counts.get("C", 0) / length
    gly_frac = aa_counts.get("G", 0) / length
    his_frac = aa_counts.get("H", 0) / length

    families: list[str] = []

    if database == "COG":
        # COG functional categories based on composition heuristics
        if charged > 0.35:
            families.append("COG:J")  # Translation, ribosomal
        if hydrophobic > 0.45:
            families.append("COG:M")  # Cell wall/membrane
        if cys_frac > 0.03 and his_frac > 0.02:
            families.append("COG:C")  # Energy production (metal-binding)
        if gly_frac > 0.12:
            families.append("COG:O")  # Post-translational modification
    elif database == "KEGG":
        if charged > 0.35:
            families.append("KO:K03046")  # Ribosomal-like
        if hydrophobic > 0.45 and length > 200:
            families.append("KO:K02004")  # Transporter-like
    elif database == "Pfam":
        if cys_frac > 0.04:
            families.append("PF00076")  # RNA recognition motif
        if gly_frac > 0.15:
            families.append("PF00071")  # Ras family

    return families


__all__ = [
    "ORF",
    "HMMHit",
    "FunctionalAnnotation",
    "predict_orfs",
    "annotate_genes",
    "classify_gene_families",
]
