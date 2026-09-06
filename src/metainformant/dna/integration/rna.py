"""DNA-RNA integration utilities for multi-omics analysis.

This module provides tools for integrating DNA sequence data with RNA
expression data, including gene annotation, regulatory element analysis,
and cross-omics correlation studies.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def find_open_reading_frames(dna_sequence: str) -> List[Tuple[int, int, str]]:
    """Find open reading frames in DNA sequence.

    Args:
        dna_sequence: DNA sequence to analyze

    Returns:
        List of tuples (start_pos, end_pos, amino_acid_sequence)

    Example:
        >>> dna = "ATGGCCATTGTAATGGGCC"
        >>> orfs = find_open_reading_frames(dna)
        >>> len(orfs) >= 0
        True
    """
    from ..sequence.core import find_orfs

    return find_orfs(dna_sequence)


def predict_transcription_start_sites(dna_sequence: str, window_size: int = 50) -> List[Tuple[int, float]]:
    """Predict transcription start sites using sequence motifs.

    Each distinct TATA box occurrence in the sequence is reported exactly
    once, keyed by its absolute start position. (The previous sliding-window
    implementation re-detected the same box from every window start and
    emitted many duplicate entries per box with varying scores.)

    Args:
        dna_sequence: DNA sequence to analyze
        window_size: Retained for API compatibility; scoring no longer
            depends on an arbitrary sliding-window offset.

    Returns:
        List of tuples (position, score) for potential TSS in ascending
        position order. ``position`` is the absolute box start plus the box
        length (as before). Scoring rule (deterministic, applied per box):
        1.0 when the box has at least 20 bp of upstream sequence context
        and at least 10 bp of downstream context after the box end, 0.5
        otherwise. Matches are case-insensitive.

    Example:
        >>> dna = "TTTTATAATTTTT"
        >>> tss = predict_transcription_start_sites(dna)
        >>> tss == [(7, 0.5)]
        True
    """
    tata_box = "TATA"
    sequence = dna_sequence.upper()

    potential_tss: List[Tuple[int, float]] = []
    search_start = 0
    while True:
        box_start = sequence.find(tata_box, search_start)
        if box_start == -1:
            break
        box_end = box_start + len(tata_box)
        has_upstream_room = box_start >= 20
        has_downstream_room = box_end + 10 <= len(sequence)
        score = 1.0 if has_upstream_room and has_downstream_room else 0.5
        potential_tss.append((box_end, score))
        search_start = box_start + 1

    return potential_tss


def find_transcription_factor_binding_sites(dna_sequence: str, tf_motifs: Dict[str, str]) -> Dict[str, List[int]]:
    """Find transcription factor binding sites in DNA sequence.

    Args:
        dna_sequence: DNA sequence to search
        tf_motifs: Dictionary mapping TF names to consensus sequences

    Returns:
        Dictionary mapping each motif sequence to the positions where it
        occurs in the sequence (keys are motif sequences, not TF names)

    Example:
        >>> dna = "GGGAATTTCCGGG"
        >>> motifs = {"SP1": "GGGCGG", "TBP": "TATAAA"}
        >>> sites = find_transcription_factor_binding_sites(dna, motifs)
        >>> isinstance(sites, dict)
        True
    """
    from ..sequence.motifs import find_motifs

    return find_motifs(dna_sequence, list(tf_motifs.values()))


def calculate_codon_usage_bias(dna_sequence: str) -> Dict[str, Any]:
    """Calculate codon usage bias for coding sequences.

    Args:
        dna_sequence: DNA coding sequence

    Returns:
        Dictionary with codon usage bias metrics

    Example:
        >>> dna = "ATGGCCATTGTAATGGGCC"
        >>> bias = calculate_codon_usage_bias(dna)
        >>> "cai" in bias
        True
    """
    from ..expression.codon import cai, codon_usage
    from ..expression.translation import translate_dna

    # Get codon usage
    usage = codon_usage(dna_sequence)

    # Calculate CAI
    cai_value = cai(dna_sequence)

    # Translate to protein
    protein = translate_dna(dna_sequence)

    return {
        "codon_usage": usage,
        "cai": cai_value,
        "protein_length": len(protein) if protein else 0,
        "coding_sequence_length": len(dna_sequence),
    }


def analyze_gene_structure(dna_sequence: str) -> Dict[str, Any]:
    """Analyze gene structure including exons, introns, and regulatory elements.

    Args:
        dna_sequence: DNA sequence containing a gene

    Returns:
        Dictionary with gene structure analysis

    Example:
        >>> dna = "ATGGCCATTGTAATGGGCC"
        >>> structure = analyze_gene_structure(dna)
        >>> "orfs" in structure
        True
    """
    # Find ORFs
    orfs = find_open_reading_frames(dna_sequence)

    # Find TSS
    tss_sites = predict_transcription_start_sites(dna_sequence)

    # Basic promoter analysis
    promoter_region = dna_sequence[:100] if len(dna_sequence) > 100 else dna_sequence
    gc_content_promoter: float = promoter_region.count("G") + promoter_region.count("C")
    gc_content_promoter = gc_content_promoter / len(promoter_region) if promoter_region else 0

    return {
        "orfs": orfs,
        "transcription_start_sites": tss_sites,
        "promoter_gc_content": gc_content_promoter,
        "sequence_length": len(dna_sequence),
        "potential_coding_regions": len(orfs),
    }


def _gene_gc_content(feature: Any) -> Optional[float]:
    """Extract per-gene GC content from a per-gene DNA feature entry.

    Accepts a DNA sequence string, a feature dict with a "gc_content" key,
    or a feature dict with a "sequence" key. Returns None when no usable
    GC value can be derived.
    """
    if isinstance(feature, str):
        sequence = feature.upper()
        if not sequence:
            return None
        return (sequence.count("G") + sequence.count("C")) / len(sequence)
    if isinstance(feature, dict):
        if "gc_content" in feature:
            return float(feature["gc_content"])
        nested = feature.get("sequence")
        if isinstance(nested, str) and nested:
            nested_upper = nested.upper()
            return (nested_upper.count("G") + nested_upper.count("C")) / len(nested_upper)
    return None


def correlate_dna_with_rna_expression(
    dna_features: Dict[str, Any], rna_expression: Dict[str, float]
) -> Dict[str, float]:
    """Correlate DNA sequence features with RNA expression levels.

    Args:
        dna_features: Dictionary mapping gene IDs to per-gene DNA features.
            Each value may be a feature dict (e.g. ``{"gc_content": 0.6}``),
            a DNA sequence string, or a feature dict with a "sequence" key;
            per-gene GC content is derived accordingly.
        rna_expression: Dictionary mapping gene IDs to expression levels

    Returns:
        Dictionary with correlation results. ``"gc_expression"`` is the
        Pearson correlation between per-gene GC content and expression
        levels over the genes present in both inputs. Returns an empty
        dict when either input is empty or fewer than two shared genes
        yield a usable GC value.

    Example:
        >>> dna_feat = {"gene1": {"gc_content": 0.6}, "gene2": {"gc_content": 0.4}}
        >>> rna_expr = {"gene1": 100.0, "gene2": 50.0}
        >>> correlation = correlate_dna_with_rna_expression(dna_feat, rna_expr)
        >>> isinstance(correlation, dict)
        True
    """
    correlations: Dict[str, float] = {}

    if not dna_features or not rna_expression:
        return correlations

    gc_values: List[float] = []
    expr_values: List[float] = []
    for gene, expression in rna_expression.items():
        if gene not in dna_features:
            continue
        gc = _gene_gc_content(dna_features[gene])
        if gc is None:
            continue
        gc_values.append(gc)
        expr_values.append(expression)

    if len(gc_values) > 1:
        correlations["gc_expression"] = calculate_correlation(gc_values, expr_values)

    return correlations


def calculate_correlation(x: List[float], y: List[float]) -> float:
    """Calculate Pearson correlation coefficient."""
    if len(x) != len(y) or len(x) < 2:
        return 0.0

    n = len(x)
    sum_x = sum(x)
    sum_y = sum(y)
    sum_xy = sum(xi * yi for xi, yi in zip(x, y))
    sum_x2 = sum(xi**2 for xi in x)
    sum_y2 = sum(yi**2 for yi in y)

    numerator = n * sum_xy - sum_x * sum_y
    denominator = ((n * sum_x2 - sum_x**2) * (n * sum_y2 - sum_y**2)) ** 0.5

    return numerator / denominator if denominator != 0 else 0.0


def predict_splice_sites(dna_sequence: str) -> Dict[str, Any]:
    """Predict splice donor and acceptor sites in DNA sequence.

    Args:
        dna_sequence: DNA sequence to analyze

    Returns:
        Dictionary with donor and acceptor site positions

    Example:
        >>> dna = "GTATCGTAG"
        >>> sites = predict_splice_sites(dna)
        >>> "donor" in sites and "acceptor" in sites
        True
    """
    # Common splice site motifs
    donor_motif = "GT"  # 5' splice site (donor)
    acceptor_motif = "AG"  # 3' splice site (acceptor)

    donor_sites = []
    acceptor_sites = []

    # Find donor sites (GT dinucleotides)
    for i in range(len(dna_sequence) - 1):
        if dna_sequence[i : i + 2].upper() == donor_motif:
            donor_sites.append(i)

    # Find acceptor sites (AG dinucleotides)
    for i in range(len(dna_sequence) - 1):
        if dna_sequence[i : i + 2].upper() == acceptor_motif:
            acceptor_sites.append(i)

    return {
        "donor": donor_sites,
        "acceptor": acceptor_sites,
        "potential_introns": min(len(donor_sites), len(acceptor_sites)),
    }


def analyze_regulatory_elements(dna_sequence: str) -> Dict[str, List[Tuple[int, str]]]:
    """Analyze regulatory elements in DNA sequence.

    Args:
        dna_sequence: DNA sequence to analyze

    Returns:
        Dictionary mapping element types to (position, sequence) tuples

    Example:
        >>> dna = "TATAATCGGGGCGGR"
        >>> elements = analyze_regulatory_elements(dna)
        >>> isinstance(elements, dict)
        True
    """
    elements: dict[str, list[tuple[int, str]]] = {"tata_box": [], "caat_box": [], "gc_box": [], "enhancer_motifs": []}

    # TATA box
    tata_positions = []
    for match in find_motif_positions(dna_sequence, "TATA"):
        tata_positions.append((match, dna_sequence[match : match + 4]))
    elements["tata_box"] = tata_positions

    # CAAT box
    caat_positions = []
    for match in find_motif_positions(dna_sequence, "CAAT"):
        caat_positions.append((match, dna_sequence[match : match + 4]))
    elements["caat_box"] = caat_positions

    # GC box (SP1 binding site)
    gc_positions = []
    for match in find_motif_positions(dna_sequence, "GGGCGG"):
        gc_positions.append((match, dna_sequence[match : match + 6]))
    elements["gc_box"] = gc_positions

    return elements


def find_motif_positions(sequence: str, motif: str) -> List[int]:
    """Find motif positions in sequence (helper function)."""
    positions = []
    motif_len = len(motif)
    for i in range(len(sequence) - motif_len + 1):
        if sequence[i : i + motif_len].upper() == motif.upper():
            positions.append(i)
    return positions


def integrate_dna_rna_data(dna_data: Dict[str, Any], rna_data: Dict[str, Any]) -> Dict[str, Any]:
    """Integrate DNA sequence data with RNA expression data.

    Args:
        dna_data: DNA sequence analysis results
        rna_data: RNA expression analysis results

    Returns:
        Integrated analysis results

    Example:
        >>> dna = {"gc_content": 0.6, "orf_count": 2}
        >>> rna = {"expression_level": 100.0, "transcript_length": 1000}
        >>> integrated = integrate_dna_rna_data(dna, rna)
        >>> isinstance(integrated, dict)
        True
    """
    integrated = {"dna_features": dna_data, "rna_features": rna_data, "integration_metrics": {}}

    # Calculate integration metrics
    if "gc_content" in dna_data and "expression_level" in rna_data:
        # GC content vs expression correlation
        integrated["integration_metrics"]["gc_expression_correlation"] = 0.0  # Would need more data

    if "orf_count" in dna_data and "transcript_length" in rna_data:
        # Coding density
        coding_density = (
            dna_data["orf_count"] / rna_data["transcript_length"] if rna_data["transcript_length"] > 0 else 0
        )
        integrated["integration_metrics"]["coding_density"] = coding_density

    return integrated


def predict_gene_function_from_sequence(dna_sequence: str) -> Dict[str, Any]:
    """Predict gene function based on DNA sequence features.

    Args:
        dna_sequence: DNA sequence of gene

    Returns:
        Dictionary with function predictions

    Example:
        >>> dna = "ATGGCCATTGTAATGGGCC"
        >>> function = predict_gene_function_from_sequence(dna)
        >>> "protein_features" in function
        True
    """
    # Neutral result for empty/whitespace-only input (same shape as below)
    if not dna_sequence.strip():
        return {
            "gc_content": 0.0,
            "orf_count": 0,
            "protein_features": {
                "codon_usage": {},
                "cai": 0.0,
                "protein_length": 0,
                "coding_sequence_length": len(dna_sequence),
            },
            "regulatory_elements": {"tata_box": [], "caat_box": [], "gc_box": [], "enhancer_motifs": []},
            "predicted_function": "unknown",
        }

    # Analyze sequence features (case-insensitive GC counting)
    upper_sequence = dna_sequence.upper()
    gc_content = (upper_sequence.count("G") + upper_sequence.count("C")) / len(dna_sequence)

    # Find ORFs
    orfs = find_open_reading_frames(dna_sequence)

    # Codon analysis
    codon_analysis = calculate_codon_usage_bias(dna_sequence)

    # Regulatory elements
    regulatory = analyze_regulatory_elements(dna_sequence)

    return {
        "gc_content": gc_content,
        "orf_count": len(orfs),
        "protein_features": codon_analysis,
        "regulatory_elements": regulatory,
        "predicted_function": "unknown",  # Would need ML model for real prediction
    }
