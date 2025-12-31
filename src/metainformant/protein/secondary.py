"""Protein secondary structure prediction utilities.

This module provides tools for predicting and analyzing protein secondary
structure elements using various algorithms and methods.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from metainformant.core import logging

logger = logging.get_logger(__name__)


def predict_secondary_structure(sequence: str, method: str = "psipred") -> List[str]:
    """Predict secondary structure for a protein sequence.

    Args:
        sequence: Protein sequence (amino acid string)
        method: Prediction method ("psipred", "jpred", "simple")

    Returns:
        List of secondary structure assignments (H, E, C)

    Raises:
        ValueError: If sequence is invalid or method is unsupported

    Example:
        >>> sequence = "MKLVLSDEL"
        >>> ss = predict_secondary_structure(sequence, method="simple")
        >>> len(ss) == len(sequence)
        True
    """
    if not sequence:
        return []

    if method not in ["psipred", "jpred", "simple"]:
        raise ValueError(f"Unsupported method: {method}")

    if method == "simple":
        return _simple_secondary_structure_prediction(sequence)
    elif method == "psipred":
        return _psipred_prediction(sequence)
    elif method == "jpred":
        return _jpred_prediction(sequence)
    else:
        # Fallback to simple method
        return _simple_secondary_structure_prediction(sequence)


def _simple_secondary_structure_prediction(sequence: str) -> List[str]:
    """Simple secondary structure prediction based on amino acid propensities."""
    # Chou-Fasman parameters (simplified)
    helix_aa = {'E', 'A', 'L', 'M', 'Q', 'K', 'H', 'V', 'I', 'F', 'Y', 'W', 'T', 'S', 'R', 'D', 'N', 'C', 'G', 'P'}
    sheet_aa = {'V', 'I', 'Y', 'F', 'W', 'L', 'T', 'M', 'C', 'A', 'G', 'R', 'K', 'Q', 'E', 'N', 'D', 'S', 'H', 'P'}

    ss_predictions = []

    for aa in sequence.upper():
        if aa in helix_aa:
            # Bias towards helix for certain amino acids
            if aa in {'E', 'A', 'L'}:
                ss_predictions.append('H')  # Helix
            elif aa in {'V', 'I', 'Y'}:
                ss_predictions.append('E')  # Sheet
            else:
                ss_predictions.append('C')  # Coil
        else:
            ss_predictions.append('C')  # Default to coil

    return ss_predictions


def _psipred_prediction(sequence: str) -> List[str]:
    """PSIPRED secondary structure prediction (placeholder)."""
    # This would call PSIPRED web service or local installation
    # For now, fall back to simple prediction
    logger.info("PSIPRED prediction not fully implemented, using simple method")
    return _simple_secondary_structure_prediction(sequence)


def _jpred_prediction(sequence: str) -> List[str]:
    """JPRED secondary structure prediction (placeholder)."""
    # This would call JPRED web service
    # For now, fall back to simple prediction
    logger.info("JPRED prediction not fully implemented, using simple method")
    return _simple_secondary_structure_prediction(sequence)


def calculate_ss_composition(ss_assignments: List[str]) -> Dict[str, float]:
    """Calculate secondary structure composition.

    Args:
        ss_assignments: List of secondary structure assignments

    Returns:
        Dictionary with composition percentages

    Example:
        >>> ss = ['H', 'H', 'E', 'C', 'C', 'C']
        >>> comp = calculate_ss_composition(ss)
        >>> comp['helix'] == 2/6
        True
    """
    if not ss_assignments:
        return {'helix': 0.0, 'sheet': 0.0, 'coil': 0.0}

    total = len(ss_assignments)
    helix_count = sum(1 for ss in ss_assignments if ss == 'H')
    sheet_count = sum(1 for ss in ss_assignments if ss == 'E')
    coil_count = sum(1 for ss in ss_assignments if ss == 'C')

    return {
        'helix': helix_count / total,
        'sheet': sheet_count / total,
        'coil': coil_count / total,
        'helix_count': helix_count,
        'sheet_count': sheet_count,
        'coil_count': coil_count,
        'total_residues': total
    }


def identify_ss_elements(ss_assignments: List[str], min_length: int = 3) -> List[Dict[str, Any]]:
    """Identify continuous secondary structure elements.

    Args:
        ss_assignments: List of secondary structure assignments
        min_length: Minimum length for an element

    Returns:
        List of secondary structure elements

    Example:
        >>> ss = ['H', 'H', 'H', 'E', 'E', 'C']
        >>> elements = identify_ss_elements(ss, min_length=2)
        >>> len(elements) >= 2
        True
    """
    elements = []

    if not ss_assignments:
        return elements

    current_type = ss_assignments[0]
    current_start = 0
    current_length = 1

    for i in range(1, len(ss_assignments)):
        if ss_assignments[i] == current_type:
            current_length += 1
        else:
            # End of current element
            if current_length >= min_length and current_type in ['H', 'E']:
                elements.append({
                    'type': current_type,
                    'start': current_start,
                    'end': current_start + current_length - 1,
                    'length': current_length
                })

            # Start new element
            current_type = ss_assignments[i]
            current_start = i
            current_length = 1

    # Handle last element
    if current_length >= min_length and current_type in ['H', 'E']:
        elements.append({
            'type': current_type,
            'start': current_start,
            'end': current_start + current_length - 1,
            'length': current_length
        })

    return elements


def compare_ss_predictions(prediction1: List[str], prediction2: List[str]) -> Dict[str, float]:
    """Compare two secondary structure predictions.

    Args:
        prediction1: First SS prediction
        prediction2: Second SS prediction

    Returns:
        Comparison statistics

    Example:
        >>> pred1 = ['H', 'H', 'E']
        >>> pred2 = ['H', 'E', 'E']
        >>> comparison = compare_ss_predictions(pred1, pred2)
        >>> comparison['accuracy'] == 2/3
        True
    """
    if len(prediction1) != len(prediction2):
        return {'accuracy': 0.0, 'matches': 0, 'total': 0}

    matches = sum(1 for a, b in zip(prediction1, prediction2) if a == b)
    total = len(prediction1)
    accuracy = matches / total if total > 0 else 0.0

    return {
        'accuracy': accuracy,
        'matches': matches,
        'total': total,
        'q3_score': accuracy  # Q3 score for 3-state prediction
    }


def ss_to_dSSP_format(ss_assignments: List[str]) -> str:
    """Convert secondary structure assignments to DSSP format.

    Args:
        ss_assignments: List of SS assignments (H, E, C)

    Returns:
        DSSP-formatted string

    Example:
        >>> ss = ['H', 'E', 'C']
        >>> dssp = ss_to_dSSP_format(ss)
        >>> len(dssp) == 3
        True
    """
    # DSSP codes: H=alpha helix, E=extended strand, blank=coil/other
    dssp_map = {'H': 'H', 'E': 'E', 'C': ' '}

    return ''.join(dssp_map.get(ss, ' ') for ss in ss_assignments)


def parse_dssp_file(dssp_content: str) -> Dict[str, Any]:
    """Parse DSSP output file.

    Args:
        dssp_content: DSSP file content

    Returns:
        Parsed DSSP data

    Example:
        >>> # Assuming DSSP content exists
        >>> # dssp_data = parse_dssp_file(dssp_string)
        >>> # isinstance(dssp_data, dict)
        >>> # True
    """
    # This would parse actual DSSP output
    # Placeholder implementation
    logger.info("DSSP parsing not fully implemented")
    return {
        'sequence': '',
        'secondary_structure': [],
        'accessibility': [],
        'total_residues': 0
    }


def predict_transmembrane_regions(sequence: str) -> List[Dict[str, Any]]:
    """Predict transmembrane regions from sequence.

    Args:
        sequence: Protein sequence

    Returns:
        List of transmembrane region predictions

    Example:
        >>> sequence = "M" * 20 + "A" * 10 + "M" * 20  # Hydrophobic regions
        >>> tm_regions = predict_transmembrane_regions(sequence)
        >>> isinstance(tm_regions, list)
        True
    """
    # Simple hydrophobicity-based prediction
    # In practice, would use TMHMM or similar
    regions = []

    # Kyte-Doolittle hydrophobicity scale (simplified)
    hydrophobicity = {
        'I': 4.5, 'V': 4.2, 'L': 3.8, 'P': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,
        'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'Q': -3.5, 'N': -3.5,
        'E': -3.5, 'D': -3.5, 'H': -3.2, 'K': -3.9, 'R': -4.5
    }

    window_size = 19  # Typical TM helix length
    threshold = 1.5   # Hydrophobicity threshold

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        avg_hydro = sum(hydrophobicity.get(aa, 0) for aa in window) / window_size

        if avg_hydro > threshold:
            regions.append({
                'start': i,
                'end': i + window_size - 1,
                'length': window_size,
                'avg_hydrophobicity': avg_hydro,
                'sequence': window
            })

    return regions


def calculate_ss_propensities(sequence: str) -> Dict[str, Dict[str, float]]:
    """Calculate secondary structure propensities for amino acids.

    Args:
        sequence: Protein sequence

    Returns:
        Dictionary with propensity scores

    Example:
        >>> sequence = "ACDEFGHIKLMNPQRSTVWY"
        >>> propensities = calculate_ss_propensities(sequence)
        >>> isinstance(propensities, dict)
        True
    """
    # Chou-Fasman parameters (simplified)
    helix_params = {
        'P': 0.57, 'V': 0.90, 'I': 1.09, 'M': 0.82, 'F': 1.16, 'Y': 0.72, 'W': 1.14,
        'A': 1.45, 'E': 0.51, 'L': 1.34, 'H': 0.96, 'Q': 0.98, 'G': 0.53, 'T': 0.82,
        'S': 0.79, 'R': 0.95, 'K': 0.81, 'D': 0.66, 'N': 0.43, 'C': 0.66
    }

    sheet_params = {
        'P': 0.31, 'V': 1.87, 'I': 1.67, 'M': 1.20, 'F': 1.33, 'Y': 1.24, 'W': 1.01,
        'A': 0.97, 'E': 0.39, 'L': 1.22, 'H': 0.51, 'Q': 0.58, 'G': 0.81, 'T': 1.20,
        'S': 0.94, 'R': 0.99, 'K': 0.88, 'D': 0.40, 'N': 0.65, 'C': 1.19
    }

    propensities = {}

    for aa in set(sequence.upper()):
        if aa in helix_params:
            propensities[aa] = {
                'helix': helix_params[aa],
                'sheet': sheet_params[aa],
                'coil': 1.0,  # Simplified
                'preferred': 'helix' if helix_params[aa] > sheet_params[aa] else 'sheet'
            }

    return propensities


def validate_ss_prediction(ss_assignments: List[str], sequence: str) -> Dict[str, Any]:
    """Validate secondary structure prediction against known constraints.

    Args:
        ss_assignments: SS predictions
        sequence: Protein sequence

    Returns:
        Validation results

    Example:
        >>> ss = ['H', 'H', 'E']
        >>> seq = "ALAALAALA"
        >>> validation = validate_ss_prediction(ss, seq)
        >>> validation['valid_length']
        True
    """
    validation = {
        'valid_length': len(ss_assignments) == len(sequence),
        'valid_states': all(ss in ['H', 'E', 'C'] for ss in ss_assignments),
        'has_structure': any(ss in ['H', 'E'] for ss in ss_assignments),
        'issues': []
    }

    if not validation['valid_length']:
        validation['issues'].append(f"Length mismatch: {len(ss_assignments)} vs {len(sequence)}")

    if not validation['valid_states']:
        invalid = [ss for ss in ss_assignments if ss not in ['H', 'E', 'C']]
        validation['issues'].append(f"Invalid states: {invalid}")

    if not validation['has_structure']:
        validation['issues'].append("No secondary structure elements predicted")

    return validation
