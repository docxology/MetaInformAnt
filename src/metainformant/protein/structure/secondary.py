"""Protein secondary structure prediction utilities.

This module provides tools for predicting and analyzing protein secondary
structure elements from sequence propensities, and for reading DSSP output.
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple, Union

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

#: Chou-Fasman conformational parameters (Pa, Pb, Pt) of the twenty amino acids: the propensity of a residue for an
#: alpha-helix, a beta-sheet and a turn (Chou & Fasman, Adv Enzymol Relat Areas Mol Biol 47:45-148, 1978). A value
#: above 1 favours the conformation. This is the one table used by every function of this module.
CHOU_FASMAN_PARAMETERS: Dict[str, Tuple[float, float, float]] = {
    "A": (1.42, 0.83, 0.66),
    "R": (0.98, 0.93, 0.95),
    "N": (0.67, 0.89, 1.56),
    "D": (1.01, 0.54, 1.46),
    "C": (0.70, 1.19, 1.19),
    "Q": (1.11, 1.10, 0.98),
    "E": (1.51, 0.37, 0.74),
    "G": (0.57, 0.75, 1.56),
    "H": (1.00, 0.87, 0.95),
    "I": (1.08, 1.60, 0.47),
    "L": (1.21, 1.30, 0.59),
    "K": (1.16, 0.74, 1.01),
    "M": (1.45, 1.05, 0.60),
    "F": (1.13, 1.38, 0.60),
    "P": (0.57, 0.55, 1.52),
    "S": (0.77, 0.75, 1.43),
    "T": (0.83, 1.19, 0.96),
    "W": (1.08, 1.37, 0.96),
    "Y": (0.69, 1.47, 1.14),
    "V": (1.06, 1.70, 0.50),
}

#: Methods of :func:`predict_secondary_structure` that run here, from the sequence alone
SEQUENCE_METHODS: Tuple[str, ...] = ("window", "simple")

#: Predictors that need an external service or program; they are recognised only to be refused
EXTERNAL_METHODS: Tuple[str, ...] = ("psipred", "jpred")

_WINDOW_SIZE = 7  # residues averaged by the "window" method
_WINDOW_THRESHOLD = 1.0  # a window mean must exceed 1, the neutral value of the parameters, to call a helix or sheet
_MIN_ELEMENT_LENGTH = 3  # shorter helices and sheets of the "window" method become coil


def predict_secondary_structure(sequence: str, method: str = "window") -> List[str]:
    """Predict secondary structure for a protein sequence from residue propensities.

    Two methods run locally and use nothing but the sequence and the Chou-Fasman parameters
    (:data:`CHOU_FASMAN_PARAMETERS`):

    - ``"window"`` (default): the mean helix and sheet parameters over a window of 7 residues centred on each
      position (shorter at the ends, counting only the standard amino acids in it). The position is helix (``"H"``)
      when the helix mean is the larger and exceeds 1, sheet (``"E"``) when the sheet mean is the larger and exceeds
      1, and coil (``"C"``) otherwise; helices and sheets shorter than 3 residues become coil.
    - ``"simple"``: each residue on its own gets the conformation of its largest parameter: helix when its helix
      parameter is the largest, sheet when its sheet parameter is (also on a tie with the turn parameter), coil
      when its turn parameter is.

    Residues that are not one of the twenty standard amino acids are coil. Neither method is a substitute for
    PSIPRED, JPred or DSSP: single-sequence propensity methods have limited accuracy. ``"psipred"`` and ``"jpred"``
    are not available: they need the PSIPRED/JPred service or a local installation of the program, this package
    calls neither, and it does not substitute another method for them.

    Args:
        sequence: Protein sequence (amino acid string)
        method: ``"window"`` or ``"simple"``

    Returns:
        List of secondary structure assignments (H, E, C), one per residue; empty for an empty sequence

    Raises:
        NotImplementedError: If ``method`` is ``"psipred"`` or ``"jpred"``
        ValueError: If ``method`` is not known

    Example:
        >>> sequence = "MKLVLSDEL"
        >>> ss = predict_secondary_structure(sequence, method="simple")
        >>> len(ss) == len(sequence)
        True
        >>> predict_secondary_structure("EEEEEEEEEE")
        ['H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']
    """
    if not sequence:
        return []

    if method in EXTERNAL_METHODS:
        raise NotImplementedError(
            f"method {method!r} is not available: it needs the {method.upper()} service or program, which this package "
            f"does not call (and it will not silently substitute another method). Use one of {SEQUENCE_METHODS}, or "
            "run the program yourself and read its output."
        )
    if method not in SEQUENCE_METHODS:
        raise ValueError(f"Unsupported method: {method}; use one of {SEQUENCE_METHODS}")

    if method == "simple":
        return _simple_secondary_structure_prediction(sequence)
    return _window_secondary_structure_prediction(sequence)


def _simple_secondary_structure_prediction(sequence: str) -> List[str]:
    """Residue-by-residue prediction from the largest Chou-Fasman parameter of each residue."""
    assignments = []
    for residue in sequence.upper():
        parameters = CHOU_FASMAN_PARAMETERS.get(residue)
        if parameters is None:
            assignments.append("C")
        elif parameters[0] > parameters[1] and parameters[0] > parameters[2]:
            assignments.append("H")
        elif parameters[1] >= parameters[2]:
            assignments.append("E")
        else:
            assignments.append("C")
    return assignments


def _window_secondary_structure_prediction(sequence: str) -> List[str]:
    """Window-averaged Chou-Fasman parameters, with short helices and sheets turned into coil."""
    sequence = sequence.upper()
    n = len(sequence)
    predictions = ["C"] * n

    for i in range(n):
        start = max(0, i - _WINDOW_SIZE // 2)
        end = min(n, i + _WINDOW_SIZE // 2 + 1)
        known = [CHOU_FASMAN_PARAMETERS[aa] for aa in sequence[start:end] if aa in CHOU_FASMAN_PARAMETERS]
        if not known:
            continue
        helix_score = sum(parameters[0] for parameters in known) / len(known)
        sheet_score = sum(parameters[1] for parameters in known) / len(known)

        if helix_score > sheet_score and helix_score > _WINDOW_THRESHOLD:
            predictions[i] = "H"
        elif sheet_score > helix_score and sheet_score > _WINDOW_THRESHOLD:
            predictions[i] = "E"

    return _apply_minimum_length_constraints(predictions)


def _apply_minimum_length_constraints(predictions: List[str]) -> List[str]:
    """Turn helices and sheets shorter than the minimum element length into coil."""
    result = predictions.copy()
    n = len(predictions)

    i = 0
    while i < n:
        current_ss = result[i]

        # Find end of current secondary structure element
        j = i + 1
        while j < n and result[j] == current_ss:
            j += 1

        length = j - i

        # If secondary structure element is too short, convert to coil
        if current_ss in ["H", "E"] and length < _MIN_ELEMENT_LENGTH:
            for k in range(i, j):
                result[k] = "C"

        i = j

    return result


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
        return {"helix": 0.0, "sheet": 0.0, "coil": 0.0}

    total = len(ss_assignments)
    helix_count = sum(1 for ss in ss_assignments if ss == "H")
    sheet_count = sum(1 for ss in ss_assignments if ss == "E")
    coil_count = sum(1 for ss in ss_assignments if ss == "C")

    return {
        "helix": helix_count / total,
        "sheet": sheet_count / total,
        "coil": coil_count / total,
        "helix_count": helix_count,
        "sheet_count": sheet_count,
        "coil_count": coil_count,
        "total_residues": total,
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
    elements: List[Dict[str, Any]] = []

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
            if current_length >= min_length and current_type in ["H", "E"]:
                elements.append(
                    {
                        "type": current_type,
                        "start": current_start,
                        "end": current_start + current_length - 1,
                        "length": current_length,
                    }
                )

            # Start new element
            current_type = ss_assignments[i]
            current_start = i
            current_length = 1

    # Handle last element
    if current_length >= min_length and current_type in ["H", "E"]:
        elements.append(
            {
                "type": current_type,
                "start": current_start,
                "end": current_start + current_length - 1,
                "length": current_length,
            }
        )

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
        return {"accuracy": 0.0, "matches": 0, "total": 0}

    matches = sum(1 for a, b in zip(prediction1, prediction2) if a == b)
    total = len(prediction1)
    accuracy = matches / total if total > 0 else 0.0

    return {
        "accuracy": accuracy,
        "matches": matches,
        "total": total,
        "q3_score": accuracy,  # Q3 score for 3-state prediction
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
    dssp_map = {"H": "H", "E": "E", "C": " "}

    return "".join(dssp_map.get(ss, " ") for ss in ss_assignments)


def parse_dssp_file(dssp_content: str) -> Dict[str, Any]:
    """Parse DSSP output file.

    Parses the standard DSSP output format, extracting per-residue secondary
    structure assignments, solvent accessibility and phi/psi angles. An
    accessibility or angle that the line does not hold is NaN (never 0), and an
    angle DSSP leaves undefined keeps its value of 360.0.

    Args:
        dssp_content: DSSP file content as string

    Returns:
        Parsed DSSP data with keys: sequence, secondary_structure,
        accessibility, phi_angles, psi_angles, residues, total_residues

    Example:
        >>> content = "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC\\n    1    1 A M  H  \\n"
        >>> data = parse_dssp_file(content)
        >>> isinstance(data, dict)
        True
    """
    lines = dssp_content.strip().split("\n")

    sequence = []
    secondary_structure = []
    accessibility = []
    phi_angles = []
    psi_angles = []
    residues = []

    # Find the start of residue data (after the header line starting with "  #  RESIDUE")
    data_start = 0
    for i, line in enumerate(lines):
        if line.strip().startswith("#  RESIDUE"):
            data_start = i + 1
            break

    if data_start == 0:
        # Try alternate: look for lines with residue data pattern
        for i, line in enumerate(lines):
            if len(line) >= 17 and line[5:10].strip().isdigit():
                data_start = i
                break

    for line in lines[data_start:]:
        if len(line) < 17:
            continue

        try:
            # DSSP format: columns are fixed-width
            # Col 1-5: residue number (DSSP numbering)
            # Col 6-10: residue number (PDB numbering)
            # Col 11: insertion code
            # Col 12: chain ID
            # Col 13: amino acid (one letter)
            # Col 16: secondary structure code
            # Col 35-38: solvent accessibility

            aa = line[13:14].strip() if len(line) > 13 else ""
            if not aa or aa == "!":  # Chain break
                continue

            ss = line[16:17].strip() if len(line) > 16 else ""
            # Normalize SS code: H=helix, E=sheet, everything else=coil
            if ss in ("H", "G", "I"):  # H=alpha, G=3-10, I=pi helix
                ss_norm = "H"
            elif ss in ("E", "B"):  # E=extended, B=bridge
                ss_norm = "E"
            else:
                ss_norm = "C"

            # Parse accessibility; a field that is missing or not a number is NaN, not zero
            acc = float("nan")
            if len(line) >= 38:
                acc_str = line[35:38].strip()
                if acc_str:
                    try:
                        acc = float(acc_str)
                    except ValueError:
                        acc = float("nan")

            # Parse phi/psi angles (DSSP writes 360.0 for an angle that is undefined); NaN where the line has none
            phi = float("nan")
            psi = float("nan")
            if len(line) >= 109:
                phi_str = line[103:109].strip()
                psi_str = line[109:115].strip() if len(line) >= 115 else ""
                try:
                    phi = float(phi_str) if phi_str else float("nan")
                except ValueError:
                    phi = float("nan")
                try:
                    psi = float(psi_str) if psi_str else float("nan")
                except ValueError:
                    psi = float("nan")

            # Parse residue number and chain
            res_num_str = line[5:10].strip()
            chain_id = line[11:12].strip()
            res_num = int(res_num_str) if res_num_str else 0

            sequence.append(aa)
            secondary_structure.append(ss_norm)
            accessibility.append(acc)
            phi_angles.append(phi)
            psi_angles.append(psi)
            residues.append(
                {
                    "residue_number": res_num,
                    "chain_id": chain_id,
                    "amino_acid": aa,
                    "secondary_structure": ss_norm,
                    "raw_ss": line[16:17].strip(),
                    "accessibility": acc,
                    "phi": phi,
                    "psi": psi,
                }
            )

        except (ValueError, IndexError) as e:
            logger.debug(f"Skipping DSSP line: {e}")
            continue

    return {
        "sequence": "".join(sequence),
        "secondary_structure": secondary_structure,
        "accessibility": accessibility,
        "phi_angles": phi_angles,
        "psi_angles": psi_angles,
        "residues": residues,
        "total_residues": len(sequence),
    }


def predict_transmembrane_regions(sequence: str) -> List[Dict[str, Any]]:
    """Predict transmembrane regions from sequence by Kyte-Doolittle hydropathy.

    The mean Kyte-Doolittle hydropathy is taken over every window of 19 residues (the typical length of a
    transmembrane helix). Windows whose mean exceeds 1.5 are hydrophobic enough to span a membrane; overlapping
    windows that follow each other are merged, so a transmembrane stretch is reported once, from the first residue of
    its first window to the last residue of its last. Letters without a hydropathy value count as 0.

    This is a hydropathy plot, not a topology predictor like TMHMM: it does not tell the orientation or tell a
    signal peptide or another hydrophobic helix from a transmembrane one.

    Args:
        sequence: Protein sequence

    Returns:
        List of regions in sequence order, each with ``start`` and ``end`` (0-based, both included), ``length``,
        ``avg_hydrophobicity`` (mean hydropathy over the region) and ``sequence``

    Example:
        >>> regions = predict_transmembrane_regions("L" * 25)  # seven overlapping windows, one region
        >>> [(region["start"], region["end"], region["length"]) for region in regions]
        [(0, 24, 25)]
    """
    # Kyte-Doolittle hydrophobicity scale (canonical table from sequence.sequences)
    from metainformant.protein.sequence.sequences import KD_HYDROPHOBICITY

    sequence = sequence.upper()
    window_size = 19  # Typical TM helix length
    threshold = 1.5  # Hydrophobicity threshold

    regions: List[Dict[str, Any]] = []
    first = last = None
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i : i + window_size]
        if sum(KD_HYDROPHOBICITY.get(aa, 0) for aa in window) / window_size > threshold:
            if first is None:
                first = i
            last = i + window_size - 1
        elif first is not None and last is not None:
            regions.append(_transmembrane_region(sequence, first, last, KD_HYDROPHOBICITY))
            first = last = None
    if first is not None and last is not None:
        regions.append(_transmembrane_region(sequence, first, last, KD_HYDROPHOBICITY))

    return regions


def _transmembrane_region(sequence: str, first: int, last: int, hydropathy: Dict[str, float]) -> Dict[str, Any]:
    """The transmembrane region ``sequence[first:last + 1]`` with its mean hydropathy."""
    segment = sequence[first : last + 1]
    return {
        "start": first,
        "end": last,
        "length": last - first + 1,
        "avg_hydrophobicity": sum(hydropathy.get(aa, 0) for aa in segment) / len(segment),
        "sequence": segment,
    }


def calculate_ss_propensities(sequence: str) -> Dict[str, Dict[str, Union[float, str]]]:
    """Look up the Chou-Fasman conformational parameters of the amino acids in a sequence.

    Args:
        sequence: Protein sequence

    Returns:
        For every distinct standard amino acid of the sequence, a dictionary with its helix (``"helix"``), sheet
        (``"sheet"``) and turn (``"turn"``) parameters from :data:`CHOU_FASMAN_PARAMETERS`, and ``"preferred"``,
        ``"helix"`` or ``"sheet"`` according to the larger of the first two. Letters that are not one of the twenty
        standard amino acids are left out.

    Example:
        >>> sequence = "ACDEFGHIKLMNPQRSTVWY"
        >>> propensities = calculate_ss_propensities(sequence)
        >>> propensities["E"]["helix"], propensities["V"]["preferred"]
        (1.51, 'sheet')
    """
    propensities: Dict[str, Dict[str, Union[float, str]]] = {}

    for aa in sorted(set(sequence.upper())):
        if aa in CHOU_FASMAN_PARAMETERS:
            helix, sheet, turn = CHOU_FASMAN_PARAMETERS[aa]
            propensities[aa] = {
                "helix": helix,
                "sheet": sheet,
                "turn": turn,
                "preferred": "helix" if helix > sheet else "sheet",
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
        >>> seq = "ALA"
        >>> validation = validate_ss_prediction(ss, seq)
        >>> validation['valid_length']
        True
        >>> validate_ss_prediction(ss, "ALAALAALA")['valid_length']
        False
    """
    validation: Dict[str, Any] = {
        "valid_length": len(ss_assignments) == len(sequence),
        "valid_states": all(ss in ["H", "E", "C"] for ss in ss_assignments),
        "has_structure": any(ss in ["H", "E"] for ss in ss_assignments),
        "issues": [],
    }

    if not validation["valid_length"]:
        validation["issues"].append(f"Length mismatch: {len(ss_assignments)} vs {len(sequence)}")

    if not validation["valid_states"]:
        invalid = [ss for ss in ss_assignments if ss not in ["H", "E", "C"]]
        validation["issues"].append(f"Invalid states: {invalid}")

    if not validation["has_structure"]:
        validation["issues"].append("No secondary structure elements predicted")

    return validation
