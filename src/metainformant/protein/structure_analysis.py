"""Advanced protein structure analysis utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from metainformant.core.logging import get_logger

logger = get_logger(__name__)


def predict_secondary_structure(
    sequence: str,
    method: str = "simple",
) -> Dict[str, Any]:
    """Predict secondary structure from protein sequence.
    
    Args:
        sequence: Protein sequence
        method: Prediction method ('simple', 'chou_fasman')
        
    Returns:
        Dictionary with:
        - 'sequence': Input sequence
        - 'helix': List of helix probabilities per residue
        - 'sheet': List of sheet probabilities per residue
        - 'coil': List of coil probabilities per residue
        - 'predicted_structure': String of predicted structure (H/E/C)
    """
    if method == "simple":
        from metainformant.protein.secondary import simple_helix_coil_propensity
        
        helix_prop = simple_helix_coil_propensity(sequence)
        # Simple: assume coil if not helix
        coil_prop = [1.0 - h for h in helix_prop]
        sheet_prop = [0.0] * len(sequence)  # Placeholder
        
        # Predict structure
        predicted = []
        for h, c in zip(helix_prop, coil_prop):
            if h > 0.5:
                predicted.append("H")
            elif c > 0.7:
                predicted.append("C")
            else:
                predicted.append("E")  # Sheet
        
        return {
            "sequence": sequence,
            "helix": helix_prop,
            "sheet": sheet_prop,
            "coil": coil_prop,
            "predicted_structure": "".join(predicted),
        }
    else:
        raise ValueError(f"Unknown method: {method}")


def identify_domains(
    sequence: str,
    domain_database: Dict[str, str] | None = None,
) -> List[Dict[str, Any]]:
    """Identify protein domains from sequence.
    
    Uses pattern matching or database lookup to identify known domains.
    
    Args:
        sequence: Protein sequence
        domain_database: Optional dictionary mapping domain names to patterns
        
    Returns:
        List of domain annotations with 'name', 'start', 'end', 'confidence'
    """
    domains = []
    
    if domain_database:
        for domain_name, pattern in domain_database.items():
            # Simple pattern matching
            if pattern in sequence:
                start = sequence.find(pattern)
                end = start + len(pattern)
                domains.append({
                    "name": domain_name,
                    "start": start,
                    "end": end,
                    "confidence": 0.8,  # Placeholder
                })
    
    return domains


def analyze_protein_stability(
    sequence: str,
    structure_coords: np.ndarray | None = None,
) -> Dict[str, float]:
    """Analyze protein stability from sequence or structure.
    
    Calculates stability-related metrics such as:
    - Hydrophobicity
    - Charge distribution
    - Flexibility
    
    Args:
        sequence: Protein sequence
        structure_coords: Optional 3D coordinates (N x 3 array)
        
    Returns:
        Dictionary with stability metrics
    """
    # Amino acid properties
    hydrophobicity = {
        "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
        "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
        "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
        "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
    }
    
    charge = {
        "D": -1, "E": -1, "K": 1, "R": 1, "H": 0.5,
    }
    
    # Calculate metrics
    seq_upper = sequence.upper()
    hydro_scores = [hydrophobicity.get(aa, 0.0) for aa in seq_upper]
    mean_hydrophobicity = np.mean(hydro_scores) if hydro_scores else 0.0
    
    charge_scores = [charge.get(aa, 0.0) for aa in seq_upper]
    net_charge = sum(charge_scores)
    
    # Flexibility (based on amino acid types)
    flexible_aas = {"G", "P", "S", "N", "D", "T"}
    flexibility = sum(1 for aa in seq_upper if aa in flexible_aas) / len(seq_upper) if seq_upper else 0.0
    
    result = {
        "mean_hydrophobicity": float(mean_hydrophobicity),
        "net_charge": float(net_charge),
        "flexibility": float(flexibility),
    }
    
    # Structure-based metrics if coordinates provided
    if structure_coords is not None:
        # Calculate radius of gyration
        center = np.mean(structure_coords, axis=0)
        distances = np.linalg.norm(structure_coords - center, axis=1)
        radius_gyration = np.sqrt(np.mean(distances ** 2))
        result["radius_of_gyration"] = float(radius_gyration)
    
    return result


def predict_protein_family(
    sequence: str,
    family_database: Dict[str, List[str]] | None = None,
) -> List[Tuple[str, float]]:
    """Predict protein family from sequence.
    
    Uses sequence similarity or pattern matching to identify protein families.
    
    Args:
        sequence: Protein sequence
        family_database: Optional dictionary mapping family names to representative sequences
        
    Returns:
        List of (family_name, similarity_score) tuples, sorted by score
    """
    if not family_database:
        return []
    
    from metainformant.protein.alignment import pairwise_identity
    
    family_scores = []
    
    for family_name, representative_seqs in family_database.items():
        max_similarity = 0.0
        for rep_seq in representative_seqs:
            similarity = pairwise_identity(sequence, rep_seq)
            max_similarity = max(max_similarity, similarity)
        
        if max_similarity > 0.3:  # Threshold
            family_scores.append((family_name, max_similarity))
    
    # Sort by similarity (descending)
    family_scores.sort(key=lambda x: x[1], reverse=True)
    
    return family_scores


def analyze_post_translational_modifications(
    sequence: str,
    modification_sites: Dict[str, List[int]] | None = None,
) -> Dict[str, Any]:
    """Analyze post-translational modification sites.
    
    Identifies potential PTM sites (phosphorylation, glycosylation, etc.)
    from sequence patterns.
    
    Args:
        sequence: Protein sequence
        modification_sites: Optional dictionary mapping modification types to site positions
        
    Returns:
        Dictionary with PTM analysis results
    """
    # Simple pattern-based prediction
    ptm_sites = {
        "phosphorylation": [],
        "glycosylation": [],
        "acetylation": [],
    }
    
    # Phosphorylation sites (S, T, Y)
    for i, aa in enumerate(sequence):
        if aa in "STY":
            ptm_sites["phosphorylation"].append(i)
    
    # N-linked glycosylation (N-X-S/T, X != P)
    for i in range(len(sequence) - 2):
        if sequence[i] == "N" and sequence[i + 1] != "P" and sequence[i + 2] in "ST":
            ptm_sites["glycosylation"].append(i)
    
    # Acetylation (K at N-terminus or specific contexts)
    if sequence and sequence[0] == "K":
        ptm_sites["acetylation"].append(0)
    
    # Add user-provided sites
    if modification_sites:
        for mod_type, sites in modification_sites.items():
            if mod_type in ptm_sites:
                ptm_sites[mod_type].extend(sites)
            else:
                ptm_sites[mod_type] = sites
    
    return {
        "sequence": sequence,
        "ptm_sites": ptm_sites,
        "total_ptm_sites": sum(len(sites) for sites in ptm_sites.values()),
    }

