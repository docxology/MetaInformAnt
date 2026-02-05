"""Protein analysis orchestration module.

Provides high-level workflows that compose multiple protein analysis steps
into complete pipelines for sequence analysis, structure analysis, and
comparative studies.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from metainformant.core import io, logging

logger = logging.get_logger(__name__)


def analyze_protein_sequence(
    sequence: str,
    name: str = "query",
    predict_ss: bool = True,
    find_motifs: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Run a complete protein sequence analysis pipeline.

    Computes physicochemical properties, secondary structure prediction,
    motif searching, and transmembrane region prediction.

    Args:
        sequence: Protein amino acid sequence
        name: Identifier for the protein
        predict_ss: Whether to predict secondary structure
        find_motifs: Optional list of motif patterns to search

    Returns:
        Dictionary with all analysis results
    """
    from .sequence.sequences import (
        amino_acid_composition,
        aromaticity,
        charge_at_ph,
        extinction_coefficient,
        gravy,
        hydropathy_score,
        instability_index,
        isoelectric_point,
        molecular_weight,
        sequence_length,
        transmembrane_regions,
        validate_protein_sequence,
    )

    logger.info(f"Running sequence analysis for '{name}' ({len(sequence)} residues)")

    results: Dict[str, Any] = {
        "name": name,
        "sequence": sequence,
        "length": sequence_length(sequence),
        "valid": validate_protein_sequence(sequence),
    }

    # Physicochemical properties
    results["physicochemical"] = {
        "molecular_weight": molecular_weight(sequence),
        "isoelectric_point": isoelectric_point(sequence),
        "gravy": gravy(sequence),
        "aromaticity": aromaticity(sequence),
        "instability_index": instability_index(sequence),
        "extinction_coefficient": extinction_coefficient(sequence),
        "charge_at_ph7": charge_at_ph(sequence, 7.0),
    }

    # Composition
    results["composition"] = amino_acid_composition(sequence)

    # Hydropathy
    if len(sequence) >= 19:
        results["hydropathy"] = hydropathy_score(sequence)
        results["transmembrane_regions"] = transmembrane_regions(sequence)
    else:
        results["hydropathy"] = []
        results["transmembrane_regions"] = []

    # Secondary structure prediction
    if predict_ss:
        from .structure.secondary import (
            calculate_ss_composition,
            identify_ss_elements,
            predict_secondary_structure,
        )

        ss = predict_secondary_structure(sequence, method="simple")
        results["secondary_structure"] = {
            "prediction": "".join(ss),
            "composition": calculate_ss_composition(ss),
            "elements": identify_ss_elements(ss),
        }

    # Motif search
    if find_motifs:
        from .sequence.sequences import find_motifs as search_motifs

        results["motifs"] = search_motifs(sequence, find_motifs)

    logger.info(f"Sequence analysis complete for '{name}'")
    return results


def analyze_protein_structure(
    pdb_path: Union[str, Path],
    compute_contacts: bool = True,
) -> Dict[str, Any]:
    """Run a complete protein structure analysis pipeline.

    Args:
        pdb_path: Path to PDB file
        compute_contacts: Whether to compute contact maps

    Returns:
        Dictionary with structural analysis results
    """
    import numpy as np

    from .structure.io import parse_pdb_file

    pdb_path = Path(pdb_path)
    logger.info(f"Running structure analysis on {pdb_path.name}")

    structure = parse_pdb_file(pdb_path)
    atoms = structure.get("atoms", [])

    if not atoms:
        return {"error": "No atoms found in structure", "path": str(pdb_path)}

    # Extract coordinates
    coords = np.array([[a["x"], a["y"], a["z"]] for a in atoms])

    results: Dict[str, Any] = {
        "path": str(pdb_path),
        "n_atoms": len(atoms),
    }

    # Structural statistics
    from .structure.general import (
        calculate_center_of_mass,
        calculate_radius_of_gyration,
        calculate_structural_statistics,
    )

    results["statistics"] = calculate_structural_statistics(coords)
    results["radius_of_gyration"] = calculate_radius_of_gyration(coords)
    results["center_of_mass"] = calculate_center_of_mass(coords).tolist()

    # Contact analysis
    if compute_contacts:
        from .structure.analysis import calculate_contact_map

        contact_map = calculate_contact_map(coords, threshold=8.0)
        results["contact_map_shape"] = contact_map.shape
        results["n_contacts"] = int(np.sum(contact_map) / 2)

    # Sequence extraction
    from .structure.pdb import get_pdb_sequence

    chains = set(a.get("chain_id", "A") for a in atoms)
    results["chains"] = {}
    for chain in sorted(chains):
        seq = get_pdb_sequence(structure, chain)
        if seq:
            results["chains"][chain] = {
                "sequence": seq,
                "length": len(seq),
            }

    logger.info(f"Structure analysis complete: {len(atoms)} atoms, {len(chains)} chains")
    return results


def batch_analyze_sequences(
    sequences: Dict[str, str],
    predict_ss: bool = False,
) -> Dict[str, Dict[str, Any]]:
    """Run sequence analysis on multiple proteins.

    Args:
        sequences: Dictionary mapping names to sequences
        predict_ss: Whether to predict secondary structure

    Returns:
        Dictionary mapping names to analysis results
    """
    logger.info(f"Batch analyzing {len(sequences)} sequences")
    results = {}
    for name, seq in sequences.items():
        results[name] = analyze_protein_sequence(seq, name=name, predict_ss=predict_ss)
    return results


def comparative_analysis(
    sequences: Dict[str, str],
    use_blosum: bool = True,
) -> Dict[str, Any]:
    """Compare multiple protein sequences with pairwise alignment and MSA.

    Args:
        sequences: Dictionary mapping names to sequences
        use_blosum: Whether to use BLOSUM62 matrix

    Returns:
        Dictionary with comparative analysis results
    """
    names = list(sequences.keys())
    seqs = list(sequences.values())
    n = len(seqs)

    logger.info(f"Comparative analysis of {n} sequences")

    if n < 2:
        return {"error": "Need at least 2 sequences for comparison"}

    # Pairwise alignments
    pairwise = {}
    if use_blosum:
        from .sequence.alignment import matrix_align
        align_fn = matrix_align
    else:
        from .sequence.alignment import global_align
        align_fn = global_align

    for i in range(n):
        for j in range(i + 1, n):
            key = f"{names[i]}_vs_{names[j]}"
            result = align_fn(seqs[i], seqs[j])
            pairwise[key] = {
                "score": result["score"],
                "identity": result["identity"],
            }

    # Distance matrix
    identity_matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        identity_matrix[i][i] = 1.0
        for j in range(i + 1, n):
            key = f"{names[i]}_vs_{names[j]}"
            ident = pairwise[key]["identity"]
            identity_matrix[i][j] = ident
            identity_matrix[j][i] = ident

    # MSA
    from .sequence.alignment import multi_sequence_alignment

    msa = multi_sequence_alignment(seqs)

    return {
        "names": names,
        "pairwise_alignments": pairwise,
        "identity_matrix": identity_matrix,
        "msa": msa,
    }


def analyze_from_fasta(
    fasta_path: Union[str, Path],
    predict_ss: bool = False,
    compare: bool = False,
) -> Dict[str, Any]:
    """Complete analysis pipeline starting from a FASTA file.

    Args:
        fasta_path: Path to input FASTA file
        predict_ss: Whether to predict secondary structure
        compare: Whether to run comparative analysis

    Returns:
        Dictionary with all analysis results
    """
    from .sequence.sequences import read_fasta

    fasta_path = Path(fasta_path)
    logger.info(f"Loading sequences from {fasta_path}")
    sequences = read_fasta(fasta_path)

    results: Dict[str, Any] = {
        "input": str(fasta_path),
        "n_sequences": len(sequences),
    }

    # Individual analyses
    results["sequences"] = batch_analyze_sequences(sequences, predict_ss=predict_ss)

    # Comparative analysis
    if compare and len(sequences) >= 2:
        results["comparative"] = comparative_analysis(sequences)

    return results


def compare_structures(
    pdb_path_a: Union[str, Path],
    pdb_path_b: Union[str, Path],
    align: bool = True,
) -> Dict[str, Any]:
    """Compare two protein structures with RMSD, alignment, and contact analysis.

    Args:
        pdb_path_a: Path to first PDB file
        pdb_path_b: Path to second PDB file
        align: Whether to perform Kabsch alignment before computing RMSD

    Returns:
        Dictionary with comparison results including RMSD, contact differences,
        and per-residue distance deviations
    """
    import numpy as np

    from .structure.io import parse_pdb_file

    pdb_path_a = Path(pdb_path_a)
    pdb_path_b = Path(pdb_path_b)

    logger.info(f"Comparing structures: {pdb_path_a.name} vs {pdb_path_b.name}")

    struct_a = parse_pdb_file(pdb_path_a)
    struct_b = parse_pdb_file(pdb_path_b)

    atoms_a = struct_a.get("atoms", [])
    atoms_b = struct_b.get("atoms", [])

    if not atoms_a or not atoms_b:
        return {"error": "One or both structures have no atoms"}

    # Extract CA coordinates
    ca_a = np.array([[a["x"], a["y"], a["z"]] for a in atoms_a if a["name"] == "CA"])
    ca_b = np.array([[a["x"], a["y"], a["z"]] for a in atoms_b if a["name"] == "CA"])

    results: Dict[str, Any] = {
        "structure_a": str(pdb_path_a),
        "structure_b": str(pdb_path_b),
        "n_ca_atoms_a": len(ca_a),
        "n_ca_atoms_b": len(ca_b),
    }

    # Truncate to common length for comparison
    n_common = min(len(ca_a), len(ca_b))
    if n_common == 0:
        results["error"] = "No CA atoms found in one or both structures"
        return results

    ca_a_common = ca_a[:n_common]
    ca_b_common = ca_b[:n_common]

    if align:
        from .structure.general import align_structures_kabsch

        aligned_b, rotation, rmsd = align_structures_kabsch(ca_a_common, ca_b_common)
        results["rmsd_aligned"] = float(rmsd)
        results["rotation_matrix"] = rotation.tolist()

        # Per-residue distance after alignment
        per_res_dist = np.linalg.norm(ca_a_common - aligned_b, axis=1)
    else:
        from .structure.general import compute_rmsd_simple

        rmsd = compute_rmsd_simple(ca_a_common, ca_b_common)
        results["rmsd_simple"] = float(rmsd)
        per_res_dist = np.linalg.norm(ca_a_common - ca_b_common, axis=1)

    results["n_compared_residues"] = n_common
    results["per_residue_distance"] = {
        "mean": float(np.mean(per_res_dist)),
        "max": float(np.max(per_res_dist)),
        "min": float(np.min(per_res_dist)),
        "std": float(np.std(per_res_dist)),
    }

    # Contact map comparison
    from .structure.analysis import calculate_contact_map

    contact_a = calculate_contact_map(ca_a_common, threshold=8.0)
    contact_b = calculate_contact_map(ca_b_common, threshold=8.0)

    contacts_only_a = int(np.sum((contact_a == 1) & (contact_b == 0)) / 2)
    contacts_only_b = int(np.sum((contact_a == 0) & (contact_b == 1)) / 2)
    shared_contacts = int(np.sum((contact_a == 1) & (contact_b == 1)) / 2)

    results["contact_comparison"] = {
        "shared_contacts": shared_contacts,
        "contacts_only_a": contacts_only_a,
        "contacts_only_b": contacts_only_b,
        "jaccard_similarity": (
            shared_contacts / max(1, shared_contacts + contacts_only_a + contacts_only_b)
        ),
    }

    logger.info(f"Structure comparison complete: RMSD={per_res_dist.mean():.2f} Å")
    return results


def assess_alphafold_quality(
    pdb_path: Union[str, Path],
) -> Dict[str, Any]:
    """Assess quality of an AlphaFold predicted structure.

    Combines confidence score analysis, structural statistics,
    and secondary structure composition into a quality report.

    Args:
        pdb_path: Path to AlphaFold PDB file

    Returns:
        Quality assessment report
    """
    import numpy as np

    pdb_path = Path(pdb_path)
    logger.info(f"Assessing AlphaFold quality for {pdb_path.name}")

    results: Dict[str, Any] = {"path": str(pdb_path)}

    # Parse confidence scores (pLDDT from B-factor column)
    from .structure.alphafold import parse_alphafold_confidence, get_alphafold_structure_quality

    quality = get_alphafold_structure_quality(pdb_path)
    results["confidence"] = quality

    # Structure statistics
    from .structure.io import parse_pdb_file
    from .structure.general import (
        calculate_radius_of_gyration,
        calculate_structural_statistics,
    )

    structure = parse_pdb_file(pdb_path)
    atoms = structure.get("atoms", [])

    if atoms:
        coords = np.array([[a["x"], a["y"], a["z"]] for a in atoms])
        results["structural_stats"] = calculate_structural_statistics(coords)

        # Extract sequence and predict secondary structure
        from .structure.pdb import get_pdb_sequence
        from .structure.secondary import predict_secondary_structure, calculate_ss_composition

        sequence = get_pdb_sequence(structure)
        if sequence:
            ss = predict_secondary_structure(sequence, method="simple")
            results["secondary_structure"] = calculate_ss_composition(ss)
            results["sequence_length"] = len(sequence)

    # Overall quality assessment
    plddt = quality.get("plddt_score", 0.0)
    if plddt >= 90:
        results["quality_rating"] = "very_high"
    elif plddt >= 70:
        results["quality_rating"] = "confident"
    elif plddt >= 50:
        results["quality_rating"] = "low"
    else:
        results["quality_rating"] = "very_low"

    logger.info(f"AlphaFold quality assessment complete: pLDDT={plddt:.1f}")
    return results


def full_protein_analysis(
    sequence: str,
    pdb_path: Optional[Union[str, Path]] = None,
    name: str = "query",
    predict_ss: bool = True,
) -> Dict[str, Any]:
    """Combined sequence and structure analysis pipeline.

    Runs sequence analysis and optionally structure analysis, merging
    results into a unified report.

    Args:
        sequence: Protein amino acid sequence
        pdb_path: Optional path to PDB structure file
        name: Protein identifier
        predict_ss: Whether to predict secondary structure

    Returns:
        Combined analysis results
    """
    logger.info(f"Full protein analysis for '{name}'")

    results: Dict[str, Any] = {"name": name}

    # Sequence analysis
    results["sequence_analysis"] = analyze_protein_sequence(
        sequence, name=name, predict_ss=predict_ss,
    )

    # Structure analysis (if PDB provided)
    if pdb_path is not None:
        results["structure_analysis"] = analyze_protein_structure(pdb_path)

    # Summary statistics
    results["summary"] = {
        "sequence_length": len(sequence),
        "has_structure": pdb_path is not None,
        "molecular_weight": results["sequence_analysis"]["physicochemical"]["molecular_weight"],
        "isoelectric_point": results["sequence_analysis"]["physicochemical"]["isoelectric_point"],
    }

    logger.info(f"Full analysis complete for '{name}'")
    return results


def batch_compare_structures(
    pdb_paths: List[Union[str, Path]],
    align: bool = True,
) -> Dict[str, Any]:
    """Pairwise structural comparison of multiple PDB files.

    Args:
        pdb_paths: List of PDB file paths
        align: Whether to use Kabsch alignment

    Returns:
        Dictionary with all pairwise comparisons and summary matrix
    """
    import numpy as np

    paths = [Path(p) for p in pdb_paths]
    n = len(paths)

    logger.info(f"Batch structure comparison of {n} structures")

    if n < 2:
        return {"error": "Need at least 2 structures for comparison"}

    pairwise = {}
    rmsd_matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            key = f"{paths[i].stem}_vs_{paths[j].stem}"
            result = compare_structures(paths[i], paths[j], align=align)
            pairwise[key] = result

            rmsd_val = result.get("rmsd_aligned", result.get("rmsd_simple", float("inf")))
            rmsd_matrix[i][j] = rmsd_val
            rmsd_matrix[j][i] = rmsd_val

    return {
        "names": [p.stem for p in paths],
        "pairwise_comparisons": pairwise,
        "rmsd_matrix": rmsd_matrix,
        "n_structures": n,
    }
