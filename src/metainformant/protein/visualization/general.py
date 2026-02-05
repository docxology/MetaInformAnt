"""Protein structure and sequence visualization functions.

This module provides comprehensive visualization capabilities for protein data,
including structure visualization, sequence analysis plots, domain architecture,
and structural analysis results.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle, Circle, FancyBboxPatch

from metainformant.core import logging, paths, validation

logger = logging.get_logger(__name__)

try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    sns = None

try:
    import plotly.graph_objects as go
    import plotly.express as px

    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False
    go = None
    px = None


def plot_sequence_logo(
    sequences: Sequence[str],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 4),
    **kwargs,
) -> Axes:
    """Create a sequence logo from aligned protein sequences.

    Args:
        sequences: List of aligned protein sequences
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(sequences, (list, tuple), "sequences")
    validation.validate_not_empty(sequences, "sequences")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Calculate position-specific amino acid frequencies
    seq_length = len(sequences[0])
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    frequencies = np.zeros((seq_length, len(amino_acids)))

    for seq in sequences:
        if len(seq) != seq_length:
            raise ValueError(f"All sequences must have same length: {seq_length}")
        for i, aa in enumerate(seq.upper()):
            if aa in amino_acids:
                idx = amino_acids.index(aa)
                frequencies[i, idx] += 1

    frequencies = frequencies / len(sequences)

    # Calculate information content
    info_content = np.zeros(seq_length)
    for i in range(seq_length):
        for j in range(len(amino_acids)):
            if frequencies[i, j] > 0:
                info_content[i] += frequencies[i, j] * np.log2(frequencies[i, j])

    info_content = -info_content + np.log2(len(amino_acids))

    # Plot sequence logo
    colors = plt.cm.tab20(np.linspace(0, 1, len(amino_acids)))
    bottom = np.zeros(seq_length)

    for i, aa in enumerate(amino_acids):
        heights = frequencies[:, i] * info_content
        ax.bar(range(seq_length), heights, bottom=bottom, color=colors[i], label=aa, alpha=0.8)
        bottom += heights

    ax.set_xlim(-0.5, seq_length - 0.5)
    ax.set_xlabel("Position")
    ax.set_ylabel("Bits")
    ax.set_title("Sequence Logo")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Sequence logo saved to {output_path}")

    return ax


def plot_domain_architecture(
    domains: List[Dict[str, Any]],
    protein_length: int,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 3),
    **kwargs,
) -> Axes:
    """Plot protein domain architecture.

    Args:
        domains: List of domain dictionaries with keys: 'start', 'end', 'name', 'type'
        protein_length: Total length of the protein
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(domains, list, "domains")
    validation.validate_type(protein_length, int, "protein_length")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot protein backbone
    ax.plot([0, protein_length], [0, 0], "k-", linewidth=2, alpha=0.7)

    # Color palette for different domain types
    colors = plt.cm.Set3(np.linspace(0, 1, 12))

    # Plot domains
    y_offset = 0.1
    domain_types = {}
    type_counter = 0

    for domain in domains:
        start = domain["start"]
        end = domain["end"]
        name = domain.get("name", "Unknown")
        domain_type = domain.get("type", "domain")

        if domain_type not in domain_types:
            domain_types[domain_type] = colors[type_counter % len(colors)]
            type_counter += 1

        color = domain_types[domain_type]

        # Draw domain rectangle
        width = end - start
        rect = Rectangle((start, -0.3), width, 0.6, facecolor=color, edgecolor="black", alpha=0.8)
        ax.add_patch(rect)

        # Add domain label
        if width > 50:  # Only label if domain is wide enough
            ax.text(start + width / 2, y_offset, name, ha="center", va="bottom", fontsize=8, rotation=45)

    ax.set_xlim(0, protein_length)
    ax.set_ylim(-0.5, 0.5)
    ax.set_xlabel("Amino Acid Position")
    ax.set_title("Protein Domain Architecture")
    ax.set_yticks([])

    # Add legend
    legend_elements = [
        plt.Rectangle((0, 0), 1, 1, facecolor=color, label=domain_type) for domain_type, color in domain_types.items()
    ]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc="upper left")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Domain architecture plot saved to {output_path}")

    return ax


def plot_secondary_structure(
    sequence: str,
    secondary_structure: str,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 4),
    **kwargs,
) -> Axes:
    """Plot protein secondary structure along sequence.

    Args:
        sequence: Protein sequence
        secondary_structure: Secondary structure string (H=helix, E=sheet, C=coil)
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(sequence, str, "sequence")
    validation.validate_type(secondary_structure, str, "secondary_structure")

    if len(sequence) != len(secondary_structure):
        raise ValueError("Sequence and secondary structure must have same length")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    positions = np.arange(len(sequence))

    # Map secondary structure to colors and labels
    color_map = {"H": "red", "E": "blue", "C": "green", "-": "gray"}
    label_map = {"H": "Helix", "E": "Sheet", "C": "Coil", "-": "Unknown"}

    # Plot each secondary structure type
    for ss_type, color in color_map.items():
        mask = np.array([c == ss_type for c in secondary_structure])
        if np.any(mask):
            ax.scatter(positions[mask], np.ones(np.sum(mask)), c=color, label=label_map[ss_type], alpha=0.7, s=20)

    ax.set_xlim(-1, len(sequence))
    ax.set_ylim(0.5, 1.5)
    ax.set_xlabel("Position")
    ax.set_title("Secondary Structure")
    ax.set_yticks([])
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Secondary structure plot saved to {output_path}")

    return ax


def plot_contact_map(
    contacts: np.ndarray,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 8),
    **kwargs,
) -> Axes:
    """Plot protein contact map.

    Args:
        contacts: Contact matrix (n_residues x n_residues)
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(contacts, np.ndarray, "contacts")

    if contacts.ndim != 2 or contacts.shape[0] != contacts.shape[1]:
        raise ValueError("Contacts must be a square matrix")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot contact map
    im = ax.imshow(contacts, cmap="Blues", aspect="equal", origin="lower", interpolation="nearest")

    ax.set_xlabel("Residue i")
    ax.set_ylabel("Residue j")
    ax.set_title("Protein Contact Map")

    # Add colorbar
    plt.colorbar(im, ax=ax, label="Contact Probability")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Contact map saved to {output_path}")

    return ax


def plot_ramachandran_plot(
    phi_angles: np.ndarray,
    psi_angles: np.ndarray,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 8),
    **kwargs,
) -> Axes:
    """Create a Ramachandran plot.

    Args:
        phi_angles: Phi dihedral angles in degrees
        psi_angles: Psi dihedral angles in degrees
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(phi_angles, np.ndarray, "phi_angles")
    validation.validate_type(psi_angles, np.ndarray, "psi_angles")

    if len(phi_angles) != len(psi_angles):
        raise ValueError("Phi and psi angles must have same length")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Create density plot
    if HAS_SEABORN:
        sns.kdeplot(x=phi_angles, y=psi_angles, ax=ax, fill=True, cmap="Blues", alpha=0.7)
    else:
        ax.scatter(phi_angles, psi_angles, alpha=0.6, s=1, c="blue")

    # Add allowed regions (simplified)
    # Alpha helix region
    alpha_patch = Circle((-60, -45), 20, facecolor="red", alpha=0.2)
    ax.add_patch(alpha_patch)
    ax.text(-60, -45, "α-helix", ha="center", va="center", fontsize=10)

    # Beta sheet region
    beta_patch = Rectangle((-120, 120), 60, 60, facecolor="green", alpha=0.2)
    ax.add_patch(beta_patch)
    ax.text(-90, 150, "β-sheet", ha="center", va="center", fontsize=10)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xlabel("φ (phi)")
    ax.set_ylabel("ψ (psi)")
    ax.set_title("Ramachandran Plot")
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Ramachandran plot saved to {output_path}")

    return ax


def plot_alignment_quality(
    alignment_scores: np.ndarray,
    positions: np.ndarray | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 4),
    **kwargs,
) -> Axes:
    """Plot alignment quality scores along sequence positions.

    Args:
        alignment_scores: Quality scores for each position
        positions: Optional position indices
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(alignment_scores, np.ndarray, "alignment_scores")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    x_data = positions if positions is not None else np.arange(len(alignment_scores))

    ax.plot(x_data, alignment_scores, "b-", linewidth=1, alpha=0.8)
    ax.fill_between(x_data, alignment_scores, alpha=0.3, color="blue")

    ax.set_xlabel("Position")
    ax.set_ylabel("Alignment Quality Score")
    ax.set_title("Sequence Alignment Quality")
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Alignment quality plot saved to {output_path}")

    return ax


def plot_conservation_scores(
    conservation_scores: np.ndarray,
    positions: np.ndarray | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 4),
    **kwargs,
) -> Axes:
    """Plot sequence conservation scores.

    Args:
        conservation_scores: Conservation scores for each position
        positions: Optional position indices
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(conservation_scores, np.ndarray, "conservation_scores")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    x_data = positions if positions is not None else np.arange(len(conservation_scores))

    ax.plot(x_data, conservation_scores, "r-", linewidth=2, alpha=0.8)
    ax.fill_between(x_data, conservation_scores, alpha=0.3, color="red")

    ax.set_xlabel("Position")
    ax.set_ylabel("Conservation Score")
    ax.set_title("Sequence Conservation")
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Conservation plot saved to {output_path}")

    return ax


def plot_structure_superposition(
    structures: List[np.ndarray],
    labels: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 8),
    **kwargs,
) -> Axes:
    """Plot superposition of protein structures.

    Args:
        structures: List of coordinate arrays (n_atoms x 3)
        labels: Optional labels for each structure
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(structures, list, "structures")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, projection="3d")

    colors = plt.cm.tab10(np.linspace(0, 1, len(structures)))

    for i, coords in enumerate(structures):
        if coords.shape[1] != 3:
            raise ValueError("Coordinates must be n_atoms x 3")

        label = labels[i] if labels else f"Structure {i+1}"
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=[colors[i]], label=label, alpha=0.6, s=1)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("Protein Structure Superposition")
    ax.legend()

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Structure superposition plot saved to {output_path}")

    return ax


def plot_protein_properties(
    sequence: str,
    *,
    properties: List[str] | None = None,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 6),
    **kwargs,
) -> Axes:
    """Plot various protein physicochemical properties along sequence.

    Args:
        sequence: Protein sequence
        properties: List of properties to plot ('hydrophobicity', 'charge', 'size')
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(sequence, str, "sequence")

    if properties is None:
        properties = ["hydrophobicity", "charge"]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    positions = np.arange(len(sequence))

    # Amino acid property scales
    hydrophobicity = {
        "A": 1.8,
        "C": 2.5,
        "D": -3.5,
        "E": -3.5,
        "F": 2.8,
        "G": -0.4,
        "H": -3.2,
        "I": 4.5,
        "K": -3.9,
        "L": 3.8,
        "M": 1.9,
        "N": -3.5,
        "P": -1.6,
        "Q": -3.5,
        "R": -4.5,
        "S": -0.8,
        "T": -0.7,
        "V": 4.2,
        "W": -0.9,
        "Y": -1.3,
    }

    charge = {
        "A": 0,
        "C": 0,
        "D": -1,
        "E": -1,
        "F": 0,
        "G": 0,
        "H": 1,
        "I": 0,
        "K": 1,
        "L": 0,
        "M": 0,
        "N": 0,
        "P": 0,
        "Q": 0,
        "R": 1,
        "S": 0,
        "T": 0,
        "V": 0,
        "W": 0,
        "Y": 0,
    }

    # Calculate properties
    if "hydrophobicity" in properties:
        hydro_scores = [hydrophobicity.get(aa.upper(), 0) for aa in sequence]
        ax.plot(positions, hydro_scores, "b-", label="Hydrophobicity", linewidth=2)

    if "charge" in properties:
        charge_scores = [charge.get(aa.upper(), 0) for aa in sequence]
        ax.plot(positions, charge_scores, "r-", label="Charge", linewidth=2)

    ax.set_xlabel("Position")
    ax.set_ylabel("Property Value")
    ax.set_title("Protein Properties")
    ax.legend()
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Protein properties plot saved to {output_path}")

    return ax


def plot_pdb_structure_quality(
    b_factors: np.ndarray,
    positions: np.ndarray | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 4),
    **kwargs,
) -> Axes:
    """Plot PDB structure quality metrics (B-factors).

    Args:
        b_factors: B-factor values for each residue
        positions: Optional residue positions
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(b_factors, np.ndarray, "b_factors")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    x_data = positions if positions is not None else np.arange(len(b_factors))

    ax.plot(x_data, b_factors, "purple", linewidth=1, alpha=0.8)
    ax.fill_between(x_data, b_factors, alpha=0.3, color="purple")

    ax.set_xlabel("Residue Position")
    ax.set_ylabel("B-factor")
    ax.set_title("PDB Structure Quality (B-factors)")
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"PDB quality plot saved to {output_path}")

    return ax


def create_interactive_structure_viewer(
    pdb_data: Dict[str, Any], *, output_path: str | Path | None = None, **kwargs
) -> Any:
    """Create an interactive 3D structure viewer using Plotly.

    Args:
        pdb_data: Dictionary containing PDB structure data
        output_path: Optional path to save the HTML file
        **kwargs: Additional arguments for Plotly customization

    Returns:
        Plotly Figure object
    """
    if not HAS_PLOTLY:
        raise ImportError("Plotly required for interactive structure viewer")

    validation.validate_type(pdb_data, dict, "pdb_data")

    # Extract coordinates
    atoms = pdb_data.get("atoms", [])
    if not atoms:
        raise ValueError("No atom data found in PDB data")

    # Create 3D scatter plot
    fig = go.Figure(
        data=[
            go.Scatter3d(
                x=[atom["x"] for atom in atoms],
                y=[atom["y"] for atom in atoms],
                z=[atom["z"] for atom in atoms],
                mode="markers",
                marker=dict(
                    size=4,
                    color=[atom.get("b_factor", 0) for atom in atoms],
                    colorscale="Viridis",
                    showscale=True,
                    colorbar=dict(title="B-factor"),
                ),
                text=[f"Residue {atom.get('residue_id', '?')}" for atom in atoms],
                hovertemplate="Atom: %{text}<br>X: %{x}<br>Y: %{y}<br>Z: %{z}<extra></extra>",
            )
        ]
    )

    fig.update_layout(
        title="Interactive Protein Structure", scene=dict(xaxis_title="X", yaxis_title="Y", zaxis_title="Z"), **kwargs
    )

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        html_path = Path(output_path).with_suffix(".html")
        fig.write_html(str(html_path))
        logger.info(f"Interactive structure viewer saved to {html_path}")

    return fig


def plot_helical_wheel(
    sequence: str,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 8),
    **kwargs,
) -> Axes:
    """Plot a helical wheel diagram for a protein segment.

    Displays residues projected onto a circle as viewed end-on along
    the helix axis (100 degrees per residue for alpha helix).

    Args:
        sequence: Protein sequence (typically 18-25 residues)
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(sequence, str, "sequence")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={"aspect": "equal"})

    # Hydrophobicity color scale (Kyte-Doolittle)
    hydro = {
        "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
        "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
        "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
        "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
    }
    h_min, h_max = -4.5, 4.5

    angle_per_residue = 100.0  # degrees for alpha helix
    n = len(sequence)
    radius = 1.0

    for i, aa in enumerate(sequence.upper()):
        angle_deg = i * angle_per_residue
        angle_rad = np.radians(angle_deg)
        x = radius * np.cos(angle_rad)
        y = radius * np.sin(angle_rad)

        # Color by hydrophobicity
        h_val = hydro.get(aa, 0.0)
        norm_val = (h_val - h_min) / (h_max - h_min)
        color = plt.cm.RdYlBu_r(norm_val)

        circle = plt.Circle((x, y), 0.12, facecolor=color, edgecolor="black", linewidth=1)
        ax.add_patch(circle)
        ax.text(x, y, f"{aa}{i + 1}", ha="center", va="center", fontsize=7, fontweight="bold")

        # Draw connection line to next residue
        if i < n - 1:
            next_angle = np.radians((i + 1) * angle_per_residue)
            nx = radius * np.cos(next_angle)
            ny = radius * np.sin(next_angle)
            ax.plot([x, nx], [y, ny], "k-", alpha=0.2, linewidth=0.5)

    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_title("Helical Wheel Diagram")
    ax.set_xticks([])
    ax.set_yticks([])

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(str(output_path), dpi=300, bbox_inches="tight")
        logger.info(f"Helical wheel saved to {output_path}")

    return ax


def plot_msa_heatmap(
    aligned_sequences: List[str],
    labels: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (14, 6),
    **kwargs,
) -> Axes:
    """Plot a heatmap visualization of a multiple sequence alignment.

    Each cell is colored by amino acid property (hydrophobic, charged,
    polar, special). Gaps are shown in white.

    Args:
        aligned_sequences: List of aligned sequences (same length, with gaps '-')
        labels: Optional sequence labels
        ax: Optional matplotlib axes
        output_path: Optional save path
        figsize: Figure size

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(aligned_sequences, list, "aligned_sequences")
    validation.validate_not_empty(aligned_sequences, "aligned_sequences")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    n_seqs = len(aligned_sequences)
    aln_len = len(aligned_sequences[0])

    # Property groups
    aa_group = {}
    for aa in "AILMFWVP":
        aa_group[aa] = 0  # hydrophobic
    for aa in "DEKRH":
        aa_group[aa] = 1  # charged
    for aa in "STNQYC":
        aa_group[aa] = 2  # polar
    for aa in "G":
        aa_group[aa] = 3  # special

    # Build matrix
    matrix = np.full((n_seqs, aln_len), np.nan)
    for i, seq in enumerate(aligned_sequences):
        for j, aa in enumerate(seq.upper()):
            if aa != "-":
                matrix[i, j] = aa_group.get(aa, 2)

    cmap = plt.cm.get_cmap("Set2", 4)
    cmap.set_bad("white")

    ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5, interpolation="nearest")

    if labels:
        ax.set_yticks(range(n_seqs))
        ax.set_yticklabels(labels[:n_seqs], fontsize=8)
    ax.set_xlabel("Alignment Position")
    ax.set_ylabel("Sequence")
    ax.set_title("Multiple Sequence Alignment")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(str(output_path), dpi=300, bbox_inches="tight")
        logger.info(f"MSA heatmap saved to {output_path}")

    return ax


def plot_property_distribution(
    properties: Dict[str, List[float]],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    kind: str = "box",
    **kwargs,
) -> Axes:
    """Plot distribution of protein properties across multiple proteins.

    Args:
        properties: Dictionary mapping property names to value lists
        ax: Optional matplotlib axes
        output_path: Optional save path
        figsize: Figure size
        kind: Plot type ('box', 'violin')

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(properties, dict, "properties")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    names = list(properties.keys())
    data = [properties[n] for n in names]

    if kind == "violin" and HAS_SEABORN:
        import pandas as pd

        rows = []
        for name, values in properties.items():
            for v in values:
                rows.append({"Property": name, "Value": v})
        df = pd.DataFrame(rows)
        sns.violinplot(data=df, x="Property", y="Value", ax=ax)
    else:
        ax.boxplot(data, labels=names)

    ax.set_ylabel("Value")
    ax.set_title("Protein Property Distributions")
    ax.grid(True, alpha=0.3, axis="y")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(str(output_path), dpi=300, bbox_inches="tight")
        logger.info(f"Property distribution saved to {output_path}")

    return ax
