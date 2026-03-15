"""Tests for protein module enhancements.

Tests new and improved functions across DSSP parsing, vectorized contacts,
new visualizations, orchestration workflows, and utility functions.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# DSSP Parser Tests
# ---------------------------------------------------------------------------


class TestDSSPParser:
    """Tests for the DSSP file parser in secondary.py."""

    def _make_dssp_content(self) -> str:
        """Build a minimal valid DSSP file for testing."""
        header = "==== Secondary Structure Definition by the program DSSP ====\n"
        header += "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC\n"
        # Simplified DSSP residue lines (fixed-width columns)
        # Col  1-5:  DSSP residue number
        # Col  6-10: PDB residue number
        # Col 12:    chain ID
        # Col 14:    amino acid
        # Col 17:    secondary structure
        # Col 36-38: accessibility
        lines = [
            "    1    1 A M  H  <         0   0  100",
            "    2    2 A K  H  >         0   0   80",
            "    3    3 A L  H  >         0   0   60",
            "    4    4 A V  E  <         0   0   50",
            "    5    5 A I  E  <         0   0   40",
            "    6    6 A S     <         0   0  120",
        ]
        return header + "\n".join(lines)

    def test_parse_dssp_basic(self):
        from metainformant.protein.structure.secondary import parse_dssp_file

        content = self._make_dssp_content()
        result = parse_dssp_file(content)

        assert isinstance(result, dict)
        assert "sequence" in result
        assert "secondary_structure" in result
        assert "accessibility" in result
        assert "residues" in result
        assert "total_residues" in result

    def test_parse_dssp_sequence_extraction(self):
        from metainformant.protein.structure.secondary import parse_dssp_file

        content = self._make_dssp_content()
        result = parse_dssp_file(content)

        assert result["total_residues"] == 6
        assert result["sequence"] == "MKLVIS"

    def test_parse_dssp_ss_normalization(self):
        from metainformant.protein.structure.secondary import parse_dssp_file

        content = self._make_dssp_content()
        result = parse_dssp_file(content)

        ss = result["secondary_structure"]
        assert len(ss) == 6
        # H stays H, E stays E, blank becomes C
        assert ss[0] == "H"
        assert ss[3] == "E"
        assert ss[5] == "C"

    def test_parse_dssp_accessibility(self):
        from metainformant.protein.structure.secondary import parse_dssp_file

        content = self._make_dssp_content()
        result = parse_dssp_file(content)

        acc = result["accessibility"]
        assert len(acc) == 6
        # All values should be numeric
        assert all(isinstance(a, float) for a in acc)

    def test_parse_dssp_residue_details(self):
        from metainformant.protein.structure.secondary import parse_dssp_file

        content = self._make_dssp_content()
        result = parse_dssp_file(content)

        residues = result["residues"]
        assert len(residues) == 6
        assert residues[0]["amino_acid"] == "M"
        assert residues[0]["chain_id"] == "A"

    def test_parse_dssp_empty_content(self):
        from metainformant.protein.structure.secondary import parse_dssp_file

        result = parse_dssp_file("")
        assert result["total_residues"] == 0
        assert result["sequence"] == ""

    def test_parse_dssp_no_data_lines(self):
        from metainformant.protein.structure.secondary import parse_dssp_file

        content = "HEADER ONLY\nNO DATA HERE\n"
        result = parse_dssp_file(content)
        assert result["total_residues"] == 0


# ---------------------------------------------------------------------------
# Vectorized Contact Pairs Tests
# ---------------------------------------------------------------------------


class TestVectorizedContacts:
    """Tests for vectorized compute_ca_contact_pairs."""

    def test_basic_contact_detection(self):
        from metainformant.protein.structure.contacts import compute_ca_contact_pairs

        coords = [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (100.0, 0.0, 0.0)]
        pairs = compute_ca_contact_pairs(coords, threshold=4.0)
        assert (0, 1) in pairs
        assert (0, 2) not in pairs
        assert (1, 2) not in pairs

    def test_all_within_threshold(self):
        from metainformant.protein.structure.contacts import compute_ca_contact_pairs

        coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, 1.0, 0.0)]
        pairs = compute_ca_contact_pairs(coords, threshold=10.0)
        assert len(pairs) == 3  # All pairs

    def test_none_within_threshold(self):
        from metainformant.protein.structure.contacts import compute_ca_contact_pairs

        coords = [(0.0, 0.0, 0.0), (100.0, 0.0, 0.0), (200.0, 0.0, 0.0)]
        pairs = compute_ca_contact_pairs(coords, threshold=1.0)
        assert len(pairs) == 0

    def test_single_atom(self):
        from metainformant.protein.structure.contacts import compute_ca_contact_pairs

        pairs = compute_ca_contact_pairs([(0.0, 0.0, 0.0)], threshold=8.0)
        assert pairs == []

    def test_empty_coords(self):
        from metainformant.protein.structure.contacts import compute_ca_contact_pairs

        pairs = compute_ca_contact_pairs([], threshold=8.0)
        assert pairs == []

    def test_exact_threshold_boundary(self):
        from metainformant.protein.structure.contacts import compute_ca_contact_pairs

        coords = [(0.0, 0.0, 0.0), (8.0, 0.0, 0.0)]
        pairs = compute_ca_contact_pairs(coords, threshold=8.0)
        assert (0, 1) in pairs  # Exactly at threshold should be included

    def test_large_set_performance(self):
        from metainformant.protein.structure.contacts import compute_ca_contact_pairs

        np.random.seed(42)
        coords = [(x, y, z) for x, y, z in np.random.rand(100, 3) * 50]
        pairs = compute_ca_contact_pairs(coords, threshold=5.0)
        assert isinstance(pairs, list)
        assert all(isinstance(p, tuple) and len(p) == 2 for p in pairs)
        # All pairs should have i < j
        assert all(p[0] < p[1] for p in pairs)

    def test_pairs_ordering(self):
        from metainformant.protein.structure.contacts import compute_ca_contact_pairs

        coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0), (3.0, 0.0, 0.0)]
        pairs = compute_ca_contact_pairs(coords, threshold=1.5)
        for i, j in pairs:
            assert i < j, f"Pair ({i},{j}) not ordered"


# ---------------------------------------------------------------------------
# Inertia Tensor Tests
# ---------------------------------------------------------------------------


class TestInertiaTensor:
    """Tests for vectorized inertia tensor calculation."""

    def test_inertia_tensor_shape(self):
        from metainformant.protein.structure.general import calculate_inertia_tensor

        coords = np.random.rand(50, 3) * 20
        tensor = calculate_inertia_tensor(coords)
        assert tensor.shape == (3, 3)

    def test_inertia_tensor_symmetric(self):
        from metainformant.protein.structure.general import calculate_inertia_tensor

        coords = np.random.rand(50, 3) * 20
        tensor = calculate_inertia_tensor(coords)
        np.testing.assert_allclose(tensor, tensor.T, atol=1e-10)

    def test_inertia_tensor_positive_diagonal(self):
        from metainformant.protein.structure.general import calculate_inertia_tensor

        coords = np.random.rand(50, 3) * 20
        tensor = calculate_inertia_tensor(coords)
        assert tensor[0, 0] >= 0
        assert tensor[1, 1] >= 0
        assert tensor[2, 2] >= 0

    def test_inertia_tensor_with_masses(self):
        from metainformant.protein.structure.general import calculate_inertia_tensor

        coords = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        masses = np.array([2.0, 3.0, 4.0])
        tensor = calculate_inertia_tensor(coords, masses)
        assert tensor.shape == (3, 3)
        assert tensor[0, 0] > 0


# ---------------------------------------------------------------------------
# New Visualization Tests
# ---------------------------------------------------------------------------


class TestHelicalWheel:
    """Tests for helical wheel visualization."""

    def test_helical_wheel_returns_axes(self):
        import matplotlib

        matplotlib.use("Agg")
        from metainformant.protein.visualization.general import plot_helical_wheel

        ax = plot_helical_wheel("AKLVFEQHIND")
        from matplotlib.axes import Axes

        assert isinstance(ax, Axes)
        import matplotlib.pyplot as plt

        plt.close("all")

    def test_helical_wheel_short_sequence(self):
        import matplotlib

        matplotlib.use("Agg")
        from metainformant.protein.visualization.general import plot_helical_wheel

        ax = plot_helical_wheel("AKL")
        assert ax is not None
        import matplotlib.pyplot as plt

        plt.close("all")

    def test_helical_wheel_save(self, tmp_path):
        import matplotlib

        matplotlib.use("Agg")
        from metainformant.protein.visualization.general import plot_helical_wheel

        out = tmp_path / "wheel.png"
        ax = plot_helical_wheel("AKLVFEQHIND", output_path=out)
        assert out.exists()
        import matplotlib.pyplot as plt

        plt.close("all")


class TestMSAHeatmap:
    """Tests for MSA heatmap visualization."""

    def test_msa_heatmap_basic(self):
        import matplotlib

        matplotlib.use("Agg")
        from metainformant.protein.visualization.general import plot_msa_heatmap

        seqs = ["AKLV-E", "AKLVFE", "AKLV-E"]
        ax = plot_msa_heatmap(seqs, labels=["seq1", "seq2", "seq3"])
        assert ax is not None
        import matplotlib.pyplot as plt

        plt.close("all")

    def test_msa_heatmap_no_labels(self):
        import matplotlib

        matplotlib.use("Agg")
        from metainformant.protein.visualization.general import plot_msa_heatmap

        seqs = ["AKLV", "EKLV"]
        ax = plot_msa_heatmap(seqs)
        assert ax is not None
        import matplotlib.pyplot as plt

        plt.close("all")

    def test_msa_heatmap_save(self, tmp_path):
        import matplotlib

        matplotlib.use("Agg")
        from metainformant.protein.visualization.general import plot_msa_heatmap

        out = tmp_path / "msa.png"
        ax = plot_msa_heatmap(["AKLV", "EKLV"], output_path=out)
        assert out.exists()
        import matplotlib.pyplot as plt

        plt.close("all")


class TestPropertyDistribution:
    """Tests for property distribution visualization."""

    def test_property_distribution_box(self):
        import matplotlib

        matplotlib.use("Agg")
        from metainformant.protein.visualization.general import plot_property_distribution

        props = {"MW": [10.0, 20.0, 30.0], "pI": [5.0, 6.0, 7.0]}
        ax = plot_property_distribution(props, kind="box")
        assert ax is not None
        import matplotlib.pyplot as plt

        plt.close("all")

    def test_property_distribution_save(self, tmp_path):
        import matplotlib

        matplotlib.use("Agg")
        from metainformant.protein.visualization.general import plot_property_distribution

        out = tmp_path / "dist.png"
        props = {"values": [1.0, 2.0, 3.0, 4.0]}
        ax = plot_property_distribution(props, output_path=out)
        assert out.exists()
        import matplotlib.pyplot as plt

        plt.close("all")


# ---------------------------------------------------------------------------
# Orchestration Tests
# ---------------------------------------------------------------------------


def _write_test_pdb(path: Path, *, n_residues: int = 10, chain: str = "A") -> None:
    """Write a minimal PDB file for testing."""
    lines = []
    lines.append(f"HEADER    TEST STRUCTURE")
    serial = 1
    for i in range(1, n_residues + 1):
        x = float(i * 3.8)
        y = float(i * 0.5)
        z = float(i * 0.3)
        line = f"ATOM  {serial:5d}  CA  ALA {chain}{i:4d}    " f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 50.00           C  "
        lines.append(line)
        serial += 1
    lines.append("END")
    path.write_text("\n".join(lines))


class TestCompareStructures:
    """Tests for structure comparison orchestration."""

    def test_compare_two_structures(self, tmp_path):
        from metainformant.protein.workflow.orchestration import compare_structures

        pdb_a = tmp_path / "a.pdb"
        pdb_b = tmp_path / "b.pdb"
        _write_test_pdb(pdb_a, n_residues=10)
        _write_test_pdb(pdb_b, n_residues=10)

        result = compare_structures(pdb_a, pdb_b, align=True)
        assert "rmsd_aligned" in result
        assert result["n_compared_residues"] == 10
        assert result["rmsd_aligned"] >= 0.0

    def test_compare_structures_no_align(self, tmp_path):
        from metainformant.protein.workflow.orchestration import compare_structures

        pdb_a = tmp_path / "a.pdb"
        pdb_b = tmp_path / "b.pdb"
        _write_test_pdb(pdb_a, n_residues=8)
        _write_test_pdb(pdb_b, n_residues=8)

        result = compare_structures(pdb_a, pdb_b, align=False)
        assert "rmsd_simple" in result

    def test_compare_different_lengths(self, tmp_path):
        from metainformant.protein.workflow.orchestration import compare_structures

        pdb_a = tmp_path / "short.pdb"
        pdb_b = tmp_path / "long.pdb"
        _write_test_pdb(pdb_a, n_residues=5)
        _write_test_pdb(pdb_b, n_residues=15)

        result = compare_structures(pdb_a, pdb_b)
        assert result["n_compared_residues"] == 5
        assert result["n_ca_atoms_a"] == 5
        assert result["n_ca_atoms_b"] == 15

    def test_compare_contact_comparison(self, tmp_path):
        from metainformant.protein.workflow.orchestration import compare_structures

        pdb_a = tmp_path / "a.pdb"
        pdb_b = tmp_path / "b.pdb"
        _write_test_pdb(pdb_a, n_residues=10)
        _write_test_pdb(pdb_b, n_residues=10)

        result = compare_structures(pdb_a, pdb_b)
        assert "contact_comparison" in result
        cc = result["contact_comparison"]
        assert "shared_contacts" in cc
        assert "jaccard_similarity" in cc
        assert 0.0 <= cc["jaccard_similarity"] <= 1.0

    def test_compare_per_residue_distance(self, tmp_path):
        from metainformant.protein.workflow.orchestration import compare_structures

        pdb_a = tmp_path / "a.pdb"
        pdb_b = tmp_path / "b.pdb"
        _write_test_pdb(pdb_a, n_residues=10)
        _write_test_pdb(pdb_b, n_residues=10)

        result = compare_structures(pdb_a, pdb_b)
        prd = result["per_residue_distance"]
        assert prd["mean"] >= 0.0
        assert prd["max"] >= prd["min"]


class TestAssessAlphaFoldQuality:
    """Tests for AlphaFold quality assessment pipeline."""

    def test_assess_alphafold_quality(self, tmp_path):
        from metainformant.protein.workflow.orchestration import assess_alphafold_quality

        pdb = tmp_path / "af_model.pdb"
        _write_test_pdb(pdb, n_residues=20)
        result = assess_alphafold_quality(pdb)

        assert "confidence" in result
        assert "quality_rating" in result
        assert result["quality_rating"] in ("very_high", "confident", "low", "very_low")

    def test_assess_quality_structural_stats(self, tmp_path):
        from metainformant.protein.workflow.orchestration import assess_alphafold_quality

        pdb = tmp_path / "af_model.pdb"
        _write_test_pdb(pdb, n_residues=15)
        result = assess_alphafold_quality(pdb)

        assert "structural_stats" in result
        stats = result["structural_stats"]
        assert "radius_of_gyration" in stats
        assert "n_atoms" in stats


class TestFullProteinAnalysis:
    """Tests for combined sequence+structure analysis."""

    def test_full_analysis_sequence_only(self):
        from metainformant.protein.workflow.orchestration import full_protein_analysis

        result = full_protein_analysis(
            "MKLVFEQHINDASSTWYCGP",
            name="test_protein",
            predict_ss=True,
        )

        assert result["name"] == "test_protein"
        assert "sequence_analysis" in result
        assert "summary" in result
        assert result["summary"]["sequence_length"] == 20
        assert result["summary"]["has_structure"] is False

    def test_full_analysis_with_structure(self, tmp_path):
        from metainformant.protein.workflow.orchestration import full_protein_analysis

        pdb = tmp_path / "structure.pdb"
        _write_test_pdb(pdb, n_residues=10)

        result = full_protein_analysis(
            "AAAAAAAAAA",
            pdb_path=pdb,
            name="with_struct",
        )

        assert result["summary"]["has_structure"] is True
        assert "structure_analysis" in result

    def test_full_analysis_summary_keys(self):
        from metainformant.protein.workflow.orchestration import full_protein_analysis

        result = full_protein_analysis("ACDEFGHIKLM")
        summary = result["summary"]
        assert "molecular_weight" in summary
        assert "isoelectric_point" in summary
        assert summary["molecular_weight"] > 0


class TestBatchCompareStructures:
    """Tests for batch structure comparison."""

    def test_batch_compare_three_structures(self, tmp_path):
        from metainformant.protein.workflow.orchestration import batch_compare_structures

        paths = []
        for name in ("a", "b", "c"):
            p = tmp_path / f"{name}.pdb"
            _write_test_pdb(p, n_residues=8)
            paths.append(p)

        result = batch_compare_structures(paths)

        assert result["n_structures"] == 3
        assert len(result["pairwise_comparisons"]) == 3  # C(3,2) = 3
        assert len(result["rmsd_matrix"]) == 3
        assert len(result["rmsd_matrix"][0]) == 3

    def test_batch_compare_single_structure(self, tmp_path):
        from metainformant.protein.workflow.orchestration import batch_compare_structures

        p = tmp_path / "only.pdb"
        _write_test_pdb(p, n_residues=5)

        result = batch_compare_structures([p])
        assert "error" in result

    def test_batch_rmsd_matrix_symmetry(self, tmp_path):
        from metainformant.protein.workflow.orchestration import batch_compare_structures

        paths = []
        for name in ("x", "y"):
            p = tmp_path / f"{name}.pdb"
            _write_test_pdb(p, n_residues=10)
            paths.append(p)

        result = batch_compare_structures(paths)
        m = result["rmsd_matrix"]
        assert abs(m[0][1] - m[1][0]) < 1e-10


# ---------------------------------------------------------------------------
# Utility Function Tests
# ---------------------------------------------------------------------------


class TestKmerFrequencies:
    """Tests for kmer_frequencies function."""

    def test_kmer_basic(self):
        from metainformant.protein.sequence.sequences import kmer_frequencies

        result = kmer_frequencies("ACAC", k=2)
        assert result["AC"] == 2
        assert result["CA"] == 1

    def test_kmer_k1(self):
        from metainformant.protein.sequence.sequences import kmer_frequencies

        result = kmer_frequencies("AABBA", k=1)
        assert result["A"] == 3
        assert result["B"] == 2

    def test_kmer_too_large(self):
        from metainformant.protein.sequence.sequences import kmer_frequencies

        result = kmer_frequencies("AC", k=5)
        assert result == {}

    def test_kmer_k0(self):
        from metainformant.protein.sequence.sequences import kmer_frequencies

        result = kmer_frequencies("ABC", k=0)
        assert result == {}

    def test_kmer_full_length(self):
        from metainformant.protein.sequence.sequences import kmer_frequencies

        result = kmer_frequencies("ABCD", k=4)
        assert result == {"ABCD": 1}


class TestCalculateAAComposition:
    """Tests for calculate_aa_composition wrapper."""

    def test_composition_fractions(self):
        from metainformant.protein.sequence.sequences import calculate_aa_composition

        comp = calculate_aa_composition("AAAA")
        assert abs(comp.get("A", 0.0) - 1.0) < 1e-6

    def test_composition_sums_to_one(self):
        from metainformant.protein.sequence.sequences import calculate_aa_composition

        comp = calculate_aa_composition("ACDEFGHIKLMNPQRSTVWY")
        total = sum(comp.values())
        assert abs(total - 1.0) < 1e-6

    def test_composition_mixed(self):
        from metainformant.protein.sequence.sequences import calculate_aa_composition

        comp = calculate_aa_composition("AACCDD")
        assert abs(comp["A"] - 1 / 3) < 1e-6
        assert abs(comp["C"] - 1 / 3) < 1e-6
        assert abs(comp["D"] - 1 / 3) < 1e-6


class TestPairwiseIdentity:
    """Tests for pairwise_identity wrapper."""

    def test_identical_sequences(self):
        from metainformant.protein.sequence.alignment import pairwise_identity

        identity = pairwise_identity("ACDEF", "ACDEF")
        assert abs(identity - 1.0) < 1e-6

    def test_completely_different(self):
        from metainformant.protein.sequence.alignment import pairwise_identity

        identity = pairwise_identity("AAAA", "DDDD")
        assert identity < 0.5

    def test_partial_match(self):
        from metainformant.protein.sequence.alignment import pairwise_identity

        identity = pairwise_identity("ACDE", "ACGE")
        assert 0.0 < identity < 1.0


class TestNeedlemanWunsch:
    """Tests for needleman_wunsch wrapper."""

    def test_returns_tuple(self):
        from metainformant.protein.sequence.alignment import needleman_wunsch

        result = needleman_wunsch("ACDE", "ACGE")
        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_score_aligned_seqs(self):
        from metainformant.protein.sequence.alignment import needleman_wunsch

        score, aligned1, aligned2 = needleman_wunsch("ACDE", "ACDE")
        assert score > 0
        assert len(aligned1) == len(aligned2)

    def test_empty_sequences(self):
        from metainformant.protein.sequence.alignment import needleman_wunsch

        score, a1, a2 = needleman_wunsch("", "")
        assert score == 0


class TestReadPdbCaCoordinates:
    """Tests for read_pdb_ca_coordinates."""

    def test_basic_reading(self, tmp_path):
        from metainformant.protein.structure.io import read_pdb_ca_coordinates

        pdb = tmp_path / "test.pdb"
        _write_test_pdb(pdb, n_residues=5)

        coords = read_pdb_ca_coordinates(pdb)
        assert len(coords) == 5
        assert all(len(c) == 3 for c in coords)

    def test_no_ca_atoms(self, tmp_path):
        from metainformant.protein.structure.io import read_pdb_ca_coordinates

        pdb = tmp_path / "no_ca.pdb"
        # Write PDB with only CB atoms
        lines = ["HEADER    TEST"]
        for i in range(1, 4):
            lines.append(
                f"ATOM  {i:5d}  CB  ALA A{i:4d}    " f"{float(i):8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00           C  "
            )
        lines.append("END")
        pdb.write_text("\n".join(lines))

        coords = read_pdb_ca_coordinates(pdb)
        assert len(coords) == 0


# ---------------------------------------------------------------------------
# Package Init Tests
# ---------------------------------------------------------------------------


class TestPackageInits:
    """Tests for package __init__.py exports."""

    def test_database_init_imports(self):
        from metainformant.protein import database

        assert hasattr(database, "interpro")
        assert hasattr(database, "uniprot")

    def test_visualization_init_imports(self):
        from metainformant.protein import visualization

        assert hasattr(visualization, "general")

    def test_database_submodule_access(self):
        from metainformant.protein.database import interpro

        assert hasattr(interpro, "fetch_interpro_domains")

    def test_visualization_submodule_access(self):
        from metainformant.protein.visualization import general

        assert hasattr(general, "plot_helical_wheel")
        assert hasattr(general, "plot_msa_heatmap")
        assert hasattr(general, "plot_property_distribution")
