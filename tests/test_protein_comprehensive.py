"""Comprehensive tests for the protein module.

Tests all submodules: sequences, alignment, structure, contacts,
secondary structure, visualization, and orchestration.
All tests use real implementations per NO_MOCKING policy.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

# ============================================================
# SEQUENCE MODULE TESTS
# ============================================================


class TestSequenceProperties:
    """Tests for protein sequence property calculations."""

    def test_molecular_weight_known_protein(self):
        """Test MW for a short known sequence."""
        from metainformant.protein.sequence.sequences import molecular_weight

        # Glycine tripeptide: G-G-G
        mw = molecular_weight("GGG")
        assert mw > 0
        assert abs(mw - 57.0519 * 3) < 1.0  # 3 Gly residues

    def test_molecular_weight_empty(self):
        from metainformant.protein.sequence.sequences import molecular_weight

        assert molecular_weight("") == 0.0

    def test_isoelectric_point_acidic(self):
        """Aspartate-rich sequence should have acidic pI."""
        from metainformant.protein.sequence.sequences import isoelectric_point

        pI = isoelectric_point("DDDDDDDDD")
        assert 2.0 < pI < 5.0

    def test_isoelectric_point_basic(self):
        """Lysine-rich sequence should have basic pI."""
        from metainformant.protein.sequence.sequences import isoelectric_point

        pI = isoelectric_point("KKKKKKKKK")
        assert 9.0 < pI < 12.0

    def test_isoelectric_point_neutral(self):
        """Balanced sequence should have neutral-ish pI."""
        from metainformant.protein.sequence.sequences import isoelectric_point

        pI = isoelectric_point("AAAAAAAAAA")
        assert 4.0 < pI < 10.0

    def test_isoelectric_point_empty(self):
        from metainformant.protein.sequence.sequences import isoelectric_point

        assert isoelectric_point("") == 0.0

    def test_gravy_hydrophobic(self):
        """All-isoleucine sequence should be highly hydrophobic."""
        from metainformant.protein.sequence.sequences import gravy

        score = gravy("IIIIII")
        assert score > 3.0

    def test_gravy_hydrophilic(self):
        """All-arginine sequence should be hydrophilic."""
        from metainformant.protein.sequence.sequences import gravy

        score = gravy("RRRRRR")
        assert score < -3.0

    def test_gravy_empty(self):
        from metainformant.protein.sequence.sequences import gravy

        assert gravy("") == 0.0

    def test_extinction_coefficient(self):
        """Test extinction coefficient calculation."""
        from metainformant.protein.sequence.sequences import extinction_coefficient

        # Sequence with 2 Trp and 3 Tyr
        result = extinction_coefficient("WWYYYA")
        assert result["reduced"] == 2 * 5500 + 3 * 1490
        assert result["n_trp"] == 2
        assert result["n_tyr"] == 3

    def test_extinction_coefficient_no_absorbers(self):
        from metainformant.protein.sequence.sequences import extinction_coefficient

        result = extinction_coefficient("AAAA")
        assert result["reduced"] == 0.0

    def test_instability_index(self):
        """Test instability index computation."""
        from metainformant.protein.sequence.sequences import instability_index

        ii = instability_index("MKTIIALSY")
        assert isinstance(ii, float)
        assert ii > 0

    def test_instability_index_short(self):
        from metainformant.protein.sequence.sequences import instability_index

        assert instability_index("A") == 0.0

    def test_aromaticity(self):
        """Test aromaticity index."""
        from metainformant.protein.sequence.sequences import aromaticity

        # 3 out of 6 are aromatic
        arom = aromaticity("FWYAAA")
        assert abs(arom - 0.5) < 0.01

    def test_aromaticity_none(self):
        from metainformant.protein.sequence.sequences import aromaticity

        assert aromaticity("AAAA") == 0.0

    def test_charge_at_ph(self):
        """Test charge calculation at different pH values."""
        from metainformant.protein.sequence.sequences import charge_at_ph

        # At very low pH, most groups are protonated -> positive charge
        charge_low = charge_at_ph("ACDE", ph=1.0)
        # At very high pH, most groups are deprotonated -> negative charge
        charge_high = charge_at_ph("ACDE", ph=14.0)
        assert charge_low > charge_high

    def test_sequence_summary(self):
        """Test comprehensive sequence summary."""
        from metainformant.protein.sequence.sequences import sequence_summary

        summary = sequence_summary("MKTIIALSY")
        assert "length" in summary
        assert "molecular_weight" in summary
        assert "isoelectric_point" in summary
        assert "gravy" in summary
        assert summary["length"] == 9

    def test_amino_acid_composition(self):
        from metainformant.protein.sequence.sequences import amino_acid_composition

        comp = amino_acid_composition("AACC")
        assert abs(comp["A"] - 50.0) < 0.1
        assert abs(comp["C"] - 50.0) < 0.1

    def test_validate_protein_sequence(self):
        from metainformant.protein.sequence.sequences import validate_protein_sequence

        assert validate_protein_sequence("ACDEFG") is True
        assert validate_protein_sequence("") is False
        assert validate_protein_sequence("ACDE12") is False

    def test_hydropathy_score(self):
        from metainformant.protein.sequence.sequences import hydropathy_score

        # Sequence shorter than window returns empty
        scores = hydropathy_score("AAAA", window_size=19)
        assert scores == []

        # Long enough sequence returns scores
        scores = hydropathy_score("A" * 25, window_size=5)
        assert len(scores) == 21

    def test_transmembrane_regions(self):
        from metainformant.protein.sequence.sequences import transmembrane_regions

        # Very hydrophobic stretch should be detected
        seq = "AAAA" * 5 + "I" * 30 + "AAAA" * 5
        regions = transmembrane_regions(seq, threshold=1.6)
        # May or may not detect depending on exact scoring
        assert isinstance(regions, list)


class TestFastaIO:
    """Tests for FASTA file reading and writing."""

    def test_write_and_read_fasta(self, tmp_path: Path):
        from metainformant.protein.sequence.sequences import read_fasta, write_fasta

        seqs = {"prot1": "MKTIIALSY", "prot2": "ACDEFGHIK"}
        fasta_path = tmp_path / "test.fasta"

        write_fasta(seqs, fasta_path)
        assert fasta_path.exists()

        loaded = read_fasta(fasta_path)
        assert loaded == seqs

    def test_write_fasta_line_wrapping(self, tmp_path: Path):
        from metainformant.protein.sequence.sequences import read_fasta, write_fasta

        long_seq = "A" * 120
        seqs = {"long": long_seq}
        fasta_path = tmp_path / "wrapped.fasta"

        write_fasta(seqs, fasta_path, line_width=60)
        content = fasta_path.read_text()
        lines = [l for l in content.strip().split("\n") if not l.startswith(">")]
        assert all(len(l) <= 60 for l in lines)

        loaded = read_fasta(fasta_path)
        assert loaded["long"] == long_seq

    def test_read_fasta_missing_file(self):
        from metainformant.protein.sequence.sequences import read_fasta

        with pytest.raises(FileNotFoundError):
            read_fasta("/nonexistent/path.fasta")

    def test_find_motifs(self):
        from metainformant.protein.sequence.sequences import find_motifs

        results = find_motifs("MKTIIALSY", ["MKT", "AL"])
        assert 0 in results["MKT"]
        assert 5 in results["AL"]


# ============================================================
# ALIGNMENT MODULE TESTS
# ============================================================


class TestBLOSUM62Alignment:
    """Tests for BLOSUM62 matrix-based alignment."""

    def test_blosum62_matrix_exists(self):
        from metainformant.protein.sequence.alignment import BLOSUM62

        assert "A" in BLOSUM62
        assert BLOSUM62["A"]["A"] == 4  # Alanine self-score
        assert BLOSUM62["W"]["W"] == 11  # Tryptophan self-score

    def test_blosum62_symmetric(self):
        from metainformant.protein.sequence.alignment import BLOSUM62

        for aa1 in BLOSUM62:
            for aa2 in BLOSUM62[aa1]:
                assert BLOSUM62[aa1][aa2] == BLOSUM62[aa2][aa1]

    def test_matrix_align_identical(self):
        from metainformant.protein.sequence.alignment import matrix_align

        result = matrix_align("ACDEFG", "ACDEFG", mode="global")
        assert result["identity"] == 1.0
        assert result["score"] > 0

    def test_matrix_align_with_gaps(self):
        from metainformant.protein.sequence.alignment import matrix_align

        result = matrix_align("ACDEFG", "ACDFG", mode="global")
        assert result["score"] > 0
        assert len(result["aligned_seq1"]) == len(result["aligned_seq2"])

    def test_matrix_align_local(self):
        from metainformant.protein.sequence.alignment import matrix_align

        result = matrix_align("XXXACDEFGYYY", "ACDEFG", mode="local")
        assert result["score"] > 0
        assert "ACDEFG" in result["aligned_seq1"] or result["identity"] > 0.5

    def test_matrix_align_empty(self):
        from metainformant.protein.sequence.alignment import matrix_align

        result = matrix_align("", "ACDE", mode="global")
        assert result["score"] == 0

    def test_matrix_align_score_higher_with_blosum(self):
        """Similar amino acids should score better with BLOSUM62."""
        from metainformant.protein.sequence.alignment import BLOSUM62

        # I and L are similar (both hydrophobic), should have positive score
        assert BLOSUM62["I"]["L"] > 0
        # I and D are different, should have negative score
        assert BLOSUM62["I"]["D"] < 0


class TestMultiSequenceAlignment:
    """Tests for multiple sequence alignment."""

    def test_msa_two_sequences(self):
        from metainformant.protein.sequence.alignment import multi_sequence_alignment

        result = multi_sequence_alignment(["ACDE", "ACDE"])
        assert result["n_sequences"] == 2
        assert len(result["aligned_sequences"]) == 2

    def test_msa_three_sequences(self):
        from metainformant.protein.sequence.alignment import multi_sequence_alignment

        result = multi_sequence_alignment(["ACDEFG", "ACDFG", "ACDEFG"])
        assert result["n_sequences"] == 3
        assert all(len(s) == result["alignment_length"] for s in result["aligned_sequences"])

    def test_msa_single_sequence(self):
        from metainformant.protein.sequence.alignment import multi_sequence_alignment

        result = multi_sequence_alignment(["ACDE"])
        assert result["n_sequences"] == 1

    def test_msa_empty(self):
        from metainformant.protein.sequence.alignment import multi_sequence_alignment

        result = multi_sequence_alignment([])
        assert result["n_sequences"] == 0

    def test_alignment_statistics(self):
        from metainformant.protein.sequence.alignment import alignment_statistics, global_align

        alignment = global_align("ACDEFG", "ACDGFG")
        stats = alignment_statistics(alignment)
        assert "identity" in stats
        assert "gaps_seq1" in stats
        assert stats["aligned_length"] == len(alignment["aligned_seq1"])


# ============================================================
# STRUCTURE MODULE TESTS
# ============================================================


class TestStructuralGeometry:
    """Tests for structural geometry calculations."""

    def test_radius_of_gyration(self):
        from metainformant.protein.structure.general import calculate_radius_of_gyration

        # Atoms spread over 10 Angstroms
        coords = np.array([[0, 0, 0], [10, 0, 0], [0, 10, 0], [0, 0, 10]], dtype=float)
        rg = calculate_radius_of_gyration(coords)
        assert rg > 0

    def test_center_of_mass_equal_masses(self):
        from metainformant.protein.structure.general import calculate_center_of_mass

        coords = np.array([[0, 0, 0], [2, 0, 0]], dtype=float)
        com = calculate_center_of_mass(coords)
        np.testing.assert_allclose(com, [1.0, 0.0, 0.0])

    def test_center_of_mass_weighted(self):
        from metainformant.protein.structure.general import calculate_center_of_mass

        coords = np.array([[0, 0, 0], [10, 0, 0]], dtype=float)
        masses = np.array([1.0, 3.0])
        com = calculate_center_of_mass(coords, masses)
        np.testing.assert_allclose(com, [7.5, 0.0, 0.0])

    def test_inertia_tensor_shape(self):
        from metainformant.protein.structure.general import calculate_inertia_tensor

        coords = np.random.rand(20, 3) * 10
        tensor = calculate_inertia_tensor(coords)
        assert tensor.shape == (3, 3)
        # Inertia tensor should be symmetric
        np.testing.assert_allclose(tensor, tensor.T, atol=1e-10)

    def test_principal_axes(self):
        from metainformant.protein.structure.general import find_principal_axes

        coords = np.random.rand(30, 3) * 10
        eigenvals, eigenvecs = find_principal_axes(coords)
        assert len(eigenvals) == 3
        assert eigenvecs.shape == (3, 3)
        # Eigenvalues should be non-negative
        assert all(v >= -1e-10 for v in eigenvals)

    def test_structural_statistics(self):
        from metainformant.protein.structure.general import calculate_structural_statistics

        coords = np.random.rand(50, 3) * 20
        stats = calculate_structural_statistics(coords)
        assert "n_atoms" in stats
        assert stats["n_atoms"] == 50
        assert "radius_of_gyration" in stats
        assert "bounding_box" in stats

    def test_rmsd_identical(self):
        from metainformant.protein.structure.general import compute_rmsd_kabsch

        coords = np.random.rand(20, 3) * 10
        rmsd = compute_rmsd_kabsch(coords, coords.copy())
        assert rmsd < 1e-8

    def test_rmsd_translated(self):
        from metainformant.protein.structure.general import compute_rmsd_kabsch

        coords = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
        translated = coords + np.array([100, 200, 300])
        rmsd = compute_rmsd_kabsch(coords, translated)
        assert rmsd < 1e-6

    def test_rmsd_shape_mismatch(self):
        from metainformant.protein.structure.general import compute_rmsd_kabsch

        a = np.random.rand(10, 3)
        b = np.random.rand(5, 3)
        with pytest.raises(ValueError):
            compute_rmsd_kabsch(a, b)

    def test_align_structures(self):
        from metainformant.protein.structure.general import align_structures_kabsch

        coords = np.random.rand(15, 3) * 10
        # Rotate around Z axis
        theta = np.pi / 4
        R = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
        rotated = (coords - coords.mean(axis=0)) @ R.T + np.array([5, 5, 5])
        aligned, rotation, rmsd = align_structures_kabsch(coords, rotated)
        assert rmsd < 1e-5
        assert aligned.shape == coords.shape

    def test_sasa_calculation(self):
        from metainformant.protein.structure.general import calculate_solvent_accessible_surface_area

        coords = np.random.rand(30, 3) * 20
        area = calculate_solvent_accessible_surface_area(coords)
        assert area > 0


class TestContactMap:
    """Tests for contact map calculations."""

    def test_contact_map_basic(self):
        from metainformant.protein.structure.analysis import calculate_contact_map

        coords = np.array([[0, 0, 0], [3, 0, 0], [100, 0, 0]], dtype=float)
        cmap = calculate_contact_map(coords, threshold=5.0)
        assert cmap[0, 1] == 1  # Close
        assert cmap[0, 2] == 0  # Far

    def test_contact_map_diagonal_zero(self):
        from metainformant.protein.structure.analysis import calculate_contact_map

        coords = np.random.rand(10, 3)
        cmap = calculate_contact_map(coords, threshold=100.0)
        for i in range(10):
            assert cmap[i, i] == 0

    def test_contact_map_symmetric(self):
        from metainformant.protein.structure.analysis import calculate_contact_map

        coords = np.random.rand(15, 3) * 10
        cmap = calculate_contact_map(coords, threshold=5.0)
        np.testing.assert_array_equal(cmap, cmap.T)


# ============================================================
# SECONDARY STRUCTURE TESTS
# ============================================================


class TestSecondaryStructure:
    """Tests for secondary structure prediction and analysis."""

    def test_simple_prediction_length(self):
        from metainformant.protein.structure.secondary import predict_secondary_structure

        seq = "MKTIIALSY"
        ss = predict_secondary_structure(seq, method="simple")
        assert len(ss) == len(seq)
        assert all(s in ("H", "E", "C") for s in ss)

    def test_prediction_empty(self):
        from metainformant.protein.structure.secondary import predict_secondary_structure

        assert predict_secondary_structure("") == []

    def test_invalid_method(self):
        from metainformant.protein.structure.secondary import predict_secondary_structure

        with pytest.raises(ValueError):
            predict_secondary_structure("ACDE", method="invalid")

    def test_ss_composition(self):
        from metainformant.protein.structure.secondary import calculate_ss_composition

        comp = calculate_ss_composition(["H", "H", "E", "C", "C", "C"])
        assert abs(comp["helix"] - 2 / 6) < 0.01
        assert abs(comp["sheet"] - 1 / 6) < 0.01
        assert abs(comp["coil"] - 3 / 6) < 0.01

    def test_ss_elements(self):
        from metainformant.protein.structure.secondary import identify_ss_elements

        ss = ["H", "H", "H", "H", "C", "C", "E", "E", "E"]
        elements = identify_ss_elements(ss, min_length=3)
        types = [e["type"] for e in elements]
        assert "H" in types
        assert "E" in types

    def test_compare_ss(self):
        from metainformant.protein.structure.secondary import compare_ss_predictions

        pred1 = ["H", "H", "E", "C"]
        pred2 = ["H", "E", "E", "C"]
        result = compare_ss_predictions(pred1, pred2)
        assert result["accuracy"] == 3 / 4

    def test_ss_to_dssp(self):
        from metainformant.protein.structure.secondary import ss_to_dSSP_format

        assert ss_to_dSSP_format(["H", "E", "C"]) == "HE "

    def test_validate_ss(self):
        from metainformant.protein.structure.secondary import validate_ss_prediction

        result = validate_ss_prediction(["H", "E", "C"], "ABC")
        assert result["valid_length"] is True
        assert result["valid_states"] is True

    def test_ss_propensities(self):
        from metainformant.protein.structure.secondary import calculate_ss_propensities

        props = calculate_ss_propensities("AEIVG")
        assert "A" in props
        assert "helix" in props["A"]
        assert "sheet" in props["A"]

    def test_tm_prediction(self):
        from metainformant.protein.structure.secondary import predict_transmembrane_regions

        # Hydrophobic stretch
        seq = "ILLLLLLLLLLLLLLLLLLLLL"
        regions = predict_transmembrane_regions(seq)
        assert isinstance(regions, list)


# ============================================================
# STRUCTURE ANALYSIS TESTS
# ============================================================


class TestStructureAnalysis:
    """Tests for structure-level analysis."""

    def test_identify_domains(self):
        from metainformant.protein.structure.analysis import identify_domains

        structure = {
            "atoms": [
                {"chain_id": "A", "res_seq": 1},
                {"chain_id": "A", "res_seq": 2},
                {"chain_id": "A", "res_seq": 3},
                {"chain_id": "A", "res_seq": 50},  # gap > 10
                {"chain_id": "A", "res_seq": 51},
            ]
        }
        domains = identify_domains(structure)
        assert len(domains) == 2

    def test_identify_domains_empty(self):
        from metainformant.protein.structure.analysis import identify_domains

        assert identify_domains({"atoms": []}) == []

    def test_surface_area(self):
        from metainformant.protein.structure.analysis import calculate_surface_area

        coords = np.random.rand(30, 3) * 20
        area = calculate_surface_area(coords)
        assert area > 0

    def test_surface_area_small(self):
        from metainformant.protein.structure.analysis import calculate_surface_area

        coords = np.array([[0, 0, 0], [1, 0, 0]], dtype=float)
        area = calculate_surface_area(coords, probe_radius=1.4)
        assert area > 0

    def test_protein_flexibility(self):
        from metainformant.protein.structure.analysis import analyze_protein_flexibility

        structure = {
            "atoms": [
                {"chain_id": "A", "res_seq": 1, "temp_factor": 10.0},
                {"chain_id": "A", "res_seq": 1, "temp_factor": 12.0},
                {"chain_id": "A", "res_seq": 2, "temp_factor": 50.0},
                {"chain_id": "A", "res_seq": 2, "temp_factor": 55.0},
            ],
            "coordinates": np.random.rand(4, 3),
        }
        flex = analyze_protein_flexibility(structure)
        assert "flexibility_score" in flex
        assert flex["flexibility_score"] > 0


# ============================================================
# ORCHESTRATION TESTS
# ============================================================


class TestOrchestration:
    """Tests for protein orchestration workflows."""

    def test_analyze_protein_sequence(self):
        from metainformant.protein.orchestration import analyze_protein_sequence

        result = analyze_protein_sequence("MKTIIALSY", name="test_protein")
        assert result["name"] == "test_protein"
        assert result["length"] == 9
        assert "physicochemical" in result
        assert "composition" in result
        assert result["physicochemical"]["molecular_weight"] > 0

    def test_analyze_sequence_with_ss(self):
        from metainformant.protein.orchestration import analyze_protein_sequence

        result = analyze_protein_sequence("MKTIIALSY", predict_ss=True)
        assert "secondary_structure" in result
        assert len(result["secondary_structure"]["prediction"]) == 9

    def test_analyze_sequence_with_motifs(self):
        from metainformant.protein.orchestration import analyze_protein_sequence

        result = analyze_protein_sequence("MKTIIALSY", find_motifs=["MKT", "AL"])
        assert "motifs" in result
        assert 0 in result["motifs"]["MKT"]

    def test_batch_analyze(self):
        from metainformant.protein.orchestration import batch_analyze_sequences

        seqs = {"p1": "MKTIIALSY", "p2": "ACDEFGHIK"}
        results = batch_analyze_sequences(seqs)
        assert "p1" in results
        assert "p2" in results
        assert results["p1"]["length"] == 9

    def test_comparative_analysis(self):
        from metainformant.protein.orchestration import comparative_analysis

        seqs = {"p1": "ACDEFG", "p2": "ACDGFG", "p3": "ACDEFG"}
        results = comparative_analysis(seqs, use_blosum=False)
        assert "pairwise_alignments" in results
        assert "identity_matrix" in results
        assert "msa" in results

    def test_comparative_analysis_blosum(self):
        from metainformant.protein.orchestration import comparative_analysis

        seqs = {"p1": "ACDEFG", "p2": "ACDGFG"}
        results = comparative_analysis(seqs, use_blosum=True)
        assert "pairwise_alignments" in results

    def test_comparative_analysis_single(self):
        from metainformant.protein.orchestration import comparative_analysis

        seqs = {"p1": "ACDEFG"}
        results = comparative_analysis(seqs)
        assert "error" in results

    def test_analyze_from_fasta(self, tmp_path: Path):
        from metainformant.protein.orchestration import analyze_from_fasta
        from metainformant.protein.sequence.sequences import write_fasta

        seqs = {"prot1": "MKTIIALSY", "prot2": "ACDEFGHIK"}
        fasta_path = tmp_path / "test.fasta"
        write_fasta(seqs, fasta_path)

        results = analyze_from_fasta(fasta_path, predict_ss=False, compare=True)
        assert results["n_sequences"] == 2
        assert "sequences" in results
        assert "comparative" in results

    def test_analyze_short_sequence(self):
        """Short sequences should not crash hydropathy."""
        from metainformant.protein.orchestration import analyze_protein_sequence

        result = analyze_protein_sequence("MK", name="short")
        assert result["hydropathy"] == []
        assert result["transmembrane_regions"] == []


# ============================================================
# PDB I/O TESTS
# ============================================================


class TestPDBIO:
    """Tests for PDB file parsing and writing."""

    def test_parse_pdb_file(self, tmp_path: Path):
        from metainformant.protein.structure.io import parse_pdb_file

        pdb_content = (
            "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N\n"
            "ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00 12.00           C\n"
            "ATOM      3  C   ALA A   1       3.000   4.000   5.000  1.00 11.00           C\n"
            "END\n"
        )
        pdb_path = tmp_path / "test.pdb"
        pdb_path.write_text(pdb_content)

        result = parse_pdb_file(pdb_path)
        assert "atoms" in result
        assert len(result["atoms"]) >= 1

    def test_validate_pdb(self, tmp_path: Path):
        from metainformant.protein.structure.io import validate_pdb_file

        pdb_content = "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N\n" "END\n"
        pdb_path = tmp_path / "test.pdb"
        pdb_path.write_text(pdb_content)

        result = validate_pdb_file(pdb_path)
        assert isinstance(result, (dict, tuple))


# ============================================================
# PROTEOME TESTS
# ============================================================


class TestProteomes:
    """Tests for proteome utilities."""

    def test_read_taxon_ids(self, tmp_path: Path):
        from metainformant.protein.sequence.proteomes import read_taxon_ids

        taxon_file = tmp_path / "taxon.txt"
        taxon_file.write_text("9606\n10090\n# Comment\n\n10116\n")
        result = read_taxon_ids(taxon_file)
        # Note: read_taxon_ids returns strings
        assert result == ["9606", "10090", "10116"]

    def test_validate_taxon_ids(self):
        from metainformant.protein.sequence.proteomes import validate_taxon_ids

        # validate_taxon_ids expects string inputs
        valid, invalid = validate_taxon_ids(["9606", "10090", "abc"])
        assert "9606" in valid
        assert "abc" in invalid
