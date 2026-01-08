"""Tests for contamination detection functionality."""

import pytest

from metainformant.quality.analysis.contamination import (
    detect_adapter_contamination,
    detect_cross_species_contamination,
    detect_mycoplasma_contamination,
    detect_rrna_contamination,
    detect_vector_contamination,
    generate_contamination_report,
)


class TestCrossSpeciesContamination:
    """Test cross-species contamination detection."""

    def test_detect_cross_species_basic(self):
        """Test basic cross-species contamination detection."""
        sequences = ["ATCGATCG", "GCTAGCTA"]
        reference_genomes = {
            "human": "ATCGATCGATCG",
            "mouse": "GCTAGCTAGCTA",
        }

        results = detect_cross_species_contamination(sequences, reference_genomes, threshold=0.8)

        assert len(results) == 2
        assert "0" in results  # First sequence matches human
        assert "1" in results  # Second sequence matches mouse

    def test_detect_cross_species_no_contamination(self):
        """Test detection when no contamination is present."""
        sequences = ["ATCGATCG", "GCTAGCTA"]
        reference_genomes = {
            "human": "TTTTTTTT",
            "mouse": "CCCCCCCC",
        }

        results = detect_cross_species_contamination(sequences, reference_genomes, threshold=0.8)

        assert len(results) == 0  # No contamination detected

    def test_detect_cross_species_empty_input(self):
        """Test handling of empty input."""
        results = detect_cross_species_contamination([], {})

        assert len(results) == 0


class TestRRNAContamination:
    """Test rRNA contamination detection."""

    def test_detect_rrna_contamination_basic(self):
        """Test basic rRNA contamination detection."""
        sequences = [
            "ATCGATCG",  # No rRNA pattern
            "GGAAGGAGCAGTG",  # Contains rRNA pattern
        ]

        results = detect_rrna_contamination(sequences)

        assert "1" in results  # Second sequence should be flagged
        assert "0" not in results  # First sequence should not be flagged

    def test_detect_rrna_contamination_custom_patterns(self):
        """Test detection with custom rRNA patterns."""
        sequences = ["ATCGATCG", "CUSTOMPATTERN"]
        custom_patterns = ["CUSTOM"]

        results = detect_rrna_contamination(sequences, custom_patterns)

        assert "1" in results  # Second sequence should be flagged
        assert "0" not in results  # First sequence should not be flagged


class TestMycoplasmaContamination:
    """Test mycoplasma contamination detection."""

    def test_detect_mycoplasma_contamination_basic(self):
        """Test basic mycoplasma contamination detection."""
        sequences = [
            "ATCGATCG",  # No mycoplasma pattern
            "TTAAATTTAAATTTAAATTT",  # Contains mycoplasma pattern
        ]

        results = detect_mycoplasma_contamination(sequences)

        assert "1" in results  # Second sequence should be flagged
        assert results["1"] is True
        assert "0" not in results  # First sequence should not be flagged

    def test_detect_mycoplasma_contamination_with_genome(self):
        """Test detection with custom mycoplasma genome."""
        # Use sequences that actually match the mycoplasma genome
        sequences = ["TTAAATTTAAATTTAAATTT", "ATCGATCG"]  # First matches, second doesn't
        mycoplasma_genome = "TTAAATTTAAATTTAAATTT"  # Exact match

        results = detect_mycoplasma_contamination(sequences, mycoplasma_genome)

        assert "0" in results  # Should detect contamination in first sequence
        assert results["0"] is True
        # Second sequence might or might not be detected depending on similarity threshold


class TestAdapterContamination:
    """Test adapter contamination detection."""

    def test_detect_adapter_contamination_basic(self):
        """Test basic adapter contamination detection."""
        sequences = [
            "ATCGATCG",  # No adapter
            "AGATCGGAAGAGATCG",  # Contains adapter
        ]

        results = detect_adapter_contamination(sequences)

        assert "1" in results  # Second sequence should be flagged
        assert "AGATCGGAAGAG" in results["1"]
        assert "0" not in results  # First sequence should not be flagged


class TestVectorContamination:
    """Test vector contamination detection."""

    def test_detect_vector_contamination_basic(self):
        """Test basic vector contamination detection."""
        sequences = [
            "ATCGATCG",  # No vector
            "GGCCGCTCTAGAACTAGTGGATC",  # Contains pUC19 sequence
        ]

        results = detect_vector_contamination(sequences)

        assert "1" in results  # Second sequence should be flagged
        assert "pUC19" in results["1"][0]
        assert "0" not in results  # First sequence should not be flagged


class TestContaminationReport:
    """Test contamination report generation."""

    def test_generate_contamination_report_basic(self):
        """Test basic contamination report generation."""
        contamination_results = {
            "cross_species": {"0": ["human"], "1": ["mouse"]},
            "rrna": {"2": [5.0]},
            "mycoplasma": {"3": [True]},
        }

        report = generate_contamination_report(contamination_results)

        assert "METAINFORMANT Contamination Analysis Report" in report
        assert "CROSS_SPECIES" in report or "cross_species" in report.lower()
        assert "RRNA" in report or "rrna" in report.lower()
        assert "MYCOPLASMA" in report or "mycoplasma" in report.lower()
        assert "Summary:" in report

    def test_generate_contamination_report_empty(self):
        """Test report generation with no contamination."""
        report = generate_contamination_report({})

        assert "METAINFORMANT Contamination Analysis Report" in report
        assert "Total samples analyzed: 0" in report
