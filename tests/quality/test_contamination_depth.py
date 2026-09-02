"""Depth tests for the contamination detector.

Uses real sequences built from the detector's own known-contaminant and
reference-genome strings so k-mer matching runs on genuine data.
"""

from __future__ import annotations

import pytest

from metainformant.quality.analysis.contamination import ContaminationDetector


@pytest.fixture()
def detector() -> ContaminationDetector:
    return ContaminationDetector()


class TestMicrobialContamination:
    def test_empty_input(self, detector: ContaminationDetector) -> None:
        assert detector.detect_microbial_contamination([]) == {"detected": False, "contaminants": []}

    def test_clean_sequences_not_flagged(self, detector: ContaminationDetector) -> None:
        # 40 low-similarity sequences from a repeat motif absent from the contaminants
        seqs = [("ACGTACGTAATTCCGG" * 10) for _ in range(40)]
        result = detector.detect_microbial_contamination(seqs)
        assert result["detected"] is False
        assert result["total_sequences_analyzed"] == 40

    def test_contaminant_sequences_detected(self, detector: ContaminationDetector) -> None:
        ecoli = detector.known_contaminants["ecoli"]
        seqs = [ecoli, ecoli] + [("ACGTTTGGCCAAGGTT" * 15) for _ in range(18)]
        result = detector.detect_microbial_contamination(seqs, threshold=0.01)
        assert result["detected"] is True
        names = {c["contaminant"] for c in result["contaminants"]}
        assert "ecoli" in names
        ecoli_entry = next(c for c in result["contaminants"] if c["contaminant"] == "ecoli")
        assert ecoli_entry["matches"] >= 2
        assert 0.0 < ecoli_entry["average_similarity"] <= 1.0

    def test_rate_below_threshold_not_reported(self, detector: ContaminationDetector) -> None:
        ecoli = detector.known_contaminants["ecoli"]
        seqs = [ecoli] + [("ACGTTTGGCCAAGGTT" * 15) for _ in range(199)]
        result = detector.detect_microbial_contamination(seqs, threshold=0.5)
        assert result["detected"] is False


class TestCrossSpeciesContamination:
    REF = {
        "mouse": "ACGTACGTACGTTTTGGGCCCAAAT" * 8,
        "rat": "TTGGCATGCATCCAAGGTTCCAAAC" * 8,
        "zebrafish": "GGCCTTAAGGCCTTAACCGGATTC" * 8,
    }

    def test_missing_target_or_empty(self, detector: ContaminationDetector) -> None:
        assert detector.detect_cross_species_contamination([], "mouse", ["rat"]) == {
            "detected": False,
            "contaminants": [],
        }
        assert detector.detect_cross_species_contamination(["ACGT"], "unknown", ["rat"]) == {
            "detected": False,
            "contaminants": [],
        }

    def test_unknown_other_species_skipped(self, detector: ContaminationDetector) -> None:
        d = ContaminationDetector(reference_genomes=self.REF)
        result = d.detect_cross_species_contamination([self.REF["mouse"]], "mouse", ["nonexistent"])
        assert result["detected"] is False

    def test_target_reads_not_flagged(self, detector: ContaminationDetector) -> None:
        d = ContaminationDetector(reference_genomes=self.REF)
        reads = [self.REF["mouse"][i : i + 60] for i in range(0, 100, 20)]
        result = d.detect_cross_species_contamination(reads, "mouse", ["rat", "zebrafish"])
        assert result["detected"] is False
        assert result["target_species"] == "mouse"

    def test_contaminant_reads_flagged(self, detector: ContaminationDetector) -> None:
        d = ContaminationDetector(reference_genomes=self.REF)
        reads = [self.REF["rat"][i : i + 60] for i in range(0, 100, 20)] * 2
        result = d.detect_cross_species_contamination(reads, "mouse", ["rat"])
        assert result["detected"] is True
        assert result["contaminants"][0]["species"] == "rat"
        assert result["contaminants"][0]["matches"] == len(reads)


class TestAdapterContamination:
    def test_empty_input(self, detector: ContaminationDetector) -> None:
        assert detector.detect_adapter_contamination([]) == {"detected": False, "adapters": []}

    def test_internal_adapter_match(self, detector: ContaminationDetector) -> None:
        adapter = "AGATCGGAAGAG"
        seqs = [("ACGTACGTACGT" + adapter + "TTGGCATTGGAA") for _ in range(20)]
        result = detector.detect_adapter_contamination(seqs)
        assert result["detected"] is True
        entry = result["adapters"][0]
        assert entry["matches"] == 20
        assert entry["positions"][0] == 12

    def test_clean_reads_not_flagged(self, detector: ContaminationDetector) -> None:
        seqs = [("ACGTACGTAATTCCGGCAAT" * 8) for _ in range(50)]
        result = detector.detect_adapter_contamination(seqs)
        assert result["detected"] is False

    def test_custom_adapter_list(self, detector: ContaminationDetector) -> None:
        seqs = ["TTTTTT" + "CACAGATCGG" for _ in range(10)]
        result = detector.detect_adapter_contamination(seqs, adapters=["CACAGATCGG"])
        assert result["detected"] is True
        assert result["adapters"][0]["adapter_sequence"] == "CACAGATCGG"


class TestDuplicationContamination:
    def test_empty_input(self, detector: ContaminationDetector) -> None:
        assert detector.detect_duplication_contamination([]) == {"detected": False, "duplicated_sequences": []}

    def test_excessive_duplicates_detected(self, detector: ContaminationDetector) -> None:
        unique = [("ACGTACGTAATTCCGGCAAT" * 6) + str(i) for i in range(10)]
        dup = "GGGGAAAATTTTCCCCGGGG" * 10
        seqs = unique + [dup] * 25
        result = detector.detect_duplication_contamination(seqs, max_duplicates=10)
        assert result["detected"] is True
        assert result["duplicated_sequences"][0]["count"] == 25
        assert result["duplication_rate"] > 0.5

    def test_diverse_library_not_flagged(self, detector: ContaminationDetector) -> None:
        seqs = [("ACGTACGTAATTCCGGCAAT" * 6) + str(i) for i in range(50)]
        result = detector.detect_duplication_contamination(seqs)
        assert result["detected"] is False
        assert result["unique_sequences"] == 50

    def test_high_rate_low_count_detected(self, detector: ContaminationDetector) -> None:
        # No single sequence over max_duplicates, but duplication rate > 0.8
        seqs = (["ACGT"] * 2 + ["TGCA"]) * 20
        result = detector.detect_duplication_contamination(seqs, max_duplicates=100)
        assert result["detected"] is True


class TestComprehensiveAndSeverity:
    def test_comprehensive_analysis_runs_all_detectors(self, detector: ContaminationDetector) -> None:
        ecoli = detector.known_contaminants["ecoli"]
        seqs = [ecoli] * 5
        result = detector.comprehensive_contamination_analysis(seqs)
        assert isinstance(result, dict)
        assert result, "comprehensive analysis should return sections"

    def test_severity_classification_thresholds(self, detector: ContaminationDetector) -> None:
        # Thresholds are percentage-scale (>=50 high, >=20 moderate, >=5 low).
        assert detector._classify_severity(0.0) == "none"
        assert detector._classify_severity(4.9) == "none"
        assert detector._classify_severity(5.0) == "low"
        assert detector._classify_severity(19.9) == "low"
        assert detector._classify_severity(20.0) == "moderate"
        assert detector._classify_severity(49.9) == "moderate"
        assert detector._classify_severity(50.0) == "high"
        assert detector._classify_severity(100.0) == "high"

    def test_comprehensive_summary_severity_level(self, detector: ContaminationDetector) -> None:
        dup = "GGGGAAAATTTTCCCCGGGG" * 10
        seqs = [("ACGTACGTAATTCCGGCAAT" * 6) + str(i) for i in range(10)] + [dup] * 25
        result = detector.comprehensive_contamination_analysis(seqs)
        assert result["summary"]["contamination_detected"] is True
        assert result["summary"]["severity_level"] == "high"
        assert "duplication_contamination" in result["summary"]["analyses_performed"]
