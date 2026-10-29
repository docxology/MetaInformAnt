"""Round-4 regression tests for ``metainformant.quality.analysis``.

- ``detect_outliers(method="modified_zscore")`` returned no outliers whenever the median absolute deviation was 0,
  which happens as soon as more than half the values equal the median, not only for constant data
  (``[5, 5, 5, 5, 5, 5, 5, 100]`` has an obvious outlier).
- ``calculate_quality_score`` on FASTQ data without any recognised metric returned ``overall_score 0`` and grade F
  without an ``error`` (the VCF and BAM scores say what is missing).
- ``compare_quality_metrics`` reported a 0% change from a baseline of 0.
- ``detect_cross_species_contamination`` accepted ``target_species`` and ``other_species`` and ignored both, so
  reads from the expected species were reported as contamination by it.
- an empty read counted as adapter, vector and mycoplasma contamination (the empty string is a prefix, suffix and
  substring of every pattern).
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.quality.analysis.contamination import (
    ContaminationDetector,
    detect_cross_species_contamination,
    detect_mycoplasma_contamination,
    detect_vector_contamination,
)
from metainformant.quality.analysis.metrics import calculate_quality_score, compare_quality_metrics, detect_outliers

# ---------------------------------------------------------------------------
# detect_outliers
# ---------------------------------------------------------------------------


class TestModifiedZScoreWithZeroMad:
    def test_an_outlier_among_mostly_identical_values_is_found(self) -> None:
        data = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 100.0]
        result = detect_outliers(data, method="modified_zscore")
        assert result["outlier_indices"] == [7]
        assert result["outliers"] == [100.0]

    def test_scale_is_the_iglewicz_hoaglin_fallback(self) -> None:
        data = np.array([10.0, 10.0, 10.0, 10.0, 10.0, 12.0, 30.0])
        median = np.median(data)
        mean_ad = np.mean(np.abs(data - median))
        z = (data - median) / (1.253314 * mean_ad)
        expected = np.where(np.abs(z) > 3.5)[0].tolist()
        assert np.median(np.abs(data - median)) == 0  # the premise: MAD = 0
        assert detect_outliers(data.tolist(), method="modified_zscore", threshold=3.5)["outlier_indices"] == expected
        assert expected == [6]

    def test_constant_data_still_has_no_outliers(self) -> None:
        assert detect_outliers([5.0, 5.0, 5.0, 5.0], method="modified_zscore")["outliers"] == []

    def test_a_positive_mad_gives_the_classic_statistic(self) -> None:
        data = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 40.0])
        median = np.median(data)
        mad = np.median(np.abs(data - median))
        classic = np.where(np.abs(0.6745 * (data - median) / mad) > 3.5)[0].tolist()
        assert detect_outliers(data.tolist(), method="modified_zscore", threshold=3.5)["outlier_indices"] == classic


# ---------------------------------------------------------------------------
# scores and comparisons
# ---------------------------------------------------------------------------


class TestFastqScoreWithoutMetrics:
    def test_empty_metrics_say_what_is_missing(self) -> None:
        result = calculate_quality_score({}, "fastq")
        assert "error" in result and "FASTQ" in result["error"]
        assert result["overall_score"] == 0.0 and result["components"] == {}

    def test_real_metrics_are_still_scored_without_an_error(self) -> None:
        result = calculate_quality_score({"basic_statistics": {"mean_quality": 35.0}}, "fastq")
        assert "error" not in result
        assert result["overall_score"] > 0.0


class TestCompareQualityMetrics:
    def test_a_change_from_zero_has_no_percentage(self) -> None:
        first = {"basic_statistics": {"mean_quality": 30.0, "mean_length": 100.0, "mean_gc": 0.0}}
        second = {"basic_statistics": {"mean_quality": 33.0, "mean_length": 100.0, "mean_gc": 12.0}}
        differences = compare_quality_metrics(first, second)["differences"]
        assert differences["mean_gc"]["difference"] == 12.0
        assert differences["mean_gc"]["percent_change"] is None
        assert differences["mean_quality"]["percent_change"] == pytest.approx(10.0)
        assert differences["mean_length"]["percent_change"] == 0.0


# ---------------------------------------------------------------------------
# contamination
# ---------------------------------------------------------------------------

HUMAN = "ACGTTGCAAGCTAGCTAGGATCCATGCATGCTAGCTAGCATCGATCGTAGC"
MOUSE = "TTAACCGGTTAACCGGAATTCCGGAATTGGCCAATTGGCCTTAAGGCCAAT"


class TestCrossSpeciesTargetAndCandidates:
    GENOMES = {"human": HUMAN, "mouse": MOUSE}

    def test_the_target_species_is_not_reported_as_contamination(self) -> None:
        reads = [HUMAN[5:35], MOUSE[3:33]]
        everything = detect_cross_species_contamination(reads, self.GENOMES, threshold=0.9)
        assert everything == {"0": "human", "1": "mouse"}  # without a target every species is a candidate

        result = detect_cross_species_contamination(reads, self.GENOMES, threshold=0.9, target_species="human")

        assert result == {"1": "mouse"}  # was {"0": "human", "1": "mouse"}: the expected species flagged itself

    def test_other_species_restricts_the_candidates(self) -> None:
        genomes = {**self.GENOMES, "yeast": "GGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCC"}
        reads = [MOUSE[3:33]]
        assert detect_cross_species_contamination(reads, genomes, threshold=0.9, other_species=["yeast"]) == {}
        assert detect_cross_species_contamination(reads, genomes, threshold=0.9, other_species=["mouse"]) == {
            "0": "mouse"
        }

    def test_a_read_that_matches_the_target_at_least_as_well_is_not_contamination(self) -> None:
        # the read lies in both genomes (a shared stretch); it is as close to the expected species as to the other
        shared = "GATTACAGATTACAGATTACA"
        genomes = {"human": HUMAN + shared, "mouse": MOUSE + shared}
        assert detect_cross_species_contamination([shared], genomes, threshold=0.9) != {}
        assert detect_cross_species_contamination([shared], genomes, threshold=0.9, target_species="human") == {}

    def test_unknown_candidate_names_are_an_error(self) -> None:
        with pytest.raises(ValueError, match="not in reference_genomes"):
            detect_cross_species_contamination([HUMAN[:20]], self.GENOMES, other_species=["mous"])

    def test_default_behaviour_without_a_target_is_unchanged(self) -> None:
        genomes = {"human": "ATCGATCGATCG", "mouse": "GCTAGCTAGCTA"}
        assert detect_cross_species_contamination(["ATCGATCG", "GCTAGCTA"], genomes, threshold=0.8) == {
            "0": "human",
            "1": "mouse",
        }


class TestEmptyReadsAreNotContamination:
    def test_adapter_detector(self) -> None:
        detector = ContaminationDetector()
        result = detector.detect_adapter_contamination(["", "", "ACGTACGTACGTACGT"])
        assert result["detected"] is False
        assert result["adapters"] == []

    def test_vector_detector(self) -> None:
        assert detect_vector_contamination(["", "GGCCGCTCTAGAACTAGTGGATC"]) == {"1": ["pUC19"]}

    def test_mycoplasma_detector(self) -> None:
        assert detect_mycoplasma_contamination(["", "TTAAATTTAAATTTAAATTT"], "TTAAATTTAAATTTAAATTT") == {"1": True}
        assert detect_mycoplasma_contamination([""]) == {}

    def test_real_reads_are_still_detected(self) -> None:
        detector = ContaminationDetector()
        adapter = "AGATCGGAAGAG"
        result = detector.detect_adapter_contamination([adapter + "TTTT", "GGGG" + adapter, "ACGTACGTACGTACGT"])
        assert result["adapters"][0]["matches"] == 2
