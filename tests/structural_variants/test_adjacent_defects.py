"""Regression tests for defects found next to the translocation, PCA, nearest-gene, genotyping and microhomology fixes.

* ``genotype_sv_population`` invented evidence for a sample that had none (a depth ratio of 1.0, 30 total reads),
  so a sample with no data was called reference with quality 50 (depth) or 60 (split): a confident call from
  nothing. It is still called reference (the documented default) but with quality 0.0, and the number of such calls
  is logged. A depth genotype for an SV without ``sv_type`` silently used the deletion thresholds.
* ``sv_ld_analysis`` documented r^2 "with all SNPs within a window" (there is no window), reported position 0 for a
  tag SNP when ``snp_positions`` was too short, and ignored ``sv_positions`` altogether.
* ``calculate_breakpoint_confidence`` put an invented mapping quality of 30 (and an invented position spread of 10)
  into the score whenever the evidence lacked one, and ``refine_breakpoints`` never supplied the MAPQ its reads carry.
* ``plot_cnv_profile`` drew the "DEL" reference line at log2(1/ploidy), which is the ratio of one copy left only for
  ploidy 2, and substituted constants for a ploidy below 1.
"""

from __future__ import annotations

import logging
import math
from typing import Any

import pytest

try:
    import matplotlib

    matplotlib.use("Agg")
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from metainformant.structural_variants.detection.breakpoints import (
    calculate_breakpoint_confidence,
    refine_breakpoints,
)
from metainformant.structural_variants.detection.sv_calling import SVEvidence, _calculate_sv_quality
from metainformant.structural_variants.population import genotype_sv_population, sv_allele_frequency, sv_ld_analysis

POPULATION_LOGGER = "metainformant.structural_variants.population.sv_population"


def _sv(sv_type: str = "DEL", **samples: dict[str, float]) -> dict[str, Any]:
    return {"chrom": "chr1", "start": 1_000, "end": 2_000, "sv_type": sv_type, "samples": dict(samples)}


class TestGenotypingWithoutEvidence:
    def test_a_sample_with_no_depth_evidence_is_reference_with_zero_quality(self) -> None:
        sv = _sv(s1={"depth_ratio": 0.1})
        result = genotype_sv_population([sv], ["s1", "s2"], method="depth")
        assert result["genotype_matrix"] == [[2, 0]]
        # s1: a clear deletion, quality (0.25 - 0.1) * 400 = 60. s2: nothing observed. (The old quality was 50.)
        assert result["quality_scores"] == [[pytest.approx(60.0), 0.0]]

    def test_a_sample_with_no_split_evidence_is_reference_with_zero_quality(self) -> None:
        sv = _sv(s1={"split_reads": 12, "discordant_pairs": 8, "total_reads": 30})
        result = genotype_sv_population([sv], ["s1", "s2"], method="split")
        # s1: 20 of 30 reads support the SV (0.67 -> heterozygous), quality min(99, 20 * 5) = 99. (The old quality
        # of s2 was (30 - 0) * 2 = 60.)
        assert result["genotype_matrix"] == [[1, 0]]
        assert result["quality_scores"] == [[pytest.approx(99.0), 0.0]]

    def test_a_sample_whose_evidence_lacks_the_measurement_is_treated_the_same(self) -> None:
        sv = _sv(s1={"split_reads": 4}, s2={"depth_ratio": 0.2})  # no total_reads for s1; no depth_ratio for s1 either
        assert genotype_sv_population([sv], ["s1"], method="split")["quality_scores"] == [[0.0]]
        assert genotype_sv_population([sv], ["s1"], method="depth")["quality_scores"] == [[0.0]]
        # An explicit zero depth ratio is a measurement (no coverage), not a missing one.
        assert genotype_sv_population([sv], ["s2"], method="depth")["genotype_matrix"] == [[2]]

    def test_an_explicit_total_read_count_of_zero_is_a_measurement(self) -> None:
        sv = _sv(s1={"split_reads": 0, "discordant_pairs": 0, "total_reads": 0})
        result = genotype_sv_population([sv], ["s1"], method="split")
        assert result["genotype_matrix"] == [[0]] and result["quality_scores"] == [[2.0]]  # (1 - 0) * 2 by the rule

    def test_the_number_of_calls_without_evidence_is_logged(self, caplog: pytest.LogCaptureFixture) -> None:
        sv = _sv(s1={"depth_ratio": 0.1})
        with caplog.at_level(logging.WARNING, logger=POPULATION_LOGGER):
            genotype_sv_population([sv, _sv()], ["s1", "s2", "s3"], method="depth")
        messages = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
        # s2 and s3 of the first SV, and all three samples of the second.
        assert any(m.startswith("5 of 6 sample genotypes had no depth_ratio evidence") for m in messages)

    def test_complete_evidence_logs_nothing(self, caplog: pytest.LogCaptureFixture) -> None:
        with caplog.at_level(logging.WARNING, logger=POPULATION_LOGGER):
            genotype_sv_population([_sv(s1={"depth_ratio": 0.5})], ["s1"], method="depth")
        assert not [r for r in caplog.records if r.levelno == logging.WARNING]

    def test_the_depth_method_needs_the_sv_type(self) -> None:
        sv = _sv(s1={"depth_ratio": 0.1})
        del sv["sv_type"]
        with pytest.raises(ValueError, match="sv_type"):
            genotype_sv_population([sv], ["s1"], method="depth")
        # The split method does not use the type.
        assert genotype_sv_population([sv], ["s1"], method="split")["genotype_matrix"] == [[0]]

    def test_an_enum_sv_type_selects_the_thresholds(self) -> None:
        from metainformant.structural_variants.detection.sv_calling import SVType

        sv = _sv(SVType.DUP, s1={"depth_ratio": 1.9})  # type: ignore[arg-type]
        assert genotype_sv_population([sv], ["s1"], method="depth")["genotype_matrix"] == [[2]]


class TestAlleleFrequencyInputs:
    def test_valid_genotypes_give_the_frequencies_by_hand(self) -> None:
        # [0, 1, 2, 1]: 4 alternate alleles of 8 -> 0.5; [0, 0, 0, 2]: 2 of 8 -> 0.25 (whole-number floats are fine).
        result = sv_allele_frequency([[0, 1, 2, 1], [0.0, 0.0, 0.0, 2.0]], sample_labels=["a", "a", "b", "b"])
        assert result["frequencies"] == [0.5, 0.25]
        assert result["population_frequencies"] == {"a": [0.25, 0.0], "b": [0.75, 0.5]}

    @pytest.mark.parametrize("bad", [-1, 3, 0.5, float("nan")])
    def test_a_value_that_is_not_a_genotype_is_refused(self, bad: float) -> None:
        # -1 (a missing call) used to be summed as an allele count and binned into the last SFS bin; 0.5 was cut to 0.
        with pytest.raises(ValueError, match="0, 1 and 2"):
            sv_allele_frequency([[0, 1, bad, 2]])

    @pytest.mark.parametrize("labels", [["a", "b"], ["a", "a", "b", "b", "b"]])
    def test_labels_must_cover_every_sample(self, labels: list[str]) -> None:
        with pytest.raises(ValueError, match="sample_labels has"):
            sv_allele_frequency([[0, 1, 2, 1]], sample_labels=labels)


class TestLdAnalysisInputs:
    SV = [[0, 1, 2, 0, 1, 2, 0, 1]]
    SNPS = [[2, 2, 2, 2, 2, 2, 2, 2], [0, 1, 2, 0, 1, 2, 0, 1]]

    def test_a_snp_a_million_bases_away_is_still_compared(self) -> None:
        # No window: r^2 is computed against every SNP passed in, wherever it is.
        result = sv_ld_analysis(self.SV, self.SNPS, [10_000], [5_000, 50_000_000])
        assert result["tag_snps"][0]["snp_index"] == 1
        assert result["tag_snps"][0]["snp_position"] == 50_000_000
        assert result["ld_matrix"][0][1] == pytest.approx(1.0)

    def test_snp_positions_must_cover_every_snp(self) -> None:
        # The tag SNP (index 1) used to be reported at position 0 when the list was short.
        with pytest.raises(ValueError, match="snp_positions has 1 entries for 2 SNPs"):
            sv_ld_analysis(self.SV, self.SNPS, [10_000], [5_000])

    def test_sv_positions_must_cover_every_sv(self) -> None:
        with pytest.raises(ValueError, match="sv_positions has 2 entries for 1 SVs"):
            sv_ld_analysis(self.SV, self.SNPS, [10_000, 20_000], [5_000, 12_000])

    def test_no_snps_means_no_tag_snp(self) -> None:
        result = sv_ld_analysis(self.SV, [], [10_000], [])  # type: ignore[arg-type]
        assert result["tag_snps"] == [{"snp_index": None, "snp_position": None, "r2": 0.0}]
        assert result["ld_matrix"] == [[]]


class TestBreakpointConfidence:
    def test_support_alone_is_scored_on_support_alone(self) -> None:
        # support 5: 1 - 1 / (1 + 5 / 5) = 0.5. The old score added invented position, fraction and MAPQ terms: 0.55.
        assert calculate_breakpoint_confidence({"support": 5}) == pytest.approx(0.5)

    def test_every_metric_given_keeps_the_documented_weights(self) -> None:
        # 0.35 * 0.5 + 0.30 * 1.0 + 0.20 * 1.0 + 0.15 * 1.0 (position std 0, all clips support, MAPQ 60)
        evidence = {"support": 5, "total_clips": 5, "position_std": 0.0, "mapq_mean": 60.0}
        assert calculate_breakpoint_confidence(evidence) == pytest.approx(0.825)
        # 0.35 * 0.5 + 0.30 * 0.5 + 0.20 * 0.5 + 0.15 * 0.5 (std 5 -> 1 / (1 + 1), 5 of 10 clips, MAPQ 30)
        half = {"support": 5, "total_clips": 10, "position_std": 5.0, "mapq_mean": 30.0}
        assert calculate_breakpoint_confidence(half) == pytest.approx(0.5)

    def test_a_missing_metric_is_left_out_and_the_weights_rescaled(self) -> None:
        # Without MAPQ: (0.35 * 0.5 + 0.30 * 1.0 + 0.20 * 1.0) / 0.85. The old score filled in a MAPQ of 30 (0.5).
        evidence = {"support": 5, "total_clips": 5, "position_std": 0.0}
        assert calculate_breakpoint_confidence(evidence) == pytest.approx(0.675 / 0.85)

    def test_none_counts_as_missing(self) -> None:
        with_none = {"support": 5, "total_clips": None, "position_std": None, "mapq_mean": None}
        assert calculate_breakpoint_confidence(with_none) == pytest.approx(0.5)

    def test_no_support_is_no_confidence(self) -> None:
        assert calculate_breakpoint_confidence({"support": 0, "mapq_mean": 60.0}) == 0.0

    def test_the_score_stays_in_the_unit_interval(self) -> None:
        extreme = {"support": 10_000, "total_clips": 1, "position_std": -3.0, "mapq_mean": 255.0}
        assert 0.0 <= calculate_breakpoint_confidence(extreme) <= 1.0


class TestRefinedConfidenceUsesTheReadsMapq:
    VARIANT = {"chrom": "chr1", "start": 5000, "end": 7000, "sv_type": "DEL"}

    def _reads(self, mapq: int | None) -> list[dict[str, Any]]:
        reads = [{"chrom": "chr1", "pos": 4980, "cigar": "20M20S", "name": f"r{i}"} for i in range(3)]
        if mapq is not None:
            for read in reads:
                read["mapq"] = mapq
        return reads

    def test_confidence_follows_the_mapping_quality_of_the_supporting_reads(self) -> None:
        # 3 reads clip at 5000: support 3 -> 0.375, std 0 -> 1, all 3 clips -> 1.
        #   MAPQ 60: 0.35 * 0.375 + 0.30 + 0.20 + 0.15 * 1.0 = 0.78125
        #   MAPQ 30: 0.35 * 0.375 + 0.30 + 0.20 + 0.15 * 0.5 = 0.70625
        high = refine_breakpoints([self.VARIANT], self._reads(60))[0].bp1.confidence
        low = refine_breakpoints([self.VARIANT], self._reads(30))[0].bp1.confidence
        assert high == pytest.approx(0.78125)
        assert low == pytest.approx(0.70625)

    def test_reads_without_a_mapq_leave_the_term_out(self) -> None:
        confidence = refine_breakpoints([self.VARIANT], self._reads(None))[0].bp1.confidence
        assert confidence == pytest.approx((0.35 * 0.375 + 0.30 + 0.20) / 0.85)

    def test_a_mapq_of_zero_is_a_measurement(self) -> None:
        confidence = refine_breakpoints([self.VARIANT], self._reads(0))[0].bp1.confidence
        assert confidence == pytest.approx(0.35 * 0.375 + 0.30 + 0.20)

    def test_the_unsupported_breakpoint_still_has_no_confidence(self) -> None:
        assert refine_breakpoints([self.VARIANT], self._reads(60))[0].bp2.confidence == 0.0


class TestSvQuality:
    def test_the_call_quality_is_a_bounded_heuristic_score_not_a_phred_value(self) -> None:
        # 60 from the support count (capped), 10 for both kinds of evidence, 20 for unique reads (capped): 90.
        evidence = SVEvidence(split_reads=5_000, discordant_pairs=5_000, evidence_reads=[f"r{i}" for i in range(500)])
        assert _calculate_sv_quality(evidence, 10_000) == 90.0
        assert _calculate_sv_quality(evidence, 0) == 0.0


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib required")
class TestCnvProfileReferenceLines:
    @staticmethod
    def _lines(ploidy: int) -> dict[str, float]:
        from metainformant.structural_variants.visualization.plots import plot_cnv_profile

        figure = plot_cnv_profile([], {"chr1": 1_000}, ploidy=ploidy)
        return {
            line.get_label(): float(line.get_ydata()[0])
            for line in figure.axes[0].lines
            if line.get_label() in ("DEL", "DUP")
        }

    @pytest.mark.parametrize(
        ("ploidy", "deletion", "duplication"),
        [
            (2, math.log2(1 / 2), math.log2(3 / 2)),  # -1.0, 0.585
            (3, math.log2(2 / 3), math.log2(4 / 3)),  # -0.585, 0.415
            (4, math.log2(3 / 4), math.log2(5 / 4)),  # -0.415, 0.322
        ],
    )
    def test_the_lines_mark_one_copy_lost_and_one_gained(
        self, ploidy: int, deletion: float, duplication: float
    ) -> None:
        lines = self._lines(ploidy)
        assert lines["DEL"] == pytest.approx(deletion, abs=1e-9)
        assert lines["DUP"] == pytest.approx(duplication, abs=1e-9)

    def test_a_haploid_profile_has_no_deletion_line(self) -> None:
        # Losing the only copy is a ratio of 0 (log2 = -infinity); the old line sat at log2(1 + 1e-10), on the axis.
        assert self._lines(1) == {"DUP": pytest.approx(1.0)}

    @pytest.mark.parametrize("ploidy", [0, -2])
    def test_a_ploidy_below_one_is_refused(self, ploidy: int) -> None:
        from metainformant.structural_variants.visualization.plots import plot_cnv_profile

        with pytest.raises(ValueError, match="ploidy"):
            plot_cnv_profile([], {"chr1": 1_000}, ploidy=ploidy)
