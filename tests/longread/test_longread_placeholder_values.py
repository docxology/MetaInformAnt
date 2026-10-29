"""Regression tests for placeholder values, silent defaults and unkept promises in the longread package.

Each class pins one defect found while fixing the base-modification, adapter and read-length defects.
All of them ran without an error and returned a plausible number:

* a signal entry without a log-likelihood ratio counted as LLR 0.0, a probability of 0.5, which the default
  threshold calls methylated;
* ``aggregate_methylation`` put a call without a position at position 0 and a call without ``methylated``
  among the unmethylated;
* ``compute_methylation_stats`` reported the minimum as the first quartile and the maximum as the third for four
  sites;
* ``calculate_phase_block_stats`` turned block quality q into a switch error rate by an invented curve
  ``1 / (1 + 10 q)``;
* ``build_phase_blocks`` scored a variant without a quality as 1.0 and reported 0 supporting reads;
* ``allele_specific_analysis`` took a read without a signal as signal 0.0;
* a single read became a consensus of Q30, and ``polish_consensus`` with no reads claimed coverage 1 and Q30 for
  zero reads; ``_simple_align`` paired base i of two unrelated sequences when it found no anchor;
* the differential step of the methylation pipeline tested the first half of a list of calls against the second;
* clustered SV calls were all genotyped 0/1; ``generate_qc_report`` counted trimmed reads by their start only;
  ``process_batch`` accepted and ignored ``timeout_per_item``.

All tests run the real functions on small hand-built inputs; there are no mocks.
"""

from __future__ import annotations

import logging
import math
from array import array
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from metainformant.longread.analysis.modified_bases import call_5mc
from metainformant.longread.analysis.phasing import PhaseBlock, calculate_phase_block_stats
from metainformant.longread.analysis.structural import detect_sv_from_long_reads, phase_structural_variants
from metainformant.longread.assembly import consensus
from metainformant.longread.methylation import calling
from metainformant.longread.phasing import haplotyping
from metainformant.longread.quality.filtering import ReadRecord
from metainformant.longread.utils.batch import process_batch
from metainformant.longread.workflow.orchestrator_core import LongReadOrchestrator, PipelineResult, PipelineStep
from metainformant.longread.workflow.pipelines import get_methylation_pipeline_config
from metainformant.longread.workflow.reporting import generate_qc_report

# ---------------------------------------------------------------------------
# methylation.calling
# ---------------------------------------------------------------------------


def _signal(positions: range, llr: float, read: str = "r") -> list[dict[str, Any]]:
    return [
        {"read_id": f"{read}{i}", "chrom": "chr1", "position": p, "strand": "+", "log_likelihood_ratio": llr}
        for i, p in enumerate(positions)
    ]


class TestMissingLogLikelihoodRatio:
    @pytest.mark.parametrize("model", ["threshold", "bayesian"])
    @pytest.mark.parametrize(
        "entry", [{"read_id": "r", "chrom": "chr1", "position": 5}, {"position": 5, "log_likelihood_ratio": None}]
    )
    def test_an_entry_without_an_llr_is_an_error_not_a_methylated_call(self, model: str, entry: dict[str, Any]) -> None:
        with pytest.raises(ValueError, match=r"signal_data\[1\] has no 'log_likelihood_ratio'"):
            calling.call_methylation_from_signal(_signal(range(1), 2.0) + [entry], model=model)

    def test_the_old_default_made_a_blank_entry_a_methylated_call_at_confidence_zero(self) -> None:
        # logistic(0.0) = 0.5 >= the 0.5 threshold; this is what the missing LLR used to turn into
        (called,) = calling.call_methylation_from_signal([{"position": 5, "log_likelihood_ratio": 0.0}])
        assert called["methylated"] is True and called["confidence"] == 0.0

    def test_a_nan_llr_is_still_an_unmethylated_call_with_nan_confidence(self) -> None:
        (called,) = calling.call_methylation_from_signal([{"position": 5, "log_likelihood_ratio": float("nan")}])
        assert called["methylated"] is False and math.isnan(called["confidence"])


class TestAggregateMethylationRequiresWhatItDocuments:
    @pytest.mark.parametrize("missing", ["position", "methylated"])
    def test_a_call_without_position_or_methylated_raises(self, missing: str) -> None:
        call = {"chrom": "chr1", "position": 10, "strand": "+", "methylated": True}
        del call[missing]
        with pytest.raises(ValueError, match=f"calls\\[1\\] has no '{missing}'"):
            calling.aggregate_methylation([{"chrom": "chr1", "position": 10, "methylated": True}, call], min_coverage=1)

    def test_a_call_without_position_no_longer_lands_at_position_zero(self) -> None:
        with pytest.raises(ValueError, match="no 'position'"):
            calling.aggregate_methylation([{"chrom": "chr1", "methylated": True}], min_coverage=1)

    def test_chrom_and_strand_keep_their_documented_defaults(self) -> None:
        (site,) = calling.aggregate_methylation([{"position": 10, "methylated": True}], min_coverage=1)["sites"]
        assert (site["chrom"], site["strand"], site["position"], site["methylation_freq"]) == ("", "+", 10, 1.0)

    def test_methylated_false_is_a_call_not_a_missing_value(self) -> None:
        (site,) = calling.aggregate_methylation([{"position": 3, "methylated": False}], min_coverage=1)["sites"]
        assert (site["n_methylated"], site["n_unmethylated"]) == (0, 1)


class TestPatternAnalysisRegion:
    CALLS = [{"read_id": "a", "chrom": "chr1", "position": 5, "methylated": True}]

    @pytest.mark.parametrize("region", [{}, {"chrom": "chr1", "start": 0}, {"start": 0, "end": 10}])
    def test_a_region_missing_chrom_start_or_end_raises_as_documented(self, region: dict[str, Any]) -> None:
        with pytest.raises(ValueError, match="region is missing"):
            calling.methylation_pattern_analysis(self.CALLS, region)

    def test_a_complete_region_still_works(self) -> None:
        result = calling.methylation_pattern_analysis(self.CALLS, {"chrom": "chr1", "start": 0, "end": 10})
        assert result["n_reads"] == 1 and result["patterns"] == ["M"]


class TestDmrDocumentedEmptyResult:
    def test_samples_without_sites_return_an_empty_list_and_warn(self, caplog: pytest.LogCaptureFixture) -> None:
        with caplog.at_level(logging.WARNING, logger=calling.logger.name):
            assert calling.detect_dmrs({"sites": []}, {"sites": [{"chrom": "chr1", "position": 1}]}) == []
        assert any("no sites" in r.getMessage() for r in caplog.records)


class TestMethylationStatsQuartiles:
    @staticmethod
    def _stats(freqs: list[float]) -> dict[str, Any]:
        return calling.compute_methylation_stats(
            [{"methylation_freq": f, "coverage": 10, "chrom": "chr1", "position": i} for i, f in enumerate(freqs)]
        )["methylation_distribution"]

    def test_four_sites_have_interpolated_quartiles_not_the_extremes(self) -> None:
        distribution = self._stats([0.1, 0.2, 0.3, 0.4])
        assert distribution["q25"] == pytest.approx(0.175)
        assert distribution["q75"] == pytest.approx(0.325)
        assert distribution["median"] == pytest.approx(0.25)
        assert (distribution["min"], distribution["max"]) == (0.1, 0.4)  # the old q25 and q75

    @pytest.mark.parametrize("count", [1, 2, 3, 5, 8, 11, 40])
    def test_quartiles_match_numpy(self, count: int) -> None:
        rng = np.random.default_rng(count)
        freqs = sorted(float(x) for x in rng.random(count))
        distribution = self._stats(freqs)
        assert distribution["q25"] == pytest.approx(np.percentile(freqs, 25))
        assert distribution["q75"] == pytest.approx(np.percentile(freqs, 75))
        assert distribution["median"] == pytest.approx(np.median(freqs))


# ---------------------------------------------------------------------------
# phasing
# ---------------------------------------------------------------------------


class TestPhaseBlockStatsSwitchErrorRate:
    def test_it_is_the_variant_weighted_link_disagreement_not_an_invented_curve(self) -> None:
        blocks = [
            PhaseBlock(block_id=0, chromosome="chr1", start=0, end=10000, num_variants=50, quality=0.95),
            PhaseBlock(block_id=1, chromosome="chr1", start=20000, end=25000, num_variants=20, quality=0.90),
        ]
        stats = calculate_phase_block_stats(blocks)
        assert stats.switch_error_rate == pytest.approx(1 - (0.95 * 50 + 0.90 * 20) / 70)  # 0.0643
        assert stats.switch_error_rate != pytest.approx(1 / (1 + 10 * 0.9357), abs=1e-3)  # the old 0.0965

    def test_blocks_whose_links_all_agree_have_rate_zero(self) -> None:
        block = PhaseBlock(block_id=0, chromosome="chr1", start=0, end=500, num_variants=4, quality=1.0)
        assert calculate_phase_block_stats([block]).switch_error_rate == 0.0  # the old curve said 0.0909


class TestBuildPhaseBlocksPlaceholders:
    VARIANTS = [{"chrom": "chr1", "position": p, "haplotype": 1} for p in (100, 200, 300)]

    def test_no_quality_is_unknown_not_a_perfect_one(self) -> None:
        (block,) = haplotyping.build_phase_blocks(self.VARIANTS)
        assert block["phase_quality"] is None  # a variant without a quality used to count as 1.0

    def test_the_mean_is_over_the_variants_that_have_a_quality(self) -> None:
        variants = [dict(v) for v in self.VARIANTS]
        variants[0]["quality"] = 20.0
        variants[2]["quality"] = 40.0
        (block,) = haplotyping.build_phase_blocks(variants)
        assert block["phase_quality"] == pytest.approx(30.0)  # not (20 + 1 + 40) / 3

    def test_the_read_count_is_unknown_not_zero(self) -> None:
        (block,) = haplotyping.build_phase_blocks(self.VARIANTS)
        assert block["n_reads"] is None


class TestAlleleSpecificAnalysisSignal:
    @staticmethod
    def _read(read_id: str, signal: Any) -> dict[str, Any]:
        read = {"read_id": read_id, "chrom": "chr1", "start": 90, "end": 230}
        if signal is not None:
            read["signal"] = signal
        return read

    def test_a_read_without_signal_is_left_out_not_averaged_in_as_zero(self, caplog: pytest.LogCaptureFixture) -> None:
        reads = [
            self._read("a", 0.8),
            self._read("b", 0.8),
            self._read("c", None),
            self._read("d", 0.2),
            self._read("e", 0.2),
        ]
        assignments = {"a": 1, "b": 1, "c": 1, "d": 2, "e": 2}
        with caplog.at_level(logging.WARNING, logger=haplotyping.logger.name):
            result = haplotyping.allele_specific_analysis(
                reads, assignments, {"chrom": "chr1", "start": 90, "end": 230}
            )
        assert result["haplotype_1_signal"] == pytest.approx(0.8)  # (0.8 + 0.8 + 0.0) / 3 = 0.533 with the placeholder
        assert result["haplotype_2_signal"] == pytest.approx(0.2)
        assert any("1 reads in the region have no signal" in r.getMessage() for r in caplog.records)


# ---------------------------------------------------------------------------
# assembly.consensus
# ---------------------------------------------------------------------------


class TestConsensusPlaceholders:
    def test_a_single_read_has_no_made_up_quality(self) -> None:
        result = consensus.generate_consensus(["ACGTACGTAC"])
        assert result.quality == [0.0] * 10 and result.mean_quality == 0.0
        assert result.coverage == [1] * 10 and result.num_reads == 1

    def test_polishing_with_no_reads_reports_no_support(self) -> None:
        result = consensus.polish_consensus("ACGTACGTAC", [])
        assert result.sequence == "ACGTACGTAC"
        assert (result.num_reads, result.mean_coverage, result.mean_quality) == (
            0,
            0.0,
            0.0,
        )  # was 0 reads, coverage 1, Q30
        assert result.coverage == [0] * 10 and result.quality == [0.0] * 10

    def test_unrelated_sequences_without_a_shared_kmer_give_an_empty_alignment(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        with caplog.at_level(logging.WARNING, logger=consensus.logger.name):
            alignment = consensus._simple_align("A" * 40, "C" * 40)
        # side by side they would pair base i with base i: forty mismatches voted into the consensus
        assert (alignment["ref_aligned"], alignment["query_aligned"], alignment["score"]) == ("", "", 0)
        assert any("No shared 11-mer" in r.getMessage() for r in caplog.records)

    def test_the_msa_warns_that_its_rows_are_not_column_aligned(self, caplog: pytest.LogCaptureFixture) -> None:
        sequences = ["ACGTACGTACGTTGCA", "ACGTACGTACGTTGCA", "ACGTACGTAAACGTTGCA"]  # one carries a 2 bp insertion
        with caplog.at_level(logging.WARNING, logger=consensus.logger.name):
            result = consensus.multiple_sequence_alignment(sequences)
        assert result.num_sequences == 3
        assert any("not column-consistent" in r.getMessage() for r in caplog.records)
        assert "UPGMA" not in (consensus.multiple_sequence_alignment.__doc__ or "")


# ---------------------------------------------------------------------------
# structural variants, reports, batch
# ---------------------------------------------------------------------------


def _deletion_alignment(name: str, tags: dict[str, Any] | None = None) -> dict[str, Any]:
    return {
        "read_name": name,
        "reference_name": "chr1",
        "reference_start": 1000,
        "mapping_quality": 60,
        "cigar_string": "500M200D500M",
        "tags": tags or {},
        "query_sequence": "",
    }


class TestSvGenotype:
    def test_clustered_calls_are_not_genotyped_from_read_support_alone(self) -> None:
        calls = detect_sv_from_long_reads([_deletion_alignment(f"r{i}") for i in range(3)], min_size=50, min_support=2)
        (deletion,) = calls
        assert deletion.supporting_reads == 3
        assert deletion.genotype == "./."  # every supported call used to be 0/1

    def test_phasing_still_assigns_a_genotype(self) -> None:
        (deletion,) = detect_sv_from_long_reads(
            [_deletion_alignment(f"r{i}") for i in range(3)], min_size=50, min_support=2
        )
        (phased,) = phase_structural_variants([deletion], {"r0": 1, "r1": 1, "r2": 1})
        assert (phased.haplotype, phased.genotype) == (1, "0/1")


class TestQcReportAdapterTrimmed:
    @staticmethod
    def _result(records: list[ReadRecord]) -> PipelineResult:
        step = PipelineStep(name="trim_adapters", function=lambda ctx: None, status="completed", result=records)
        return PipelineResult(pipeline_name="qc", steps=[step], success=True, total_duration=0.0, output_dir=Path("."))

    def test_a_read_trimmed_only_at_the_end_counts(self) -> None:
        records = [
            ReadRecord("start", "ACGT", metadata={"trimmed_start": 28, "trimmed_end": 0}),
            ReadRecord("end", "ACGT", metadata={"trimmed_start": 0, "trimmed_end": 22}),
            ReadRecord("both", "ACGT", metadata={"trimmed_start": 28, "trimmed_end": 22}),
            ReadRecord("none", "ACGT", metadata={"trimmed_start": 0, "trimmed_end": 0}),
        ]
        assert generate_qc_report(self._result(records)).adapter_trimmed == 3  # 2 before: end-only was missed

    def test_dict_records_count_too(self) -> None:
        result = self._result([{"metadata": {"trimmed_end": 5}}, {"metadata": {}}])  # type: ignore[list-item]
        assert generate_qc_report(result).adapter_trimmed == 1


class TestBatchTimeout:
    def test_a_timeout_that_cannot_be_enforced_is_refused(self) -> None:
        with pytest.raises(ValueError, match="timeout_per_item is not supported"):
            process_batch({"a": 1}, lambda item_id, item: item, timeout_per_item=5.0)

    def test_no_timeout_still_runs(self) -> None:
        result = process_batch({"a": 1, "b": 2}, lambda item_id, item: item * 2, timeout_per_item=None)
        assert result.results == {"a": 2, "b": 4}


# ---------------------------------------------------------------------------
# the methylation pipeline's differential step
# ---------------------------------------------------------------------------

SEQUENCE = "ACGTACGTACGTACGT"  # C at 1, 5, 9, 13


def _tagged_read(read_id: str, sample: str | None, probability: int) -> dict[str, Any]:
    read: dict[str, Any] = {
        "read_id": read_id,
        "sequence": SEQUENCE,
        "reference_name": "chr1",
        "reference_start": 1000,
        "cigar_string": "16M",
        "is_reverse": False,
        "tags": {"MM": "C+m,0;", "ML": array("B", [probability])},
    }
    if sample is not None:
        read["sample"] = sample
    return read


def _run_methylation(tmp_path: Path, reads: list[dict[str, Any]]) -> dict[str, Any]:
    config = get_methylation_pipeline_config(modification_types=["5mC"], min_coverage=5)
    result = LongReadOrchestrator(config, tmp_path).run_methylation_pipeline(reads)
    step = next(s for s in result.steps if s.name == "differential_analysis")
    assert step.status == "completed", step.error
    return step.result  # type: ignore[no-any-return]


class TestDifferentialStepOfTheMethylationPipeline:
    def test_two_labelled_samples_are_compared_with_each_other(self, tmp_path: Path) -> None:
        reads = [_tagged_read(f"a{i}", "tumour", 250) for i in range(6)] + [
            _tagged_read(f"b{i}", "normal", 5) for i in range(6)
        ]
        outcome = _run_methylation(tmp_path, reads)
        assert outcome["status"] == "completed"
        assert outcome["samples"] == ["normal", "tumour"]
        assert (outcome["num_tested"], outcome["num_significant"]) == (1, 1)
        (site,) = outcome["results"]
        assert (site.chromosome, site.position) == ("chr1", 1001)
        assert site.sample1_methylation < 0.1 < 0.9 < site.sample2_methylation  # normal, then tumour

    def test_one_sample_is_not_split_in_two_and_tested_against_itself(self, tmp_path: Path) -> None:
        # twelve calls of a single sample: the old half split tested six of them against the other six
        reads = [_tagged_read(f"a{i}", "tumour", 250 if i < 6 else 5) for i in range(12)]
        outcome = _run_methylation(tmp_path, reads)
        assert outcome["status"] == "skipped"
        assert outcome["num_tested"] == 0 and outcome["results"] == []
        assert "reads carry 1 sample label" in outcome["reason"]

    def test_unlabelled_reads_are_not_compared_either(self, tmp_path: Path) -> None:
        reads = [_tagged_read(f"a{i}", None, 250 if i < 6 else 5) for i in range(12)]
        outcome = _run_methylation(tmp_path, reads)
        assert outcome["status"] == "skipped" and "0 sample label" in outcome["reason"]

    def test_three_samples_are_refused_rather_than_guessed(self, tmp_path: Path) -> None:
        reads = [_tagged_read(f"{label}{i}", label, 100) for label in ("x", "y", "z") for i in range(5)]
        assert _run_methylation(tmp_path, reads)["status"] == "skipped"

    def test_calls_carry_the_read_they_came_from(self) -> None:
        (call,) = call_5mc([{"position": 7, "probability": 0.5, "read_id": "r9"}])
        assert call.read_id == "r9"
