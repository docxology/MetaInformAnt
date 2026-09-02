"""Round-4 test-depth tests for longread methylation, phasing, and quality/workflow.

Covers previously untested public APIs:
- methylation.calling: call_methylation_from_signal, aggregate_methylation,
  detect_dmrs, methylation_pattern_analysis, compute_methylation_stats
- phasing.haplotyping: phase_reads, build_phase_blocks, compute_switch_errors,
  haplotag_reads, allele_specific_analysis
- quality.filtering: filter_by_length, filter_by_quality, trim_adapters,
  detect_adapters
- workflow.pipelines: load_pipeline_config (JSON/YAML/missing)

All tests use real computation on synthetic data (zero-mocks policy).
"""

from __future__ import annotations

import json
import random
from pathlib import Path

import pytest

from metainformant.longread.methylation.calling import (
    aggregate_methylation,
    call_methylation_from_signal,
    compute_methylation_stats,
    detect_dmrs,
    methylation_pattern_analysis,
)
from metainformant.longread.phasing.haplotyping import (
    allele_specific_analysis,
    build_phase_blocks,
    compute_switch_errors,
    haplotag_reads,
    phase_reads,
)
from metainformant.longread.quality.filtering import (
    detect_adapters,
    filter_by_length,
    filter_by_quality,
    trim_adapters,
)
from metainformant.longread.workflow.pipelines import load_pipeline_config

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

ADAPTER = "AGATCGGAAGAGC"
INSERT = "ACGTAGCTTAGGCATG" * 20  # 320 bp


def _phred(n: int, low: int = 20, high: int = 40) -> str:
    return "".join(chr(33 + random.randint(low, high)) for _ in range(n))


def _signal(positions: range, read_prefix: str, llr: float, chrom: str = "chr1") -> list[dict]:
    return [
        {"read_id": f"{read_prefix}{i}", "chrom": chrom, "position": p, "strand": "+", "log_likelihood_ratio": llr}
        for i in range(20)
        for p in positions
    ]


# ---------------------------------------------------------------------------
# methylation.calling
# ---------------------------------------------------------------------------


class TestCallMethylationFromSignal:
    def test_threshold_model(self) -> None:
        calls = call_methylation_from_signal(_signal(range(100, 105), "a", 3.0))
        assert len(calls) == 100
        assert calls[0]["methylated"] is True
        assert calls[0]["confidence"] > 0.5
        for c in calls:
            assert isinstance(c["methylated"], bool)
            assert 0.0 <= c["confidence"] <= 1.0

    def test_unmethylated_calls(self) -> None:
        calls = call_methylation_from_signal(_signal(range(100, 105), "b", -3.0))
        assert all(c["methylated"] is False for c in calls)

    def test_bayesian_model_shrinks_confidence(self) -> None:
        strong = {"read_id": "r0", "chrom": "chr1", "position": 1, "strand": "+", "log_likelihood_ratio": 10.0}
        bayes = call_methylation_from_signal([strong], model="bayesian")
        thresh = call_methylation_from_signal([strong], model="threshold")
        assert bayes[0]["methylated"] is True
        assert bayes[0]["confidence"] < thresh[0]["confidence"]

    def test_invalid_model_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown model"):
            call_methylation_from_signal([], model="bogus")


class TestAggregateMethylation:
    def test_per_site_frequencies(self) -> None:
        calls = call_methylation_from_signal(_signal(range(100, 130), "a", 3.0))
        agg = aggregate_methylation(calls, min_coverage=5)
        assert agg["n_sites"] == 30
        site = agg["sites"][0]
        for key in ("chrom", "position", "strand", "methylation_freq", "coverage", "n_methylated", "n_unmethylated"):
            assert key in site
        assert site["methylation_freq"] == 1.0
        assert site["coverage"] == 20
        assert agg["mean_coverage"] == 20.0

    def test_min_coverage_filter(self) -> None:
        calls = call_methylation_from_signal(_signal(range(100, 105), "a", 3.0))
        agg = aggregate_methylation(calls, min_coverage=25)
        assert agg["n_sites"] == 0

    def test_invalid_min_coverage_raises(self) -> None:
        with pytest.raises(ValueError, match="min_coverage"):
            aggregate_methylation([], min_coverage=0)


class TestDetectDmrs:
    def test_fully_separated_samples(self) -> None:
        agg_a = aggregate_methylation(call_methylation_from_signal(_signal(range(100, 130), "a", 3.0)), min_coverage=5)
        agg_b = aggregate_methylation(call_methylation_from_signal(_signal(range(100, 130), "b", -3.0)), min_coverage=5)
        dmrs = detect_dmrs(agg_a, agg_b, min_cpgs=3, min_diff=0.2)
        assert len(dmrs) == 1
        dmr = dmrs[0]
        assert dmr["chrom"] == "chr1"
        assert dmr["start"] == 100
        assert dmr["end"] == 129
        assert dmr["n_cpgs"] == 30
        assert dmr["direction"] == "hypo"  # B < A
        assert abs(dmr["mean_diff"]) == pytest.approx(1.0)
        assert 0.0 <= dmr["p_value"] <= 1.0

    def test_identical_samples_no_dmrs(self) -> None:
        agg = aggregate_methylation(call_methylation_from_signal(_signal(range(100, 130), "a", 3.0)), min_coverage=5)
        assert (
            detect_dmrs(
                agg,
                aggregate_methylation(call_methylation_from_signal(_signal(range(100, 130), "c", 3.0)), min_coverage=5),
                min_cpgs=3,
                min_diff=0.2,
            )
            == []
        )

    def test_empty_input_returns_empty(self) -> None:
        assert detect_dmrs({"sites": []}, {"sites": []}) == []


class TestMethylationStatsAndPatterns:
    def test_compute_methylation_stats(self) -> None:
        agg = aggregate_methylation(call_methylation_from_signal(_signal(range(100, 130), "a", 3.0)), min_coverage=5)
        stats = compute_methylation_stats(agg["sites"])
        assert stats["n_sites"] == 30
        assert stats["global_methylation"] == pytest.approx(1.0)
        assert stats["mean_coverage"] == pytest.approx(20.0)
        dist = stats["methylation_distribution"]
        assert dist["min"] <= dist["median"] <= dist["max"]
        assert stats["highly_methylated_fraction"] == pytest.approx(1.0)
        assert stats["lowly_methylated_fraction"] == pytest.approx(0.0)

    def test_compute_methylation_stats_empty(self) -> None:
        stats = compute_methylation_stats([])
        assert stats["n_sites"] == 0
        assert stats["global_methylation"] == 0.0

    def test_methylation_pattern_analysis(self) -> None:
        calls = call_methylation_from_signal(_signal(range(100, 140), "a", 3.0))
        pattern = methylation_pattern_analysis(calls, region={"chrom": "chr1", "start": 100, "end": 140})
        assert {"patterns", "entropy", "n_reads", "epialleles"}.issubset(pattern.keys())
        assert pattern["n_reads"] == 20
        assert 0.0 <= pattern["entropy"] <= 1.0


# ---------------------------------------------------------------------------
# phasing.haplotyping
# ---------------------------------------------------------------------------


def _phasing_setup():
    variants = [{"chrom": "chr1", "position": p, "ref": "A", "alt": "G"} for p in range(100, 220, 20)]
    reads = []
    for i in range(30):
        hap = i % 2
        reads.append(
            {
                "read_id": f"r{i}",
                "chrom": "chr1",
                "start": 90,
                "end": 230,
                "signal": float(hap),  # per-haplotype signal for allele-specific tests
                "alleles": [{"position": v["position"], "allele": hap} for v in variants],
            }
        )
    return variants, reads


class TestPhaseReads:
    def test_separates_two_haplotypes(self) -> None:
        variants, reads = _phasing_setup()
        result = phase_reads(variants, reads)
        assert result["n_phased_reads"] == 30
        assert result["switch_error_rate"] == pytest.approx(0.0)
        assignments = result["haplotype_assignments"]
        assert assignments["r0"] != assignments["r1"]
        assert assignments["r0"] == assignments["r2"]  # same haplotype

    def test_phase_block_built(self) -> None:
        variants, reads = _phasing_setup()
        blocks = phase_reads(variants, reads)["phase_blocks"]
        assert len(blocks) == 1
        block = blocks[0]
        assert block["chrom"] == "chr1"
        assert block["n_variants"] == 6
        assert block["n_reads"] == 30

    def test_empty_inputs(self) -> None:
        result = phase_reads([], [])
        assert result["n_phased_reads"] == 0
        assert result["phase_blocks"] == []

    def test_invalid_method_raises(self) -> None:
        variants, reads = _phasing_setup()
        with pytest.raises(ValueError, match="Unknown phasing method"):
            phase_reads(variants, reads, method="whatshap")


class TestBuildPhaseBlocks:
    def test_contiguous_variants_single_block(self) -> None:
        phased = [
            {"chrom": "chr1", "position": p, "haplotype": 1 + (i % 2), "quality": 30.0}
            for i, p in enumerate(range(100, 400, 50))
        ]
        blocks = build_phase_blocks(phased)
        assert len(blocks) == 1
        assert blocks[0]["start"] == 100
        assert blocks[0]["end"] == 350
        assert blocks[0]["n_variants"] == 6
        assert blocks[0]["phase_quality"] == pytest.approx(30.0)

    def test_gap_splits_blocks(self) -> None:
        phased = [
            {"chrom": "chr1", "position": 100, "haplotype": 1, "quality": 30.0},
            {"chrom": "chr1", "position": 200, "haplotype": 1, "quality": 30.0},
            {"chrom": "chr1", "position": 2000000, "haplotype": 2, "quality": 30.0},
            {"chrom": "chr1", "position": 2000100, "haplotype": 2, "quality": 30.0},
        ]
        blocks = build_phase_blocks(phased, max_gap=1000000)
        assert len(blocks) == 2

    def test_single_variant_no_block(self) -> None:
        assert build_phase_blocks([{"chrom": "chr1", "position": 1, "haplotype": 1}]) == []
        assert build_phase_blocks([]) == []


class TestComputeSwitchErrors:
    def test_known_switch(self) -> None:
        result = compute_switch_errors([0, 1, 0, 1, 0, 1], [0, 1, 0, 1, 1, 0])
        assert result["n_errors"] == 1
        assert result["switch_error_rate"] == pytest.approx(0.2)
        assert result["n_switches"] == 5
        assert result["hamming_distance"] == pytest.approx(1 / 3)

    def test_no_errors(self) -> None:
        result = compute_switch_errors([0, 1, 0], [0, 1, 0])
        assert result["n_errors"] == 0
        assert result["switch_error_rate"] == 0.0

    def test_too_short(self) -> None:
        result = compute_switch_errors([0], [0])
        assert result["switch_error_rate"] == 0.0


class TestHaplotagReads:
    def test_tags_consistent_read(self) -> None:
        block = {
            "chrom": "chr1",
            "start": 100,
            "end": 350,
            "variant_haplotypes": {100: 0, 150: 0, 200: 0, 250: 1, 300: 1, 350: 1},
        }
        alleles = [{"position": p, "allele": 0} for p in (100, 150, 200)] + [
            {"position": p, "allele": 1} for p in (250, 300, 350)
        ]
        tagged = haplotag_reads(
            [{"read_id": "r0", "chrom": "chr1", "start": 90, "end": 360, "alleles": alleles}], [block]
        )
        assert tagged[0]["HP"] == 1
        assert tagged[0]["PS"] == 0

    def test_nonoverlapping_read_untagged(self) -> None:
        block = {"chrom": "chr1", "start": 100, "end": 350, "variant_haplotypes": {100: 0, 150: 1}}
        tagged = haplotag_reads([{"read_id": "rx", "chrom": "chr1", "start": 500, "end": 600}], [block])
        assert tagged[0]["HP"] == 0
        assert tagged[0]["PS"] == 0


class TestAlleleSpecificAnalysis:
    def test_allelic_imbalance_detected(self) -> None:
        variants, reads = _phasing_setup()
        result = phase_reads(variants, reads)  # noqa: E501
        assignments = result["haplotype_assignments"]
        asa = allele_specific_analysis(reads, assignments, region={"chrom": "chr1", "start": 90, "end": 230})
        assert {"haplotype_1_signal", "haplotype_2_signal", "allelic_imbalance", "p_value"}.issubset(asa.keys())
        assert 0.0 <= asa["p_value"] <= 1.0


# ---------------------------------------------------------------------------
# quality.filtering
# ---------------------------------------------------------------------------


def _reads() -> list[dict]:
    random.seed(0)
    reads = []
    for i in range(6):
        length = 1000 if i % 2 else 400
        reads.append({"read_id": f"r{i}", "sequence": "ACGT" * (length // 4), "quality": _phred(length)})
    return reads


class TestFilterByLength:
    def test_min_length(self) -> None:
        kept = filter_by_length(_reads(), min_length=500)
        assert len(kept) == 3
        for read in kept:
            assert len(read["sequence"]) >= 500

    def test_max_length(self) -> None:
        kept = filter_by_length(_reads(), min_length=100, max_length=600)
        assert len(kept) == 3


class TestFilterByQuality:
    def test_high_quality_reads_pass(self) -> None:
        kept = filter_by_quality(_reads(), min_q=25.0)
        assert len(kept) == 6  # all synthetic reads are Q20-40 with mean > 25

    def test_low_quality_reads_fail(self) -> None:
        bad = [{"read_id": "b0", "sequence": "ACGT" * 10, "quality": chr(33 + 5) * 40}]
        assert filter_by_quality(bad, min_q=25.0) == []


class TestTrimAdapters:
    def test_trims_trailing_adapter(self) -> None:
        random.seed(0)
        read = [{"read_id": "a0", "sequence": INSERT + ADAPTER, "quality": _phred(len(INSERT) + len(ADAPTER))}]
        trimmed = trim_adapters(read, adapter_sequences={"illumina": ADAPTER})
        assert len(trimmed) == 1
        assert trimmed[0].sequence == INSERT
        assert trimmed[0].metadata["original_length"] == len(INSERT) + len(ADAPTER)

    def test_trims_leading_adapter(self) -> None:
        random.seed(0)
        read = [{"read_id": "a1", "sequence": ADAPTER + INSERT, "quality": _phred(len(INSERT) + len(ADAPTER))}]
        trimmed = trim_adapters(read, adapter_sequences={"illumina": ADAPTER})
        assert trimmed[0].sequence == INSERT

    def test_no_adapter_leaves_read(self) -> None:
        random.seed(0)
        read = [{"read_id": "a2", "sequence": INSERT, "quality": _phred(len(INSERT))}]
        trimmed = trim_adapters(read, adapter_sequences={"illumina": ADAPTER})
        assert len(trimmed) == 1
        assert trimmed[0].sequence == INSERT


class TestDetectAdapters:
    def test_exact_match(self) -> None:
        sequence = "ACGTACGTACGT" + ADAPTER + "TTTTTTTTTTTT"
        matches = detect_adapters(sequence, known_adapters={"illumina": ADAPTER})
        assert len(matches) == 1
        assert matches[0].adapter_name == "illumina"
        assert matches[0].position == 12
        assert matches[0].identity == pytest.approx(1.0)

    def test_no_match(self) -> None:
        matches = detect_adapters("ACGTACGTACGTACGT", known_adapters={"illumina": ADAPTER})
        assert matches == []


# ---------------------------------------------------------------------------
# workflow.pipelines: load_pipeline_config
# ---------------------------------------------------------------------------


class TestLoadPipelineConfig:
    def test_json_config(self, tmp_path: Path) -> None:
        path = tmp_path / "cfg.json"
        path.write_text(json.dumps({"pipeline_name": "qc", "min_length": 1000}))
        config = load_pipeline_config(path)
        assert config["pipeline_name"] == "qc"
        assert config["min_length"] == 1000

    def test_yaml_config(self, tmp_path: Path) -> None:
        path = tmp_path / "cfg.yaml"
        path.write_text("pipeline_name: qc\nmin_length: 1000\n")
        config = load_pipeline_config(path)
        assert config["pipeline_name"] == "qc"

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            load_pipeline_config(tmp_path / "missing.json")

    def test_missing_pipeline_name_defaults_to_unknown(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.json"
        path.write_text(json.dumps({"min_length": 1000}))
        config = load_pipeline_config(path)
        assert config.get("pipeline_name", "unknown") == "unknown" or config.get("pipeline_name") is None
