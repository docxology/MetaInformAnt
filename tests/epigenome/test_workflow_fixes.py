"""Epigenome workflow: BED input, the q-value filter, the integration inputs and the ChIP/ATAC analysis.

* A plain ``.bed`` peak file was read as narrowPeak (ten columns needed), so every line was skipped and the workflow
  reported zero peaks without an error.
* ``chipseq_qvalue_threshold`` was also passed to ``filter_peaks_by_score`` as the minimum peak *score*, and the ATAC
  workflow used ``atacseq_qvalue_threshold`` as an accessibility-score floor for peaks without a q-value: peaks with a
  score of 0 (every unscored BED peak) were dropped for it.
* The ``run_*_workflow`` functions never returned the ``sites`` / ``peaks`` lists that ``integrate_epigenome_results``
  reads, so integrating their output produced only "no sites provided" warnings.
* The ChIP-seq/ATAC-seq integration returned a single finding of type ``placeholder``.

Every test runs the real workflow on real files in ``tmp_path``. REAL IMPLEMENTATION.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Dict, List

import pytest

from metainformant.epigenome.assays.chipseq import ChIPPeak, load_chip_peaks
from metainformant.epigenome.workflow import workflow as workflow_module
from metainformant.epigenome.workflow.workflow import (
    EpigenomeConfig,
    _detect_peak_format,
    integrate_epigenome_results,
    run_atacseq_workflow,
    run_chipseq_workflow,
    run_methylation_workflow,
)

# Runs a test once for the ChIP-seq workflow and once for the ATAC-seq workflow
CHIP_AND_ATAC = pytest.mark.parametrize(
    "runner", [run_chipseq_workflow, run_atacseq_workflow], ids=["chipseq", "atacseq"]
)


def _write(directory: Path, name: str, text: str) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    path.write_text(text)
    return path


def _narrowpeak_line(
    chrom: str,
    start: int,
    end: int,
    score: float,
    signal: float,
    p_log10: float | str,
    q_log10: float | str,
    offset: int,
) -> str:
    """One narrowPeak line. ``p_log10`` and ``q_log10`` are the values stored in columns 8 and 9, which are
    ``-log10`` of the probabilities: a q-value of 0.01 is written 2."""
    return f"{chrom}\t{start}\t{end}\tp_{start}\t{score}\t.\t{signal}\t{p_log10}\t{q_log10}\t{offset}\n"


# A MACS2-style line: -log10 p-value 12.3 and -log10 q-value 9.1, so q is about 7.9e-10
MACS2_LINE = "chr1\t100\t500\tpeak1\t50\t.\t7.25\t12.3\t9.1\t200\n"


# ---------------------------------------------------------------------------
# Format detection
# ---------------------------------------------------------------------------


class TestDetectPeakFormatFromContent:
    def test_a_bed_file_is_read_as_bed(self, tmp_path: Path) -> None:
        path = _write(tmp_path, "peaks.bed", "chr1\t100\t200\tp1\t5\t+\n")

        assert _detect_peak_format(path) == "bed"

    def test_a_narrowpeak_named_file_with_nine_columns_is_broadpeak(self, tmp_path: Path) -> None:
        path = _write(tmp_path, "sample.narrowPeak", "chr1\t100\t200\tp1\t5\t.\t1.0\t2.0\t3.0\n")

        assert _detect_peak_format(path) == "broadpeak"

    def test_a_narrowpeak_named_file_with_six_columns_is_bed(self, tmp_path: Path) -> None:
        path = _write(tmp_path, "sample.narrowPeak", "chr1\t100\t200\tp1\t5\t.\n")

        assert _detect_peak_format(path) == "bed"

    def test_a_broadpeak_named_file_with_six_columns_is_bed(self, tmp_path: Path) -> None:
        path = _write(tmp_path, "sample.broadPeak", "chr1\t100\t200\tp1\t5\t.\n")

        assert _detect_peak_format(path) == "bed"

    def test_a_true_narrowpeak_file_stays_narrowpeak(self, tmp_path: Path) -> None:
        path = _write(tmp_path, "sample.narrowPeak", _narrowpeak_line("chr1", 100, 200, 5, 1.0, 2.0, 3.0, 50))

        assert _detect_peak_format(path) == "narrowpeak"

    def test_header_lines_are_not_counted(self, tmp_path: Path) -> None:
        path = _write(
            tmp_path,
            "sample.narrowPeak",
            "# generated\ntrack name=x\nbrowser position chr1\n"
            + _narrowpeak_line("chr1", 100, 200, 5, 1.0, 2.0, 3.0, 50),
        )

        assert _detect_peak_format(path) == "narrowpeak"

    def test_an_unreadable_or_empty_file_falls_back_to_its_name(self, tmp_path: Path) -> None:
        empty = _write(tmp_path, "empty.bed", "")

        assert _detect_peak_format(empty) == "bed"
        assert _detect_peak_format(tmp_path / "missing.broadPeak") == "broadpeak"


# ---------------------------------------------------------------------------
# BED input and the q-value filter
# ---------------------------------------------------------------------------


class TestBedInput:
    def test_a_bed_directory_yields_its_peaks(self, tmp_path: Path) -> None:
        # Three 100 bp peaks; two have score 0
        _write(
            tmp_path / "in",
            "h3k4me3.bed",
            "chr1\t1000\t1100\tp1\t500\t.\nchr1\t1200\t1300\tp2\t0\t.\nchr1\t5000\t5100\tp3\t0\t.\n",
        )

        results = run_chipseq_workflow(tmp_path / "in", tmp_path / "out", EpigenomeConfig())

        assert results["errors"] == []
        assert results["statistics"]["total_peaks"] == 3
        assert [(p["chrom"], p["start"], p["end"]) for p in results["peaks"]] == [
            ("chr1", 1000, 1100),
            ("chr1", 1200, 1300),
            ("chr1", 5000, 5100),
        ]

    def test_atac_bed_peaks_with_no_score_are_kept(self, tmp_path: Path) -> None:
        # The accessibility-score floor (0.05, the q-value threshold) used to remove every peak scored 0
        _write(tmp_path / "in", "open.bed", "chr1\t1000\t1100\tp1\t0\t.\nchr1\t2000\t2100\tp2\t0\t.\n")

        results = run_atacseq_workflow(tmp_path / "in", tmp_path / "out", EpigenomeConfig())

        assert results["statistics"]["total_peaks"] == 2

    def test_a_bed3_file_is_read(self, tmp_path: Path) -> None:
        _write(tmp_path / "in", "regions.bed", "chr1\t1000\t1100\nchr2\t300\t450\n")

        results = run_chipseq_workflow(tmp_path / "in", tmp_path / "out", EpigenomeConfig())

        assert results["statistics"]["total_peaks"] == 2


class TestQValueFilterIsNotAScoreFilter:
    # The q column holds -log10(q): 3 is q = 0.001, 0.30103 is q = 0.5 and 1.39794 is q = 0.04
    CHIP_LINES = (
        _narrowpeak_line("chr1", 1000, 1200, 0, 5.0, 3.0, 3, 100)  # score 0, q 0.001: significant
        + _narrowpeak_line("chr1", 2000, 2200, 900, 5.0, 3.0, 0.30103, 100)  # high score but q 0.5: not significant
        + _narrowpeak_line("chr1", 3000, 3200, 0.01, 5.0, 3.0, 1.39794, 100)  # score 0.01, q 0.04: significant
    )

    def test_chipseq_filters_on_q_value_only(self, tmp_path: Path) -> None:
        # With the q-value threshold also used as a score floor (0.05) the first and third peak were dropped as well
        _write(tmp_path / "in", "sample.narrowPeak", self.CHIP_LINES)

        results = run_chipseq_workflow(tmp_path / "in", tmp_path / "out", EpigenomeConfig())

        assert [p["start"] for p in results["peaks"]] == [1000, 3000]
        assert results["statistics"]["total_peaks"] == 2

    def test_the_threshold_comes_from_the_configuration(self, tmp_path: Path) -> None:
        _write(tmp_path / "in", "sample.narrowPeak", self.CHIP_LINES)

        results = run_chipseq_workflow(tmp_path / "in", tmp_path / "out", EpigenomeConfig(chipseq_qvalue_threshold=0.9))

        assert [p["start"] for p in results["peaks"]] == [1000, 2000, 3000]

    def test_kept_peaks_stay_in_file_order(self, tmp_path: Path) -> None:
        # Peaks used to come back sorted by descending score as a side effect of the score filter
        text = (
            _narrowpeak_line("chr1", 1000, 1200, 10, 5.0, 3.0, 3, 100)
            + _narrowpeak_line("chr1", 2000, 2200, 800, 5.0, 3.0, 3, 100)
            + _narrowpeak_line("chr1", 3000, 3200, 400, 5.0, 3.0, 3, 100)
        )
        _write(tmp_path / "in", "sample.narrowPeak", text)

        results = run_chipseq_workflow(tmp_path / "in", tmp_path / "out", EpigenomeConfig())

        assert [p["start"] for p in results["peaks"]] == [1000, 2000, 3000]

    def test_atacseq_filters_on_q_value_only(self, tmp_path: Path) -> None:
        _write(tmp_path / "in", "open.narrowPeak", self.CHIP_LINES)

        results = run_atacseq_workflow(tmp_path / "in", tmp_path / "out", EpigenomeConfig())

        assert [p["start"] for p in results["peaks"]] == [1000, 3000]

    def test_peaks_without_a_q_value_are_kept_whatever_their_score(self, tmp_path: Path) -> None:
        # BED peaks carry no q-value, so the FDR threshold cannot be applied to them
        _write(tmp_path / "in", "unscored.bed", "chr1\t1000\t1200\tp1\t0\t.\nchr1\t2000\t2200\tp2\t1000\t.\n")

        results = run_chipseq_workflow(tmp_path / "in", tmp_path / "out", EpigenomeConfig())

        assert results["statistics"]["total_peaks"] == 2

    def test_the_length_limits_still_apply(self, tmp_path: Path) -> None:
        too_short = _narrowpeak_line("chr1", 1000, 1010, 5, 5.0, 3.0, 3, 5)  # 10 bp: below the 50 bp minimum
        long_enough = _narrowpeak_line("chr1", 2000, 2200, 5, 5.0, 3.0, 3, 100)
        _write(tmp_path / "in", "sample.narrowPeak", too_short + long_enough)

        results = run_chipseq_workflow(tmp_path / "in", tmp_path / "out", EpigenomeConfig())

        assert [p["start"] for p in results["peaks"]] == [2000]


class TestQColumnIsMinusLog10:
    """A file stores -log10(q); the threshold is an FDR limit: a peak is kept when its column is >= -log10(limit)."""

    @CHIP_AND_ATAC
    def test_a_macs2_style_peak_survives_the_default_fdr_limit(self, tmp_path: Path, runner: Any) -> None:
        # Read verbatim, the q column 9.1 would be "q = 9.1" and fail q <= 0.05; as -log10 it is q of about 7.9e-10
        _write(tmp_path / "in", "sample.narrowPeak", MACS2_LINE)

        results = runner(tmp_path / "in", tmp_path / "out", EpigenomeConfig())

        (peak,) = results["peaks"]
        assert (peak["chrom"], peak["start"], peak["end"]) == ("chr1", 100, 500)
        assert peak["q_value"] == pytest.approx(10**-9.1, rel=1e-12)
        assert peak["q_value"] <= 0.05

    @CHIP_AND_ATAC
    def test_the_workflow_output_is_a_macs2_style_file_again(self, tmp_path: Path, runner: Any) -> None:
        _write(tmp_path / "in", "sample.narrowPeak", MACS2_LINE)

        runner(tmp_path / "in", tmp_path / "out", EpigenomeConfig())

        (line,) = [
            row.split("\t")
            for row in (tmp_path / "out" / "sample_peaks.narrowPeak").read_text().splitlines()
            if not row.startswith("#")
        ]
        assert line[6] == "7.25"
        assert float(line[7]) == pytest.approx(12.3, rel=1e-12)
        assert float(line[8]) == pytest.approx(9.1, rel=1e-12)
        assert line[9] == "200"

    @CHIP_AND_ATAC
    def test_weak_columns_are_dropped_and_strong_ones_kept(self, tmp_path: Path, runner: Any) -> None:
        # -log10(0.05) = 1.30103: a column of 1.31 is q = 0.0490 (kept), 1.30 is q = 0.0501 (dropped), 0 is q = 1
        columns = {1000: 9.1, 2000: 1.31, 3000: 1.30, 4000: 1.0, 5000: 0}
        text = "".join(
            _narrowpeak_line("chr1", start, start + 200, 50, 5.0, 3.0, q, 100) for start, q in columns.items()
        )
        _write(tmp_path / "in", "sample.narrowPeak", text)

        results = runner(tmp_path / "in", tmp_path / "out", EpigenomeConfig())

        assert [p["start"] for p in results["peaks"]] == [1000, 2000]

    @pytest.mark.parametrize(
        ("runner", "option"),
        [(run_chipseq_workflow, "chipseq_qvalue_threshold"), (run_atacseq_workflow, "atacseq_qvalue_threshold")],
        ids=["chipseq", "atacseq"],
    )
    @pytest.mark.parametrize("threshold", [0.2, 0.05, 0.01, 0.1, 0.9])
    def test_a_peak_exactly_on_the_threshold_is_kept(
        self, tmp_path: Path, runner: Any, option: str, threshold: float
    ) -> None:
        # The column -log10(threshold) reads back as threshold plus or minus rounding (for 0.2: 0.20000000000000004),
        # which must not cost the peak. One a millionth of a unit weaker is a real difference and is dropped.
        on_the_line = repr(-math.log10(threshold))
        just_weaker = repr(-math.log10(threshold) - 1e-6)
        text = _narrowpeak_line("chr1", 1000, 1200, 50, 5.0, 3.0, on_the_line, 100)
        text += _narrowpeak_line("chr1", 2000, 2200, 50, 5.0, 3.0, just_weaker, 100)
        _write(tmp_path / "in", "sample.narrowPeak", text)

        results = runner(tmp_path / "in", tmp_path / "out", EpigenomeConfig(**{option: threshold}))

        assert [p["start"] for p in results["peaks"]] == [1000]


class TestQValuesAboveOneAreFlagged:
    """A q-value above 1 is not a probability. The loaders never produce one (they convert the -log10 column of a
    file), so it can only come from peaks built in memory with -log10 values."""

    def test_in_memory_peaks_holding_log10_values_trigger_a_warning(self, caplog: pytest.LogCaptureFixture) -> None:
        peaks = [
            ChIPPeak("chr1", 0, 100, q_value=5.2),  # a -log10 value where a probability belongs
            ChIPPeak("chr1", 200, 300, q_value=0.01),
            ChIPPeak("chr1", 400, 500),
        ]

        with caplog.at_level(logging.WARNING):
            # looked up here, not imported at the top, so that only this test depends on the helper
            kept = workflow_module._filter_by_qvalue(Path("peaks.bed"), peaks, 0.05)

        messages = [record.getMessage() for record in caplog.records]
        assert any("not probabilities" in message and "10 ** -value" in message for message in messages)
        # the 5.2 is judged as the probability it claims to be and falls above the limit; no q-value is kept
        assert [p.start for p in kept] == [200, 400]

    @CHIP_AND_ATAC
    def test_a_file_column_above_one_is_a_probability_and_warns_nothing(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture, runner: Any
    ) -> None:
        # The same 5.2, read from a file, is -log10 of a q-value of 6.3e-6: a strong peak, nothing to warn about
        _write(tmp_path / "in", "sample.narrowPeak", _narrowpeak_line("chr1", 1000, 1200, 500, 5.0, 12.3, 5.2, 100))

        with caplog.at_level(logging.WARNING):
            results = runner(tmp_path / "in", tmp_path / "out", EpigenomeConfig())

        assert not any("not probabilities" in record.getMessage() for record in caplog.records)
        (peak,) = results["peaks"]
        assert peak["q_value"] == pytest.approx(10**-5.2, rel=1e-12)

    @CHIP_AND_ATAC
    def test_ordinary_columns_warn_nothing(self, tmp_path: Path, caplog: pytest.LogCaptureFixture, runner: Any) -> None:
        _write(tmp_path / "in", "sample.narrowPeak", _narrowpeak_line("chr1", 1000, 1200, 500, 5.0, 3, 2, 100))

        with caplog.at_level(logging.WARNING):
            runner(tmp_path / "in", tmp_path / "out", EpigenomeConfig())

        assert not any("not probabilities" in record.getMessage() for record in caplog.records)

    def test_peaks_without_a_q_value_are_counted_in_a_warning(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        _write(tmp_path / "in", "unscored.bed", "chr1\t1000\t1200\tp1\t0\t.\nchr1\t2000\t2200\tp2\t0\t.\n")

        with caplog.at_level(logging.WARNING):
            run_chipseq_workflow(tmp_path / "in", tmp_path / "out", EpigenomeConfig())

        assert any("2 of 2 peaks have no q-value" in record.getMessage() for record in caplog.records)


class TestWorkflowOutputFilesRoundTrip:
    def test_the_saved_peak_file_loads_back_as_the_kept_peaks(self, tmp_path: Path) -> None:
        # The workflow saves with save_chip_peaks; the file used to come back empty from load_chip_peaks.
        # -log10 columns: p 3 and 4 (1e-3, 1e-4), q 3 and 2.7 (1e-3, about 2e-3)
        _write(
            tmp_path / "in",
            "sample.narrowPeak",
            _narrowpeak_line("chr1", 1000, 1200, 55, 5.5, 3.0, 3, 100)
            + _narrowpeak_line("chr2", 2000, 2400, 60, 6.5, 4.0, 2.7, 150),
        )

        results = run_chipseq_workflow(tmp_path / "in", tmp_path / "out", EpigenomeConfig())

        saved = load_chip_peaks(tmp_path / "out" / "sample_peaks.narrowPeak")
        assert [(p.chromosome, p.start, p.end, p.summit, p.signal_value) for p in saved] == [
            ("chr1", 1000, 1200, 1100, 5.5),
            ("chr2", 2000, 2400, 2150, 6.5),
        ]
        assert [p.p_value for p in saved] == pytest.approx([1e-3, 1e-4], rel=1e-12)
        assert [p.q_value for p in saved] == pytest.approx([1e-3, 10**-2.7], rel=1e-12)
        assert len(saved) == results["statistics"]["total_peaks"]


# ---------------------------------------------------------------------------
# What the workflows hand to the integration step
# ---------------------------------------------------------------------------


def _assay_directories(tmp_path: Path) -> Dict[str, Path]:
    """Twenty methylation sites 10 bp apart from 1000, three ChIP-seq BED peaks and one ATAC-seq narrowPeak."""
    methylation = tmp_path / "methylation"
    _write(
        methylation,
        "s1.bedgraph",
        "".join(f"chr1\t{1000 + 10 * i}\t{1001 + 10 * i}\t{0.1 + 0.04 * i:.2f}\n" for i in range(20)),
    )
    chip = tmp_path / "chip"
    _write(
        chip, "h3k4me3.bed", "chr1\t1000\t1100\tp1\t500\t.\nchr1\t1200\t1300\tp2\t0\t.\nchr1\t5000\t5100\tp3\t0\t.\n"
    )
    atac = tmp_path / "atac"
    _write(atac, "open.narrowPeak", _narrowpeak_line("chr1", 1050, 1250, 100, 5.0, 3.0, 3, 100))  # q = 10^-3
    return {"methylation": methylation, "chip": chip, "atac": atac}


class TestWorkflowResultsFeedTheIntegration:
    def test_the_methylation_workflow_returns_its_sites(self, tmp_path: Path) -> None:
        directories = _assay_directories(tmp_path)

        results = run_methylation_workflow(directories["methylation"], tmp_path / "o1", EpigenomeConfig())

        assert len(results["sites"]) == results["statistics"]["total_sites"] == 20
        first = results["sites"][0]
        assert first["sample"] == "s1"
        assert (first["chrom"], first["position"]) == ("chr1", 1000)
        assert first["methylation_level"] == pytest.approx(0.1)
        assert first["coverage"] == 100

    def test_the_peak_workflows_return_their_peaks(self, tmp_path: Path) -> None:
        directories = _assay_directories(tmp_path)

        chip = run_chipseq_workflow(directories["chip"], tmp_path / "o2", EpigenomeConfig())
        atac = run_atacseq_workflow(directories["atac"], tmp_path / "o3", EpigenomeConfig())

        assert [(p["sample"], p["chrom"], p["start"], p["end"]) for p in chip["peaks"]] == [
            ("h3k4me3", "chr1", 1000, 1100),
            ("h3k4me3", "chr1", 1200, 1300),
            ("h3k4me3", "chr1", 5000, 5100),
        ]
        assert [(p["chrom"], p["start"], p["end"]) for p in atac["peaks"]] == [("chr1", 1050, 1250)]
        assert [p["q_value"] for p in atac["peaks"]] == pytest.approx([0.001], rel=1e-12)  # read from the column 3

    def test_integrating_the_workflow_output_gives_real_statistics(self, tmp_path: Path) -> None:
        directories = _assay_directories(tmp_path)
        config = EpigenomeConfig()
        methylation = run_methylation_workflow(directories["methylation"], tmp_path / "o1", config)
        chip = run_chipseq_workflow(directories["chip"], tmp_path / "o2", config)
        atac = run_atacseq_workflow(directories["atac"], tmp_path / "o3", config)

        integrated = integrate_epigenome_results(methylation, chip, atac, tmp_path / "integrated", config)

        analyses = integrated["integrated_analyses"]
        assert integrated["errors"] == []
        for analysis in analyses.values():
            assert all(finding["type"] != "warning" for finding in analysis["findings"])

        # Sites at 1000, 1010, ..., 1190. ChIP-seq peaks [1000, 1100], [1200, 1300], [5000, 5100]: the eleven sites
        # from 1000 to 1100 lie inside a peak and the other nine are within 1000 bp of [1200, 1300].
        chip_stats = analyses["methylation_chipseq"]["statistics"]
        assert (chip_stats["overlapping"], chip_stats["proximal"], chip_stats["distal"]) == (11, 9, 0)

        # ATAC-seq region [1050, 1250]: fifteen sites inside it (1050 to 1190), five (1000 to 1040) within 1000 bp
        atac_stats = analyses["methylation_atacseq"]["statistics"]
        assert (
            atac_stats["sites_in_accessible"],
            atac_stats["sites_near_accessible"],
            atac_stats["sites_in_closed"],
        ) == (15, 5, 0)
        assert atac_stats["mean_methylation_accessible"] == pytest.approx(0.58)
        assert atac_stats["mean_methylation_near_accessible"] == pytest.approx(0.18)

        # The ChIP-seq/ATAC-seq analysis reports real overlap numbers, not a placeholder
        chip_atac = analyses["chipseq_atacseq"]
        assert [finding["type"] for finding in chip_atac["findings"]] == ["overlap"]
        assert chip_atac["statistics"]["chip_peaks_overlapping_atac"] == 2


# ---------------------------------------------------------------------------
# ChIP-seq / ATAC-seq integration
# ---------------------------------------------------------------------------


def _integrate_chip_atac(
    tmp_path: Path, chip_peaks: List[Dict[str, Any]], atac_peaks: List[Dict[str, Any]], atac_key: str = "peaks"
) -> Dict[str, Any]:
    result = integrate_epigenome_results(
        {"statistics": {"total_samples": 0}},
        {"statistics": {"total_samples": 1}, "peaks": chip_peaks},
        {"statistics": {"total_samples": 1}, atac_key: atac_peaks},
        tmp_path / "integrated",
        EpigenomeConfig(),
    )
    analysis: Dict[str, Any] = result["integrated_analyses"]["chipseq_atacseq"]
    return analysis


class TestChipAtacIntegration:
    def test_overlap_counts_and_base_pair_jaccard(self, tmp_path: Path) -> None:
        # ChIP-seq [1000,1100) [1200,1300) [5000,5100); ATAC-seq [1050,1250)
        # shared bases 50 + 50 = 100; ChIP-seq bases 300, ATAC-seq bases 200; union 300 + 200 - 100 = 400
        chip = [{"chrom": "chr1", "start": s, "end": s + 100} for s in (1000, 1200, 5000)]
        atac = [{"chrom": "chr1", "start": 1050, "end": 1250}]

        analysis = _integrate_chip_atac(tmp_path, chip, atac)

        stats = analysis["statistics"]
        assert stats["total_chip_peaks"] == 3
        assert stats["total_atac_peaks"] == 1
        assert stats["chip_peaks_overlapping_atac"] == 2
        assert stats["atac_peaks_overlapping_chip"] == 1
        assert stats["chip_overlap_rate"] == pytest.approx(2 / 3)
        assert stats["atac_overlap_rate"] == pytest.approx(1.0)
        assert (stats["chip_bases"], stats["atac_bases"], stats["shared_bases"], stats["union_bases"]) == (
            300,
            200,
            100,
            400,
        )
        assert stats["jaccard_index"] == pytest.approx(0.25)
        (finding,) = analysis["findings"]
        assert finding["type"] == "overlap"
        assert finding["count"] == 2
        assert finding["jaccard_index"] == pytest.approx(0.25)

    def test_no_placeholder_finding_remains(self, tmp_path: Path) -> None:
        analysis = _integrate_chip_atac(
            tmp_path, [{"chrom": "chr1", "start": 0, "end": 10}], [{"chrom": "chr1", "start": 5, "end": 15}]
        )

        assert all(finding["type"] != "placeholder" for finding in analysis["findings"])
        assert "note" not in analysis["findings"][0]

    def test_disjoint_peak_sets_have_zero_overlap(self, tmp_path: Path) -> None:
        analysis = _integrate_chip_atac(
            tmp_path, [{"chrom": "chr1", "start": 0, "end": 100}], [{"chrom": "chr1", "start": 100, "end": 200}]
        )

        stats = analysis["statistics"]
        # touching intervals share no base
        assert (stats["chip_peaks_overlapping_atac"], stats["atac_peaks_overlapping_chip"], stats["shared_bases"]) == (
            0,
            0,
            0,
        )
        assert stats["jaccard_index"] == 0.0

    def test_peaks_on_different_chromosomes_do_not_overlap(self, tmp_path: Path) -> None:
        analysis = _integrate_chip_atac(
            tmp_path, [{"chrom": "chr1", "start": 0, "end": 100}], [{"chrom": "chr2", "start": 0, "end": 100}]
        )

        assert analysis["statistics"]["shared_bases"] == 0
        assert analysis["statistics"]["union_bases"] == 200

    def test_identical_peak_sets_have_a_jaccard_index_of_one(self, tmp_path: Path) -> None:
        peaks = [{"chrom": "chr1", "start": 0, "end": 100}, {"chrom": "chr2", "start": 50, "end": 80}]

        analysis = _integrate_chip_atac(tmp_path, peaks, [dict(p) for p in peaks])

        assert analysis["statistics"]["jaccard_index"] == pytest.approx(1.0)
        assert analysis["statistics"]["chip_peaks_overlapping_atac"] == 2

    def test_overlapping_peaks_inside_an_assay_count_their_bases_once(self, tmp_path: Path) -> None:
        # ChIP-seq [0,100) and [50,150) cover 150 bases together, not 200. ATAC-seq [100,200): shared bases 50.
        chip = [{"chrom": "chr1", "start": 0, "end": 100}, {"chrom": "chr1", "start": 50, "end": 150}]
        atac = [{"chrom": "chr1", "start": 100, "end": 200}]

        stats = _integrate_chip_atac(tmp_path, chip, atac)["statistics"]

        assert (stats["chip_bases"], stats["atac_bases"], stats["shared_bases"], stats["union_bases"]) == (
            150,
            100,
            50,
            200,
        )
        assert stats["jaccard_index"] == pytest.approx(0.25)
        # only the second ChIP-seq peak reaches the ATAC-seq region
        assert stats["chip_peaks_overlapping_atac"] == 1
        assert stats["atac_peaks_overlapping_chip"] == 1

    def test_a_long_peak_hides_nothing_behind_a_short_one(self, tmp_path: Path) -> None:
        # ChIP-seq [0,1000) contains [10,20); the ATAC-seq peak [500,600) overlaps only the long one
        chip = [{"chrom": "chr1", "start": 0, "end": 1000}, {"chrom": "chr1", "start": 10, "end": 20}]
        atac = [{"chrom": "chr1", "start": 500, "end": 600}]

        stats = _integrate_chip_atac(tmp_path, chip, atac)["statistics"]

        assert stats["chip_peaks_overlapping_atac"] == 1
        assert stats["atac_peaks_overlapping_chip"] == 1

    def test_accessible_regions_is_accepted_for_atacseq(self, tmp_path: Path) -> None:
        analysis = _integrate_chip_atac(
            tmp_path,
            [{"chrom": "chr1", "start": 0, "end": 100}],
            [{"chrom": "chr1", "start": 50, "end": 150}],
            atac_key="accessible_regions",
        )

        assert analysis["statistics"]["shared_bases"] == 50

    def test_chromosome_and_position_aliases_are_read(self, tmp_path: Path) -> None:
        analysis = _integrate_chip_atac(
            tmp_path,
            [{"chromosome": "chr1", "start": 0, "end": 100}],
            [{"chromosome": "chr1", "position": 50, "end": 150}],
        )

        assert analysis["statistics"]["shared_bases"] == 50

    def test_missing_peaks_are_reported_as_warnings(self, tmp_path: Path) -> None:
        no_chip = _integrate_chip_atac(tmp_path, [], [{"chrom": "chr1", "start": 0, "end": 10}])
        no_atac = _integrate_chip_atac(tmp_path, [{"chrom": "chr1", "start": 0, "end": 10}], [])

        for analysis in (no_chip, no_atac):
            assert [finding["type"] for finding in analysis["findings"]] == ["warning"]
            assert analysis["statistics"] == {}

    def test_malformed_peaks_are_skipped(self, tmp_path: Path) -> None:
        chip = [
            {"chrom": "chr1", "start": 0, "end": 100},
            {"chrom": "chr1", "start": 5},
            {"start": 0, "end": 10},
            {"chrom": "chr1", "start": 50, "end": 50},
        ]
        atac = [{"chrom": "chr1", "start": 50, "end": 150}]

        stats = _integrate_chip_atac(tmp_path, chip, atac)["statistics"]

        assert stats["total_chip_peaks"] == 1
        assert stats["shared_bases"] == 50

    def test_only_malformed_peaks_give_a_warning(self, tmp_path: Path) -> None:
        analysis = _integrate_chip_atac(
            tmp_path, [{"chrom": "chr1", "start": 5}], [{"chrom": "chr1", "start": 0, "end": 10}]
        )

        assert [finding["type"] for finding in analysis["findings"]] == ["warning"]

    def test_the_peak_objects_of_the_assay_modules_round_trip_into_the_analysis(self, tmp_path: Path) -> None:
        # ChIPPeak.to_dict() uses "chromosome", "start" and "end": the analysis reads them as they are
        peaks = [ChIPPeak("chr1", 0, 100), ChIPPeak("chr1", 200, 300)]

        analysis = _integrate_chip_atac(
            tmp_path, [p.to_dict() for p in peaks], [{"chrom": "chr1", "start": 50, "end": 250}]
        )

        assert analysis["statistics"]["chip_peaks_overlapping_atac"] == 2
        assert analysis["statistics"]["shared_bases"] == 100
