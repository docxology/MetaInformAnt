"""Regression tests: ``batch_quality_analysis`` refuses data types it cannot analyze.

Defect: for any ``data_type`` other than ``"fastq"`` the per-file analysis was skipped ("Placeholder for other data
types": ``quality_data = {}``) and the empty data was scored anyway. With ``"vcf"`` or ``"bam"`` every file, even one
that does not exist, came back as ``status == "success"`` with ``overall_score == 0.0`` and grade F, and the summary
counted them as successful analyses. An unknown type was reported once per file as a failure instead of once as a
bad argument. Unsupported types now raise before any file is touched.

Oracles: the per-file scores of the FASTQ path are recomputed directly from ``analyze_fastq_quality`` and
``calculate_quality_score``, and the summary is recomputed from them with the standard library's ``statistics``.
"""

from __future__ import annotations

import statistics
from pathlib import Path

import pytest

from metainformant.quality.analysis.metrics import batch_quality_analysis, calculate_quality_score
from metainformant.quality.io.fastq import analyze_fastq_quality


def _write_fastq(path: Path, quality_char: str, n_reads: int = 20) -> Path:
    path.write_text("".join(f"@read{i}\n{'ACGTACGTAC' * 3}\n+\n{quality_char * 30}\n" for i in range(n_reads)))
    return path


@pytest.mark.parametrize("data_type", ["vcf", "bam"])
def test_data_types_without_a_file_analyzer_raise_not_implemented(tmp_path: Path, data_type: str) -> None:
    missing = tmp_path / f"does_not_exist.{data_type}"  # the old code called this a successful analysis

    with pytest.raises(NotImplementedError, match="FASTQ files only") as excinfo:
        batch_quality_analysis([missing], data_type=data_type)

    assert f"data_type={data_type!r}" in str(excinfo.value)
    assert "calculate_quality_score" in str(excinfo.value)  # the supported way to score such metrics


@pytest.mark.parametrize("data_type", ["xyz", "", "FASTQ", None])
def test_an_unknown_data_type_is_a_value_error_raised_once(tmp_path: Path, data_type: object) -> None:
    files = [tmp_path / "a", tmp_path / "b"]

    with pytest.raises(ValueError, match="Unsupported data type"):
        batch_quality_analysis(files, data_type=data_type)  # type: ignore[arg-type]


@pytest.mark.parametrize("data_type", ["vcf", "bam", "xyz"])
def test_the_data_type_is_checked_even_for_an_empty_file_list(data_type: str) -> None:
    with pytest.raises((NotImplementedError, ValueError)):
        batch_quality_analysis([], data_type=data_type)


def test_the_message_matches_calculate_quality_score_for_an_unknown_type() -> None:
    with pytest.raises(ValueError) as from_batch:
        batch_quality_analysis([], data_type="xyz")
    with pytest.raises(ValueError) as from_score:
        calculate_quality_score({}, "xyz")

    assert str(from_batch.value) == str(from_score.value)


def test_vcf_and_bam_metrics_can_still_be_scored_with_calculate_quality_score() -> None:
    """The alternative named in the error message works, and says honestly that it has nothing to score."""
    score = calculate_quality_score({}, "vcf")

    assert score["overall_score"] == 0.0
    assert "No VCF quality metrics provided" in score["error"]


def test_fastq_files_are_analysed_and_summarised(tmp_path: Path) -> None:
    good = _write_fastq(tmp_path / "good.fastq", "I")  # Phred 40 everywhere
    poor = _write_fastq(tmp_path / "poor.fastq", "+")  # Phred 10 everywhere

    result = batch_quality_analysis([good, poor], data_type="fastq")

    scores = []
    for path in (good, poor):
        entry = result["results"][str(path)]
        direct = calculate_quality_score(analyze_fastq_quality(path), "fastq")
        assert entry["status"] == "success"
        assert entry["quality_score"] == direct
        scores.append(direct["overall_score"])
    assert scores[0] > scores[1]  # the high-quality file scores higher
    assert result["summary"] == {
        "total_files": 2,
        "successful_analyses": 2,
        "failed_analyses": 0,
        "mean_score": statistics.mean(scores),
        "median_score": statistics.median(scores),
        "min_score": min(scores),
        "max_score": max(scores),
    }


def test_a_missing_fastq_file_is_a_failed_analysis_not_a_success(tmp_path: Path) -> None:
    good = _write_fastq(tmp_path / "good.fastq", "I")
    missing = tmp_path / "missing.fastq"

    result = batch_quality_analysis([good, missing], data_type="fastq")

    assert result["results"][str(good)]["status"] == "success"
    assert result["results"][str(missing)]["status"] == "failed"
    assert "does not exist" in result["results"][str(missing)]["error"]
    assert result["summary"]["successful_analyses"] == 1
    assert result["summary"]["failed_analyses"] == 1
    assert result["summary"]["total_files"] == 2


def test_the_default_data_type_is_fastq(tmp_path: Path) -> None:
    good = _write_fastq(tmp_path / "good.fastq", "I")

    assert batch_quality_analysis([good])["results"][str(good)]["status"] == "success"


def test_n_reads_limits_the_reads_analysed_per_file(tmp_path: Path) -> None:
    path = _write_fastq(tmp_path / "many.fastq", "I", n_reads=20)

    result = batch_quality_analysis([path], data_type="fastq", n_reads=5)

    assert result["results"][str(path)]["quality_data"]["basic_statistics"]["total_reads"] == 5


def test_no_files_gives_an_empty_summary() -> None:
    result = batch_quality_analysis([], data_type="fastq")

    assert result["results"] == {}
    assert result["summary"] == {"total_files": 0, "successful_analyses": 0, "failed_analyses": 0}
