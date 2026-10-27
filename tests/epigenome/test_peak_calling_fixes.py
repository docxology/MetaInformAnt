"""Peak merging across chromosomes, merged q-values and the chromosome label of called peaks.

* ``merge_peaks`` compared peaks on ``start`` and ``end`` only, so peaks on different chromosomes that overlapped in
  coordinates were merged into one (and a peak sorted in between two overlapping peaks of another chromosome could
  vanish). The merged peak also kept the *first* peak's q-value while taking the best p-value, signal and summit from
  the strongest member, so p-value and q-value described different peaks.
* ``call_peaks_simple`` and ``call_peaks_broad`` hard-coded ``chrom = "chr_unknown"``.

REAL IMPLEMENTATION.
"""

from __future__ import annotations

from typing import Any, Dict, List

import pytest

from metainformant.epigenome.peak_calling.peak_detection import call_peaks_broad, call_peaks_simple, merge_peaks


def _peak(chrom: str, start: int, end: int, **fields: Any) -> Dict[str, Any]:
    return {"chrom": chrom, "start": start, "end": end, **fields}


class TestMergePeaksKeepsChromosomesApart:
    def test_overlapping_coordinates_on_different_chromosomes_are_not_merged(self) -> None:
        peaks = [
            _peak("chr1", 100, 200, p_value=1e-3, q_value=0.5, signal_value=5.0, score=3.0, fold_enrichment=2.0),
            _peak("chr2", 150, 300, p_value=1e-9, q_value=1e-6, signal_value=9.0, score=9.0, fold_enrichment=6.0),
        ]

        merged = merge_peaks(peaks, distance=0)

        assert merged == peaks

    def test_a_peak_between_two_overlapping_peaks_of_another_chromosome_survives(self) -> None:
        # Sorted by start alone the chr2 peak sat between the two chr1 peaks and was swallowed by them
        peaks = [_peak("chr1", 100, 200), _peak("chr2", 120, 130), _peak("chr1", 150, 250)]

        merged = merge_peaks(peaks, distance=0)

        assert [(p["chrom"], p["start"], p["end"]) for p in merged] == [("chr1", 100, 250), ("chr2", 120, 130)]

    def test_the_merge_distance_does_not_reach_across_chromosomes(self) -> None:
        peaks = [_peak("chr1", 100, 200), _peak("chr2", 205, 300)]

        assert len(merge_peaks(peaks, distance=50)) == 2

    def test_chromosomes_come_out_in_order_of_first_appearance_each_sorted_by_start(self) -> None:
        peaks = [
            _peak("chr2", 500, 600),
            _peak("chr1", 300, 400),
            _peak("chr2", 100, 200),
            _peak("chr1", 10, 20),
        ]

        merged = merge_peaks(peaks, distance=0)

        assert [(p["chrom"], p["start"]) for p in merged] == [("chr2", 100), ("chr2", 500), ("chr1", 10), ("chr1", 300)]

    def test_peaks_without_a_chrom_key_are_one_group(self) -> None:
        # The helper functions of this module treat a missing chrom as "chr_unknown"
        merged = merge_peaks([{"start": 100, "end": 200}, {"start": 150, "end": 250}], distance=0)

        assert [(p["start"], p["end"]) for p in merged] == [(100, 250)]

    def test_input_peaks_are_not_modified(self) -> None:
        peaks = [_peak("chr1", 100, 200, q_value=0.5), _peak("chr1", 150, 250, q_value=0.01)]
        snapshot = [dict(p) for p in peaks]

        merge_peaks(peaks, distance=0)

        assert peaks == snapshot


class TestMergedQValue:
    def test_the_merged_q_value_goes_with_the_retained_p_value(self) -> None:
        # The strong second peak supplies the p-value, signal and summit; the q-value stayed that of the weak first peak
        peaks = [
            _peak("chr1", 100, 200, summit=150, signal_value=5.0, p_value=1e-3, q_value=0.5),
            _peak("chr1", 150, 300, summit=200, signal_value=9.0, p_value=1e-9, q_value=1e-6),
        ]

        (merged,) = merge_peaks(peaks, distance=0)

        assert merged["p_value"] == 1e-9
        assert merged["q_value"] == 1e-6
        assert (merged["summit"], merged["signal_value"]) == (200, 9.0)

    def test_the_best_q_value_wins_whichever_peak_comes_first(self) -> None:
        peaks = [
            _peak("chr1", 100, 200, p_value=1e-9, q_value=1e-6),
            _peak("chr1", 150, 300, p_value=1e-3, q_value=0.5),
        ]

        (merged,) = merge_peaks(peaks, distance=0)

        assert merged["q_value"] == 1e-6

    def test_a_q_value_is_taken_from_a_member_when_the_first_peak_has_none(self) -> None:
        peaks = [_peak("chr1", 100, 200), _peak("chr1", 150, 300, q_value=0.02)]

        (merged,) = merge_peaks(peaks, distance=0)

        assert merged["q_value"] == 0.02

    def test_no_q_value_appears_when_no_member_has_one(self) -> None:
        (merged,) = merge_peaks([_peak("chr1", 100, 200), _peak("chr1", 150, 300)], distance=0)

        assert "q_value" not in merged

    def test_three_peaks_use_the_lowest_q_value_of_all(self) -> None:
        peaks = [
            _peak("chr1", 100, 200, q_value=0.4),
            _peak("chr1", 190, 300, q_value=0.001),
            _peak("chr1", 290, 400, q_value=0.2),
        ]

        (merged,) = merge_peaks(peaks, distance=0)

        assert (merged["start"], merged["end"], merged["q_value"]) == (100, 400, 0.001)


class TestChromosomeOfCalledPeaks:
    @staticmethod
    def _signal() -> List[float]:
        # 200 background bins, 60 strongly enriched, 200 background bins
        return [1.0] * 200 + [60.0] * 60 + [1.0] * 200

    def test_narrow_peaks_carry_the_requested_chromosome(self) -> None:
        peaks = call_peaks_simple(self._signal(), threshold=3.0, min_length=20, chrom="chr7")

        assert len(peaks) == 1
        assert {p["chrom"] for p in peaks} == {"chr7"}
        assert (peaks[0]["start"], peaks[0]["end"]) == (200, 260)

    def test_broad_peaks_carry_the_requested_chromosome(self) -> None:
        peaks = call_peaks_broad(self._signal(), p_threshold=0.05, broad_cutoff=0.1, min_length=20, chrom="chrX")

        assert len(peaks) == 1
        assert {p["chrom"] for p in peaks} == {"chrX"}

    def test_the_placeholder_remains_the_default(self) -> None:
        assert {p["chrom"] for p in call_peaks_simple(self._signal(), threshold=3.0, min_length=20)} == {"chr_unknown"}
        assert {p["chrom"] for p in call_peaks_broad(self._signal(), p_threshold=0.05, min_length=20)} == {
            "chr_unknown"
        }

    def test_an_empty_chromosome_name_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="chrom"):
            call_peaks_simple(self._signal(), chrom="")
        with pytest.raises(ValueError, match="chrom"):
            call_peaks_broad(self._signal(), chrom="")

    def test_peaks_called_on_two_chromosomes_stay_two_peaks_when_merged(self) -> None:
        # Identical coordinates on two chromosomes: every peak was "chr_unknown" before, so merging them gave one
        peaks = call_peaks_simple(self._signal(), threshold=3.0, min_length=20, chrom="chr1") + call_peaks_simple(
            self._signal(), threshold=3.0, min_length=20, chrom="chr2"
        )

        merged = merge_peaks(peaks, distance=0)

        assert [(p["chrom"], p["start"], p["end"]) for p in merged] == [("chr1", 200, 260), ("chr2", 200, 260)]
