"""Results that were wrong or dropped without any notice, found while fixing the open epigenome items.

* ``calculate_atac_specific_metrics`` counted a peak of exactly 150 or 250 bp in two length classes, so the three
  nucleosome fractions could add up to more than 1.
* ``merge_overlapping_peaks`` kept the *first* peak's summit while taking the higher signal value of the two, so the
  merged peak's summit was not where its signal value was measured.
* ``find_differentially_methylated_regions`` chained consecutive significant sites whatever the direction of change and
  labelled the region by its first site, so a run of hypermethylated sites followed by hypomethylated ones was one
  "hyper" region.
* ``interpret_states`` matched the mark names to the columns of the emission means by position without checking that
  there were as many names as columns: too few names left the extra columns out of the labels without a word.
* ``filter_peaks`` removed peaks that had no fold enrichment, or neither a q-value nor a p-value, without saying so.
* ``GenomicTrack.normalize_values`` ignored a method name it did not know, changed nothing, and logged that the values
  had been normalized.

REAL IMPLEMENTATION.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List

import pytest

from metainformant.epigenome.analysis.tracks import GenomicTrack
from metainformant.epigenome.assays.atacseq import ATACPeak, calculate_atac_specific_metrics
from metainformant.epigenome.assays.chipseq import ChIPPeak, merge_overlapping_peaks
from metainformant.epigenome.assays.methylation import MethylationSite, find_differentially_methylated_regions
from metainformant.epigenome.chromatin_state.state_learning import interpret_states
from metainformant.epigenome.peak_calling.peak_detection import filter_peaks


class TestNucleosomeLengthClasses:
    def test_a_peak_on_a_class_boundary_is_in_one_class_only(self) -> None:
        lengths = [50, 100, 150, 200, 250, 300, 350, 400]
        peaks = [ATACPeak("chr1", 1000 * i, 1000 * i + length) for i, length in enumerate(lengths)]

        metrics = calculate_atac_specific_metrics(peaks)

        # [50, 150): 50 and 100; [150, 250): 150 and 200; [250, 350): 250 and 300; 350 and 400 are in none
        assert metrics["nfr_peak_fraction"] == 2 / 8
        assert metrics["mononucleosome_peak_fraction"] == 2 / 8
        assert metrics["dinucleosome_peak_fraction"] == 2 / 8

    def test_the_fractions_never_add_up_to_more_than_one(self) -> None:
        peaks = [ATACPeak("chr1", 1000 * i, 1000 * i + length) for i, length in enumerate([150] * 4 + [250] * 4)]

        metrics = calculate_atac_specific_metrics(peaks)

        total = (
            metrics["nfr_peak_fraction"]
            + metrics["mononucleosome_peak_fraction"]
            + metrics["dinucleosome_peak_fraction"]
        )
        assert total == pytest.approx(1.0)


class TestMergedPeakSummit:
    def test_the_merged_summit_is_that_of_the_stronger_peak(self) -> None:
        weak = ChIPPeak("chr1", 100, 200, summit=150, signal_value=5.0, score=10.0)
        strong = ChIPPeak("chr1", 190, 300, summit=250, signal_value=20.0, score=30.0)

        (merged,) = merge_overlapping_peaks([weak, strong])

        assert (merged.start, merged.end) == (100, 300)
        assert merged.signal_value == 20.0
        assert merged.summit == 250  # it used to be 150, where the signal was only 5

    def test_the_first_peak_keeps_its_summit_when_it_is_the_stronger_or_ties(self) -> None:
        first = ChIPPeak("chr1", 100, 200, summit=150, signal_value=20.0)
        second = ChIPPeak("chr1", 190, 300, summit=250, signal_value=5.0)
        tie = ChIPPeak("chr1", 190, 300, summit=260, signal_value=20.0)

        assert merge_overlapping_peaks([first, second])[0].summit == 150
        assert merge_overlapping_peaks([first, tie])[0].summit == 150

    def test_three_peaks_end_with_the_summit_of_the_strongest(self) -> None:
        peaks = [
            ChIPPeak("chr1", 100, 200, summit=120, signal_value=3.0),
            ChIPPeak("chr1", 180, 300, summit=210, signal_value=9.0),
            ChIPPeak("chr1", 290, 400, summit=350, signal_value=6.0),
        ]

        (merged,) = merge_overlapping_peaks(peaks)

        assert merged.summit == 210 and merged.signal_value == 9.0


def _sites(levels: List[float], first: int = 1000) -> List[MethylationSite]:
    return [MethylationSite("chr1", first + 10 * i, round(level * 100), 100) for i, level in enumerate(levels)]


class TestDmrDirection:
    def test_a_run_that_switches_direction_is_cut_where_it_switches(self) -> None:
        condition1 = _sites([0.9, 0.9, 0.9, 0.1, 0.1, 0.1])
        condition2 = _sites([0.1, 0.1, 0.1, 0.9, 0.9, 0.9])

        dmrs = find_differentially_methylated_regions({"chr1": condition1}, {"chr1": condition2}, min_sites=3)

        assert [(d["direction"], d["start"], d["end"], len(d["sites"])) for d in dmrs] == [
            ("hyper", 1000, 1020, 3),
            ("hypo", 1030, 1050, 3),
        ]
        assert all(d["mean_delta"] == pytest.approx(0.8) for d in dmrs)

    def test_a_short_stretch_of_the_other_direction_does_not_join_or_make_a_region(self) -> None:
        condition1 = _sites([0.9, 0.9, 0.9, 0.1, 0.9, 0.9, 0.9])
        condition2 = _sites([0.1, 0.1, 0.1, 0.9, 0.1, 0.1, 0.1])

        dmrs = find_differentially_methylated_regions({"chr1": condition1}, {"chr1": condition2}, min_sites=3)

        assert [(d["direction"], len(d["sites"])) for d in dmrs] == [("hyper", 3), ("hyper", 3)]

    def test_a_uniform_run_is_still_one_region(self) -> None:
        condition1 = _sites([0.9] * 6)
        condition2 = _sites([0.1] * 6)

        (dmr,) = find_differentially_methylated_regions({"chr1": condition1}, {"chr1": condition2}, min_sites=3)

        assert dmr["direction"] == "hyper" and len(dmr["sites"]) == 6


class TestInterpretStatesMarkNames:
    MEANS = {"means": [[10.0, 1.0, 1.0], [1.0, 10.0, 1.0]]}

    def test_fewer_names_than_columns_is_an_error_not_a_silently_ignored_column(self) -> None:
        with pytest.raises(ValueError, match="2 marks but the emission means have 3 columns"):
            interpret_states(self.MEANS, ["H3K4me3", "H3K27me3"])

    def test_more_names_than_columns_is_an_error(self) -> None:
        with pytest.raises(ValueError, match="4 marks but the emission means have 3 columns"):
            interpret_states(self.MEANS, ["H3K4me3", "H3K27me3", "H3K9me3", "H3K36me3"])

    def test_one_name_per_column_is_interpreted(self) -> None:
        interpretations = interpret_states(self.MEANS, ["H3K4me3", "H3K27me3", "H3K9me3"])

        assert [i["state_id"] for i in interpretations] == [0, 1]
        assert interpretations[0]["dominant_marks"] == ["H3K4me3"]

    def test_no_states_is_not_an_error(self) -> None:
        assert interpret_states({"means": []}, ["H3K4me3"]) == []

    def test_a_mark_that_is_the_same_in_every_state_is_dominant_in_none(self) -> None:
        # H3K9me3 is 1.0 in both states: with ">=" against its own mean both states were "high" for it, and a state
        # with nothing else high was labelled Heterochromatin
        interpretations = interpret_states({"means": [[1.0, 1.0], [1.0, 1.0]]}, ["H3K4me3", "H3K9me3"])

        assert [i["dominant_marks"] for i in interpretations] == [[], []]
        assert [i["label"] for i in interpretations] == ["Quiescent", "Quiescent"]


class TestFilterPeaksSaysWhatItRemoved:
    PEAKS: List[Dict[str, Any]] = [
        {"fold_enrichment": 5.0, "q_value": 0.01},
        {"q_value": 0.01},  # no fold enrichment
        {"fold_enrichment": 5.0},  # neither q-value nor p-value
        {"fold_enrichment": 5.0, "p_value": 0.2},  # fails the threshold on its p-value
    ]

    def test_the_peaks_that_cannot_be_judged_are_counted_in_a_warning(self, caplog: pytest.LogCaptureFixture) -> None:
        with caplog.at_level(logging.WARNING):
            kept = filter_peaks(self.PEAKS, min_fold=2.0, max_qvalue=0.05)

        assert kept == [self.PEAKS[0]]
        messages = [record.getMessage() for record in caplog.records if record.levelno == logging.WARNING]
        assert any("1 of 4 peaks have no fold_enrichment" in message for message in messages)
        assert any("1 of 3 peaks have neither q_value nor p_value" in message for message in messages)

    def test_complete_peaks_draw_no_warning(self, caplog: pytest.LogCaptureFixture) -> None:
        with caplog.at_level(logging.WARNING):
            filter_peaks([self.PEAKS[0], self.PEAKS[3]], min_fold=2.0, max_qvalue=0.05)

        assert [r for r in caplog.records if r.levelno >= logging.WARNING] == []


class TestNormalizeValuesMethod:
    @pytest.mark.parametrize("method", ["minmax", "zscore", "robust"])
    def test_the_documented_methods_change_the_values(self, method: str) -> None:
        track = GenomicTrack()
        for i, value in enumerate([1.0, 2.0, 4.0, 8.0, 100.0]):
            track.add_feature("chr1", 10 * i, 10 * i + 5, value=value)

        track.normalize_values(method=method)

        assert [f["value"] for f in track.get_features("chr1")] != [1.0, 2.0, 4.0, 8.0, 100.0]

    @pytest.mark.parametrize("method", ["MinMax", "z-score", "", "median"])
    def test_an_unknown_method_is_an_error_not_a_silent_no_op(self, method: str) -> None:
        track = GenomicTrack()
        track.add_feature("chr1", 0, 10, value=1.0)
        track.add_feature("chr1", 10, 20, value=5.0)

        with pytest.raises(ValueError, match="Unsupported normalization method"):
            track.normalize_values(method=method)

        assert [f["value"] for f in track.get_features("chr1")] == [1.0, 5.0]

    def test_an_unknown_method_is_an_error_even_without_features(self) -> None:
        with pytest.raises(ValueError, match="Unsupported normalization method"):
            GenomicTrack().normalize_values(method="bogus")
