"""Bin refinement applies the contamination criterion to every bin it keeps, and binning input is checked.

``refine_bins`` filtered split sub-bins on completeness only and returned a contaminated bin that it
could not split (two contigs or fewer) unchanged; sub-bins also came back with a GC content of 0.0.
``bin_contigs`` clustered on empty feature vectors when asked for coverage binning without coverage,
and invented a coverage of 1 for contigs missing from the coverage table. Real implementations.
"""

from __future__ import annotations

import pytest

from metainformant.metagenomics.shotgun.binning import GenomeBin, bin_contigs, refine_bins

MARKERS = [f"m{i}" for i in range(1, 11)]


def _gc(sequences: list[str]) -> float:
    """Independent GC fraction."""
    joined = "".join(sequences).upper()
    return sum(joined.count(base) for base in "GC") / sum(joined.count(base) for base in "ACGT")


class TestContaminationCriterion:
    """A bin is kept only if it is complete enough AND clean enough."""

    CONTIGS = {"c1": "ATCG" * 500, "c2": "GCTA" * 500}

    def test_contaminated_two_contig_bin_is_removed(self) -> None:
        """Two contigs cannot be split by TNF clustering, so the bin keeps its contamination of 1.0."""
        bins = [GenomeBin("b", ["c1", "c2"], total_length=4000, num_contigs=2)]
        markers = {"c1": MARKERS, "c2": MARKERS}  # every marker twice: completeness 1.0, contamination 1.0
        refined = refine_bins(bins, self.CONTIGS, 0.5, 0.1, marker_genes=markers, marker_set=MARKERS)
        assert refined == []

    def test_clean_two_contig_bin_is_kept(self) -> None:
        bins = [GenomeBin("b", ["c1", "c2"], total_length=4000, num_contigs=2)]
        markers = {"c1": MARKERS[:5], "c2": MARKERS[5:]}
        refined = refine_bins(bins, self.CONTIGS, 0.5, 0.1, marker_genes=markers, marker_set=MARKERS)
        assert [(b.completeness, b.contamination) for b in refined] == [(1.0, 0.0)]

    def test_hand_built_bins_do_not_need_num_contigs(self) -> None:
        """``num_contigs`` defaults to 0 for a hand-built bin; the contig list is what counts."""
        contigs = {
            "x1": "ATCG" * 500,
            "x2": "ATCG" * 400,
            "y1": "GGGGCCCCAAAATTTT" * 125,
            "y2": "GGGGCCCCAAAATTTT" * 100,
        }
        markers = {"x1": MARKERS[:5], "x2": MARKERS[5:], "y1": MARKERS[:5], "y2": MARKERS[5:]}
        merged = GenomeBin("m", list(contigs), total_length=sum(map(len, contigs.values())))
        assert merged.num_contigs == 0
        refined = refine_bins([merged], contigs, 0.9, 0.1, marker_genes=markers, marker_set=MARKERS)
        assert sorted(sorted(b.contig_ids) for b in refined) == [["x1", "x2"], ["y1", "y2"]]

    def test_subbin_that_is_still_contaminated_is_removed(self) -> None:
        """The x group carries two full marker sets, so even after the split it fails the contamination limit."""
        contigs = {
            "x1": "ATCG" * 500,
            "x2": "ATCG" * 400,
            "y1": "GGGGCCCCAAAATTTT" * 125,
            "y2": "GGGGCCCCAAAATTTT" * 100,
        }
        markers = {"x1": MARKERS, "x2": MARKERS, "y1": MARKERS[:5], "y2": MARKERS[5:]}
        merged = GenomeBin("m", list(contigs), total_length=sum(map(len, contigs.values())), num_contigs=4)
        refined = refine_bins([merged], contigs, 0.9, 0.1, marker_genes=markers, marker_set=MARKERS)
        assert [sorted(b.contig_ids) for b in refined] == [["y1", "y2"]]
        assert (refined[0].completeness, refined[0].contamination) == (1.0, 0.0)
        # every bin that is returned meets both criteria
        assert all(b.completeness >= 0.9 and b.contamination <= 0.1 for b in refined)

    def test_thresholds_are_inclusive_limits(self) -> None:
        bins = [GenomeBin("b", ["c1", "c2"], total_length=4000)]
        markers = {"c1": MARKERS, "c2": MARKERS[:1]}  # one extra copy of one marker: contamination 0.1
        kept = refine_bins(bins, self.CONTIGS, 1.0, 0.1, marker_genes=markers, marker_set=MARKERS)
        assert len(kept) == 1
        bins = [GenomeBin("b", ["c1", "c2"], total_length=4000)]
        dropped = refine_bins(bins, self.CONTIGS, 1.0, 0.09, marker_genes=markers, marker_set=MARKERS)
        assert dropped == []


class TestSplitBinsDescribeThemselves:
    """Sub-bins used to come back with a GC content of 0.0 and no coverage."""

    CONTIGS = {
        "x1": "AAAC" * 400,
        "x2": "AAAC" * 300,
        "y1": "GGAC" * 400,
        "y2": "GGAC" * 300,
    }
    MARKERS = {"x1": MARKERS[:5], "x2": MARKERS[5:], "y1": MARKERS[:5], "y2": MARKERS[5:]}

    def _refine(self, coverage: dict[str, float | list[float]] | None = None, seed: int = 42) -> list[GenomeBin]:
        merged = GenomeBin("m", list(self.CONTIGS), total_length=sum(map(len, self.CONTIGS.values())))
        return refine_bins(
            [merged],
            self.CONTIGS,
            0.9,
            0.1,
            marker_genes=self.MARKERS,
            marker_set=MARKERS,
            coverage=coverage,  # type: ignore[arg-type]
            seed=seed,
        )

    def test_gc_content_is_computed_from_the_members(self) -> None:
        refined = self._refine()
        assert sorted(sorted(b.contig_ids) for b in refined) == [["x1", "x2"], ["y1", "y2"]]
        for genome_bin in refined:
            assert genome_bin.gc_content == pytest.approx(_gc([self.CONTIGS[c] for c in genome_bin.contig_ids]))
            assert genome_bin.total_length == sum(len(self.CONTIGS[c]) for c in genome_bin.contig_ids)
            assert genome_bin.num_contigs == 2
        assert sorted(round(b.gc_content, 2) for b in refined) == [0.25, 0.75]

    def test_mean_coverage_comes_from_the_supplied_coverage(self) -> None:
        coverage = {"x1": 10.0, "x2": 20.0, "y1": [4.0, 6.0], "y2": 8.0}  # y1 has two samples: mean 5
        by_contig = {tuple(sorted(b.contig_ids)): b for b in self._refine(coverage=coverage)}
        assert by_contig[("x1", "x2")].mean_coverage == pytest.approx(15.0)
        assert by_contig[("y1", "y2")].mean_coverage == pytest.approx(6.5)

    def test_without_coverage_the_value_is_unknown_zero(self) -> None:
        assert all(b.mean_coverage == 0.0 for b in self._refine())

    def test_a_contig_without_a_sequence_cannot_be_split(self) -> None:
        contigs = {k: v for k, v in self.CONTIGS.items() if k != "y2"}
        merged = GenomeBin("m", list(self.CONTIGS), total_length=1)
        with pytest.raises(ValueError, match="without a sequence.*y2"):
            refine_bins([merged], contigs, 0.9, 0.1, marker_genes=self.MARKERS, marker_set=MARKERS)

    def test_the_split_seed_is_a_parameter(self) -> None:
        # a different k-means seed still separates two well-separated compositions
        refined = self._refine(seed=7)
        assert sorted(sorted(b.contig_ids) for b in refined) == [["x1", "x2"], ["y1", "y2"]]


class TestBinContigsInput:
    CONTIGS = {"c1": "ATCG" * 500, "c2": "GCTA" * 500, "c3": "AAAA" * 500}

    def test_coverage_binning_needs_coverage(self) -> None:
        for coverage in (None, {}):
            with pytest.raises(ValueError, match="requires coverage"):
                bin_contigs(self.CONTIGS, coverage, method="coverage", n_bins=2, min_contig_length=100)

    def test_a_contig_without_a_coverage_value_is_an_error(self) -> None:
        coverage = {"c1": 10.0, "c2": 11.0}  # c3 missing: the old code gave it log-coverage 0 (coverage 1)
        for method in ("coverage", "combined"):
            with pytest.raises(ValueError, match="no value for contig 'c3'"):
                bin_contigs(self.CONTIGS, coverage, method=method, n_bins=2, min_contig_length=100)

    def test_unequal_numbers_of_samples_are_an_error(self) -> None:
        coverage = {"c1": [10.0, 12.0], "c2": [11.0, 9.0], "c3": [50.0]}
        with pytest.raises(ValueError, match="same number of samples"):
            bin_contigs(self.CONTIGS, coverage, method="coverage", n_bins=2, min_contig_length=100)

    def test_coverage_binning_separates_contigs_by_coverage(self) -> None:
        coverage = {"c1": 10.0, "c2": 10.5, "c3": 200.0}
        result = bin_contigs(self.CONTIGS, coverage, method="coverage", n_bins=2, min_contig_length=100)
        groups = sorted(sorted(b.contig_ids) for b in result.bins)
        assert groups == [["c1", "c2"], ["c3"]]

    def test_bin_statistics_are_computed_from_the_members(self) -> None:
        coverage = {"c1": [10.0, 14.0], "c2": [20.0, 30.0], "c3": [50.0, 50.0]}  # per-contig means 12, 25, 50
        result = bin_contigs(self.CONTIGS, coverage, method="combined", n_bins=1, min_contig_length=100)
        (genome_bin,) = result.bins
        assert genome_bin.gc_content == pytest.approx(_gc(list(self.CONTIGS.values())))
        assert genome_bin.mean_coverage == pytest.approx((12.0 + 25.0 + 50.0) / 3.0)

    def test_combined_without_coverage_is_composition_only(self) -> None:
        with_none = bin_contigs(self.CONTIGS, None, method="combined", n_bins=2, min_contig_length=100)
        composition = bin_contigs(self.CONTIGS, None, method="composition", n_bins=2, min_contig_length=100)
        assert [b.contig_ids for b in with_none.bins] == [b.contig_ids for b in composition.bins]
