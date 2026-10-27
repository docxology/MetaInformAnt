"""Methylation entropy, DMR p-values and BEDgraph level scales.

* ``calculate_methylation_entropy`` binned levels with ``lo <= level < hi`` over [0, 1], so a fully methylated site
  (1.0) fell into no bin yet still counted toward the window total, understating the entropy.
* ``find_differentially_methylated_regions`` built its standard error from each sample's own proportion. That is
  exactly zero whenever both proportions are 0 or 1, however different they are, and the p-value was then set to 1:
  the most different sites (0% against 100%) could never be significant. ``mean_p_value`` also held the smallest
  p-value rather than a mean.
* ``load_methylation_bedgraph`` rejected the whole file when one level was outside 0-1, so a percentage file (0-100)
  could not be loaded at all.

Expected p-values come from an exact integer computation (hypergeometric weights from ``math.comb``), independent of
the implementation. REAL IMPLEMENTATION.
"""

from __future__ import annotations

import math
from fractions import Fraction
from pathlib import Path
from statistics import NormalDist
from typing import Dict, List, Tuple

import pytest

from metainformant.epigenome.assays import methylation as methylation_module
from metainformant.epigenome.assays.methylation import (
    MethylationSite,
    calculate_methylation_entropy,
    find_differentially_methylated_regions,
    load_methylation_bedgraph,
)


def _two_proportion_p_value(methylated1: int, total1: int, methylated2: int, total2: int) -> float:
    """The p-value helper under test, looked up when called so that only the tests of the helper depend on it."""
    return float(methylation_module._two_proportion_p_value(methylated1, total1, methylated2, total2))


def _site(position: int, methylated: int, total: int, chrom: str = "chr1") -> MethylationSite:
    return MethylationSite(chromosome=chrom, position=position, methylated_reads=methylated, total_reads=total)


def _fisher_reference(methylated1: int, total1: int, methylated2: int, total2: int) -> float:
    """Two-sided Fisher exact p-value by exact integer arithmetic: the tables no more likely than the observed one."""
    methylated = methylated1 + methylated2
    weights = {
        x: math.comb(total1, x) * math.comb(total2, methylated - x)
        for x in range(max(0, methylated - total2), min(total1, methylated) + 1)
    }
    observed = weights[methylated1]
    return float(Fraction(sum(w for w in weights.values() if w <= observed), math.comb(total1 + total2, methylated)))


def _z_reference(methylated1: int, total1: int, methylated2: int, total2: int) -> float:
    """Pooled two-proportion z-test with continuity correction, computed with the normal distribution."""
    pooled = (methylated1 + methylated2) / (total1 + total2)
    standard_error = math.sqrt(pooled * (1 - pooled) * (1 / total1 + 1 / total2))
    difference = abs(methylated1 / total1 - methylated2 / total2) - 0.5 * (1 / total1 + 1 / total2)
    return 2 * (1 - NormalDist().cdf(max(0.0, difference) / standard_error))


def _chromosome(sites_spec: List[tuple], start: int = 1000, step: int = 10) -> Dict[str, List[MethylationSite]]:
    """One chromosome of sites from (methylated, total) pairs placed ``step`` bp apart."""
    return {"chr1": [_site(start + i * step, m, t) for i, (m, t) in enumerate(sites_spec)]}


# ---------------------------------------------------------------------------
# Methylation entropy
# ---------------------------------------------------------------------------


class TestEntropyCountsFullyMethylatedSites:
    @staticmethod
    def _profile(levels_as_reads: List[int], count: int = 60) -> List[float]:
        """Entropies of 200 bp windows over sites 10 bp apart cycling through methylated-read counts (of 10 reads)."""
        sites = [_site(1000 + i * 10, levels_as_reads[i % len(levels_as_reads)], 10) for i in range(count)]
        profile = calculate_methylation_entropy({"chr1": sites}, window_size=200)
        return [entropy for _, entropy in profile["chr1"]]

    def test_half_fully_methylated_and_half_low_is_one_bit(self) -> None:
        # Every window holds ten sites at 1.0 and ten at 0.1: two bins of one half each. The 1.0 sites fell in no bin,
        # so the old code saw p = 0.5 for the low bin alone and returned 0.5 bits.
        entropies = self._profile([10, 1])

        assert len(entropies) > 1
        assert entropies == pytest.approx([1.0] * len(entropies))

    def test_all_five_bins_including_the_last_are_balanced(self) -> None:
        # Levels 0.1, 0.3, 0.5, 0.7, 1.0 fill the five bins evenly: log2(5) bits (the old code gave 1.86)
        entropies = self._profile([1, 3, 5, 7, 10])

        assert entropies == pytest.approx([math.log2(5)] * len(entropies))

    def test_only_fully_methylated_sites_have_no_entropy(self) -> None:
        assert self._profile([10]) == pytest.approx([0.0] * len(self._profile([10])))

    def test_a_site_just_below_one_shares_the_last_bin_with_one(self) -> None:
        # 0.9 and 1.0 are both in the last bin [0.8, 1.0], so a window of them has one occupied bin: 0 bits. The 1.0
        # sites used to be dropped, leaving the 0.9 sites at p = 0.5 and an entropy of 0.5 bits.
        entropies = self._profile([9, 10])

        assert entropies == pytest.approx([0.0] * len(entropies))


# ---------------------------------------------------------------------------
# DMR p-values
# ---------------------------------------------------------------------------


class TestDifferentialMethylationTest:
    def test_fully_methylated_against_unmethylated_sites_form_a_dmr(self) -> None:
        # 0 of 10 reads against 10 of 10 at five consecutive sites. Both proportions are 0 or 1, so the old standard
        # error was 0, the p-value 1, and no region was reported.
        condition1 = _chromosome([(0, 10)] * 5)
        condition2 = _chromosome([(10, 10)] * 5)

        dmrs = find_differentially_methylated_regions(condition1, condition2, min_sites=3)

        assert len(dmrs) == 1
        (dmr,) = dmrs
        assert (dmr["start"], dmr["end"], dmr["sites"]) == (1000, 1040, [1000, 1010, 1020, 1030, 1040])
        assert dmr["direction"] == "hypo"
        assert dmr["mean_delta"] == pytest.approx(1.0)
        # exactly 2 / C(20, 10) per site
        assert dmr["mean_p_value"] == pytest.approx(2 / math.comb(20, 10), rel=1e-9)
        assert dmr["min_p_value"] == pytest.approx(2 / math.comb(20, 10), rel=1e-9)

    def test_deep_coverage_gives_a_vanishing_p_value_not_one(self) -> None:
        condition1 = _chromosome([(0, 200)] * 4)
        condition2 = _chromosome([(200, 200)] * 4)

        (dmr,) = find_differentially_methylated_regions(condition1, condition2, min_sites=3)

        assert 0.0 < dmr["mean_p_value"] < 1e-20

    def test_a_few_reads_at_the_extremes_are_not_significant(self) -> None:
        # 0 of 2 reads against 2 of 2: Fisher's exact p is 1/3. (A z-test would call this significant, p = 0.046.)
        condition1 = _chromosome([(0, 2)] * 5)
        condition2 = _chromosome([(2, 2)] * 5)

        assert find_differentially_methylated_regions(condition1, condition2, min_sites=3) == []
        assert _two_proportion_p_value(0, 2, 2, 2) == pytest.approx(1 / 3)

    def test_mean_p_value_is_the_mean_and_min_p_value_the_minimum(self) -> None:
        # Three consecutive significant sites with different p-values; mean_p_value held only the smallest of them
        tables = [(0, 10, 10, 10), (1, 10, 9, 10), (2, 10, 8, 10)]
        condition1 = _chromosome([(m1, n1) for m1, n1, _, _ in tables])
        condition2 = _chromosome([(m2, n2) for _, _, m2, n2 in tables])
        expected = [_fisher_reference(*table) for table in tables]
        assert expected[0] < expected[1] < expected[2] < 0.05

        (dmr,) = find_differentially_methylated_regions(condition1, condition2, min_sites=3)

        assert dmr["mean_p_value"] == pytest.approx(sum(expected) / 3, rel=1e-9)
        assert dmr["min_p_value"] == pytest.approx(min(expected), rel=1e-9)
        assert dmr["mean_p_value"] > dmr["min_p_value"]

    def test_the_running_p_value_list_is_not_leaked(self) -> None:
        condition1 = _chromosome([(0, 10)] * 4)
        condition2 = _chromosome([(10, 10)] * 4)

        (dmr,) = find_differentially_methylated_regions(condition1, condition2, min_sites=3)

        assert set(dmr) == {
            "chromosome",
            "start",
            "end",
            "sites",
            "mean_delta",
            "mean_p_value",
            "min_p_value",
            "direction",
        }

    def test_regions_split_by_max_gap_each_get_their_own_p_values(self) -> None:
        near = [(0, 10)] * 3
        far = [(2, 10)] * 3
        condition1 = {
            "chr1": [_site(1000 + i * 10, m, t) for i, (m, t) in enumerate(near)]
            + [_site(900000 + i * 10, m, t) for i, (m, t) in enumerate(far)]
        }
        condition2 = {
            "chr1": [_site(1000 + i * 10, 10, 10) for i in range(3)] + [_site(900000 + i * 10, 8, 10) for i in range(3)]
        }

        first, second = find_differentially_methylated_regions(condition1, condition2, min_sites=3, max_gap=1000)

        assert first["mean_p_value"] == pytest.approx(_fisher_reference(0, 10, 10, 10), rel=1e-9)
        assert second["mean_p_value"] == pytest.approx(_fisher_reference(2, 10, 8, 10), rel=1e-9)


class TestTwoProportionPValue:
    @pytest.mark.parametrize("sizes", [(4, 6), (5, 5), (8, 3), (1, 7)])
    def test_small_tables_equal_the_exact_test(self, sizes: tuple) -> None:
        total1, total2 = sizes
        for methylated1 in range(total1 + 1):
            for methylated2 in range(total2 + 1):
                got = _two_proportion_p_value(methylated1, total1, methylated2, total2)
                assert got == pytest.approx(
                    _fisher_reference(methylated1, total1, methylated2, total2), rel=1e-9, abs=1e-12
                )

    def test_tea_tasting_table(self) -> None:
        # The classic 2x2 example, [[3, 1], [1, 3]]: two-sided p = 34/70
        assert _two_proportion_p_value(3, 4, 1, 4) == pytest.approx(34 / 70, rel=1e-9)

    @pytest.mark.parametrize("table", [(40, 100, 55, 100), (60, 200, 90, 200), (30, 60, 45, 70)])
    def test_large_tables_use_the_continuity_corrected_z_test(self, table: tuple) -> None:
        assert _two_proportion_p_value(*table) == pytest.approx(_z_reference(*table), rel=1e-9)
        # and agree with the exact test to within the usual difference between the two
        assert 0.8 <= _two_proportion_p_value(*table) / _fisher_reference(*table) <= 1.6

    def test_symmetric_in_the_two_samples(self) -> None:
        for table, swapped in [((3, 12, 9, 15), (9, 15, 3, 12)), ((40, 100, 55, 100), (55, 100, 40, 100))]:
            assert _two_proportion_p_value(*table) == pytest.approx(_two_proportion_p_value(*swapped), rel=1e-9)

    def test_identical_reads_cannot_differ(self) -> None:
        # Nothing methylated at all, everything methylated, and equal proportions: p = 1, with no division by zero
        assert _two_proportion_p_value(0, 10, 0, 30) == 1.0
        assert _two_proportion_p_value(10, 10, 30, 30) == 1.0
        assert _two_proportion_p_value(50, 100, 50, 100) == pytest.approx(1.0)

    def test_a_sample_without_reads_gives_one(self) -> None:
        assert _two_proportion_p_value(0, 0, 5, 10) == 1.0
        assert _two_proportion_p_value(5, 10, 0, 0) == 1.0

    def test_p_value_stays_in_range_for_deep_coverage(self) -> None:
        # The true p-value is far below the smallest float; it underflows to 0.0 instead of failing
        got = _two_proportion_p_value(0, 100000, 100000, 100000)

        assert 0.0 <= got < 1e-300


# ---------------------------------------------------------------------------
# BEDgraph level scale
# ---------------------------------------------------------------------------


def _bedgraph(tmp_path: Path, rows: List[Tuple[int, str]], name: str = "levels.bedgraph") -> Path:
    path = tmp_path / name
    path.write_text("".join(f"chr1\t{pos}\t{pos + 1}\t{level}\n" for pos, level in rows))
    return path


class TestBedgraphLevelScale:
    def test_a_percentage_file_loads_and_is_scaled(self, tmp_path: Path) -> None:
        # 80.0 as a fraction would need 8000 methylated reads out of 100, which aborted the whole load
        path = _bedgraph(tmp_path, [(100, "80.0"), (200, "20.0"), (300, "55")])

        data = load_methylation_bedgraph(path)

        assert [(s.position, s.methylation_level) for s in data["chr1"]] == pytest.approx(
            [(100, 0.8), (200, 0.2), (300, 0.55)]
        )

    def test_a_fraction_file_is_left_alone(self, tmp_path: Path) -> None:
        path = _bedgraph(tmp_path, [(100, "0.8"), (200, "0.2"), (300, "1.0")])

        data = load_methylation_bedgraph(path)

        assert [s.methylation_level for s in data["chr1"]] == pytest.approx([0.8, 0.2, 1.0])

    def test_values_outside_the_scale_are_skipped_not_fatal(self, tmp_path: Path) -> None:
        # Percent scale (80 > 1): 150 and -3 are outside 0-100; the other sites still load
        path = _bedgraph(tmp_path, [(100, "80"), (200, "150"), (300, "-3"), (400, "40")])

        data = load_methylation_bedgraph(path)

        assert [s.position for s in data["chr1"]] == [100, 400]

    def test_a_negative_fraction_is_skipped_not_fatal(self, tmp_path: Path) -> None:
        path = _bedgraph(tmp_path, [(100, "0.8"), (200, "-0.2"), (300, "0.5")])

        data = load_methylation_bedgraph(path)

        assert [s.position for s in data["chr1"]] == [100, 300]

    def test_a_non_numeric_level_is_skipped(self, tmp_path: Path) -> None:
        path = _bedgraph(tmp_path, [(100, "0.8"), (200, "nan"), (300, "0.5")])

        data = load_methylation_bedgraph(path)

        assert [s.position for s in data["chr1"]] == [100, 300]

    def test_fraction_scale_can_be_forced(self, tmp_path: Path) -> None:
        # A fraction file with one stray 1.5: "auto" would call the whole file percentages (every level 100 times too
        # small), "fraction" skips only the stray value
        path = _bedgraph(tmp_path, [(100, "0.8"), (200, "1.5"), (300, "0.5")])

        forced = load_methylation_bedgraph(path, level_scale="fraction")
        automatic = load_methylation_bedgraph(path)

        assert [(s.position, s.methylation_level) for s in forced["chr1"]] == pytest.approx([(100, 0.8), (300, 0.5)])
        assert [s.methylation_level for s in automatic["chr1"]] == pytest.approx([0.01, 0.02, 0.0])

    def test_percent_scale_can_be_forced(self, tmp_path: Path) -> None:
        # Levels of exactly 1 and 0 are fractions to "auto" but 1% and 0% when the file is known to be percentages
        path = _bedgraph(tmp_path, [(100, "1.0"), (200, "0.0")])

        automatic = load_methylation_bedgraph(path)
        forced = load_methylation_bedgraph(path, level_scale="percent")

        assert [s.methylation_level for s in automatic["chr1"]] == pytest.approx([1.0, 0.0])
        assert [s.methylation_level for s in forced["chr1"]] == pytest.approx([0.01, 0.0])

    def test_an_unknown_scale_is_rejected(self, tmp_path: Path) -> None:
        path = _bedgraph(tmp_path, [(100, "0.5")])

        with pytest.raises(ValueError, match="level_scale"):
            load_methylation_bedgraph(path, level_scale="ratio")
