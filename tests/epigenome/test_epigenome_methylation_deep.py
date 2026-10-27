"""Tests for untested epigenome methylation loaders, DMR, CpG islands, entropy, export.

Covers:
    - load_methylation_bedgraph / export_methylation_bedgraph (round trip)
    - load_methylation_cov (Bismark format)
    - find_differentially_methylated_regions
    - identify_cpg_islands
    - calculate_methylation_entropy
    - generate_methylation_report

Uses real files in tmp_path and real implementations only
(real-implementation policy).
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from metainformant.core.utils.errors import ValidationError
from metainformant.epigenome.assays.methylation import (
    MethylationSite,
    calculate_methylation_entropy,
    export_methylation_bedgraph,
    find_differentially_methylated_regions,
    generate_methylation_report,
    identify_cpg_islands,
    load_methylation_bedgraph,
    load_methylation_cov,
)


def _site(chrom: str, pos: int, methylated: int, total: int) -> MethylationSite:
    return MethylationSite(chromosome=chrom, position=pos, methylated_reads=methylated, total_reads=total)


@pytest.fixture
def bedgraph_file(tmp_path: Path) -> Path:
    """Small BEDgraph file with two chromosomes, one malformed line, one skipped line."""
    content = (
        'track type=bedGraph name="demo"\n'
        "# a comment line\n"
        "chr1\t100\t101\t0.8\n"
        "chr1\t200\t201\t0.2\n"
        "chr1\t300\t301\t0.5\n"
        "chr2\t50\t51\t0.9\n"
        "chr2\t150\t151\t0.1\n"
        "badline\n"
        "chr2\t250\t252\t0.6\n"
    )
    p = tmp_path / "demo.bedgraph"
    p.write_text(content)
    return p


@pytest.fixture
def cov_file(tmp_path: Path) -> Path:
    """Bismark .cov: chrom, start, end, percent methylated, count methylated, count unmethylated."""
    content = (
        "# a comment line\n"
        "chr1\t100\t100\t80.0\t8\t2\n"
        "chr1\t200\t200\t20.0\t2\t8\n"
        "chr1\t300\t300\t50.0\t5\t5\n"
        "chr2\t50\t50\t90.0\t9\t1\n"
        "short\tline\n"
        "chr2\t150\t150\t10.0\t1\t9\n"
    )
    p = tmp_path / "demo.cov"
    p.write_text(content)
    return p


class TestLoaders:
    def test_bedgraph_parses_sites_sorted(self, bedgraph_file: Path) -> None:
        data = load_methylation_bedgraph(bedgraph_file)
        assert set(data.keys()) == {"chr1", "chr2"}
        assert [s.position for s in data["chr1"]] == [100, 200, 300]
        # Malformed 'badline' skipped without aborting the load.
        assert len(data["chr2"]) == 3

    def test_bedgraph_level_consistency(self, bedgraph_file: Path) -> None:
        data = load_methylation_bedgraph(bedgraph_file)
        site = data["chr1"][0]
        assert site.methylation_level == pytest.approx(0.8)

    def test_bedgraph_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(Exception):
            load_methylation_bedgraph(tmp_path / "missing.bedgraph")

    def test_cov_parses_bismark_columns(self, cov_file: Path) -> None:
        data = load_methylation_cov(cov_file)
        assert set(data.keys()) == {"chr1", "chr2"}
        first = data["chr1"][0]
        assert first.position == 100
        assert first.methylated_reads == 8
        assert first.total_reads == 10  # 8 methylated + 2 unmethylated
        assert first.methylation_level == pytest.approx(0.8)
        # 'short\tline' skipped; chr2 has 2 valid sites.
        assert len(data["chr2"]) == 2

    def test_cov_min_coverage_filter(self, cov_file: Path) -> None:
        data = load_methylation_cov(cov_file, min_coverage=20)
        assert data == {}


class TestBismarkCovLayout:
    """``load_methylation_cov`` must read the real Bismark layout, not an invented one.

    A Bismark coverage file is tab-separated with no header and the columns
    ``chromosome, start, end, methylation percentage, count methylated, count unmethylated``. The loader used to
    read column 3 as a strand and column 6 as the *total*, so ``chr1 100 100 75.0 15 5`` aborted the whole load
    ("total_reads must be >= 15, got 5") and a minority-methylated site such as ``chr1 200 200 25.0 5 15`` was
    silently read as 15 reads instead of 20.
    """

    @staticmethod
    def _write(tmp_path: Path, text: str, name: str = "sample.cov") -> Path:
        path = tmp_path / name
        path.write_text(text)
        return path

    def test_total_is_methylated_plus_unmethylated(self, tmp_path: Path) -> None:
        data = load_methylation_cov(self._write(tmp_path, "chr1\t100\t100\t75.0\t15\t5\n"))

        (site,) = data["chr1"]
        assert site.position == 100
        assert site.methylated_reads == 15
        assert site.total_reads == 20
        assert site.coverage == 20
        assert site.methylation_level == pytest.approx(0.75)

    def test_minority_methylated_site_keeps_its_unmethylated_reads(self, tmp_path: Path) -> None:
        data = load_methylation_cov(self._write(tmp_path, "chr1\t200\t200\t25.0\t5\t15\n"))

        (site,) = data["chr1"]
        assert (site.methylated_reads, site.total_reads) == (5, 20)
        assert site.methylation_level == pytest.approx(0.25)

    def test_reads_a_multi_chromosome_file_in_bismark_form(self, tmp_path: Path) -> None:
        # real-looking records: per-cytosine lines (end == start), a merged CpG (end == start + 1), a fractional percent
        text = (
            "chr1\t10470\t10470\t66.666667\t2\t1\n"
            "chr1\t10468\t10468\t100\t3\t0\n"
            "chr1\t10471\t10472\t0\t0\t7\n"
            "chr2\t5\t5\t50\t1\t1\n"
        )
        data = load_methylation_cov(self._write(tmp_path, text))

        assert sorted(data) == ["chr1", "chr2"]
        # sites are sorted by position within a chromosome
        assert [(s.position, s.methylated_reads, s.total_reads) for s in data["chr1"]] == [
            (10468, 3, 3),
            (10470, 2, 3),
            (10471, 0, 7),
        ]
        assert [s.methylation_level for s in data["chr1"]] == pytest.approx([1.0, 2 / 3, 0.0])
        assert [(s.position, s.methylated_reads, s.total_reads) for s in data["chr2"]] == [(5, 1, 2)]

    def test_gzipped_bismark_file(self, tmp_path: Path) -> None:
        path = tmp_path / "sample.bismark.cov.gz"
        with gzip.open(path, "wt") as handle:
            handle.write("chr1\t100\t100\t75.0\t15\t5\nchr1\t200\t200\t25.0\t5\t15\n")

        data = load_methylation_cov(path)

        assert [(s.position, s.total_reads) for s in data["chr1"]] == [(100, 20), (200, 20)]

    def test_min_coverage_counts_methylated_and_unmethylated_reads(self, tmp_path: Path) -> None:
        path = self._write(tmp_path, "chr1\t100\t100\t75.0\t15\t5\n")

        assert load_methylation_cov(path, min_coverage=20)["chr1"][0].total_reads == 20
        assert load_methylation_cov(path, min_coverage=21) == {}

    def test_a_coverage_file_has_no_strand(self, tmp_path: Path) -> None:
        # the strand used to be read from column 3, which in this format is the end coordinate
        data = load_methylation_cov(self._write(tmp_path, "chr1\t100\t100\t75.0\t15\t5\n"))

        assert data["chr1"][0].strand == MethylationSite("chr1", 100, 15, 20).strand

    def test_older_strand_layout_is_skipped_rather_than_misread(self, tmp_path: Path) -> None:
        # chromosome, position, strand, percentage, methylated, total: the third column is not an end coordinate
        data = load_methylation_cov(self._write(tmp_path, "chr1\t100\t+\t80.0\t8\t10\nchr1\t200\t-\t20.0\t2\t10\n"))

        assert data == {}

    @pytest.mark.parametrize(
        "bad_line",
        [
            "chr1\t100\t100\tabc\t15\t5",  # percentage is not a number
            "chr1\t100\t100\t75.0\t15.5\t5",  # counts are integers
            "chr1\t100\t100\t75.0\t-1\t5",  # negative methylated count
            "chr1\t100\t100\t75.0\t15\t-5",  # negative unmethylated count
            "chr1\t100\t99\t75.0\t15\t5",  # end before start
            "chr1\t-4\t-4\t75.0\t15\t5",  # negative position
            "chr1\t100\t100\t75.0\t15",  # a column short
        ],
    )
    def test_malformed_lines_are_skipped_and_the_rest_still_loads(self, tmp_path: Path, bad_line: str) -> None:
        data = load_methylation_cov(self._write(tmp_path, f"{bad_line}\nchr1\t300\t300\t50.0\t5\t5\n"))

        assert [(s.position, s.total_reads) for s in data["chr1"]] == [(300, 10)]


class TestBedgraphRoundTrip:
    def test_export_then_reload_preserves_positions_and_levels(self, tmp_path: Path, bedgraph_file: Path) -> None:
        original = load_methylation_bedgraph(bedgraph_file)
        out = tmp_path / "roundtrip.bedgraph"
        export_methylation_bedgraph(original, out)
        assert out.exists()

        text = out.read_text()
        assert text.startswith("track type=bedGraph")

        reloaded = load_methylation_bedgraph(out)
        for chrom, sites in original.items():
            assert chrom in reloaded
            assert [s.position for s in reloaded[chrom]] == [s.position for s in sites]
            for orig_site, new_site in zip(sites, reloaded[chrom]):
                # Export writes levels; loader derives reads from levels, so the
                # level must round-trip while read counts may differ.
                assert new_site.methylation_level == pytest.approx(orig_site.methylation_level, abs=0.05)


def _uniform_region(chrom: str, start: int, count: int, level: float) -> list[MethylationSite]:
    return [_site(chrom, start + i * 10, int(level * 10), 10) for i in range(count)]


class TestFindDMRs:
    def test_hyper_and_hypo_dmrs_detected(self) -> None:
        # Direction convention: relative to condition 1 (site_data1).
        # chr1: cond1 0.2 vs cond2 0.8 -> cond1 hypomethylated ("hypo").
        # chr2: cond1 0.8 vs cond2 0.2 -> cond1 hypermethylated ("hyper").
        control = {
            "chr1": _uniform_region("chr1", 1000, 4, 0.2),
            "chr2": _uniform_region("chr2", 2000, 4, 0.8),
        }
        treatment = {
            "chr1": _uniform_region("chr1", 1000, 4, 0.8),
            "chr2": _uniform_region("chr2", 2000, 4, 0.2),
        }
        dmrs = find_differentially_methylated_regions(control, treatment, min_sites=3)
        assert len(dmrs) == 2
        by_chrom = {d["chromosome"]: d for d in dmrs}
        assert by_chrom["chr1"]["direction"] == "hypo"
        assert by_chrom["chr2"]["direction"] == "hyper"
        for dmr in dmrs:
            assert len(dmr["sites"]) == 4
            assert dmr["mean_delta"] == pytest.approx(0.6)

    def test_no_difference_no_dmrs(self) -> None:
        sites = {"chr1": _uniform_region("chr1", 1000, 5, 0.5)}
        assert find_differentially_methylated_regions(sites, {"chr1": _uniform_region("chr1", 1000, 5, 0.5)}) == []

    def test_disjoint_chromosomes_yield_nothing(self) -> None:
        a = {"chr1": _uniform_region("chr1", 1000, 4, 0.2)}
        b = {"chrX": _uniform_region("chrX", 1000, 4, 0.8)}
        assert find_differentially_methylated_regions(a, b) == []

    def test_min_sites_gate(self) -> None:
        control = {"chr1": _uniform_region("chr1", 1000, 4, 0.1)}
        treatment = {"chr1": _uniform_region("chr1", 1000, 4, 0.9)}
        # Requiring 5+ sites filters the 4-site region out.
        assert find_differentially_methylated_regions(control, treatment, min_sites=5) == []


def _genome_with_island() -> dict[str, str]:
    """chr1 = 500 bp of AT, a 400 bp CpG island (CG x 200), 500 bp of AT; chr3 is AT only.

    The island occupies [500, 900). In it G + C = 100%, the CpG count is 200 and C = G = 200, so the observed/expected
    CpG ratio is 200 * 400 / (200 * 200) = 2.0.
    """
    at_flank = "AT" * 250
    return {"chr1": at_flank + "CG" * 200 + at_flank, "chr3": "ATATATAT" * 100}


class TestCpGIslands:
    """CpG islands are defined by the sequence (Gardiner-Garden & Frommer), not by methylation levels."""

    def test_island_is_called_from_the_sequence(self) -> None:
        islands = identify_cpg_islands({}, step_size=1, sequences=_genome_with_island())

        assert [island["chromosome"] for island in islands] == ["chr1"]  # the AT-only chr3 has none
        (island,) = islands
        # The stretch of qualifying windows contains the planted island and reaches less than a window beyond it
        assert island["start"] <= 500 and island["end"] >= 900
        assert island["start"] >= 500 - 200 and island["end"] <= 900 + 200
        assert island["length"] == island["end"] - island["start"] >= 200
        assert island["gc_content"] >= 0.5
        assert island["cpg_ratio"] >= 0.6
        assert island["cpg_count"] == 200

    def test_methylation_levels_cannot_create_islands(self) -> None:
        # Dense, highly methylated sites over a sequence with no CpG at all used to be reported as an island
        sites = {"chr1": [_site("chr1", 1000 + i * 20, 9, 10) for i in range(20)]}
        sequences = {"chr1": "AT" * 1000}

        assert identify_cpg_islands(sites, window_size=200, step_size=50, sequences=sequences) == []

    def test_sequence_is_required(self) -> None:
        sites = {"chr1": [_site("chr1", 1000 + i * 20, 9, 10) for i in range(20)]}
        with pytest.raises(ValueError, match="sequence"):
            identify_cpg_islands(sites, window_size=200, step_size=50)

    def test_gc_rich_but_cpg_depleted_sequence_is_not_an_island(self) -> None:
        # G + C is 89% but no C is ever followed by G: observed/expected CpG is 0
        assert identify_cpg_islands({}, sequences={"chr1": "GGGGCCCCA" * 150}) == []

    def test_cpg_dense_but_gc_poor_sequence_is_not_an_island(self) -> None:
        # Every CpG is there, but G + C is only 20%
        assert identify_cpg_islands({}, sequences={"chr1": "CGAAAAAAAA" * 100}) == []

    def test_min_length_is_enforced(self) -> None:
        # The planted island is called as a 600 bp stretch (see above); 700 bp is not reached
        assert identify_cpg_islands({}, step_size=1, min_length=700, sequences=_genome_with_island()) == []

    def test_sites_inside_the_island_are_summarised(self) -> None:
        sites = {
            "chr1": [
                _site("chr1", 600, 9, 10),
                _site("chr1", 650, 1, 10),
                _site("chr1", 5, 5, 10),  # outside the island
                _site("chr1", 1900, 5, 10),  # outside the island
            ]
        }

        (island,) = identify_cpg_islands(sites, step_size=1, sequences=_genome_with_island())

        assert island["sites"] == 2
        assert island["avg_methylation"] == pytest.approx(0.5)

    def test_island_without_sites_has_no_average(self) -> None:
        (island,) = identify_cpg_islands(
            {"chr2": [_site("chr2", 600, 9, 10)]}, step_size=1, sequences=_genome_with_island()
        )

        assert island["sites"] == 0
        assert island["avg_methylation"] is None

    def test_overlapping_windows_collapse_into_one_island(self) -> None:
        # One long CpG-dense block: the many overlapping qualifying windows are reported as a single island
        islands = identify_cpg_islands({}, window_size=100, step_size=25, sequences={"chr1": "CG" * 400})

        assert len(islands) == 1
        assert (islands[0]["start"], islands[0]["end"]) == (0, 800)

    def test_bad_window_parameters_are_rejected(self) -> None:
        with pytest.raises(ValueError):
            identify_cpg_islands({}, window_size=0, sequences={"chr1": "CG" * 400})
        with pytest.raises(ValueError):
            identify_cpg_islands({}, step_size=0, sequences={"chr1": "CG" * 400})


class TestMethylationEntropy:
    def test_uniform_distribution_maximizes_entropy(self) -> None:
        # Sites spread evenly across bins [0,0.2,...,1.0); 25 bp spacing so the
        # 20-site block spans > window_size (200 bp).
        levels = [0.0, 0.2, 0.4, 0.6, 0.8]
        sites = [_site("chr1", 1000 + i * 25, int(level * 100), 100) for i, level in enumerate(levels * 4)]
        profile = calculate_methylation_entropy({"chr1": sites}, window_size=200)
        assert "chr1" in profile
        positions, entropies = zip(*profile["chr1"])
        for h in entropies:
            # Max entropy with 5 populated bins is log2(5) ~ 2.32.
            assert 0.0 < h <= 2.322 + 1e-9

    def test_single_bin_minimizes_entropy(self) -> None:
        # All sites in the same bin => entropy near zero.
        sites = [_site("chr1", 1000 + i * 10, 95, 100) for i in range(20)]
        profile = calculate_methylation_entropy({"chr1": sites}, window_size=200)
        for _, h in profile["chr1"]:
            assert h == pytest.approx(0.0, abs=1e-9)

    def test_too_few_sites_skips_chromosome(self) -> None:
        sites = [_site("chr1", 1000 + i * 10, 50, 100) for i in range(3)]
        assert calculate_methylation_entropy({"chr1": sites}) == {}


class TestMethylationReport:
    def test_report_contains_key_sections_and_saves(self, tmp_path: Path) -> None:
        data = {
            "chr1": _uniform_region("chr1", 1000, 6, 0.7),
            "chr2": _uniform_region("chr2", 2000, 6, 0.3),
        }
        out = tmp_path / "report.txt"
        report = generate_methylation_report(data, out)
        assert out.exists()
        assert out.read_text() == report
        assert "DNA METHYLATION ANALYSIS REPORT" in report
        assert "Total Sites" in report
        assert "chr1" in report


class TestFindDMRsMaxGap:
    """max_gap opt-in splits distant significant sites into separate DMRs."""

    def _two_distant_blocks(self) -> tuple[dict, dict]:
        near1 = [_site("chr1", 1000 + i * 20, 2, 10) for i in range(4)]
        far1 = [_site("chr1", 500000 + i * 20, 2, 10) for i in range(4)]
        near2 = [_site("chr1", 1000 + i * 20, 8, 10) for i in range(4)]
        far2 = [_site("chr1", 500000 + i * 20, 8, 10) for i in range(4)]
        return {"chr1": near1 + far1}, {"chr1": near2 + far2}

    def test_default_chains_distant_sites(self) -> None:
        cond1, cond2 = self._two_distant_blocks()
        dmrs = find_differentially_methylated_regions(cond1, cond2, min_sites=3)
        assert len(dmrs) == 1  # historical behaviour

    def test_max_gap_splits_blocks(self) -> None:
        cond1, cond2 = self._two_distant_blocks()
        dmrs = find_differentially_methylated_regions(cond1, cond2, min_sites=3, max_gap=1000)
        assert len(dmrs) == 2
        assert dmrs[0]["end"] < dmrs[1]["start"]


class TestMethylationSiteValidation:
    """MethylationSite rejects structurally invalid data."""

    def test_rejects_negative_position(self) -> None:
        with pytest.raises(ValidationError):
            MethylationSite(chromosome="chr1", position=-1, methylated_reads=1, total_reads=10)

    def test_rejects_total_below_methylated(self) -> None:
        with pytest.raises(ValidationError):
            MethylationSite(chromosome="chr1", position=100, methylated_reads=11, total_reads=10)

    def test_rejects_negative_methylated_reads(self) -> None:
        with pytest.raises(ValidationError):
            MethylationSite(chromosome="chr1", position=100, methylated_reads=-1, total_reads=10)
