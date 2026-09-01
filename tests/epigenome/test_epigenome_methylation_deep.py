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

from pathlib import Path

import pytest

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
    """Bismark .cov: chrom, start, strand, percent methylated, count methylated, count total."""
    content = (
        "# Bismark coverage\n"
        "chr1\t100\t+\t80.0\t8\t10\n"
        "chr1\t200\t-\t20.0\t2\t10\n"
        "chr1\t300\t+\t50.0\t5\t10\n"
        "chr2\t50\t+\t90.0\t9\t10\n"
        "short\tline\n"
        "chr2\t150\t-\t10.0\t1\t10\n"
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

    def test_cov_parses_percent_and_strand(self, cov_file: Path) -> None:
        data = load_methylation_cov(cov_file)
        assert set(data.keys()) == {"chr1", "chr2"}
        first = data["chr1"][0]
        assert first.position == 100
        assert first.strand == "+"
        assert first.methylated_reads == 8
        assert first.total_reads == 10
        assert first.methylation_level == pytest.approx(0.8)
        # 'short\tline' skipped; chr2 has 2 valid sites.
        assert len(data["chr2"]) == 2

    def test_cov_min_coverage_filter(self, cov_file: Path) -> None:
        data = load_methylation_cov(cov_file, min_coverage=20)
        assert data == {}


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


class TestCpGIslands:
    def test_high_methylation_windows_flagged(self) -> None:
        # 20 dense sites at 0.9 methylation spanning 400 bp -> windows qualify.
        sites = {"chr1": [_site("chr1", 1000 + i * 20, 9, 10) for i in range(20)]}
        islands = identify_cpg_islands(sites, window_size=200, step_size=50)
        assert len(islands) >= 1
        first = islands[0]
        assert first["chromosome"] == "chr1"
        assert 0.0 <= first["avg_methylation"] <= 1.0
        assert first["cpg_ratio"] >= 0.6

    def test_low_methylation_windows_excluded(self) -> None:
        sites = {"chr1": [_site("chr1", 1000 + i * 20, 0, 10) for i in range(20)]}
        islands = identify_cpg_islands(sites, window_size=200, step_size=50)
        assert islands == []

    def test_too_few_sites_skips_chromosome(self) -> None:
        sites = {"chr1": [_site("chr1", 1000 + i * 20, 9, 10) for i in range(5)]}
        assert identify_cpg_islands(sites) == []

    def test_overlapping_windows_deduplicated(self) -> None:
        # One long dense block: overlapping qualifying windows must collapse.
        sites = {"chr1": [_site("chr1", 1000 + i * 10, 9, 10) for i in range(40)]}
        islands = identify_cpg_islands(sites, window_size=100, step_size=25)
        # No two reported islands on the same chromosome may overlap.
        for i in range(len(islands)):
            for j in range(i + 1, len(islands)):
                a, b = islands[i], islands[j]
                assert a["end"] <= b["start"] or b["end"] <= a["start"]


class TestMethylationEntropy:
    def test_uniform_distribution_maximizes_entropy(self) -> None:
        # Sites spread evenly across bins [0,0.2,...,1.0); 25 bp spacing so the
        # 20-site block spans > window_size (200 bp).
        levels = [0.0, 0.2, 0.4, 0.6, 0.8]
        sites = [_site("chr1", 1000 + i * 25, int(l * 100), 100) for i, l in enumerate(levels * 4)]
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
