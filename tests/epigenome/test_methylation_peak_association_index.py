"""Methylation-versus-peak association: an interval index instead of a scan of every peak for every site.

Defect: ``_analyze_methylation_chip_associations`` and ``_analyze_methylation_atac_associations`` compared every
methylation site with every peak, ``O(sites x peaks)``: 20,000 sites against 5,000 peaks took 4.9 s, and a
whole-genome bisulfite sample (tens of millions of sites) against a peak set (tens of thousands of peaks) would take
the better part of a day. The peaks are now merged once into sorted disjoint intervals per chromosome and each site is
placed with two binary searches (0.01 s for the same 20,000 sites).

The counts must not change, so the tests compare the module with the scan it replaced, written out here as the
specification (a site is inside a peak when ``start <= position <= end``, ends included; its distance is to the
nearest peak of its chromosome), on random data with nested and overlapping peaks, zero-length peaks, peaks given as a
single position, alias keys, chromosomes without peaks and sites without a location.

The old scan also crashed on a peak given as a single position at 0 (``start + 1 if start else None`` is ``None`` for
a start of 0, and the comparison with it raised ``TypeError``, so the whole analysis was lost), and silently skipped
sites and peaks it could not place.

REAL IMPLEMENTATION.
"""

from __future__ import annotations

import json
import os
import random
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pytest

import metainformant
from metainformant.epigenome.workflow.workflow import (
    EpigenomeConfig,
    _analyze_methylation_atac_associations,
    _analyze_methylation_chip_associations,
    integrate_epigenome_results,
)

Site = Dict[str, Any]
Peak = Dict[str, Any]


# ---------------------------------------------------------------------------
# The scan over all peaks, as the specification
# ---------------------------------------------------------------------------


def _scan_distance(chromosome: Any, position: Any, peaks: List[Peak]) -> Optional[Any]:
    """Distance from a position to the nearest peak of its chromosome by looking at every peak (``None``: no peak)."""
    nearest = None
    for peak in peaks:
        if peak.get("chrom", peak.get("chromosome")) != chromosome:
            continue
        start = peak.get("start", peak.get("position"))
        if start is None:
            continue
        end = peak.get("end")
        if end is None:
            end = start + 1  # a single position
        if end < start:
            continue  # unplaceable
        if start <= position <= end:
            distance = 0
        elif position < start:
            distance = start - position
        else:
            distance = position - end
        nearest = distance if nearest is None else min(nearest, distance)
    return nearest


def _scan_chip(sites: List[Site], peaks: List[Peak], max_distance: int) -> Dict[str, Any]:
    overlapping = proximal = distal = 0
    for site in sites:
        chromosome = site.get("chrom", site.get("chromosome"))
        position = site.get("position", site.get("pos"))
        if chromosome is None or position is None:
            continue
        distance = _scan_distance(chromosome, position, peaks)
        if distance is None:
            continue
        if distance == 0:
            overlapping += 1
        elif distance <= max_distance:
            proximal += 1
        else:
            distal += 1
    total = len(sites)
    return {
        "total_methylation_sites": total,
        "total_chip_peaks": len(peaks),
        "overlapping": overlapping,
        "proximal": proximal,
        "distal": distal,
        "overlap_rate": overlapping / total,
        "proximal_rate": proximal / total,
    }


def _scan_atac(sites: List[Site], peaks: List[Peak], max_distance: int) -> Dict[str, Any]:
    levels: Dict[str, List[float]] = {"accessible": [], "near": [], "closed": []}
    counts = {"accessible": 0, "near": 0, "closed": 0}
    for site in sites:
        chromosome = site.get("chrom", site.get("chromosome"))
        position = site.get("position", site.get("pos"))
        level = site.get("methylation_level", site.get("beta", None))
        if chromosome is None or position is None:
            continue
        distance = _scan_distance(chromosome, position, peaks)
        if distance == 0:
            group = "accessible"
        elif distance is not None and distance <= max_distance:
            group = "near"
        else:
            group = "closed"
        counts[group] += 1
        if level is not None:
            levels[group].append(level)
    total = len(sites)
    statistics: Dict[str, Any] = {
        "total_methylation_sites": total,
        "total_accessible_regions": len(peaks),
        "sites_in_accessible": counts["accessible"],
        "sites_near_accessible": counts["near"],
        "sites_in_closed": counts["closed"],
        "accessible_rate": (counts["accessible"] + counts["near"]) / total,
    }
    for group, key in (
        ("accessible", "mean_methylation_accessible"),
        ("near", "mean_methylation_near_accessible"),
        ("closed", "mean_methylation_closed"),
    ):
        if levels[group]:
            statistics[key] = sum(levels[group]) / len(levels[group])
    return statistics


def _random_case(rng: random.Random, n_sites: int, n_peaks: int, span: int) -> Tuple[List[Site], List[Peak]]:
    chromosomes = ["chr1", "chr2", "chr3"][: rng.randint(1, 3)]
    peaks: List[Peak] = []
    for _ in range(n_peaks):
        chromosome = rng.choice(chromosomes)
        start = rng.randint(0, span)
        length = rng.choice([0, 1, 5, 50, 300, 2000])  # zero-length peaks, short ones, and long ones that nest others
        roll = rng.random()
        if roll < 0.1:
            peaks.append({"chromosome": chromosome, "position": start})  # a single position, start 0 included
        elif roll < 0.15:
            peaks.append({"chrom": chromosome, "start": start, "end": start + length, "score": 5})
        else:
            peaks.append({"chrom": chromosome, "start": start, "end": start + length})
    sites: List[Site] = []
    for _ in range(n_sites):
        chromosome = rng.choice(chromosomes + ["chrZ"])  # chrZ has no peaks
        position = rng.randint(0, span + 3000)
        roll = rng.random()
        if roll < 0.1:
            sites.append({"chromosome": chromosome, "pos": position, "beta": rng.random()})
        elif roll < 0.13:
            sites.append({"chrom": None, "position": position})
        elif roll < 0.16:
            sites.append({"chrom": chromosome})
        else:
            sites.append({"chrom": chromosome, "position": position, "methylation_level": rng.random()})
    return sites, peaks


class TestSameCountsAsScanningEveryPeak:
    @pytest.mark.parametrize("seed", range(60))
    def test_chip_seq_association(self, seed: int) -> None:
        rng = random.Random(seed)
        sites, peaks = _random_case(rng, rng.randint(1, 80), rng.randint(1, 30), rng.choice([40, 500, 5000]))
        max_distance = rng.choice([0, 1, 10, 100, 1000])

        result = _analyze_methylation_chip_associations(
            {"sites": sites},
            {"peaks": peaks},
            EpigenomeConfig(max_distance_for_integration=max_distance),
        )

        assert result["statistics"] == _scan_chip(sites, peaks, max_distance)
        findings = {f["type"]: f["count"] for f in result["findings"] if f["type"] in ("overlap", "proximal")}
        expected = result["statistics"]
        assert findings == {
            key: expected[name]
            for key, name in (("overlap", "overlapping"), ("proximal", "proximal"))
            if expected[name]
        }

    @pytest.mark.parametrize("seed", range(60))
    def test_atac_seq_association(self, seed: int) -> None:
        rng = random.Random(500 + seed)
        sites, peaks = _random_case(rng, rng.randint(1, 80), rng.randint(1, 30), rng.choice([40, 500, 5000]))
        max_distance = rng.choice([0, 1, 10, 100, 1000])

        result = _analyze_methylation_atac_associations(
            {"sites": sites},
            {"peaks": peaks},
            EpigenomeConfig(max_distance_for_integration=max_distance),
        )

        assert result["statistics"] == _scan_atac(sites, peaks, max_distance)

    def test_closed_intervals_include_both_ends_and_distances_count_from_the_nearest_end(self) -> None:
        sites = [{"chrom": "chr1", "position": p} for p in (99, 100, 150, 200, 201, 260, 261, 1300)]
        peaks = [{"chrom": "chr1", "start": 100, "end": 200}]

        stats = _analyze_methylation_chip_associations(
            {"sites": sites}, {"peaks": peaks}, EpigenomeConfig(max_distance_for_integration=60)
        )["statistics"]

        # inside: 100, 150, 200; within 60 of an end: 99, 201, 260; farther: 261, 1300
        assert (stats["overlapping"], stats["proximal"], stats["distal"]) == (3, 3, 2)

    def test_overlapping_and_nested_peaks_act_as_one_region(self) -> None:
        sites = [{"chrom": "chr1", "position": p} for p in (10, 55, 95, 205, 300)]
        peaks = [{"chrom": "chr1", "start": s, "end": e} for s, e in ((50, 100), (0, 20), (60, 70), (90, 120))]

        stats = _analyze_methylation_chip_associations(
            {"sites": sites}, {"peaks": peaks}, EpigenomeConfig(max_distance_for_integration=100)
        )["statistics"]

        assert (stats["overlapping"], stats["proximal"], stats["distal"]) == (3, 1, 1)

    def test_sites_on_a_chromosome_without_peaks_keep_their_old_places(self) -> None:
        sites = [{"chrom": "chrZ", "position": 5, "methylation_level": 0.25}, {"chrom": "chr1", "position": 5}]
        peaks = [{"chrom": "chr1", "start": 0, "end": 10}]

        chip = _analyze_methylation_chip_associations({"sites": sites}, {"peaks": peaks}, EpigenomeConfig())
        atac = _analyze_methylation_atac_associations({"sites": sites}, {"peaks": peaks}, EpigenomeConfig())

        # ChIP-seq: in none of the three categories; ATAC-seq: in closed chromatin
        assert (chip["statistics"]["overlapping"], chip["statistics"]["proximal"], chip["statistics"]["distal"]) == (
            1,
            0,
            0,
        )
        assert atac["statistics"]["sites_in_closed"] == 1
        assert atac["statistics"]["mean_methylation_closed"] == 0.25


class TestPeaksThatUsedToBreakTheScan:
    def test_a_single_position_peak_at_zero_spans_zero_and_one(self) -> None:
        sites = [{"chrom": "chr1", "position": p} for p in (0, 1, 5)]
        peaks = [{"chrom": "chr1", "start": 0}]

        stats = _analyze_methylation_chip_associations(
            {"sites": sites}, {"peaks": peaks}, EpigenomeConfig(max_distance_for_integration=10)
        )["statistics"]

        assert (stats["overlapping"], stats["proximal"], stats["distal"]) == (2, 1, 0)

    def test_the_integration_no_longer_loses_the_analysis_to_such_a_peak(self, tmp_path: Path) -> None:
        methylation = {"statistics": {"total_samples": 1}, "sites": [{"chrom": "chr1", "position": 0}]}
        peaks = {"statistics": {"total_samples": 1}, "peaks": [{"chrom": "chr1", "position": 0}]}

        result = integrate_epigenome_results(methylation, peaks, peaks, tmp_path, EpigenomeConfig())

        assert result["errors"] == []
        assert result["integrated_analyses"]["methylation_chipseq"]["statistics"]["overlapping"] == 1
        assert result["integrated_analyses"]["methylation_atacseq"]["statistics"]["sites_in_accessible"] == 1

    def test_a_peak_with_an_explicit_empty_end_is_a_single_position(self) -> None:
        sites = [{"chrom": "chr1", "position": 7}]
        peaks = [{"chrom": "chr1", "start": 7, "end": None}]

        stats = _analyze_methylation_chip_associations({"sites": sites}, {"peaks": peaks}, EpigenomeConfig())[
            "statistics"
        ]

        assert stats["overlapping"] == 1

    def test_peaks_that_cannot_be_placed_are_left_out_and_reported(self) -> None:
        sites = [{"chrom": "chr1", "position": 5}, {"chrom": "chr1", "position": 500}]
        peaks = [
            {"chrom": "chr1", "start": 0, "end": 10},
            {"chrom": "chr1", "start": 400, "end": 100},  # ends before it starts
            {"chrom": "chr1", "start": "x", "end": 100},  # not a number
            {"chrom": None, "start": 0, "end": 1000},  # no chromosome
            {"start": 0, "end": 1000},  # no chromosome
            {"chrom": "chr1", "start": float("nan"), "end": 5},
        ]

        result = _analyze_methylation_chip_associations(
            {"sites": sites}, {"peaks": peaks}, EpigenomeConfig(max_distance_for_integration=100)
        )

        stats = result["statistics"]
        assert (stats["overlapping"], stats["proximal"], stats["distal"]) == (1, 0, 1)
        assert stats["total_chip_peaks"] == 6
        (warning,) = [f for f in result["findings"] if f["type"] == "warning"]
        assert warning["count"] == 5 and "5 of 6" in warning["description"]

    def test_sites_without_a_location_are_reported_not_silently_skipped(self) -> None:
        sites = [{"chrom": "chr1", "position": 5}, {"chrom": None, "position": 5}, {"position": 5}, {"chrom": "chr1"}]
        peaks = [{"chrom": "chr1", "start": 0, "end": 10}]

        for analyse in (_analyze_methylation_chip_associations, _analyze_methylation_atac_associations):
            result = analyse({"sites": sites}, {"peaks": peaks}, EpigenomeConfig())

            assert result["statistics"]["total_methylation_sites"] == 4
            (warning,) = [f for f in result["findings"] if f["type"] == "warning"]
            assert warning["count"] == 3

    def test_well_formed_input_draws_no_warning(self) -> None:
        sites = [{"chrom": "chr1", "position": 5}]
        peaks = [{"chrom": "chr1", "start": 0, "end": 10}]

        for analyse in (_analyze_methylation_chip_associations, _analyze_methylation_atac_associations):
            result = analyse({"sites": sites}, {"peaks": peaks}, EpigenomeConfig())

            assert all(f["type"] != "warning" for f in result["findings"])


_SCALE_SCENARIO = """
import json
import random
import sys
import time
from bisect import bisect_right
from pathlib import Path

import metainformant
from metainformant.epigenome.workflow.workflow import EpigenomeConfig, integrate_epigenome_results

rng = random.Random(1)
n_sites, n_peaks, span = 150_000, 40_000, 50_000_000
chromosomes = [f"chr{i}" for i in range(1, 6)]
sites = [
    {"chrom": rng.choice(chromosomes), "position": rng.randrange(span), "methylation_level": rng.random()}
    for _ in range(n_sites)
]
peaks = []
for _ in range(n_peaks):
    start = rng.randrange(span)
    peaks.append({"chrom": rng.choice(chromosomes), "start": start, "end": start + rng.randrange(100, 2000)})

started = time.perf_counter()
result = integrate_epigenome_results(
    {"statistics": {"total_samples": 1}, "sites": sites},
    {"statistics": {"total_samples": 1}, "peaks": peaks},
    {"statistics": {"total_samples": 1}, "peaks": peaks},
    Path(sys.argv[1]),
    EpigenomeConfig(generate_reports=False),
)
elapsed = time.perf_counter() - started

# An independent count of the sites inside a peak: merge the peaks of each chromosome and bisect
merged = {}
for chrom in chromosomes:
    starts, ends = [], []
    for start, end in sorted((p["start"], p["end"]) for p in peaks if p["chrom"] == chrom):
        if ends and start <= ends[-1]:
            ends[-1] = max(ends[-1], end)
        else:
            starts.append(start)
            ends.append(end)
    merged[chrom] = (starts, ends)
inside = 0
for site in sites:
    starts, ends = merged[site["chrom"]]
    k = bisect_right(starts, site["position"]) - 1
    inside += k >= 0 and site["position"] <= ends[k]

chip = result["integrated_analyses"]["methylation_chipseq"]["statistics"]
atac = result["integrated_analyses"]["methylation_atacseq"]["statistics"]
print(
    "RESULT"
    + json.dumps(
        {
            "package": str(Path(metainformant.__file__).resolve().parent),
            "errors": result["errors"],
            "elapsed": elapsed,
            "inside": inside,
            "chip": chip,
            "atac": atac,
        }
    )
)
"""


class TestScale:
    @pytest.mark.slow
    def test_a_genome_scale_sample_is_integrated_in_seconds(self, tmp_path: Path) -> None:
        # 150,000 sites against 40,000 peaks is 6e9 pairs for the scan (about five minutes); the index takes a second
        package_root = Path(metainformant.__file__).resolve().parents[1]
        env = {
            "PATH": os.environ.get("PATH", ""),
            "HOME": os.environ.get("HOME", ""),
            "TMPDIR": os.environ.get("TMPDIR", tempfile.gettempdir()),
            "PYTHONPATH": str(package_root),  # the child must import the package this test is checking
            "PYTHONDONTWRITEBYTECODE": "1",
        }
        try:
            completed = subprocess.run(
                [sys.executable, "-c", _SCALE_SCENARIO, str(tmp_path)],
                capture_output=True,
                text=True,
                check=True,
                env=env,
                timeout=60,
            )
        except subprocess.TimeoutExpired:
            pytest.fail("integrating 150,000 sites with 40,000 peaks did not finish in 60 s: it is not indexed")

        line = next(line for line in completed.stdout.splitlines() if line.startswith("RESULT"))
        report = json.loads(line[len("RESULT") :])

        assert Path(report["package"]) == package_root / "metainformant"
        assert report["errors"] == []
        assert report["elapsed"] < 60
        chip, atac = report["chip"], report["atac"]
        assert chip["overlapping"] == report["inside"] == atac["sites_in_accessible"]
        assert chip["overlapping"] + chip["proximal"] + chip["distal"] == chip["total_methylation_sites"] == 150_000
        assert atac["sites_in_accessible"] + atac["sites_near_accessible"] + atac["sites_in_closed"] == 150_000
