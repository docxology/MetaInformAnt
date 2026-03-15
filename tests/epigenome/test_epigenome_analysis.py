"""Tests for epigenome analysis module (tracks).

Tests GenomicTrack, load_bedgraph_track, load_bed_track,
save_bed_track, save_bedgraph_track, merge_tracks,
calculate_track_statistics, extract_track_region, compare_tracks,
and generate_track_report using real implementations. NO MOCKING.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.epigenome.analysis.tracks import (
    GenomicTrack,
    calculate_track_statistics,
    compare_tracks,
    extract_track_region,
    generate_track_report,
    load_bed_track,
    load_bedgraph_track,
    merge_feature_group,
    merge_tracks,
    save_bed_track,
    save_bedgraph_track,
)

# ---------------------------------------------------------------------------
# GenomicTrack
# ---------------------------------------------------------------------------


class TestGenomicTrack:
    """Tests for the GenomicTrack class."""

    def test_creation(self) -> None:
        track = GenomicTrack(name="test", description="test track")
        assert track.name == "test"
        assert track.description == "test track"
        assert track.get_total_features() == 0

    def test_add_feature(self) -> None:
        track = GenomicTrack()
        track.add_feature("chr1", 100, 200, value=5.0)
        assert track.get_total_features() == 1
        features = track.get_features("chr1")
        assert len(features) == 1
        assert features[0]["start"] == 100
        assert features[0]["end"] == 200
        assert features[0]["value"] == 5.0

    def test_get_features_by_region(self) -> None:
        track = GenomicTrack()
        track.add_feature("chr1", 100, 200, value=1.0)
        track.add_feature("chr1", 300, 400, value=2.0)
        track.add_feature("chr1", 500, 600, value=3.0)

        # Should get only the middle feature
        features = track.get_features("chr1", start=250, end=450)
        assert len(features) == 1
        assert features[0]["value"] == 2.0

    def test_get_features_empty_chrom(self) -> None:
        track = GenomicTrack()
        assert track.get_features("chr99") == []

    def test_get_chromosomes(self) -> None:
        track = GenomicTrack()
        track.add_feature("chr1", 0, 100, value=1.0)
        track.add_feature("chr2", 0, 100, value=2.0)
        chroms = track.get_chromosomes()
        assert set(chroms) == {"chr1", "chr2"}

    def test_merge_overlapping_features(self) -> None:
        track = GenomicTrack()
        track.add_feature("chr1", 100, 300, value=5.0)
        track.add_feature("chr1", 200, 400, value=10.0)
        track.add_feature("chr1", 500, 600, value=3.0)

        track.merge_overlapping_features("chr1")
        features = track.get_features("chr1")
        assert len(features) == 2
        assert features[0]["end"] == 400  # Merged
        assert features[0]["value"] == 10.0  # Max value

    def test_normalize_minmax(self) -> None:
        track = GenomicTrack()
        track.add_feature("chr1", 0, 100, value=0.0)
        track.add_feature("chr1", 100, 200, value=10.0)
        track.add_feature("chr1", 200, 300, value=5.0)

        track.normalize_values(method="minmax")
        features = track.get_features("chr1")
        values = [f["value"] for f in features]
        assert min(values) == pytest.approx(0.0)
        assert max(values) == pytest.approx(1.0)
        assert values[2] == pytest.approx(0.5)

    def test_normalize_zscore(self) -> None:
        track = GenomicTrack()
        for i in range(10):
            track.add_feature("chr1", i * 100, (i + 1) * 100, value=float(i))

        track.normalize_values(method="zscore")
        features = track.get_features("chr1")
        values = [f["value"] for f in features]
        # Mean should be approximately 0 after z-score normalization
        mean_val = sum(values) / len(values)
        assert mean_val == pytest.approx(0.0, abs=0.01)

    def test_normalize_empty_track(self) -> None:
        track = GenomicTrack()
        track.normalize_values(method="minmax")  # Should not raise

    def test_add_feature_with_extra_kwargs(self) -> None:
        track = GenomicTrack()
        track.add_feature("chr1", 100, 200, value=5.0, name="peak_1", strand="+")
        features = track.get_features("chr1")
        assert features[0]["name"] == "peak_1"
        assert features[0]["strand"] == "+"


# ---------------------------------------------------------------------------
# load_bedgraph_track and save_bedgraph_track
# ---------------------------------------------------------------------------


class TestBedgraphTrackIO:
    """Tests for loading and saving BEDgraph tracks."""

    def test_load_bedgraph(self) -> None:
        repo_root = Path(__file__).resolve().parents[1]
        path = repo_root / "tests/data/epigenome/example.bedgraph"
        track = load_bedgraph_track(path, name="test_bg")
        assert track.get_total_features() == 3
        assert "chr1" in track.get_chromosomes()
        assert "chr2" in track.get_chromosomes()

    def test_save_and_reload_bedgraph(self, tmp_path: Path) -> None:
        track = GenomicTrack(name="save_test")
        track.add_feature("chr1", 0, 100, value=1.5)
        track.add_feature("chr1", 100, 200, value=3.0)
        track.add_feature("chr2", 0, 50, value=2.0)

        out_path = tmp_path / "output.bedgraph"
        save_bedgraph_track(track, out_path)
        assert out_path.exists()

        reloaded = load_bedgraph_track(out_path)
        assert reloaded.get_total_features() == 3


# ---------------------------------------------------------------------------
# load_bed_track and save_bed_track
# ---------------------------------------------------------------------------


class TestBedTrackIO:
    """Tests for loading and saving BED tracks."""

    def test_save_and_reload_bed(self, tmp_path: Path) -> None:
        track = GenomicTrack(name="bed_test")
        track.add_feature("chr1", 100, 500, value=50.0, name="peak1", strand="+")
        track.add_feature("chr2", 200, 400, value=30.0, name="peak2", strand="-")

        out_path = tmp_path / "output.bed"
        save_bed_track(track, out_path)
        assert out_path.exists()

        reloaded = load_bed_track(out_path)
        assert reloaded.get_total_features() == 2

    def test_load_bed_skips_header_lines(self, tmp_path: Path) -> None:
        bed_content = (
            "# comment\n"
            'track name="test" description="test"\n'
            "chr1\t100\t200\tpeak1\t50.0\t+\n"
            "chr1\t300\t400\tpeak2\t30.0\t-\n"
        )
        bed_path = tmp_path / "header_test.bed"
        bed_path.write_text(bed_content)

        track = load_bed_track(bed_path)
        assert track.get_total_features() == 2


# ---------------------------------------------------------------------------
# calculate_track_statistics
# ---------------------------------------------------------------------------


class TestCalculateTrackStatistics:
    """Tests for calculate_track_statistics."""

    def test_basic_statistics(self) -> None:
        track = GenomicTrack()
        track.add_feature("chr1", 0, 100, value=5.0)
        track.add_feature("chr1", 100, 300, value=10.0)
        track.add_feature("chr2", 0, 50, value=2.0)

        stats = calculate_track_statistics(track)
        assert stats["total_features"] == 3
        assert stats["total_chromosomes"] == 2
        assert stats["min_value"] == 2.0
        assert stats["max_value"] == 10.0

    def test_empty_track(self) -> None:
        track = GenomicTrack()
        assert calculate_track_statistics(track) == {}

    def test_length_statistics(self) -> None:
        track = GenomicTrack()
        track.add_feature("chr1", 0, 100, value=1.0)
        track.add_feature("chr1", 100, 300, value=2.0)

        stats = calculate_track_statistics(track)
        assert stats["min_length"] == 100
        assert stats["max_length"] == 200
        assert stats["mean_length"] == 150.0


# ---------------------------------------------------------------------------
# extract_track_region
# ---------------------------------------------------------------------------


class TestExtractTrackRegion:
    """Tests for extract_track_region."""

    def test_extract_subset(self) -> None:
        track = GenomicTrack()
        track.add_feature("chr1", 0, 100, value=1.0)
        track.add_feature("chr1", 100, 200, value=2.0)
        track.add_feature("chr1", 200, 300, value=3.0)
        track.add_feature("chr1", 300, 400, value=4.0)

        region = extract_track_region(track, "chr1", 150, 350)
        features = region.get_features("chr1")
        # Should capture features overlapping [150, 350)
        assert len(features) >= 2
        # Features should be clipped to region boundaries
        for f in features:
            assert f["start"] >= 150
            assert f["end"] <= 350

    def test_extract_empty_region(self) -> None:
        track = GenomicTrack()
        track.add_feature("chr1", 0, 100, value=1.0)

        region = extract_track_region(track, "chr1", 500, 600)
        assert region.get_total_features() == 0


# ---------------------------------------------------------------------------
# merge_tracks
# ---------------------------------------------------------------------------


class TestMergeTracks:
    """Tests for merge_tracks."""

    def test_union_merge(self) -> None:
        t1 = GenomicTrack()
        t1.add_feature("chr1", 0, 100, value=1.0)
        t2 = GenomicTrack()
        t2.add_feature("chr1", 200, 300, value=2.0)
        t2.add_feature("chr2", 0, 100, value=3.0)

        merged = merge_tracks([t1, t2], operation="union")
        assert merged.get_total_features() >= 3

    def test_empty_input(self) -> None:
        merged = merge_tracks([])
        assert merged.get_total_features() == 0

    def test_mean_merge(self) -> None:
        t1 = GenomicTrack()
        t1.add_feature("chr1", 0, 100, value=10.0)
        t2 = GenomicTrack()
        t2.add_feature("chr1", 50, 150, value=20.0)

        merged = merge_tracks([t1, t2], operation="mean")
        # Should have features on chr1
        assert merged.get_total_features() >= 1


# ---------------------------------------------------------------------------
# merge_feature_group
# ---------------------------------------------------------------------------


class TestMergeFeatureGroup:
    """Tests for merge_feature_group helper."""

    def test_mean_merge(self) -> None:
        features = [
            {"start": 100, "end": 200, "value": 10.0},
            {"start": 150, "end": 250, "value": 20.0},
        ]
        merged = merge_feature_group(features, "mean")
        assert merged["start"] == 100
        assert merged["end"] == 250
        assert merged["value"] == pytest.approx(15.0)

    def test_max_merge(self) -> None:
        features = [
            {"start": 100, "end": 200, "value": 10.0},
            {"start": 150, "end": 250, "value": 20.0},
        ]
        merged = merge_feature_group(features, "max")
        assert merged["value"] == 20.0


# ---------------------------------------------------------------------------
# compare_tracks
# ---------------------------------------------------------------------------


class TestCompareTracks:
    """Tests for compare_tracks."""

    def test_compare_same_tracks(self) -> None:
        track = GenomicTrack()
        track.add_feature("chr1", 0, 100, value=5.0)
        track.add_feature("chr1", 100, 200, value=10.0)

        result = compare_tracks(track, track)
        assert result["track1_features"] == 2
        assert result["track2_features"] == 2
        assert result["shared_chromosomes"] == 1

    def test_compare_different_chroms(self) -> None:
        t1 = GenomicTrack()
        t1.add_feature("chr1", 0, 100, value=5.0)
        t2 = GenomicTrack()
        t2.add_feature("chr2", 0, 100, value=5.0)

        result = compare_tracks(t1, t2)
        assert result["shared_chromosomes"] == 0
        assert result["unique_to_track1"] == 1
        assert result["unique_to_track2"] == 1


# ---------------------------------------------------------------------------
# generate_track_report
# ---------------------------------------------------------------------------


class TestGenerateTrackReport:
    """Tests for generate_track_report."""

    def test_report_generation(self) -> None:
        track = GenomicTrack(name="test_report")
        track.add_feature("chr1", 0, 100, value=5.0)
        track.add_feature("chr1", 100, 200, value=10.0)

        report = generate_track_report(track)
        assert "GENOMIC TRACK ANALYSIS REPORT" in report
        assert "test_report" in report

    def test_report_save_to_file(self, tmp_path: Path) -> None:
        track = GenomicTrack(name="saved_report")
        track.add_feature("chr1", 0, 100, value=5.0)

        out_path = tmp_path / "report.txt"
        report = generate_track_report(track, output_path=out_path)
        assert out_path.exists()
        assert "saved_report" in out_path.read_text()

    def test_empty_track_report(self) -> None:
        track = GenomicTrack()
        report = generate_track_report(track)
        assert "GENOMIC TRACK ANALYSIS REPORT" in report
