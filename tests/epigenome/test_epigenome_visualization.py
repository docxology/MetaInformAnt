"""Tests for epigenome visualization functions.

Tests plot_methylation_profile, plot_chipseq_peaks, plot_atacseq_signal,
plot_histone_modification_heatmap, plot_differential_methylation,
plot_chromatin_states, plot_epigenetic_correlation_heatmap,
plot_genome_browser_tracks, and plot_dna_methylation_clusters using
real matplotlib rendering with the Agg backend. NO MOCKING.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.axes import Axes

from metainformant.epigenome.visualization.visualization import (
    plot_atacseq_signal,
    plot_chipseq_peaks,
    plot_chromatin_states,
    plot_differential_methylation,
    plot_dna_methylation_clusters,
    plot_epigenetic_correlation_heatmap,
    plot_genome_browser_tracks,
    plot_histone_modification_heatmap,
    plot_methylation_profile,
)


@pytest.fixture(autouse=True)
def _close_figures():
    """Close all matplotlib figures after each test to prevent memory leaks."""
    yield
    plt.close("all")


# ---------------------------------------------------------------------------
# plot_methylation_profile
# ---------------------------------------------------------------------------


class TestPlotMethylationProfile:
    """Tests for plot_methylation_profile."""

    def test_returns_axes(self) -> None:
        data = np.random.RandomState(42).random(100)
        ax = plot_methylation_profile(data)
        assert isinstance(ax, Axes)

    def test_with_positions(self) -> None:
        data = np.random.RandomState(42).random(50)
        positions = np.arange(1000, 1050)
        ax = plot_methylation_profile(data, positions=positions)
        assert isinstance(ax, Axes)

    def test_with_provided_axes(self) -> None:
        fig, ax = plt.subplots()
        data = np.random.RandomState(42).random(20)
        returned_ax = plot_methylation_profile(data, ax=ax)
        assert returned_ax is ax

    def test_large_data_with_smoothing(self) -> None:
        # Triggers the smoothed trend line (len > 10)
        data = np.random.RandomState(42).random(200)
        ax = plot_methylation_profile(data)
        assert isinstance(ax, Axes)


# ---------------------------------------------------------------------------
# plot_chipseq_peaks
# ---------------------------------------------------------------------------


class TestPlotChipseqPeaks:
    """Tests for plot_chipseq_peaks."""

    def test_returns_axes(self) -> None:
        peaks = [
            {"start": 100, "end": 200, "score": 50.0},
            {"start": 300, "end": 500, "score": 100.0},
            {"start": 700, "end": 800, "score": 75.0},
        ]
        ax = plot_chipseq_peaks(peaks, "chr1")
        assert isinstance(ax, Axes)

    def test_with_provided_axes(self) -> None:
        fig, ax = plt.subplots()
        peaks = [{"start": 0, "end": 100, "score": 10.0}]
        returned_ax = plot_chipseq_peaks(peaks, "chr1", ax=ax)
        assert returned_ax is ax

    def test_empty_peaks(self) -> None:
        ax = plot_chipseq_peaks([], "chr1")
        assert isinstance(ax, Axes)


# ---------------------------------------------------------------------------
# plot_atacseq_signal
# ---------------------------------------------------------------------------


class TestPlotAtacseqSignal:
    """Tests for plot_atacseq_signal."""

    def test_returns_axes(self) -> None:
        rng = np.random.RandomState(42)
        signal = rng.random(100) * 50
        positions = np.arange(1000, 1100)
        ax = plot_atacseq_signal(signal, positions)
        assert isinstance(ax, Axes)

    def test_with_provided_axes(self) -> None:
        fig, ax = plt.subplots()
        signal = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        positions = np.array([100, 200, 300, 400, 500])
        returned_ax = plot_atacseq_signal(signal, positions, ax=ax)
        assert returned_ax is ax


# ---------------------------------------------------------------------------
# plot_histone_modification_heatmap
# ---------------------------------------------------------------------------


class TestPlotHistoneModificationHeatmap:
    """Tests for plot_histone_modification_heatmap."""

    def test_returns_axes(self) -> None:
        rng = np.random.RandomState(42)
        histone_data = {
            "H3K4me3": rng.random(50) * 10,
            "H3K27ac": rng.random(50) * 8,
            "H3K27me3": rng.random(50) * 5,
        }
        ax = plot_histone_modification_heatmap(histone_data)
        assert isinstance(ax, Axes)

    def test_with_positions(self) -> None:
        rng = np.random.RandomState(42)
        histone_data = {
            "H3K4me3": rng.random(30),
            "H3K4me1": rng.random(30),
        }
        positions = np.arange(30)
        ax = plot_histone_modification_heatmap(histone_data, positions=positions)
        assert isinstance(ax, Axes)


# ---------------------------------------------------------------------------
# plot_differential_methylation
# ---------------------------------------------------------------------------


class TestPlotDifferentialMethylation:
    """Tests for plot_differential_methylation (volcano plot)."""

    def test_returns_axes(self) -> None:
        rng = np.random.RandomState(42)
        results = [
            {"delta": rng.uniform(-0.5, 0.5), "pvalue": rng.uniform(0.001, 0.5), "position": i} for i in range(50)
        ]
        ax = plot_differential_methylation(results)
        assert isinstance(ax, Axes)

    def test_with_custom_threshold(self) -> None:
        results = [
            {"delta": 0.3, "pvalue": 0.01, "position": 100},
            {"delta": -0.1, "pvalue": 0.5, "position": 200},
        ]
        ax = plot_differential_methylation(results, sig_threshold=0.01)
        assert isinstance(ax, Axes)

    def test_empty_results(self) -> None:
        ax = plot_differential_methylation([])
        assert isinstance(ax, Axes)


# ---------------------------------------------------------------------------
# plot_chromatin_states
# ---------------------------------------------------------------------------


class TestPlotChromatinStates:
    """Tests for plot_chromatin_states."""

    def test_returns_axes(self) -> None:
        states = np.array([0, 0, 1, 1, 2, 2, 0, 0, 1, 1])
        positions = np.arange(10) * 200
        ax = plot_chromatin_states(states, positions)
        assert isinstance(ax, Axes)

    def test_with_state_labels(self) -> None:
        states = np.array([0, 0, 1, 1, 2, 2])
        positions = np.arange(6) * 100
        labels = {0: "Promoter", 1: "Enhancer", 2: "Repressed"}
        ax = plot_chromatin_states(states, positions, state_labels=labels)
        assert isinstance(ax, Axes)


# ---------------------------------------------------------------------------
# plot_epigenetic_correlation_heatmap
# ---------------------------------------------------------------------------


class TestPlotEpigeneticCorrelationHeatmap:
    """Tests for plot_epigenetic_correlation_heatmap."""

    def test_returns_axes(self) -> None:
        rng = np.random.RandomState(42)
        data = {
            "H3K4me3": rng.random(100),
            "H3K27ac": rng.random(100),
            "H3K27me3": rng.random(100),
        }
        ax = plot_epigenetic_correlation_heatmap(data)
        assert isinstance(ax, Axes)

    def test_with_correlated_data(self) -> None:
        rng = np.random.RandomState(42)
        base = rng.random(50)
        data = {
            "mark_A": base + rng.normal(0, 0.1, 50),
            "mark_B": base * 0.8 + rng.normal(0, 0.1, 50),
        }
        ax = plot_epigenetic_correlation_heatmap(data)
        assert isinstance(ax, Axes)


# ---------------------------------------------------------------------------
# plot_genome_browser_tracks
# ---------------------------------------------------------------------------


class TestPlotGenomeBrowserTracks:
    """Tests for plot_genome_browser_tracks."""

    def test_returns_axes(self) -> None:
        rng = np.random.RandomState(42)
        tracks = {
            "ChIP-seq": {"data": rng.random(100).tolist(), "color": "red", "label": "H3K4me3"},
            "ATAC-seq": {"data": rng.random(100).tolist(), "color": "blue", "label": "Accessibility"},
        }
        ax = plot_genome_browser_tracks(tracks, region_start=0, region_end=10000)
        assert isinstance(ax, Axes)

    def test_empty_tracks(self) -> None:
        ax = plot_genome_browser_tracks({}, region_start=0, region_end=1000)
        assert isinstance(ax, Axes)


# ---------------------------------------------------------------------------
# plot_dna_methylation_clusters
# ---------------------------------------------------------------------------


class TestPlotDnaMethylationClusters:
    """Tests for plot_dna_methylation_clusters."""

    def test_returns_axes(self) -> None:
        rng = np.random.RandomState(42)
        matrix = rng.random((20, 50))
        ax = plot_dna_methylation_clusters(matrix)
        assert isinstance(ax, Axes)

    def test_with_cluster_labels(self) -> None:
        rng = np.random.RandomState(42)
        matrix = rng.random((20, 50))
        labels = np.array([0] * 10 + [1] * 10)
        ax = plot_dna_methylation_clusters(matrix, cluster_labels=labels)
        assert isinstance(ax, Axes)

    def test_with_provided_axes(self) -> None:
        fig, ax = plt.subplots()
        rng = np.random.RandomState(42)
        matrix = rng.random((10, 30))
        returned_ax = plot_dna_methylation_clusters(matrix, ax=ax)
        assert returned_ax is ax
