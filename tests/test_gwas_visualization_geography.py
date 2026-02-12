"""Tests for GWAS geographic visualization functions."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import numpy as np
import pytest

from metainformant.gwas.visualization.population.geography import (
    allele_frequency_map,
    population_count_map,
    sample_map,
)


def _make_geo_metadata(
    n_per_pop: int = 10,
    populations: list[str] | None = None,
    centers: list[tuple[float, float]] | None = None,
    rng: np.random.RandomState | None = None,
) -> Dict[str, Dict[str, Any]]:
    """Build synthetic geographic metadata for testing.

    Args:
        n_per_pop: Number of samples per population.
        populations: Population names.  Defaults to three populations.
        centers: (latitude, longitude) centre for each population.
        rng: Optional random state for reproducibility.

    Returns:
        Metadata dictionary keyed by sample id.
    """
    if rng is None:
        rng = np.random.RandomState(42)
    if populations is None:
        populations = ["Pop_A", "Pop_B", "Pop_C"]
    if centers is None:
        centers = [(35.0, -120.0), (48.0, 2.0), (-33.0, 151.0)]

    metadata: Dict[str, Dict[str, Any]] = {}
    idx = 0
    for pop, (lat_c, lon_c) in zip(populations, centers):
        for _ in range(n_per_pop):
            sid = f"sample_{idx}"
            metadata[sid] = {
                "latitude": lat_c + rng.normal(0, 1.5),
                "longitude": lon_c + rng.normal(0, 1.5),
                "population": pop,
                "coverage": rng.uniform(10, 50),
            }
            idx += 1
    return metadata


# ---------------------------------------------------------------------------
# sample_map tests
# ---------------------------------------------------------------------------


def test_sample_map_basic(tmp_path: Path) -> None:
    """Test sample_map with 30 samples across 3 populations."""
    metadata = _make_geo_metadata(n_per_pop=10)
    output_file = tmp_path / "sample_map.png"

    result = sample_map(metadata, output_file=output_file)

    assert result["status"] == "success"
    assert result["n_samples_plotted"] == 30
    assert result["n_skipped"] == 0
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_sample_map_with_missing_coords(tmp_path: Path) -> None:
    """Test sample_map skips samples missing latitude or longitude."""
    metadata = _make_geo_metadata(n_per_pop=10)

    # Add 5 samples with missing coordinates
    metadata["no_lat"] = {"longitude": 10.0, "population": "Pop_X"}
    metadata["no_lon"] = {"latitude": 20.0, "population": "Pop_X"}
    metadata["no_both"] = {"population": "Pop_X"}
    metadata["none_lat"] = {"latitude": None, "longitude": 5.0, "population": "Pop_X"}
    metadata["none_lon"] = {"latitude": 5.0, "longitude": None, "population": "Pop_X"}

    output_file = tmp_path / "sample_map_skipped.png"
    result = sample_map(metadata, output_file=output_file)

    assert result["status"] == "success"
    assert result["n_samples_plotted"] == 30
    assert result["n_skipped"] == 5
    assert output_file.exists()


def test_sample_map_with_size_by(tmp_path: Path) -> None:
    """Test sample_map with size_by parameter scaling points."""
    metadata = _make_geo_metadata(n_per_pop=10)
    output_file = tmp_path / "sample_map_sized.png"

    result = sample_map(metadata, output_file=output_file, size_by="coverage")

    assert result["status"] == "success"
    assert result["n_samples_plotted"] == 30
    assert result["n_skipped"] == 0
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_sample_map_empty_metadata() -> None:
    """Test sample_map with empty metadata returns failed."""
    result = sample_map({})

    assert result["status"] == "failed"
    assert result["n_samples_plotted"] == 0


def test_sample_map_no_output_file() -> None:
    """Test sample_map without saving to file still succeeds."""
    metadata = _make_geo_metadata(n_per_pop=5)

    result = sample_map(metadata)

    assert result["status"] == "success"
    assert result["output_path"] is None
    assert result["n_samples_plotted"] == 15


# ---------------------------------------------------------------------------
# allele_frequency_map tests
# ---------------------------------------------------------------------------


def test_allele_frequency_map_basic(tmp_path: Path) -> None:
    """Test allele_frequency_map with 3 populations."""
    metadata = _make_geo_metadata(n_per_pop=10)
    allele_freqs = {
        "Pop_A": {"ref_freq": 0.8, "alt_freq": 0.2},
        "Pop_B": {"ref_freq": 0.5, "alt_freq": 0.5},
        "Pop_C": {"ref_freq": 0.1, "alt_freq": 0.9},
    }
    output_file = tmp_path / "allele_freq_map.png"

    result = allele_frequency_map(
        metadata,
        allele_freqs,
        variant_id="rs12345",
        output_file=output_file,
    )

    assert result["status"] == "success"
    assert result["n_populations"] == 3
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_allele_frequency_map_custom_title(tmp_path: Path) -> None:
    """Test allele_frequency_map with a custom title."""
    metadata = _make_geo_metadata(n_per_pop=5, populations=["EU", "AS"], centers=[(50.0, 10.0), (35.0, 105.0)])
    allele_freqs = {
        "EU": {"ref_freq": 0.7, "alt_freq": 0.3},
        "AS": {"ref_freq": 0.4, "alt_freq": 0.6},
    }
    output_file = tmp_path / "allele_freq_custom.png"

    result = allele_frequency_map(
        metadata,
        allele_freqs,
        variant_id="rs99999",
        output_file=output_file,
        title="Custom Allele Map Title",
    )

    assert result["status"] == "success"
    assert result["n_populations"] == 2
    assert output_file.exists()


def test_allele_frequency_map_empty_inputs() -> None:
    """Test allele_frequency_map with empty metadata returns failed."""
    result = allele_frequency_map({}, {}, variant_id="rs0")

    assert result["status"] == "failed"
    assert result["n_populations"] == 0


# ---------------------------------------------------------------------------
# population_count_map tests
# ---------------------------------------------------------------------------


def test_population_count_map_basic(tmp_path: Path) -> None:
    """Test population_count_map with 3 populations."""
    metadata = _make_geo_metadata(n_per_pop=10)
    output_file = tmp_path / "pop_count_map.png"

    result = population_count_map(metadata, output_file=output_file)

    assert result["status"] == "success"
    assert result["n_populations"] == 3
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_population_count_map_custom_title(tmp_path: Path) -> None:
    """Test population_count_map with custom title."""
    metadata = _make_geo_metadata(n_per_pop=8)
    output_file = tmp_path / "pop_count_titled.png"

    result = population_count_map(
        metadata,
        output_file=output_file,
        title="My Population Counts",
    )

    assert result["status"] == "success"
    assert result["n_populations"] == 3
    assert output_file.exists()


def test_population_count_map_empty_metadata() -> None:
    """Test population_count_map with empty metadata returns failed."""
    result = population_count_map({})

    assert result["status"] == "failed"
    assert result["n_populations"] == 0


def test_population_count_map_no_output_file() -> None:
    """Test population_count_map without saving to file still succeeds."""
    metadata = _make_geo_metadata(n_per_pop=5)

    result = population_count_map(metadata)

    assert result["status"] == "success"
    assert result["output_path"] is None
    assert result["n_populations"] == 3
