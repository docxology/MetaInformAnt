"""Tests for metabolomics I/O formats."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metainformant.metabolomics.io.formats import (
    MassSpectrum,
    MetabolomicsDataset,
    extract_chromatogram,
    filter_spectra,
    read_csv,
    read_mgf,
    write_csv,
    write_mgf,
)


class TestCSVReadWrite:
    """Tests for CSV metabolomics data I/O."""

    def test_roundtrip(self, tmp_path: Path) -> None:
        """Write and read back should preserve data."""
        dataset = MetabolomicsDataset(
            metabolites=["met_A", "met_B", "met_C"],
            samples=["S1", "S2", "S3"],
            intensities=np.array([
                [100.0, 200.0, 300.0],
                [50.0, 150.0, 250.0],
                [10.0, 20.0, 30.0],
            ]),
        )
        filepath = tmp_path / "test.csv"
        write_csv(dataset, filepath)
        loaded = read_csv(filepath)

        assert loaded.metabolites == dataset.metabolites
        assert loaded.samples == dataset.samples
        np.testing.assert_allclose(loaded.intensities, dataset.intensities)

    def test_read_with_header(self, tmp_path: Path) -> None:
        """CSV with header row should parse sample names."""
        filepath = tmp_path / "data.csv"
        filepath.write_text("metabolite,sample1,sample2\nglucose,100,200\nfructose,50,75\n")
        dataset = read_csv(filepath)
        assert dataset.samples == ["sample1", "sample2"]
        assert dataset.metabolites == ["glucose", "fructose"]
        assert dataset.intensities.shape == (2, 2)

    def test_read_handles_invalid_values(self, tmp_path: Path) -> None:
        """Invalid numeric values should be replaced with 0.0."""
        filepath = tmp_path / "bad.csv"
        filepath.write_text("met,S1,S2\nA,100,NA\nB,missing,50\n")
        dataset = read_csv(filepath)
        assert dataset.intensities[0, 1] == 0.0
        assert dataset.intensities[1, 0] == 0.0
        assert dataset.intensities[1, 1] == 50.0

    def test_custom_delimiter(self, tmp_path: Path) -> None:
        """Tab-delimited files should work."""
        filepath = tmp_path / "tab.tsv"
        filepath.write_text("met\tS1\tS2\nA\t10\t20\n")
        dataset = read_csv(filepath, delimiter="\t")
        assert dataset.samples == ["S1", "S2"]
        assert dataset.intensities[0, 0] == 10.0


class TestMGFReadWrite:
    """Tests for MGF format I/O."""

    def test_roundtrip(self, tmp_path: Path) -> None:
        """Write and read MGF should preserve spectral data."""
        spectra = [
            MassSpectrum(
                scan_number=1,
                retention_time=60.0,
                ms_level=2,
                mz_array=np.array([100.0, 200.0, 300.0]),
                intensity_array=np.array([1000.0, 500.0, 250.0]),
                precursor_mz=350.5,
                total_ion_current=1750.0,
            ),
            MassSpectrum(
                scan_number=2,
                retention_time=120.0,
                ms_level=2,
                mz_array=np.array([150.0, 250.0]),
                intensity_array=np.array([800.0, 400.0]),
                precursor_mz=275.0,
                total_ion_current=1200.0,
            ),
        ]
        filepath = tmp_path / "test.mgf"
        write_mgf(spectra, filepath)
        loaded = read_mgf(filepath)

        assert len(loaded) == 2
        np.testing.assert_allclose(loaded[0].mz_array, spectra[0].mz_array, atol=0.001)
        np.testing.assert_allclose(loaded[0].intensity_array, spectra[0].intensity_array, atol=1.0)
        assert loaded[0].precursor_mz == pytest.approx(350.5)

    def test_read_mgf_retention_time(self, tmp_path: Path) -> None:
        """MGF parser should extract retention time."""
        mgf_content = (
            "BEGIN IONS\n"
            "PEPMASS=500.0\n"
            "RTINSECONDS=300.5\n"
            "100.0 1000\n"
            "200.0 500\n"
            "END IONS\n"
        )
        filepath = tmp_path / "rt.mgf"
        filepath.write_text(mgf_content)
        spectra = read_mgf(filepath)
        assert len(spectra) == 1
        assert spectra[0].retention_time == pytest.approx(300.5)

    def test_empty_mgf(self, tmp_path: Path) -> None:
        filepath = tmp_path / "empty.mgf"
        filepath.write_text("")
        spectra = read_mgf(filepath)
        assert spectra == []


class TestFilterSpectra:
    """Tests for spectral filtering."""

    def setup_method(self) -> None:
        self.spectra = [
            MassSpectrum(0, 10.0, 1, np.arange(10.0), np.arange(10.0), total_ion_current=45.0),
            MassSpectrum(1, 20.0, 2, np.arange(3.0), np.ones(3), precursor_mz=100.0, total_ion_current=3.0),
            MassSpectrum(2, 30.0, 1, np.arange(20.0), np.arange(20.0), total_ion_current=190.0),
        ]

    def test_filter_by_min_peaks(self) -> None:
        result = filter_spectra(self.spectra, min_peaks=5)
        assert len(result) == 2  # first and third have >= 5 peaks

    def test_filter_by_ms_level(self) -> None:
        result = filter_spectra(self.spectra, ms_level=1, min_peaks=0)
        assert len(result) == 2
        assert all(s.ms_level == 1 for s in result)

    def test_filter_by_rt_range(self) -> None:
        result = filter_spectra(self.spectra, rt_range=(15.0, 25.0), min_peaks=0)
        assert len(result) == 1
        assert result[0].retention_time == 20.0

    def test_filter_by_min_tic(self) -> None:
        result = filter_spectra(self.spectra, min_tic=100.0, min_peaks=0)
        assert len(result) == 1
        assert result[0].total_ion_current >= 100.0


class TestExtractChromatogram:
    """Tests for EIC extraction."""

    def test_basic_extraction(self) -> None:
        spectra = [
            MassSpectrum(i, float(i * 10), 1, np.array([100.0, 200.0, 300.0]),
                         np.array([float(i * 100), 50.0, 25.0]))
            for i in range(5)
        ]
        rts, intensities = extract_chromatogram(spectra, target_mz=100.0, ppm_tolerance=50.0)
        assert len(rts) == 5
        assert len(intensities) == 5
        assert intensities[0] == 0.0  # scan 0: 0*100 = 0
        assert intensities[4] == 400.0  # scan 4: 4*100 = 400

    def test_no_ms1_spectra(self) -> None:
        """MS2-only spectra should return empty arrays."""
        spectra = [MassSpectrum(0, 10.0, 2, np.array([100.0]), np.array([50.0]))]
        rts, intensities = extract_chromatogram(spectra, target_mz=100.0)
        assert len(rts) == 0

    def test_target_not_present(self) -> None:
        spectra = [
            MassSpectrum(0, 10.0, 1, np.array([200.0, 300.0]), np.array([100.0, 50.0])),
        ]
        rts, intensities = extract_chromatogram(spectra, target_mz=500.0, ppm_tolerance=5.0)
        assert len(rts) == 1
        assert intensities[0] == 0.0
