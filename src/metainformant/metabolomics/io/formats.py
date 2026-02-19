"""Mass spectrometry file format reading and writing.

Supports CSV-based metabolomics tables with metabolite × sample matrices,
plus stub readers for mzML and mzXML formats.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


@dataclass
class MetabolomicsDataset:
    """Container for metabolomics data.

    Attributes:
        metabolites: List of metabolite names/IDs.
        samples: List of sample names.
        intensities: 2D array (metabolites × samples) of intensity values.
        metadata: Optional dict of per-metabolite metadata (m/z, RT, etc.).
    """

    metabolites: list[str]
    samples: list[str]
    intensities: np.ndarray
    metadata: dict[str, Any] = field(default_factory=dict)


def read_csv(
    filepath: str | Path,
    delimiter: str = ",",
    metabolite_col: int = 0,
    data_start_col: int = 1,
) -> MetabolomicsDataset:
    """Read a metabolomics intensity matrix from CSV.

    Expected format: first column = metabolite IDs, remaining columns = samples.
    First row = header with sample names.

    Args:
        filepath: Path to CSV file.
        delimiter: Column delimiter.
        metabolite_col: Column index for metabolite names.
        data_start_col: First column of intensity data.

    Returns:
        MetabolomicsDataset with parsed data.
    """
    path = Path(filepath)
    lines = path.read_text().strip().split("\n")
    header = lines[0].split(delimiter)
    samples = [h.strip() for h in header[data_start_col:]]

    metabolites = []
    rows = []
    for line in lines[1:]:
        parts = line.split(delimiter)
        metabolites.append(parts[metabolite_col].strip())
        row = []
        for val in parts[data_start_col:]:
            val = val.strip()
            try:
                row.append(float(val))
            except ValueError:
                row.append(0.0)
        rows.append(row)

    intensities = np.array(rows, dtype=np.float64)
    return MetabolomicsDataset(
        metabolites=metabolites,
        samples=samples,
        intensities=intensities,
    )


def write_csv(
    dataset: MetabolomicsDataset,
    filepath: str | Path,
    delimiter: str = ",",
) -> None:
    """Write a metabolomics dataset to CSV.

    Args:
        dataset: MetabolomicsDataset to write.
        filepath: Output file path.
        delimiter: Column delimiter.
    """
    path = Path(filepath)
    header = delimiter.join(["metabolite"] + dataset.samples)
    lines = [header]
    for i, met in enumerate(dataset.metabolites):
        vals = delimiter.join(str(v) for v in dataset.intensities[i])
        lines.append(f"{met}{delimiter}{vals}")
    path.write_text("\n".join(lines) + "\n")


@dataclass
class MassSpectrum:
    """A single mass spectrum.

    Attributes:
        scan_number: Scan index.
        retention_time: RT in seconds.
        ms_level: 1 for MS1, 2 for MS2, etc.
        mz_array: 1D array of m/z values.
        intensity_array: 1D array of intensities.
        precursor_mz: Precursor m/z for MS2+ spectra (None for MS1).
        total_ion_current: Sum of all intensities.
    """

    scan_number: int
    retention_time: float
    ms_level: int
    mz_array: np.ndarray
    intensity_array: np.ndarray
    precursor_mz: float | None = None
    total_ion_current: float = 0.0


def read_mgf(filepath: str | Path) -> list[MassSpectrum]:
    """Read spectra from MGF (Mascot Generic Format) file.

    MGF is a simple text format commonly used for MS2 spectra in metabolomics
    and proteomics.

    Args:
        filepath: Path to MGF file.

    Returns:
        List of MassSpectrum objects.
    """
    path = Path(filepath)
    text = path.read_text()
    spectra: list[MassSpectrum] = []
    scan_idx = 0

    blocks = text.split("BEGIN IONS")
    for block in blocks[1:]:  # skip text before first BEGIN IONS
        end_idx = block.find("END IONS")
        if end_idx == -1:
            continue
        content = block[:end_idx].strip()
        lines = content.split("\n")

        rt = 0.0
        precursor = None
        mz_list: list[float] = []
        intensity_list: list[float] = []

        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith("PEPMASS="):
                parts = line.split("=")[1].split()
                precursor = float(parts[0])
            elif line.startswith("RTINSECONDS="):
                rt = float(line.split("=")[1])
            elif line.startswith("SCANS="):
                scan_idx = int(line.split("=")[1])
            elif line[0].isdigit():
                parts = line.split()
                if len(parts) >= 2:
                    mz_list.append(float(parts[0]))
                    intensity_list.append(float(parts[1]))

        mz_arr = np.array(mz_list, dtype=np.float64)
        int_arr = np.array(intensity_list, dtype=np.float64)

        spectra.append(
            MassSpectrum(
                scan_number=scan_idx,
                retention_time=rt,
                ms_level=2,
                mz_array=mz_arr,
                intensity_array=int_arr,
                precursor_mz=precursor,
                total_ion_current=float(int_arr.sum()),
            )
        )
        scan_idx += 1

    return spectra


def write_mgf(
    spectra: list[MassSpectrum],
    filepath: str | Path,
) -> None:
    """Write spectra to MGF format.

    Args:
        spectra: List of MassSpectrum objects.
        filepath: Output file path.
    """
    path = Path(filepath)
    lines: list[str] = []

    for spec in spectra:
        lines.append("BEGIN IONS")
        if spec.precursor_mz is not None:
            lines.append(f"PEPMASS={spec.precursor_mz}")
        lines.append(f"RTINSECONDS={spec.retention_time}")
        lines.append(f"SCANS={spec.scan_number}")
        for mz, intensity in zip(spec.mz_array, spec.intensity_array):
            lines.append(f"{mz:.6f} {intensity:.1f}")
        lines.append("END IONS")
        lines.append("")

    path.write_text("\n".join(lines))


def filter_spectra(
    spectra: list[MassSpectrum],
    min_peaks: int = 5,
    min_tic: float = 0.0,
    ms_level: int | None = None,
    rt_range: tuple[float, float] | None = None,
) -> list[MassSpectrum]:
    """Filter spectra based on quality criteria.

    Args:
        spectra: Input spectra list.
        min_peaks: Minimum number of peaks required.
        min_tic: Minimum total ion current.
        ms_level: Keep only this MS level (None = all).
        rt_range: (min_rt, max_rt) in seconds (None = all).

    Returns:
        Filtered list of MassSpectrum.
    """
    result: list[MassSpectrum] = []
    for spec in spectra:
        if len(spec.mz_array) < min_peaks:
            continue
        if spec.total_ion_current < min_tic:
            continue
        if ms_level is not None and spec.ms_level != ms_level:
            continue
        if rt_range is not None and not (rt_range[0] <= spec.retention_time <= rt_range[1]):
            continue
        result.append(spec)
    return result


def extract_chromatogram(
    spectra: list[MassSpectrum],
    target_mz: float,
    ppm_tolerance: float = 10.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract ion chromatogram (EIC) for a specific m/z across scans.

    Args:
        spectra: List of MS1 spectra ordered by retention time.
        target_mz: Target m/z value.
        ppm_tolerance: Mass tolerance in ppm.

    Returns:
        Tuple of (retention_times, intensities) arrays.
    """
    rts: list[float] = []
    intensities: list[float] = []

    mz_tol = target_mz * ppm_tolerance / 1e6

    for spec in spectra:
        if spec.ms_level != 1:
            continue
        rts.append(spec.retention_time)
        mask = np.abs(spec.mz_array - target_mz) <= mz_tol
        if np.any(mask):
            intensities.append(float(spec.intensity_array[mask].max()))
        else:
            intensities.append(0.0)

    return np.array(rts), np.array(intensities)
