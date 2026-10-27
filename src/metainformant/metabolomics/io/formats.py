"""Mass spectrometry file format reading and writing.

Reads and writes two formats: CSV metabolite × sample intensity tables
(:func:`read_csv`, :func:`write_csv`) and MGF (Mascot Generic Format) spectra
(:func:`read_mgf`, :func:`write_mgf`). It also filters spectra and extracts ion
chromatograms from them. mzML and mzXML are not supported; convert such files to
MGF first.

The readers never invent values. A missing intensity stays missing (NaN), and a CSV
cell or an MGF header or peak line that cannot be interpreted raises :class:`ValueError`
naming where it is.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

# Cell texts (compared case-insensitively, after stripping) that mean "no value" in an intensity table.
_MISSING_TOKENS = frozenset({"", "na", "n/a", "#n/a", "nan", "null", "none"})
# A European decimal comma ("1,5"): worth naming in the error because it is the usual cause of a non-numeric cell.
_DECIMAL_COMMA = re.compile(r"^[+-]?\d+,\d+$")
# MGF lines starting with one of these characters are comments (Matrix Science MGF specification).
_MGF_COMMENT_PREFIXES = ("#", ";", "!", "/")


@dataclass
class MetabolomicsDataset:
    """Container for metabolomics data.

    Attributes:
        metabolites: List of metabolite names/IDs.
        samples: List of sample names.
        intensities: 2D array (metabolites × samples) of intensity values.
        metadata: Optional dict of per-metabolite metadata (m/z, RT, etc.), keyed by metabolite
            name: ``{metabolite: {column_name: value}}``. :func:`read_csv` fills it from the
            annotation columns that precede the intensity columns.
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
    First row = header with sample names. Blank lines are ignored.

    Missing values stay missing. A blank cell or one of ``NA``, ``N/A``, ``#N/A``, ``NaN``,
    ``null`` or ``None`` (any case) is returned as ``NaN``, never as ``0.0``: a true zero and
    an absent measurement mean different things, and
    :func:`~metainformant.metabolomics.analysis.identification.missing_value_imputation`
    treats NaN and zero as missing alike. Any other cell that is not a number (text, or a
    decimal comma such as ``1,5``) raises :class:`ValueError` naming the line, metabolite and
    sample instead of being replaced by a value.

    Columns before ``data_start_col`` other than ``metabolite_col`` are per-metabolite
    annotations (for example ``mz`` and ``rt``). They are returned in
    ``MetabolomicsDataset.metadata`` as ``{metabolite: {column_name: value}}``; numeric
    values become floats, missing-value tokens become ``NaN`` and anything else stays text.
    With the default layout there are no such columns and ``metadata`` is empty.

    Args:
        filepath: Path to CSV file.
        delimiter: Column delimiter.
        metabolite_col: Column index for metabolite names.
        data_start_col: First column of intensity data.

    Returns:
        MetabolomicsDataset with parsed data.

    Raises:
        ValueError: If a row does not have as many cells as the header, a cell is neither a
            number nor a missing-value token, ``metabolite_col``/``data_start_col`` do not
            fit the header, or annotation columns are present and a metabolite or
            annotation column name occurs twice (the metadata would be ambiguous).
    """
    path = Path(filepath)
    # (1-based line number, text) for every non-blank line; the header keeps its leading delimiters
    numbered = [(n, line.rstrip("\r")) for n, line in enumerate(path.read_text().split("\n"), start=1) if line.strip()]
    if not numbered:
        return MetabolomicsDataset(metabolites=[], samples=[], intensities=np.array([], dtype=np.float64))

    # a UTF-8 byte-order mark (Excel's "CSV UTF-8") would otherwise become part of the first header cell
    header = numbered[0][1].lstrip("\ufeff").split(delimiter)
    if not 0 <= metabolite_col < len(header) or not 0 <= data_start_col <= len(header):
        raise ValueError(
            f"{path}: metabolite_col={metabolite_col} and data_start_col={data_start_col} do not fit "
            f"the {len(header)} header columns"
        )
    samples = [h.strip() for h in header[data_start_col:]]
    annotation_cols = [i for i in range(data_start_col) if i != metabolite_col]
    annotation_names = [header[i].strip() or f"column_{i}" for i in annotation_cols]
    if len(set(annotation_names)) != len(annotation_names):
        raise ValueError(f"{path}: duplicate annotation column names {annotation_names}")

    metabolites: list[str] = []
    rows: list[list[float]] = []
    annotations: list[dict[str, float | str]] = []
    for line_no, line in numbered[1:]:
        parts = line.split(delimiter)
        if len(parts) != len(header):
            raise ValueError(
                f"{path}: line {line_no} has {len(parts)} cells but the header has {len(header)}; "
                "check the delimiter (a decimal comma in a comma-delimited file splits one number "
                "into two cells)"
            )
        name = parts[metabolite_col].strip()
        row: list[float] = []
        for sample, cell in zip(samples, parts[data_start_col:]):
            text = cell.strip()
            try:
                row.append(float(text))
            except ValueError:
                if text.lower() not in _MISSING_TOKENS:
                    hint = (
                        " (a decimal comma is not supported; use a decimal point)" if _DECIMAL_COMMA.match(text) else ""
                    )
                    raise ValueError(
                        f"{path}: line {line_no}, metabolite {name!r}, sample {sample!r}: cannot read "
                        f"{text!r} as a number{hint}; write a missing value as a blank cell or NA"
                    ) from None
                row.append(math.nan)
        metabolites.append(name)
        rows.append(row)
        annotations.append({n: _annotation_value(parts[i]) for n, i in zip(annotation_names, annotation_cols)})

    metadata: dict[str, Any] = {}
    if annotation_cols:
        for name, annotation in zip(metabolites, annotations):
            if name in metadata:
                raise ValueError(
                    f"{path}: metabolite {name!r} occurs more than once, so its annotation columns "
                    f"{annotation_names} would be ambiguous"
                )
            metadata[name] = annotation

    intensities = np.array(rows, dtype=np.float64)
    return MetabolomicsDataset(
        metabolites=metabolites,
        samples=samples,
        intensities=intensities,
        metadata=metadata,
    )


def _annotation_value(cell: str) -> float | str:
    """Interpret an annotation cell: a float, ``NaN`` for a missing-value token, otherwise the text."""
    text = cell.strip()
    try:
        return float(text)
    except ValueError:
        return math.nan if text.lower() in _MISSING_TOKENS else text


def write_csv(
    dataset: MetabolomicsDataset,
    filepath: str | Path,
    delimiter: str = ",",
) -> None:
    """Write a metabolomics dataset to CSV.

    Missing values (NaN) are written as ``nan``, which :func:`read_csv` reads back as NaN.
    ``dataset.metadata`` is not written.

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

    Inside a ``BEGIN IONS`` / ``END IONS`` block a line is interpreted by its shape:

    * blank lines and comment lines (starting with ``#``, ``;``, ``!`` or ``/``) are skipped;
    * a line containing ``=`` is a header (``PEPMASS``, ``RTINSECONDS`` and ``SCANS`` are read,
      every other header such as ``TITLE`` or ``CHARGE`` is ignored);
    * every other line is a peak and must start with two finite numbers, ``m/z intensity``
      (further columns, such as a charge, are ignored). A peak line is never skipped silently, so
      ``.5 100`` and ``+100.0 50`` are read and a malformed peak line is an error.

    Text before the first ``BEGIN IONS`` is ignored, and so is a ``BEGIN IONS`` block that is never
    closed by ``END IONS`` (for example the last spectrum of a truncated file).

    Args:
        filepath: Path to MGF file.

    Returns:
        List of MassSpectrum objects.

    Raises:
        ValueError: If ``PEPMASS`` has no value, a ``PEPMASS``, ``RTINSECONDS`` or ``SCANS``
            value is not a number (``SCANS`` an integer), or a peak line is malformed. The message
            names the file, the spectrum (counted from 1 in file order) and the offending line.
    """
    path = Path(filepath)
    text = path.read_text()
    spectra: list[MassSpectrum] = []
    scan_idx = 0

    blocks = text.split("BEGIN IONS")
    for block_no, block in enumerate(blocks[1:], start=1):  # skip text before first BEGIN IONS
        end_idx = block.find("END IONS")
        if end_idx == -1:
            continue
        where = f"{path}: spectrum {block_no}"
        content = block[:end_idx].strip()
        lines = content.split("\n")

        rt = 0.0
        precursor = None
        mz_list: list[float] = []
        intensity_list: list[float] = []

        for line in lines:
            line = line.strip()
            if not line or line.startswith(_MGF_COMMENT_PREFIXES):
                continue
            if line.startswith("PEPMASS="):
                fields = line.partition("=")[2].split()
                if not fields:
                    raise ValueError(f"{where}: PEPMASS has no value ({line!r})")
                precursor = _mgf_float(fields[0], "PEPMASS", where, line)
            elif line.startswith("RTINSECONDS="):
                rt = _mgf_float(line.partition("=")[2].strip(), "RTINSECONDS", where, line)
            elif line.startswith("SCANS="):
                scan_idx = _mgf_int(line.partition("=")[2].strip(), "SCANS", where, line)
            elif "=" in line:
                continue  # any other header (TITLE=, CHARGE=, ...)
            else:
                fields = line.split()
                if len(fields) < 2:
                    raise ValueError(f"{where}: malformed peak line {line!r}: expected 'm/z intensity'")
                try:
                    mz, intensity = float(fields[0]), float(fields[1])
                except ValueError:
                    raise ValueError(f"{where}: malformed peak line {line!r}: expected 'm/z intensity'") from None
                if not (math.isfinite(mz) and math.isfinite(intensity)):
                    raise ValueError(f"{where}: peak line {line!r} has a non-finite m/z or intensity")
                mz_list.append(mz)
                intensity_list.append(intensity)

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


def _mgf_float(text: str, key: str, where: str, line: str) -> float:
    """Parse the numeric value of an MGF header line, naming the spectrum and line on failure."""
    try:
        return float(text)
    except ValueError:
        raise ValueError(f"{where}: {key} must be a number, got {text!r} ({line!r})") from None


def _mgf_int(text: str, key: str, where: str, line: str) -> int:
    """Parse the integer value of an MGF header line, naming the spectrum and line on failure."""
    try:
        return int(text)
    except ValueError:
        raise ValueError(f"{where}: {key} must be an integer, got {text!r} ({line!r})") from None


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
