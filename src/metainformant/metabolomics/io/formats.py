"""Mass spectrometry file format reading and writing.

Reads and writes two formats: CSV metabolite × sample intensity tables
(:func:`read_csv`, :func:`write_csv`) and MGF (Mascot Generic Format) spectra
(:func:`read_mgf`, :func:`write_mgf`). It also filters spectra and extracts ion
chromatograms from them. mzML and mzXML are not supported; convert such files to
MGF first.

The readers never invent values and never drop data silently. A missing intensity stays
missing (NaN), a missing MGF retention time stays NaN, and a CSV cell, an MGF header or
peak line that cannot be interpreted, or an MGF ``BEGIN IONS`` block that is not closed
by ``END IONS`` (a truncated file), raises :class:`ValueError` naming where it is.
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


def _as_real_number(value: Any) -> float:
    """``value`` as a float if it is a real number (not text, not None), else NaN."""
    if isinstance(value, (str, bytes)):
        return math.nan
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def _nonnegative_ppm_tolerance(ppm_tolerance: float) -> float:
    """Return ``ppm_tolerance`` as a float, rejecting values that cannot be a mass tolerance.

    A tolerance must be a finite number of ppm that is zero or positive. Zero is meaningful (only a peak
    at exactly the target m/z is inside the window); a negative or NaN tolerance would match nothing and an
    infinite one everything, each without saying why, so they raise instead.

    Raises:
        ValueError: If ``ppm_tolerance`` is not a number (text included), or is NaN, infinite or negative.
    """
    tolerance = _as_real_number(ppm_tolerance)
    if not math.isfinite(tolerance) or tolerance < 0:
        raise ValueError(f"ppm_tolerance must be a finite number >= 0, got {ppm_tolerance!r}")
    return tolerance


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
        retention_time: RT in seconds; NaN when the source gives none (an MGF block without
            ``RTINSECONDS``), never a made-up ``0.0``.
        ms_level: 1 for MS1, 2 for MS2, etc.
        mz_array: 1D array of m/z values.
        intensity_array: 1D array of intensities.
        precursor_mz: Precursor m/z for MS2+ spectra (None for MS1).
        total_ion_current: Sum of all intensities. When it is left out (or NaN) it is computed from
            ``intensity_array``; a value given explicitly is kept as given.
    """

    scan_number: int
    retention_time: float
    ms_level: int
    mz_array: np.ndarray
    intensity_array: np.ndarray
    precursor_mz: float | None = None
    total_ion_current: float = math.nan

    def __post_init__(self) -> None:
        if math.isnan(self.total_ion_current):
            self.total_ion_current = float(np.sum(self.intensity_array))


def read_mgf(filepath: str | Path) -> list[MassSpectrum]:
    """Read spectra from MGF (Mascot Generic Format) file.

    MGF is a simple text format commonly used for MS2 spectra in metabolomics
    and proteomics. Every spectrum is returned with ``ms_level=2``.

    Inside a ``BEGIN IONS`` / ``END IONS`` block a line is interpreted by its shape:

    * blank lines and comment lines (starting with ``#``, ``;``, ``!`` or ``/``) are skipped;
    * a line containing ``=`` is a header (``PEPMASS``, ``RTINSECONDS`` and ``SCANS`` are read,
      every other header such as ``TITLE`` or ``CHARGE`` is ignored);
    * every other line is a peak and must start with two finite numbers, ``m/z intensity``
      (further columns, such as a charge, are ignored). A peak line is never skipped silently, so
      ``.5 100`` and ``+100.0 50`` are read and a malformed peak line is an error.

    ``BEGIN IONS`` and ``END IONS`` are recognised as whole lines. Text outside a block (a preamble of
    global parameters, or anything between two blocks) is ignored, but a block must be closed: a
    ``BEGIN IONS`` that is not followed by ``END IONS`` before the next ``BEGIN IONS`` or the end of the
    file (the signature of a truncated file) is an error, as is an ``END IONS`` outside a block, because
    either would otherwise lose a spectrum without a trace. Nothing is returned for a file that raises.

    A spectrum without ``RTINSECONDS`` has a NaN retention time, not ``0.0``; one without ``PEPMASS`` has
    ``precursor_mz=None``. A spectrum without ``SCANS`` gets the previous spectrum's scan number plus one
    (``0`` for the first spectrum).

    Args:
        filepath: Path to MGF file.

    Returns:
        List of MassSpectrum objects.

    Raises:
        ValueError: If a block is not closed by ``END IONS`` (truncated file), an ``END IONS`` has no
            ``BEGIN IONS``, ``PEPMASS`` has no value, a ``PEPMASS``, ``RTINSECONDS`` or ``SCANS`` value is
            not a number (``SCANS`` an integer), or a peak line is malformed. The message names the file,
            the spectrum (counted from 1 in file order, with its ``TITLE`` when it has one) and the
            offending line.
    """
    path = Path(filepath)
    spectra: list[MassSpectrum] = []
    scan_idx = 0

    block_no = 0  # BEGIN IONS lines seen so far: the number of the spectrum being read
    block_line = 0  # line number of the open block's BEGIN IONS (0 while outside a block)
    where = ""
    title = ""
    rt = math.nan
    precursor: float | None = None
    mz_list: list[float] = []
    intensity_list: list[float] = []

    for line_no, raw_line in enumerate(path.read_text().split("\n"), start=1):
        line = raw_line.strip()

        if line == "BEGIN IONS":
            if block_line:
                raise ValueError(
                    f"{where}: BEGIN IONS (line {block_line}){_mgf_title_note(title)} is not closed by END IONS "
                    f"before the next BEGIN IONS (line {line_no}); the file is truncated or damaged"
                )
            block_no += 1
            block_line = line_no
            where = f"{path}: spectrum {block_no}"
            title = ""
            rt = math.nan
            precursor = None
            mz_list = []
            intensity_list = []
            continue

        if line == "END IONS":
            if not block_line:
                raise ValueError(f"{path}: line {line_no}: END IONS without a matching BEGIN IONS")
            int_arr = np.array(intensity_list, dtype=np.float64)
            spectra.append(
                MassSpectrum(
                    scan_number=scan_idx,
                    retention_time=rt,
                    ms_level=2,
                    mz_array=np.array(mz_list, dtype=np.float64),
                    intensity_array=int_arr,
                    precursor_mz=precursor,
                    total_ion_current=float(int_arr.sum()),
                )
            )
            scan_idx += 1
            block_line = 0
            continue

        if not block_line or not line or line.startswith(_MGF_COMMENT_PREFIXES):
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
            if line.startswith("TITLE="):
                title = line.partition("=")[2].strip()  # kept only to name the spectrum in an error
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

    if block_line:
        raise ValueError(
            f"{where}: BEGIN IONS (line {block_line}){_mgf_title_note(title)} is not closed by END IONS before "
            "the end of the file; the file is truncated or damaged"
        )

    return spectra


def _mgf_title_note(title: str) -> str:
    """`` with TITLE='...'`` for an error message about a spectrum that has a title, else nothing."""
    return f" with TITLE={title!r}" if title else ""


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

    Peaks are written with the shortest decimal text that reads back as the same number (m/z to six
    decimals), so an intensity such as ``0.0123`` survives a write and read. ``PEPMASS`` is left out when
    ``precursor_mz`` is None and ``RTINSECONDS`` when the retention time is NaN, so :func:`read_mgf` gives
    back the same None and NaN.

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
        if not math.isnan(spec.retention_time):
            lines.append(f"RTINSECONDS={spec.retention_time}")
        lines.append(f"SCANS={spec.scan_number}")
        for mz, intensity in zip(spec.mz_array, spec.intensity_array):
            lines.append(f"{float(mz):.6f} {float(intensity)!r}")
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
        rt_range: (min_rt, max_rt) in seconds (None = all). A spectrum whose retention time is unknown
            (NaN) cannot be shown to lie inside the range, so it is dropped when a range is given.

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

    Each MS1 scan contributes one point: the intensity of its most intense peak inside the window
    ``target_mz +/- target_mz * ppm_tolerance / 1e6`` (the maximum, not the sum), or ``0.0`` when no peak
    falls inside it.

    Args:
        spectra: List of MS1 spectra ordered by retention time.
        target_mz: Target m/z value; a positive, finite number.
        ppm_tolerance: Half-width of the m/z window in ppm of ``target_mz``; a finite number >= 0. Zero is
            allowed and keeps only peaks at exactly ``target_mz``.

    Returns:
        Tuple of (retention_times, intensities) arrays.

    Raises:
        ValueError: If ``ppm_tolerance`` is not a finite number >= 0 or ``target_mz`` is not a positive,
            finite number. A negative or NaN window would otherwise match nothing and return an all-zero
            chromatogram that looks like an ion that is simply absent.
    """
    ppm_tolerance = _nonnegative_ppm_tolerance(ppm_tolerance)
    target = _as_real_number(target_mz)
    if not (math.isfinite(target) and target > 0):
        raise ValueError(f"target_mz must be a positive finite number, got {target_mz!r}")
    target_mz = target

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
