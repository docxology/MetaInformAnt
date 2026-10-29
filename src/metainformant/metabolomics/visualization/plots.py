"""Metabolomics visualization utilities.

Provides volcano plots, PCA ordination, and concentration heatmaps for
metabolomics datasets.
"""

from __future__ import annotations

import operator
from dataclasses import dataclass

import numpy as np

from metainformant.metabolomics.io.formats import _nonnegative_ppm_tolerance


@dataclass
class VolcanoPoint:
    """A single point on a volcano plot.

    Attributes:
        name: Metabolite name.
        log2_fc: Log2 fold change.
        neg_log10_p: -log10(p-value).
        significant: Whether the point passes threshold.
    """

    name: str
    log2_fc: float
    neg_log10_p: float
    significant: bool


def volcano_plot_data(
    metabolite_names: list[str],
    log2_fold_changes: np.ndarray,
    p_values: np.ndarray,
    fc_threshold: float = 1.0,
    p_threshold: float = 0.05,
) -> list[VolcanoPoint]:
    """Prepare data for a volcano plot.

    Classifies each metabolite as significant or not based on fold change
    and p-value thresholds.

    Args:
        metabolite_names: List of metabolite names.
        log2_fold_changes: 1D array of log2 fold changes.
        p_values: 1D array of p-values, each in [0, 1] or NaN. A NaN p-value (a test that could not be
            run) gives a NaN ``neg_log10_p`` and is never significant. A p-value of exactly 0 is plotted
            at 300 (``-log10(1e-300)``) so that it stays finite.
        fc_threshold: Minimum absolute log2 FC for significance.
        p_threshold: Maximum p-value for significance.

    Returns:
        List of VolcanoPoint for plotting.

    Raises:
        ValueError: If the names and the two arrays differ in length, or a p-value lies outside [0, 1].
            Such a value used to be clipped silently into range, so a negative p-value was plotted as
            the most significant metabolite.
    """
    p_array = np.asarray(p_values, dtype=float)
    fc_array = np.asarray(log2_fold_changes, dtype=float)
    if p_array.ndim != 1 or fc_array.ndim != 1 or not len(metabolite_names) == len(p_array) == len(fc_array):
        raise ValueError(
            f"metabolite_names ({len(metabolite_names)}), log2_fold_changes ({fc_array.size}) and "
            f"p_values ({p_array.size}) must be 1D and of the same length"
        )
    out_of_range = ~np.isnan(p_array) & ((p_array < 0.0) | (p_array > 1.0))
    if out_of_range.any():
        first = int(np.argmax(out_of_range))
        raise ValueError(
            f"p_values must lie in [0, 1]; {metabolite_names[first]!r} has {float(p_array[first])!r} "
            f"({int(out_of_range.sum())} value(s) out of range)"
        )

    neg_log10_p = -np.log10(np.clip(p_array, 1e-300, 1.0))
    points = []
    for i, name in enumerate(metabolite_names):
        lfc = float(fc_array[i])
        nlp = float(neg_log10_p[i])
        sig = bool(abs(lfc) >= fc_threshold and p_array[i] <= p_threshold)
        points.append(VolcanoPoint(name=name, log2_fc=lfc, neg_log10_p=nlp, significant=sig))
    return points


def pca_metabolomics(
    intensities: np.ndarray,
    n_components: int = 2,
) -> tuple[np.ndarray, np.ndarray]:
    """Perform PCA on metabolomics intensity data.

    Uses SVD-based approach without external ML dependencies. Centers data
    before decomposition.

    Args:
        intensities: 2D array (metabolites × samples), finite. Transposed internally
            so samples are rows for PCA.
        n_components: Number of principal components to return, from 1 to
            min(n_samples, n_metabolites).

    Returns:
        Tuple of (scores, explained_variance_ratio).
        - scores: (n_samples, n_components) projected coordinates.
        - explained_variance_ratio: fraction of variance per component (zeros for data with no variance).

    Raises:
        ValueError: If ``n_components`` is outside 1..min(n_samples, n_metabolites) (fewer components than
            requested used to be returned silently, so the shapes did not match the documentation) or the
            data contain NaN or infinite values (impute them first, for example with
            :func:`~metainformant.metabolomics.analysis.identification.missing_value_imputation`).
    """
    X = np.asarray(intensities, dtype=float).T
    if not np.all(np.isfinite(X)):
        raise ValueError(
            "intensities must be finite for PCA; impute missing values first "
            "(see metainformant.metabolomics.analysis.identification.missing_value_imputation)"
        )
    max_components = min(X.shape)
    if not 1 <= n_components <= max_components:
        raise ValueError(
            f"n_components must be between 1 and min(n_samples, n_metabolites) = {max_components}, got {n_components}"
        )
    X_centered = X - X.mean(axis=0)

    U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)

    scores = U[:, :n_components] * S[:n_components]

    total_var = np.sum(S**2)
    explained = (S[:n_components] ** 2) / total_var if total_var > 0 else np.zeros(n_components)

    return scores, explained


def intensity_heatmap_data(
    intensities: np.ndarray,
    metabolite_names: list[str],
    sample_names: list[str],
    normalize: bool = True,
) -> dict:
    """Prepare data for an intensity heatmap.

    Optionally z-score normalizes rows (metabolites) for better visualization.

    Args:
        intensities: 2D array (metabolites × samples).
        metabolite_names: Row labels.
        sample_names: Column labels.
        normalize: If True, z-score normalize each row.

    Returns:
        Dict with keys 'matrix', 'row_labels', 'col_labels'.
    """
    data = intensities.copy().astype(float)
    if normalize:
        row_mean = data.mean(axis=1, keepdims=True)
        row_std = data.std(axis=1, keepdims=True)
        row_std = np.where(row_std > 0, row_std, 1.0)
        data = (data - row_mean) / row_std

    return {
        "matrix": data,
        "row_labels": metabolite_names,
        "col_labels": sample_names,
    }


@dataclass
class ChromatogramPeak:
    """A detected peak in a chromatogram.

    Attributes:
        peak_idx: Index of the peak apex in the RT array.
        retention_time: RT at apex (seconds).
        intensity: Intensity at apex.
        left_idx: Left boundary index.
        right_idx: Right boundary index.
        area: Integrated area under the peak.
    """

    peak_idx: int
    retention_time: float
    intensity: float
    left_idx: int
    right_idx: int
    area: float


def _trapezoid_area(y: np.ndarray, x: np.ndarray) -> float:
    """Area under ``y(x)`` by the trapezoidal rule.

    Spelled out because NumPy renamed ``np.trapz`` to ``np.trapezoid`` in 2.0 and removed the old name
    in later releases, while this package supports NumPy >= 1.24: neither name exists on every
    supported version. The arithmetic is exactly NumPy's own.
    """
    return float(np.sum(np.diff(x) * (y[1:] + y[:-1]) / 2.0))


def detect_chromatographic_peaks(
    retention_times: np.ndarray,
    intensities: np.ndarray,
    min_intensity: float = 0.0,
    min_peak_width: int = 3,
) -> list[ChromatogramPeak]:
    """Detect peaks in an extracted ion chromatogram.

    Uses a simple local-maximum approach with boundary detection.

    Args:
        retention_times: 1D array of retention times.
        intensities: 1D array of intensities.
        min_intensity: Minimum intensity for a peak apex.
        min_peak_width: Minimum number of points for a valid peak.

    Returns:
        List of ChromatogramPeak objects.
    """
    n = len(intensities)
    if n < 3:
        return []

    peaks: list[ChromatogramPeak] = []

    for i in range(1, n - 1):
        if intensities[i] <= min_intensity:
            continue
        if intensities[i] <= intensities[i - 1] or intensities[i] <= intensities[i + 1]:
            continue

        # Find left boundary
        left = i - 1
        while left > 0 and intensities[left] > intensities[left - 1]:
            left -= 1

        # Find right boundary
        right = i + 1
        while right < n - 1 and intensities[right] > intensities[right + 1]:
            right += 1

        width = right - left + 1
        if width < min_peak_width:
            continue

        # Trapezoidal integration for area
        area = _trapezoid_area(intensities[left : right + 1], retention_times[left : right + 1])

        peaks.append(
            ChromatogramPeak(
                peak_idx=i,
                retention_time=float(retention_times[i]),
                intensity=float(intensities[i]),
                left_idx=left,
                right_idx=right,
                area=area,
            )
        )

    return peaks


@dataclass
class MassSpectrumPlotData:
    """Data prepared for plotting a mass spectrum.

    Attributes:
        mz_values: m/z values for stick plot.
        intensities: Relative intensities (0-100 scale).
        annotations: Dict mapping m/z to annotation text.
        base_peak_mz: m/z of the most intense peak; NaN for an empty spectrum, which has none.
        base_peak_intensity: Raw intensity of the base peak; NaN for an empty spectrum.
    """

    mz_values: np.ndarray
    intensities: np.ndarray
    annotations: dict[float, str]
    base_peak_mz: float
    base_peak_intensity: float


def prepare_spectrum_plot(
    mz_array: np.ndarray,
    intensity_array: np.ndarray,
    top_n: int = 20,
    annotation_db: dict[float, str] | None = None,
    ppm_tolerance: float = 10.0,
) -> MassSpectrumPlotData:
    """Prepare mass spectrum data for visualization.

    Normalizes intensities to 0-100 relative to base peak and optionally
    annotates top peaks.

    Args:
        mz_array: 1D array of m/z values.
        intensity_array: 1D array of intensities, the same length as ``mz_array``.
        top_n: Number of top peaks to include; a positive integer. Of equally intense peaks the one
            with the higher m/z is kept first.
        annotation_db: Optional dict of m/z -> annotation label. A peak takes the label of the first entry
            (in the dict's order) that lies within ``ppm_tolerance`` of it.
        ppm_tolerance: Tolerance for annotation matching in ppm of the peak's m/z; a finite number >= 0.
            Zero annotates only a peak whose m/z equals a database m/z exactly.

    Returns:
        MassSpectrumPlotData ready for plotting.

    Raises:
        ValueError: If ``ppm_tolerance`` is not a finite number >= 0 (it is checked even when no
            ``annotation_db`` is given, and for an empty spectrum), ``top_n`` is not a positive integer
            (``top_n=0`` used to return every peak and a negative one silently dropped the weakest), or the
            two arrays are not 1D arrays of the same length.
    """
    ppm_tolerance = _nonnegative_ppm_tolerance(ppm_tolerance)
    try:
        top_n = operator.index(top_n)
    except TypeError:
        raise ValueError(f"top_n must be a positive integer, got {top_n!r}") from None
    if top_n < 1:
        raise ValueError(f"top_n must be a positive integer, got {top_n!r}")
    mz_array = np.asarray(mz_array, dtype=float)
    intensity_array = np.asarray(intensity_array, dtype=float)
    if mz_array.ndim != 1 or mz_array.shape != intensity_array.shape:
        raise ValueError(
            f"mz_array and intensity_array must be 1D arrays of the same length, got shapes "
            f"{mz_array.shape} and {intensity_array.shape}"
        )

    if len(mz_array) == 0:
        return MassSpectrumPlotData(
            mz_values=np.array([]),
            intensities=np.array([]),
            annotations={},
            base_peak_mz=float("nan"),
            base_peak_intensity=float("nan"),
        )

    base_peak_idx = int(np.argmax(intensity_array))
    base_peak_int = float(intensity_array[base_peak_idx])
    base_peak_mz = float(mz_array[base_peak_idx])

    # Relative intensities
    rel_int = intensity_array / base_peak_int * 100.0 if base_peak_int > 0 else intensity_array

    # Select top N peaks (a stable sort, so equally intense peaks are chosen the same way on every platform)
    if len(mz_array) > top_n:
        top_indices = np.argsort(rel_int, kind="stable")[-top_n:]
        top_indices = np.sort(top_indices)  # sort by m/z order
    else:
        top_indices = np.arange(len(mz_array))

    mz_out = mz_array[top_indices]
    int_out = rel_int[top_indices]

    # Annotations
    annotations: dict[float, str] = {}
    if annotation_db is not None:
        for mz in mz_out:
            window = mz * ppm_tolerance / 1e6
            for db_mz, label in annotation_db.items():
                if abs(mz - db_mz) <= window:
                    annotations[float(mz)] = label
                    break

    return MassSpectrumPlotData(
        mz_values=mz_out,
        intensities=int_out,
        annotations=annotations,
        base_peak_mz=base_peak_mz,
        base_peak_intensity=base_peak_int,
    )


def retention_time_alignment_data(
    reference_rts: np.ndarray,
    sample_rts: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Compute retention time deviation between sample and reference.

    Uses a simple polynomial fit for RT correction assessment.

    Args:
        reference_rts: Reference retention times.
        sample_rts: Corresponding sample retention times.

    Returns:
        Tuple of (deviations, corrected_rts, rmse).
        - deviations: RT differences (sample - reference).
        - corrected_rts: Sample RTs after polynomial correction.
        - rmse: Root mean square error after correction.
    """
    deviations = sample_rts - reference_rts

    # Fit quadratic polynomial for RT correction
    if len(reference_rts) >= 3:
        coeffs = np.polyfit(sample_rts, deviations, deg=min(2, len(reference_rts) - 1))
        predicted_dev = np.polyval(coeffs, sample_rts)
        corrected = sample_rts - predicted_dev
    else:
        # Simple shift correction
        mean_dev = np.mean(deviations)
        corrected = sample_rts - mean_dev

    residuals = corrected - reference_rts
    rmse = float(np.sqrt(np.mean(residuals**2)))

    return deviations, corrected, rmse
