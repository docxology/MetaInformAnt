"""Statistical plotting functions for data analysis and visualization.

This module provides statistical visualization functions including histograms,
box plots, violin plots, Q-Q plots, correlation heatmaps, density plots,
ridge plots, ROC curves, precision-recall curves, residual plots, and leverage plots.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, cast

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes

from metainformant.core.data import validation
from metainformant.core.io import paths
from metainformant.core.utils import logging
from metainformant.visualization.config.conventions import save_figure_deterministic
from metainformant.visualization.config.figures import root_figure

logger = logging.get_logger(__name__)

# Optional imports with graceful fallbacks
try:
    from scipy import stats

    HAS_SCIPY = True
except ImportError:
    stats = cast("Any", None)
    HAS_SCIPY = False

try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    sns = None
    HAS_SEABORN = False

try:
    import sklearn.metrics as metrics

    HAS_SKLEARN = True
except ImportError:
    metrics = None
    HAS_SKLEARN = False


def _save_plot(ax: Axes, output_path: str | Path, label: str) -> str:
    """Ensure the output directory exists and save the plotted figure deterministically.

    Consolidates the repeated ensure_directory / save_figure_deterministic / logger
    triple used by every plot function in this module. Saves ``ax.figure`` so the
    plotted figure is written even when a different pyplot figure is current.
    """
    paths.ensure_directory(Path(output_path).parent)
    save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
    logger.info(f"{label} saved to {output_path}")
    return str(output_path)


def histogram(
    data: np.ndarray, *, bins: int = 30, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
) -> Axes:
    """Create a histogram.

    Args:
        data: Data to plot
        bins: Number of bins
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib hist().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data is empty
    """
    validation.validate_type(data, np.ndarray, "data")
    if len(data) == 0:
        raise ValueError("Data array cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    ax.hist(data, bins=bins, **kwargs)

    if output_path:
        _save_plot(ax, output_path, "Histogram")

    return ax


def box_plot(
    data: list[np.ndarray], *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
) -> Axes:
    """Create a box plot.

    Args:
        data: List of arrays to plot as box plots
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib boxplot().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data is empty or contains invalid arrays
    """
    validation.validate_type(data, list, "data")
    if not data:
        raise ValueError("Data list cannot be empty")

    for i, arr in enumerate(data):
        validation.validate_type(arr, np.ndarray, f"data[{i}]")
        if len(arr) == 0:
            raise ValueError(f"Data array at index {i} cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    ax.boxplot(data, **kwargs)

    if output_path:
        _save_plot(ax, output_path, "Box plot")

    return ax


def violin_plot(
    data: list[np.ndarray], *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
) -> Axes:
    """Create a violin plot.

    Args:
        data: List of arrays to plot as violin plots
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib violinplot().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data is empty or contains invalid arrays
    """
    validation.validate_type(data, list, "data")
    if not data:
        raise ValueError("Data list cannot be empty")

    for i, arr in enumerate(data):
        validation.validate_type(arr, np.ndarray, f"data[{i}]")
        if len(arr) == 0:
            raise ValueError(f"Data array at index {i} cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    # Use seaborn if available, otherwise matplotlib
    if HAS_SEABORN:
        # Convert to DataFrame for seaborn
        df_data = []
        for i, arr in enumerate(data):
            df_data.extend([(i + 1, val) for val in arr])
        df = pd.DataFrame(df_data, columns=["group", "value"])
        sns.violinplot(data=df, x="group", y="value", ax=ax, **kwargs)
    else:
        # Fallback to matplotlib boxplot with warning
        logger.warning("Seaborn not available, using boxplot instead of violinplot")
        ax.boxplot(data, **kwargs)

    if output_path:
        _save_plot(ax, output_path, "Violin plot")

    return ax


def qq_plot(
    data: np.ndarray,
    *,
    distribution: str = "norm",
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Create a Q-Q plot against a normal or uniform distribution fitted to the data.

    The i-th smallest of n observations is plotted against the quantile of the fitted distribution at the
    plotting position ``(i - 0.5) / n``, so a sample from that distribution lies along the dashed ``y = x`` line.
    The fitted distribution has the sample mean and standard deviation ("norm") or the sample range ("uniform").

    Args:
        data: Data to plot
        distribution: Distribution to compare against: "norm" or "uniform"
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data is empty or constant, or distribution is unsupported
    """
    validation.validate_type(data, np.ndarray, "data")
    if len(data) == 0:
        raise ValueError("Data array cannot be empty")

    if distribution not in ("norm", "uniform"):
        raise ValueError(f"Unsupported distribution: {distribution}")

    if not HAS_SCIPY:
        raise ImportError("scipy required for Q-Q plot")

    if np.max(data) == np.min(data):
        raise ValueError("Q-Q plot needs data that is not constant: a constant sample has no fitted distribution")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    # Plotting positions (i - 0.5) / n: strictly inside (0, 1) and symmetric, so the tails are not squeezed
    data_sorted = np.sort(data)
    n = len(data)
    probabilities = (np.arange(1, n + 1) - 0.5) / n

    # Get distribution quantiles
    if distribution == "norm":
        sample_quantiles = stats.norm.ppf(probabilities, loc=np.mean(data), scale=np.std(data))
    else:
        sample_quantiles = stats.uniform.ppf(probabilities, loc=np.min(data), scale=np.max(data) - np.min(data))

    ax.scatter(sample_quantiles, data_sorted, **kwargs)
    ax.plot(
        [np.min(sample_quantiles), np.max(sample_quantiles)],
        [np.min(sample_quantiles), np.max(sample_quantiles)],
        "r--",
        alpha=0.7,
    )
    ax.set_xlabel(f"Theoretical {distribution.title()} Quantiles")
    ax.set_ylabel("Sample Quantiles")

    if output_path:
        _save_plot(ax, output_path, "Q-Q plot")

    return ax


def correlation_heatmap(
    data: pd.DataFrame, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
) -> Axes:
    """Create a correlation heatmap.

    Args:
        data: DataFrame with numeric columns
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to seaborn heatmap().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data contains no numeric columns
    """
    validation.validate_type(data, pd.DataFrame, "data")

    # Select only numeric columns
    numeric_data = data.select_dtypes(include=[np.number])
    if numeric_data.empty:
        raise ValueError("DataFrame must contain numeric columns for correlation heatmap")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 8)))

    # Calculate correlation matrix
    corr_matrix = numeric_data.corr()

    # Use seaborn if available, otherwise matplotlib
    if HAS_SEABORN:
        # Extract specific kwargs to avoid duplicate keyword argument error
        annot = kwargs.pop("annot", True)
        cmap = kwargs.pop("cmap", "coolwarm")
        sns.heatmap(corr_matrix, ax=ax, annot=annot, cmap=cmap, **kwargs)
    else:
        # Fallback to matplotlib imshow
        logger.warning("Seaborn not available, using basic heatmap")
        # Remove seaborn-specific kwargs that imshow doesn't understand
        cmap = kwargs.pop("cmap", "coolwarm")
        kwargs.pop("annot", None)  # Remove annot if present
        im = ax.imshow(corr_matrix.values, cmap=cmap, **kwargs)
        plt.colorbar(im, ax=ax)

        # Add text annotations
        for i in range(len(corr_matrix)):
            for j in range(len(corr_matrix)):
                ax.text(j, i, f"{corr_matrix.iloc[i, j]:.2f}", ha="center", va="center", color="w")

    if output_path:
        _save_plot(ax, output_path, "Correlation heatmap")

    return ax


def density_plot(
    data: np.ndarray,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Create a kernel density estimate plot.

    Args:
        data: Data to plot
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib or seaborn.

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data is empty
    """
    validation.validate_type(data, np.ndarray, "data")
    if len(data) == 0:
        raise ValueError("Data array cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    # Use seaborn if available, otherwise matplotlib histogram
    if HAS_SEABORN:
        sns.kdeplot(data=data, ax=ax, **kwargs)
    else:
        logger.warning("Seaborn not available, using histogram instead of density plot")
        ax.hist(data, density=True, alpha=0.7, **kwargs)

    if output_path:
        _save_plot(ax, output_path, "Density plot")

    return ax


def ridge_plot(
    data: list[np.ndarray], *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
) -> Axes:
    """Create a ridge plot (overlapping density plots).

    Args:
        data: List of arrays to plot as ridge plot
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to plotting functions.

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data is empty or contains invalid arrays
    """
    validation.validate_type(data, list, "data")
    if not data:
        raise ValueError("Data list cannot be empty")

    for i, arr in enumerate(data):
        validation.validate_type(arr, np.ndarray, f"data[{i}]")
        if len(arr) == 0:
            raise ValueError(f"Data array at index {i} cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 6)))

    # Simple ridge plot implementation

    for i, arr in enumerate(data):
        if HAS_SEABORN:
            # Use seaborn for better KDE
            sns.kdeplot(data=arr, ax=ax, fill=True, alpha=0.7, label=f"Group {i+1}", **kwargs)
        else:
            # Fallback to histogram
            ax.hist(arr, density=True, alpha=0.5, bins=20, label=f"Group {i+1}", **kwargs)

    ax.legend()

    if output_path:
        _save_plot(ax, output_path, "Ridge plot")

    return ax


def roc_curve(
    y_true: np.ndarray,
    y_scores: np.ndarray,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Create a ROC curve plot.

    Args:
        y_true: True binary labels
        y_scores: Target scores (probabilities or confidence values)
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib plot().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If arrays have mismatched lengths or invalid values
    """
    validation.validate_type(y_true, np.ndarray, "y_true")
    validation.validate_type(y_scores, np.ndarray, "y_scores")

    if len(y_true) != len(y_scores):
        raise ValueError("y_true and y_scores must have same length")

    if metrics is None:
        raise ImportError("scikit-learn required for ROC curve")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    fpr, tpr, _ = metrics.roc_curve(y_true, y_scores)
    roc_auc = metrics.roc_auc_score(y_true, y_scores)

    ax.plot(fpr, tpr, label=f"ROC curve (AUC = {roc_auc:.2f})", **kwargs)
    ax.plot([0, 1], [0, 1], "k--", alpha=0.7)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC Curve")
    ax.legend()

    if output_path:
        _save_plot(ax, output_path, "ROC curve")

    return ax


def precision_recall_curve(
    y_true: np.ndarray,
    y_scores: np.ndarray,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Create a precision-recall curve plot.

    Args:
        y_true: True binary labels
        y_scores: Target scores (probabilities or confidence values)
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib plot().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If arrays have mismatched lengths or invalid values
    """
    validation.validate_type(y_true, np.ndarray, "y_true")
    validation.validate_type(y_scores, np.ndarray, "y_scores")

    if len(y_true) != len(y_scores):
        raise ValueError("y_true and y_scores must have same length")

    if metrics is None:
        raise ImportError("scikit-learn required for precision-recall curve")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    precision, recall, _ = metrics.precision_recall_curve(y_true, y_scores)
    pr_auc = metrics.average_precision_score(y_true, y_scores)

    ax.plot(recall, precision, label=f"PR curve (AP = {pr_auc:.2f})", **kwargs)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Precision-Recall Curve")
    ax.legend()

    if output_path:
        _save_plot(ax, output_path, "Precision-recall curve")

    return ax


def residual_plot(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Create a residual plot.

    Args:
        y_true: True values
        y_pred: Predicted values
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If arrays have mismatched lengths
    """
    validation.validate_type(y_true, np.ndarray, "y_true")
    validation.validate_type(y_pred, np.ndarray, "y_pred")

    if len(y_true) != len(y_pred):
        raise ValueError("y_true and y_pred must have same length")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    residuals = y_true - y_pred

    ax.scatter(y_pred, residuals, **kwargs)
    ax.axhline(y=0, color="r", linestyle="--", alpha=0.7)
    ax.set_xlabel("Predicted Values")
    ax.set_ylabel("Residuals")
    ax.set_title("Residual Plot")

    if output_path:
        _save_plot(ax, output_path, "Residual plot")

    return ax


def leverage_plot(
    X: np.ndarray, y: np.ndarray, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
) -> Axes:
    """Create a residuals-versus-leverage plot for an ordinary least squares fit with an intercept.

    The leverage of observation i is the i-th diagonal element ``h_ii`` of the hat matrix
    ``H = D (D'D)^+ D'`` of the design ``D = [1, X]`` (the values sum to the number of fitted parameters, and a
    simple regression has ``h_ii = 1/n + (x_i - mean(x))^2 / Sxx``). The vertical axis is the internally
    studentized residual ``e_i / sqrt(s^2 (1 - h_ii))`` with ``s^2 = RSS / (n - p)``.

    Args:
        X: Feature matrix of shape (n_samples, n_features), or a vector for one feature
        y: Target values
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If X and y have different numbers of samples, or there are not more samples than fitted
            parameters
    """
    validation.validate_type(X, np.ndarray, "X")
    validation.validate_type(y, np.ndarray, "y")

    if len(X) != len(y):
        raise ValueError("X and y must have same number of samples")

    features = np.asarray(X, dtype=float)
    if features.ndim == 1:
        features = features.reshape(-1, 1)
    target = np.asarray(y, dtype=float)
    n_samples = len(target)

    design = np.column_stack([np.ones(n_samples), features])
    left, singular_values, _ = np.linalg.svd(design, full_matrices=False)
    rank = int(np.sum(singular_values > singular_values.max() * max(design.shape) * np.finfo(float).eps))
    if n_samples <= rank:
        raise ValueError(f"Leverage plot needs more samples ({n_samples}) than fitted parameters ({rank})")

    # Hat values: squared row norms of the orthonormal basis of the design's column space
    leverage = np.sum(left[:, :rank] ** 2, axis=1)
    residuals = target - left[:, :rank] @ (left[:, :rank].T @ target)
    residual_variance = float(np.sum(residuals**2)) / (n_samples - rank)

    with np.errstate(divide="ignore", invalid="ignore"):
        standardized_residuals = residuals / np.sqrt(residual_variance * (1.0 - leverage))
    # An observation with leverage 1 fits itself exactly and has no studentized residual; it is not drawn
    standardized_residuals[~np.isfinite(standardized_residuals)] = np.nan

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    ax.scatter(leverage, standardized_residuals, **kwargs)
    ax.axhline(y=0, color="r", linestyle="--", alpha=0.7)
    ax.set_xlabel("Leverage")
    ax.set_ylabel("Standardized Residuals")
    ax.set_title("Leverage Plot")

    if output_path:
        _save_plot(ax, output_path, "Leverage plot")

    return ax
