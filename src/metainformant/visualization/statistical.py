"""Statistical plotting functions for data analysis.

This module provides statistical visualization functions including histograms,
box plots, violin plots, Q-Q plots, correlation heatmaps, density plots,
and statistical diagnostic plots.
"""

from __future__ import annotations

from typing import Sequence

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Use non-interactive backend by default for tests/headless
matplotlib.use("Agg", force=True)

try:
    import seaborn as sns

    SEABORN_AVAILABLE = True
except ImportError:
    SEABORN_AVAILABLE = False
    sns = None


def histogram(
    data: Sequence[float],
    *,
    bins: int | str | Sequence[float] = 30,
    ax: plt.Axes | None = None,
    density: bool = False,
    alpha: float = 0.7,
    color: str = "blue",
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
) -> plt.Axes:
    """Create a histogram.

    Args:
        data: Data to plot
        bins: Number of bins or bin edges
        ax: Matplotlib axes (creates new if None)
        density: Whether to normalize to density
        alpha: Bar transparency
        color: Bar color
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import histogram
        >>> import numpy as np
        >>> data = np.random.normal(0, 1, 1000)
        >>> ax = histogram(data, bins=30)
    """
    if ax is None:
        _, ax = plt.subplots()

    ax.hist(data, bins=bins, density=density, alpha=alpha, color=color)

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    return ax


def box_plot(
    data: Sequence[float] | Sequence[Sequence[float]],
    *,
    ax: plt.Axes | None = None,
    positions: Sequence[float] | None = None,
    labels: Sequence[str] | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
) -> plt.Axes:
    """Create a box plot.

    Args:
        data: Data to plot (single array or list of arrays)
        ax: Matplotlib axes (creates new if None)
        positions: Box positions
        labels: Box labels
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import box_plot
        >>> import numpy as np
        >>> data = [np.random.normal(0, 1, 100) for _ in range(3)]
        >>> ax = box_plot(data, labels=["A", "B", "C"])
    """
    if ax is None:
        _, ax = plt.subplots()

    bp = ax.boxplot(data, positions=positions, labels=labels)

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    return ax


def violin_plot(
    data: Sequence[float] | Sequence[Sequence[float]],
    *,
    ax: plt.Axes | None = None,
    positions: Sequence[float] | None = None,
    labels: Sequence[str] | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
) -> plt.Axes:
    """Create a violin plot.

    Args:
        data: Data to plot (single array or list of arrays)
        ax: Matplotlib axes (creates new if None)
        positions: Violin positions
        labels: Violin labels
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import violin_plot
        >>> import numpy as np
        >>> data = [np.random.normal(0, 1, 100) for _ in range(3)]
        >>> ax = violin_plot(data, labels=["A", "B", "C"])
    """
    if ax is None:
        _, ax = plt.subplots()

    if SEABORN_AVAILABLE:
        if isinstance(data[0], (list, tuple)) if data else False:
            # Multiple datasets
            for i, dataset in enumerate(data):
                pos = positions[i] if positions else i + 1
                sns.violinplot(data=dataset, ax=ax, position=pos)
        else:
            # Single dataset
            sns.violinplot(data=data, ax=ax)
    else:
        # Fallback to boxplot
        return box_plot(data, ax=ax, positions=positions, labels=labels,
                       xlabel=xlabel, ylabel=ylabel, title=title)

    if labels:
        ax.set_xticklabels(labels)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    return ax


def qq_plot(
    p_values: list[float],
    *,
    ax: plt.Axes | None = None,
    title: str = "Q-Q Plot",
    **kwargs
) -> plt.Axes:
    """Create a Q-Q plot for p-value distribution analysis.

    Args:
        p_values: List of p-values to plot
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments for scatter

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import qq_plot
        >>> import numpy as np
        >>> pvals = np.random.uniform(0, 1, 1000)
        >>> ax = qq_plot(pvals.tolist())
    """
    if ax is None:
        _, ax = plt.subplots()

    # Remove p-values of 0 (which would cause issues with log)
    p_values = [p for p in p_values if p > 0]

    if len(p_values) == 0:
        ax.text(0.5, 0.5, "No valid p-values", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title)
        return ax

    # Sort p-values
    p_sorted = np.sort(p_values)

    # Generate expected p-values under null hypothesis
    n = len(p_sorted)
    expected = np.array([(i + 0.5) / n for i in range(n)])

    # Create Q-Q plot
    ax.scatter(-np.log10(expected), -np.log10(p_sorted), alpha=0.6, **kwargs)

    # Add diagonal line (y = x)
    max_val = max(-np.log10(expected).max(), -np.log10(p_sorted).max())
    ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.8, label='Expected under null')

    ax.set_xlabel("Expected -log₁₀(P)")
    ax.set_ylabel("Observed -log₁₀(P)")
    ax.set_title(title)
    ax.legend()

    return ax


def correlation_heatmap(
    data: pd.DataFrame,
    *,
    method: str = "pearson",
    cmap: str = "RdBu_r",
    annot: bool = True,
    ax: plt.Axes | None = None,
    **kwargs
) -> plt.Axes:
    """Create a correlation heatmap.

    Args:
        data: DataFrame with numeric columns
        method: Correlation method ('pearson', 'spearman', 'kendall')
        cmap: Colormap for the heatmap
        annot: Whether to annotate cells with correlation values
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments for seaborn.heatmap

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import correlation_heatmap
        >>> import pandas as pd
        >>> import numpy as np
        >>> df = pd.DataFrame(np.random.random((10, 5)))
        >>> ax = correlation_heatmap(df)
    """
    if ax is None:
        _, ax = plt.subplots()

    if SEABORN_AVAILABLE:
        corr_matrix = data.corr(method=method)
        sns.heatmap(corr_matrix, cmap=cmap, annot=annot, ax=ax, **kwargs)
    else:
        # Fallback to matplotlib
        corr_matrix = data.corr(method=method)
        im = ax.imshow(corr_matrix, cmap=cmap)
        if annot:
            for i in range(len(corr_matrix)):
                for j in range(len(corr_matrix)):
                    ax.text(j, i, f"{corr_matrix.iloc[i, j]:.2f}",
                           ha="center", va="center", color="white" if abs(corr_matrix.iloc[i, j]) > 0.5 else "black")

    return ax


def density_plot(
    data: Sequence[float],
    *,
    ax: plt.Axes | None = None,
    fill: bool = True,
    alpha: float = 0.5,
    color: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
) -> plt.Axes:
    """Create a density plot (kernel density estimation).

    Args:
        data: Data to plot
        ax: Matplotlib axes (creates new if None)
        fill: Whether to fill under the curve
        alpha: Transparency
        color: Line/fill color
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import density_plot
        >>> import numpy as np
        >>> data = np.random.normal(0, 1, 1000)
        >>> ax = density_plot(data)
    """
    if ax is None:
        _, ax = plt.subplots()

    if SEABORN_AVAILABLE:
        sns.kdeplot(data=data, ax=ax, fill=fill, alpha=alpha, color=color)
    else:
        # Fallback to histogram
        from scipy import stats
        kde = stats.gaussian_kde(data)
        x_range = np.linspace(min(data), max(data), 200)
        y = kde(x_range)
        ax.plot(x_range, y, color=color)
        if fill:
            ax.fill_between(x_range, y, alpha=alpha, color=color)

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    return ax


def ridge_plot(
    data: list[Sequence[float]],
    labels: list[str] | None = None,
    *,
    ax: plt.Axes | None = None,
    overlap: float = 0.5,
    **kwargs
) -> plt.Axes:
    """Create a ridge plot (overlapping density plots).

    Args:
        data: List of data arrays
        labels: Labels for each distribution
        ax: Matplotlib axes (creates new if None)
        overlap: Overlap factor between plots (0-1)
        **kwargs: Additional arguments for density plots

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import ridge_plot
        >>> import numpy as np
        >>> data = [np.random.normal(i, 1, 100) for i in range(3)]
        >>> ax = ridge_plot(data, labels=["A", "B", "C"])
    """
    if ax is None:
        _, ax = plt.subplots()

    n = len(data)
    colors = plt.cm.viridis(np.linspace(0, 1, n))

    for i, (d, color) in enumerate(zip(data, colors)):
        y_offset = i * (1 - overlap)
        if SEABORN_AVAILABLE:
            kde = sns.kdeplot(data=d, ax=ax, fill=True, alpha=0.6, color=color)
            # Manually adjust y position
            for line in kde.lines:
                y_data = line.get_ydata()
                line.set_ydata(y_data + y_offset)
        else:
            from scipy import stats
            kde = stats.gaussian_kde(d)
            x_range = np.linspace(min(d), max(d), 200)
            y = kde(x_range) + y_offset
            ax.fill_between(x_range, y, y_offset, alpha=0.6, color=color)

        if labels and i < len(labels):
            ax.text(max(d), y_offset + 0.1, labels[i], ha="right", va="bottom")

    return ax


def roc_curve(
    y_true: Sequence[int],
    y_scores: Sequence[float],
    *,
    ax: plt.Axes | None = None,
    title: str = "ROC Curve",
    **kwargs
) -> plt.Axes:
    """Create a ROC curve plot.

    Args:
        y_true: True binary labels
        y_scores: Predicted scores/probabilities
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments for plot

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import roc_curve
        >>> import numpy as np
        >>> y_true = [0, 1, 0, 1, 1]
        >>> y_scores = [0.1, 0.9, 0.2, 0.8, 0.7]
        >>> ax = roc_curve(y_true, y_scores)
    """
    if ax is None:
        _, ax = plt.subplots()

    try:
        from sklearn.metrics import roc_curve, auc
        fpr, tpr, _ = roc_curve(y_true, y_scores)
        roc_auc = auc(fpr, tpr)

        ax.plot(fpr, tpr, label=f'ROC curve (AUC = {roc_auc:.2f})', **kwargs)
        ax.plot([0, 1], [0, 1], 'k--', label='Random classifier')
        ax.set_xlabel('False Positive Rate')
        ax.set_ylabel('True Positive Rate')
        ax.set_title(title)
        ax.legend()
    except ImportError:
        ax.text(0.5, 0.5, "scikit-learn required for ROC curve", ha="center", va="center", transform=ax.transAxes)

    return ax


def precision_recall_curve(
    y_true: Sequence[int],
    y_scores: Sequence[float],
    *,
    ax: plt.Axes | None = None,
    title: str = "Precision-Recall Curve",
    **kwargs
) -> plt.Axes:
    """Create a precision-recall curve plot.

    Args:
        y_true: True binary labels
        y_scores: Predicted scores/probabilities
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments for plot

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import precision_recall_curve
        >>> import numpy as np
        >>> y_true = [0, 1, 0, 1, 1]
        >>> y_scores = [0.1, 0.9, 0.2, 0.8, 0.7]
        >>> ax = precision_recall_curve(y_true, y_scores)
    """
    if ax is None:
        _, ax = plt.subplots()

    try:
        from sklearn.metrics import precision_recall_curve, auc
        precision, recall, _ = precision_recall_curve(y_true, y_scores)
        pr_auc = auc(recall, precision)

        ax.plot(recall, precision, label=f'PR curve (AUC = {pr_auc:.2f})', **kwargs)
        ax.set_xlabel('Recall')
        ax.set_ylabel('Precision')
        ax.set_title(title)
        ax.legend()
    except ImportError:
        ax.text(0.5, 0.5, "scikit-learn required for PR curve", ha="center", va="center", transform=ax.transAxes)

    return ax


def residual_plot(
    y_true: Sequence[float],
    y_pred: Sequence[float],
    *,
    ax: plt.Axes | None = None,
    title: str = "Residual Plot",
    **kwargs
) -> plt.Axes:
    """Create a residual plot for regression diagnostics.

    Args:
        y_true: True values
        y_pred: Predicted values
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments for scatter

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import residual_plot
        >>> import numpy as np
        >>> y_true = np.array([1, 2, 3, 4, 5])
        >>> y_pred = np.array([1.1, 1.9, 3.2, 3.8, 5.1])
        >>> ax = residual_plot(y_true, y_pred)
    """
    if ax is None:
        _, ax = plt.subplots()

    residuals = np.array(y_true) - np.array(y_pred)

    ax.scatter(y_pred, residuals, **kwargs)
    ax.axhline(y=0, color='r', linestyle='--')
    ax.set_xlabel('Predicted Values')
    ax.set_ylabel('Residuals')
    ax.set_title(title)

    return ax


def leverage_plot(
    X: Sequence[Sequence[float]],
    y: Sequence[float],
    *,
    ax: plt.Axes | None = None,
    title: str = "Leverage Plot",
    **kwargs
) -> plt.Axes:
    """Create a leverage plot for regression diagnostics.

    Args:
        X: Feature matrix
        y: Target values
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments for scatter

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import leverage_plot
        >>> import numpy as np
        >>> X = np.random.random((100, 2))
        >>> y = np.random.random(100)
        >>> ax = leverage_plot(X, y)
    """
    if ax is None:
        _, ax = plt.subplots()

    try:
        from sklearn.linear_model import LinearRegression
        from sklearn.preprocessing import StandardScaler

        X = np.array(X)
        y = np.array(y)

        # Fit model
        model = LinearRegression()
        model.fit(X, y)
        y_pred = model.predict(X)

        # Calculate leverage
        X_scaled = StandardScaler().fit_transform(X)
        H = X_scaled @ np.linalg.pinv(X_scaled.T @ X_scaled) @ X_scaled.T
        leverage = np.diag(H)

        # Calculate residuals
        residuals = y - y_pred

        ax.scatter(leverage, residuals, **kwargs)
        ax.set_xlabel('Leverage')
        ax.set_ylabel('Residuals')
        ax.set_title(title)
    except ImportError:
        ax.text(0.5, 0.5, "scikit-learn required for leverage plot", ha="center", va="center", transform=ax.transAxes)

    return ax

