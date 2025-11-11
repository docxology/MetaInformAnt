"""Multi-dimensional visualization functions.

This module provides visualization functions for multi-dimensional data including
pair plots, parallel coordinates, radar charts, scatter plot matrices, and 3D plots.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Sequence, Union

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from matplotlib.figure import Figure
    try:
        from seaborn import PairGrid
    except ImportError:
        PairGrid = None  # type: ignore

# Use non-interactive backend by default for tests/headless
matplotlib.use("Agg", force=True)

try:
    import seaborn as sns

    SEABORN_AVAILABLE = True
except ImportError:
    SEABORN_AVAILABLE = False
    sns = None


def pairplot_dataframe(df: pd.DataFrame, *, hue: str | None = None) -> Union[Figure, "PairGrid"]:
    """Pairplot for a tidy DataFrame, returns seaborn PairGrid or matplotlib figure.

    Args:
        df: DataFrame with numeric columns
        hue: Optional column name for coloring

    Returns:
        seaborn PairGrid or matplotlib Figure

    Example:
        >>> from metainformant.visualization import pairplot_dataframe
        >>> import pandas as pd
        >>> import numpy as np
        >>> df = pd.DataFrame(np.random.random((100, 3)), columns=['A', 'B', 'C'])
        >>> grid = pairplot_dataframe(df)
    """
    if SEABORN_AVAILABLE:
        grid = sns.pairplot(df, hue=hue)
        return grid
    else:
        # Fallback: simple scatter plot matrix with matplotlib
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        n_cols = len(numeric_cols)

        fig, axes = plt.subplots(n_cols, n_cols, figsize=(n_cols * 3, n_cols * 3))
        if n_cols == 1:
            axes = np.array([[axes]])
        elif n_cols == 2:
            axes = axes.reshape(2, 2)

        for i, col1 in enumerate(numeric_cols):
            for j, col2 in enumerate(numeric_cols):
                ax = axes[i, j]
                if i == j:
                    # Diagonal: histogram
                    ax.hist(df[col1], alpha=0.7)
                    ax.set_xlabel(col1)
                else:
                    # Off-diagonal: scatter plot
                    ax.scatter(df[col2], df[col1], alpha=0.6)
                    ax.set_xlabel(col2)
                    ax.set_ylabel(col1)

        plt.tight_layout()
        return fig


def parallel_coordinates_plot(
    data: pd.DataFrame,
    class_column: str | None = None,
    *,
    ax: plt.Axes | None = None,
    title: str = "Parallel Coordinates Plot",
    **kwargs
) -> plt.Axes:
    """Create a parallel coordinates plot.

    Args:
        data: DataFrame with numeric columns
        class_column: Optional column name for class coloring
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import parallel_coordinates_plot
        >>> import pandas as pd
        >>> import numpy as np
        >>> df = pd.DataFrame(np.random.random((20, 4)), columns=['A', 'B', 'C', 'D'])
        >>> df['class'] = ['X' if i < 10 else 'Y' for i in range(20)]
        >>> ax = parallel_coordinates_plot(df, class_column='class')
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    numeric_cols = data.select_dtypes(include=[np.number]).columns.tolist()
    if class_column:
        numeric_cols = [c for c in numeric_cols if c != class_column]

    n_cols = len(numeric_cols)
    if n_cols < 2:
        raise ValueError("Need at least 2 numeric columns for parallel coordinates plot")

    # Normalize each column to 0-1
    normalized_data = data[numeric_cols].copy()
    for col in numeric_cols:
        col_min = normalized_data[col].min()
        col_max = normalized_data[col].max()
        if col_max > col_min:
            normalized_data[col] = (normalized_data[col] - col_min) / (col_max - col_min)

    # Plot lines
    x_positions = np.arange(n_cols)
    colors = None
    if class_column and class_column in data.columns:
        unique_classes = data[class_column].unique()
        color_map = plt.cm.tab10(np.linspace(0, 1, len(unique_classes)))
        class_colors = {cls: color_map[i] for i, cls in enumerate(unique_classes)}
        colors = [class_colors[cls] for cls in data[class_column]]

    for i in range(len(data)):
        y_values = normalized_data.iloc[i].values
        color = colors[i] if colors else None
        alpha = 0.3 if class_column else 0.5
        ax.plot(x_positions, y_values, color=color, alpha=alpha, **kwargs)

    ax.set_xticks(x_positions)
    ax.set_xticklabels(numeric_cols, rotation=45, ha='right')
    ax.set_ylabel("Normalized Value")
    ax.set_title(title)
    ax.set_ylim(-0.1, 1.1)

    return ax


def radar_chart(
    categories: Sequence[str],
    values: Sequence[float],
    *,
    ax: plt.Axes | None = None,
    title: str = "Radar Chart",
    **kwargs
) -> plt.Axes:
    """Create a radar chart (spider chart).

    Args:
        categories: Category names
        values: Values for each category
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import radar_chart
        >>> categories = ['A', 'B', 'C', 'D', 'E']
        >>> values = [0.8, 0.6, 0.9, 0.7, 0.5]
        >>> ax = radar_chart(categories, values)
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection='polar'))
    else:
        # Ensure polar projection
        if not hasattr(ax, 'set_theta_zero_location'):
            fig = ax.figure
            ax.remove()
            ax = fig.add_subplot(111, projection='polar')

    n = len(categories)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False).tolist()
    values = list(values)
    
    # Close the plot
    angles += angles[:1]
    values += values[:1]

    ax.plot(angles, values, 'o-', linewidth=2, **kwargs)
    ax.fill(angles, values, alpha=0.25, **kwargs)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories)
    ax.set_ylim(0, 1)
    ax.set_title(title, pad=20)

    return ax


def splom_plot(
    data: pd.DataFrame,
    *,
    hue: str | None = None,
    figsize: tuple[int, int] | None = None,
    **kwargs
):
    """Create a scatter plot matrix (SPLOM).

    Args:
        data: DataFrame with numeric columns
        hue: Optional column name for coloring
        figsize: Figure size tuple
        **kwargs: Additional arguments

    Returns:
        matplotlib Figure

    Example:
        >>> from metainformant.visualization import splom_plot
        >>> import pandas as pd
        >>> import numpy as np
        >>> df = pd.DataFrame(np.random.random((50, 4)), columns=['A', 'B', 'C', 'D'])
        >>> fig = splom_plot(df)
    """
    numeric_cols = data.select_dtypes(include=[np.number]).columns.tolist()
    n_cols = len(numeric_cols)

    if figsize is None:
        figsize = (n_cols * 2.5, n_cols * 2.5)

    fig, axes = plt.subplots(n_cols, n_cols, figsize=figsize)

    if n_cols == 1:
        axes = np.array([[axes]])
    elif n_cols == 2:
        axes = axes.reshape(2, 2)
    else:
        axes = axes.flatten().reshape(n_cols, n_cols)

    colors = None
    if hue and hue in data.columns:
        unique_classes = data[hue].unique()
        color_map = plt.cm.tab10(np.linspace(0, 1, len(unique_classes)))
        class_colors = {cls: color_map[i] for i, cls in enumerate(unique_classes)}
        colors = [class_colors[cls] for cls in data[hue]]

    for i, col1 in enumerate(numeric_cols):
        for j, col2 in enumerate(numeric_cols):
            ax = axes[i, j]
            if i == j:
                # Diagonal: histogram
                ax.hist(data[col1], alpha=0.7, bins=20)
                ax.set_xlabel(col1)
            else:
                # Off-diagonal: scatter plot
                if colors:
                    for cls, color in class_colors.items():
                        mask = data[hue] == cls
                        ax.scatter(data.loc[mask, col2], data.loc[mask, col1],
                                  c=[color], alpha=0.6, label=cls, **kwargs)
                else:
                    ax.scatter(data[col2], data[col1], alpha=0.6, **kwargs)
                ax.set_xlabel(col2)
                ax.set_ylabel(col1)

    if hue and n_cols > 1:
        axes[0, n_cols-1].legend(loc='upper right')

    plt.tight_layout()
    return fig


def scatter_3d(
    x: Sequence[float],
    y: Sequence[float],
    z: Sequence[float],
    *,
    color: Sequence[float] | Sequence[str] | None = None,
    ax: plt.Axes | None = None,
    title: str = "3D Scatter Plot",
    **kwargs
) -> plt.Axes:
    """Create a 3D scatter plot.

    Args:
        x: X-axis values
        y: Y-axis values
        z: Z-axis values
        color: Optional color values
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object (3D)

    Example:
        >>> from metainformant.visualization import scatter_3d
        >>> import numpy as np
        >>> x = np.random.random(100)
        >>> y = np.random.random(100)
        >>> z = np.random.random(100)
        >>> ax = scatter_3d(x, y, z)
    """
    if ax is None:
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
    else:
        # Ensure 3D projection
        if not hasattr(ax, 'zaxis'):
            fig = ax.figure
            ax.remove()
            ax = fig.add_subplot(111, projection='3d')

    ax.scatter(x, y, z, c=color, **kwargs)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title(title)

    return ax

