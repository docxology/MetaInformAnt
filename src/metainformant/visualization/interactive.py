"""Interactive visualization utilities.

This module provides utilities for creating interactive plots using Plotly
(optional dependency) with hover tooltips and selection capabilities.
"""

from __future__ import annotations

from typing import Any

try:
    import plotly.graph_objects as go
    import plotly.express as px
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    go = None
    px = None


def create_interactive_scatter(
    x: list[float],
    y: list[float],
    *,
    labels: list[str] | None = None,
    colors: list[str] | None = None,
    title: str = "Interactive Scatter Plot",
    **kwargs
) -> Any:
    """Create an interactive scatter plot.

    Args:
        x: X-axis values
        y: Y-axis values
        labels: Optional labels for hover tooltips
        colors: Optional color values
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Plotly figure object

    Example:
        >>> from metainformant.visualization.interactive import create_interactive_scatter
        >>> import numpy as np
        >>> x = np.random.random(100)
        >>> y = np.random.random(100)
        >>> fig = create_interactive_scatter(x, y)
        >>> fig.show()
    """
    if not PLOTLY_AVAILABLE:
        raise ImportError("Plotly not available. Install with: uv pip install plotly")

    fig = go.Figure()

    hover_text = labels if labels else [f"Point {i}" for i in range(len(x))]

    fig.add_trace(go.Scatter(
        x=x,
        y=y,
        mode='markers',
        marker=dict(color=colors, size=10) if colors else dict(size=10),
        text=hover_text,
        hovertemplate='<b>%{text}</b><br>X: %{x}<br>Y: %{y}<extra></extra>',
        **kwargs
    ))

    fig.update_layout(
        title=title,
        xaxis_title="X",
        yaxis_title="Y",
    )

    return fig


def create_interactive_heatmap(
    data: list[list[float]],
    *,
    x_labels: list[str] | None = None,
    y_labels: list[str] | None = None,
    title: str = "Interactive Heatmap",
    **kwargs
) -> Any:
    """Create an interactive heatmap.

    Args:
        data: 2D array of values
        x_labels: Optional x-axis labels
        y_labels: Optional y-axis labels
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Plotly figure object

    Example:
        >>> from metainformant.visualization.interactive import create_interactive_heatmap
        >>> import numpy as np
        >>> data = np.random.random((10, 10)).tolist()
        >>> fig = create_interactive_heatmap(data)
        >>> fig.show()
    """
    if not PLOTLY_AVAILABLE:
        raise ImportError("Plotly not available. Install with: uv pip install plotly")

    fig = go.Figure(data=go.Heatmap(
        z=data,
        x=x_labels,
        y=y_labels,
        hoverongaps=False,
        **kwargs
    ))

    fig.update_layout(
        title=title,
    )

    return fig


def convert_matplotlib_to_plotly(fig: Any) -> Any:
    """Convert a matplotlib figure to Plotly format.

    Args:
        fig: Matplotlib figure

    Returns:
        Plotly figure object

    Example:
        >>> from metainformant.visualization.interactive import convert_matplotlib_to_plotly
        >>> import matplotlib.pyplot as plt
        >>> fig, ax = plt.subplots()
        >>> ax.plot([1, 2, 3], [4, 5, 6])
        >>> plotly_fig = convert_matplotlib_to_plotly(fig)
    """
    if not PLOTLY_AVAILABLE:
        raise ImportError("Plotly not available. Install with: uv pip install plotly")

    try:
        import plotly.tools as tls
        return tls.mpl_to_plotly(fig)
    except ImportError:
        # Fallback: return None if conversion not available
        return None


def is_plotly_available() -> bool:
    """Check if Plotly is available.

    Returns:
        True if Plotly is installed

    Example:
        >>> from metainformant.visualization.interactive import is_plotly_available
        >>> if is_plotly_available():
        ...     fig = create_interactive_scatter([1, 2, 3], [4, 5, 6])
    """
    return PLOTLY_AVAILABLE

