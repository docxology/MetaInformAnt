"""Real plot geometry and invalid-input controls for sensitivity intervals."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from metainformant.visualization.plots.cross_species import (
    _draw_stability_intervals,
    plot_divergence_stability,
)


def test_interval_geometry_preserves_endpoints_outside_point() -> None:
    # Given an estimate outside the feature-resampling percentile interval.
    frame = pd.DataFrame(
        {
            "point_estimate": [0.1],
            "sensitivity_lower": [0.4],
            "sensitivity_upper": [0.8],
        }
    )
    fig, ax = plt.subplots()
    try:
        # When drawn onto a real matplotlib axis.
        _draw_stability_intervals(ax, frame)
        # Then interval geometry is exactly the recorded endpoints, with a separate point.
        np.testing.assert_allclose(
            ax.collections[0].get_segments()[0], [[0.4, 0.0], [0.8, 0.0]]
        )
        np.testing.assert_allclose(ax.lines[0].get_xdata(), [0.1])
    finally:
        plt.close(fig)


@pytest.mark.parametrize(
    "lower,upper,point",
    [(0.8, 0.4, 0.5), (np.nan, 0.8, 0.5), (0.4, 2.5, 0.5), (0.4, 0.8, -0.1)],
)
def test_invalid_stability_is_rejected_before_output(
    tmp_path: Path, lower: float, upper: float, point: float
) -> None:
    # Given invalid stability evidence.
    frame = pd.DataFrame(
        {
            "species_a": ["a"],
            "species_b": ["b"],
            "point_estimate": [point],
            "sensitivity_lower": [lower],
            "sensitivity_upper": [upper],
            "sensitivity_iqr": [0.1],
        }
    )
    output = tmp_path / "invalid.png"
    # When plotted; then no figure misrepresents that evidence.
    with pytest.raises(ValueError):
        plot_divergence_stability(frame, output)
    assert not output.exists()
