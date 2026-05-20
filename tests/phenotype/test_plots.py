import matplotlib.pyplot as plt
import pandas as pd
import pytest

from metainformant.phenotype.visualization.plots import (
    plot_boxplot_with_swarm,
    plot_categorical_proportions,
    plot_correlation_heatmap,
    plot_violin,
)


@pytest.fixture
def sample_df():
    data = {
        "group": ["A", "A", "A", "B", "B", "B"],
        "category": ["X", "Y", "X", "Y", "Y", "X"],
        "value": [1.1, 2.2, 1.5, 5.5, 6.1, 5.8],
        "value2": [10, 20, 15, 50, 60, 55],
    }
    return pd.DataFrame(data)


def test_plot_boxplot_with_swarm(sample_df):
    fig = plot_boxplot_with_swarm(sample_df, "group", "value", title="Test Boxplot")
    assert isinstance(fig, plt.Figure)
    assert len(fig.axes) > 0
    assert fig.axes[0].get_title() == "Test Boxplot"


def test_plot_violin(sample_df):
    fig = plot_violin(sample_df, "group", "value", hue="category", title="Test Violin")
    assert isinstance(fig, plt.Figure)
    assert len(fig.axes) > 0


def test_plot_categorical_proportions(sample_df):
    fig = plot_categorical_proportions(sample_df, "group", "category", title="Test Bars")
    assert isinstance(fig, plt.Figure)
    assert len(fig.axes) > 0


def test_plot_correlation_heatmap(sample_df):
    corr = sample_df[["value", "value2"]].corr()
    fig = plot_correlation_heatmap(corr, title="Test Heatmap")
    assert isinstance(fig, plt.Figure)
    assert len(fig.axes) > 0
