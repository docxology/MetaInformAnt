"""Tests for GWAS visualization configuration and theming."""

from __future__ import annotations

import pytest

from metainformant.gwas.visualization.config import (
    THEMES,
    PlotStyle,
    apply_style,
    get_style,
    style_from_config,
)


class TestPlotStyleDefaults:
    """Verify PlotStyle dataclass default values."""

    def test_default_figsize(self) -> None:
        style = PlotStyle()
        assert style.figsize == (12, 6)

    def test_default_dpi(self) -> None:
        style = PlotStyle()
        assert style.dpi == 300

    def test_default_font_size(self) -> None:
        style = PlotStyle()
        assert style.font_size == 12

    def test_default_title_size(self) -> None:
        style = PlotStyle()
        assert style.title_size == 14

    def test_default_colormap(self) -> None:
        style = PlotStyle()
        assert style.colormap == "viridis"

    def test_default_significance_color(self) -> None:
        style = PlotStyle()
        assert style.significance_color == "red"

    def test_default_alpha_and_point_size(self) -> None:
        style = PlotStyle()
        assert style.alpha == 0.7
        assert style.point_size == 8.0

    def test_default_grid_settings(self) -> None:
        style = PlotStyle()
        assert style.grid is True
        assert style.grid_alpha == 0.3


class TestThemes:
    """Verify the THEMES dictionary."""

    def test_themes_has_three_entries(self) -> None:
        assert len(THEMES) == 3

    def test_themes_keys(self) -> None:
        assert set(THEMES.keys()) == {"publication", "presentation", "poster"}

    def test_publication_theme_values(self) -> None:
        pub = THEMES["publication"]
        assert pub.figsize == (10, 6)
        assert pub.dpi == 300
        assert pub.font_size == 10
        assert pub.title_size == 12

    def test_presentation_theme_values(self) -> None:
        pres = THEMES["presentation"]
        assert pres.figsize == (14, 8)
        assert pres.dpi == 150
        assert pres.font_size == 14
        assert pres.title_size == 18

    def test_poster_theme_values(self) -> None:
        poster = THEMES["poster"]
        assert poster.figsize == (20, 12)
        assert poster.dpi == 300
        assert poster.font_size == 18
        assert poster.title_size == 24
        assert poster.point_size == 12.0


class TestGetStyle:
    """Verify get_style theme selection and overrides."""

    def test_valid_theme(self) -> None:
        style = get_style("presentation")
        assert style.figsize == (14, 8)
        assert style.dpi == 150

    def test_invalid_theme_falls_back(self) -> None:
        style = get_style("nonexistent")
        # Should fall back to publication defaults
        assert style.figsize == (10, 6)
        assert style.dpi == 300

    def test_overrides_applied(self) -> None:
        style = get_style("publication", dpi=600, alpha=0.5)
        assert style.dpi == 600
        assert style.alpha == 0.5
        # Non-overridden fields keep theme values
        assert style.figsize == (10, 6)

    def test_unknown_override_ignored(self) -> None:
        # Should not raise, just warn
        style = get_style("publication", nonexistent_field=42)
        assert style.figsize == (10, 6)


class TestApplyStyle:
    """Verify apply_style sets matplotlib rcParams."""

    def test_apply_sets_rcparams(self) -> None:
        import matplotlib as mpl

        style = PlotStyle(
            font_size=16,
            title_size=20,
            dpi=150,
            figsize=(8, 4),
            font_family="serif",
            grid=False,
            grid_alpha=0.5,
        )
        apply_style(style)

        assert mpl.rcParams["font.size"] == 16
        assert mpl.rcParams["axes.titlesize"] == 20
        assert mpl.rcParams["figure.dpi"] == 150
        assert list(mpl.rcParams["figure.figsize"]) == [8, 4]
        assert mpl.rcParams["font.family"] == ["serif"]
        assert mpl.rcParams["axes.grid"] is False
        assert mpl.rcParams["grid.alpha"] == 0.5


class TestStyleFromConfig:
    """Verify style_from_config parses nested dicts."""

    def test_nested_config(self) -> None:
        config = {
            "output": {
                "visualization": {
                    "figsize": [16, 10],
                    "dpi": 200,
                    "font_size": 14,
                    "colormap": "plasma",
                }
            }
        }
        style = style_from_config(config)
        assert style.figsize == (16, 10)
        assert style.dpi == 200
        assert style.font_size == 14
        assert style.colormap == "plasma"
        # Defaults for unset fields
        assert style.significance_color == "red"

    def test_empty_config_uses_defaults(self) -> None:
        style = style_from_config({})
        defaults = PlotStyle()
        assert style.figsize == defaults.figsize
        assert style.dpi == defaults.dpi
        assert style.font_size == defaults.font_size

    def test_partial_config(self) -> None:
        config = {"output": {"visualization": {"dpi": 72}}}
        style = style_from_config(config)
        assert style.dpi == 72
        assert style.figsize == (12, 6)  # default
