"""Regression tests: Okabe-Ito palette must have a single canonical definition.

The shared chart-convention contract (commit f99767686) defines OKABE_ITO in
``metainformant.visualization.config.conventions``. All other public palette
names must be derived from that definition — duplicated hex literals were the
drift risk this guards against.
"""

from __future__ import annotations

from metainformant.visualization import OKABE_ITO as PACKAGE_OKABE_ITO, WONG
from metainformant.visualization.config import palettes, themes
from metainformant.visualization.config.conventions import OKABE_ITO as CANONICAL_OKABE_ITO


def test_package_reexport_is_the_canonical_object() -> None:
    """metainformant.visualization.OKABE_ITO is a re-export, not a copy."""
    assert PACKAGE_OKABE_ITO is CANONICAL_OKABE_ITO


def test_wong_is_a_derivation_of_okabe_ito() -> None:
    """WONG (palettes + package level) must contain exactly the canonical hues."""
    assert sorted(WONG) == sorted(CANONICAL_OKABE_ITO)
    assert sorted(palettes.WONG) == sorted(CANONICAL_OKABE_ITO)
    # Wong ordering convention: black first.
    assert WONG[0] == "#000000"


def test_no_duplicate_okabe_hex_literals_outside_conventions() -> None:
    """No Okabe-Ito hex literal may appear in visualization config files that
    do not source it from conventions.py."""
    from pathlib import Path

    config_dir = Path(palettes.__file__).parent
    offenders: list[str] = []
    okabe_hexes = set(CANONICAL_OKABE_ITO)
    for py in config_dir.glob("*.py"):
        if py.name in {"conventions.py", "palettes.py", "themes.py"}:
            continue
        text = py.read_text()
        for hexcode in okabe_hexes:
            if hexcode in text:
                offenders.append(f"{py.name}: {hexcode}")
    assert not offenders, f"duplicated Okabe-Ito literals: {offenders}"


def test_colorblind_theme_cycle_uses_okabe_ito_colors() -> None:
    """The 'colorblind' theme's prop cycle must draw only canonical hues."""
    theme = themes.get_theme("colorblind")
    colors = theme["axes.prop_cycle"].by_key()["color"]
    assert colors, "colorblind theme must define a color cycle"
    assert set(colors) <= set(CANONICAL_OKABE_ITO)
