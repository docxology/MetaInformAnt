"""Tests for ``visualization.interactive_dashboards.dashboards``, which had none.

The module builds JSON-serialisable plot structures and the HTML that shows them. Findings fixed here, each tested
against the code before the fix:

* ``create_dashboard`` printed every panel as a JSON dump cut off mid-token at 2000 characters (no plot at all),
  accepted any ``layout`` string, ignored a panel's documented ``width`` and put titles and data into the page
  unescaped;
* ``export_to_html`` put titles and labels into the page unescaped and the data into ``<script>`` unescaped, so a gene
  name containing ``</script>`` ended the script and injected markup; it drew nothing for heatmap or track data and
  wrote the file in the platform's default encoding although the page says UTF-8;
* ``create_interactive_heatmap`` called a greedy nearest-neighbour chain "hierarchical clustering" and returned only
  an order as its "dendrograms";
* ``create_genome_browser_track`` invented position 0, ``end = start + 1``, the name "unnamed", chromosome "chr1" and
  the range 0-1000 for an empty track; ``create_interactive_volcano`` read a missing p-value as 1.0 and a missing fold
  change as 0.0.

The HTML is parsed with the standard-library parser, so the assertions are about elements, not substrings.
"""

from __future__ import annotations

import json
import locale
from html.parser import HTMLParser
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pytest
from scipy.cluster.hierarchy import linkage

from metainformant.visualization.interactive_dashboards import dashboards as d


class _Elements(HTMLParser):
    """Collects start tags (with attributes) and the text of the document's own ``<title>`` (the first one: SVG
    hover titles come later and are not the page title)."""

    def __init__(self) -> None:
        super().__init__()
        self.tags: List[Tuple[str, Dict[str, str | None]]] = []
        self.title = ""
        self._in_title = False
        self._title_seen = False

    def handle_starttag(self, tag: str, attrs: List[Tuple[str, str | None]]) -> None:
        self.tags.append((tag, dict(attrs)))
        self._in_title = tag == "title" and not self._title_seen

    def handle_endtag(self, tag: str) -> None:
        if tag == "title":
            self._in_title = False
            self._title_seen = True

    def handle_data(self, data: str) -> None:
        if self._in_title:
            self.title += data

    def count(self, tag: str) -> int:
        return sum(1 for name, _ in self.tags if name == tag)


def _parse(markup: str) -> _Elements:
    parser = _Elements()
    parser.feed(markup)
    return parser


def _scatter() -> dict:
    return d.create_interactive_scatter([1.0, 2.0, 3.0], [3.0, 1.0, 2.0], labels=["a", "b", "c"])


def _heatmap() -> dict:
    return d.create_interactive_heatmap([[1.0, 2.0], [3.0, 4.0], [5.0, 9.0]], ["r1", "r2", "r3"], ["c1", "c2"])


def _track() -> dict:
    return d.create_genome_browser_track(
        [{"start": 10, "end": 50, "name": "geneA"}, {"start": 70, "end": 90, "name": "geneB"}]
    )


# ---------------------------------------------------------------------------
# create_dashboard
# ---------------------------------------------------------------------------


def test_dashboard_draws_every_supported_panel_instead_of_dumping_json() -> None:
    volcano = d.create_interactive_volcano(
        [{"gene": "G1", "log2fc": 2.0, "pvalue": 0.001}, {"gene": "G2", "log2fc": -0.2, "pvalue": 0.5}]
    )

    dashboard = d.create_dashboard(
        [
            {"title": "Scatter", "plot_data": _scatter()},
            {"title": "Heatmap", "plot_data": _heatmap()},
            {"title": "Track", "plot_data": _track()},
            {"title": "Volcano", "plot_data": volcano},
        ]
    )

    page = _parse(dashboard["html"])
    assert page.count("svg") == 4 and page.count("pre") == 0  # one drawing per panel, nothing printed as JSON
    circles = [attrs for tag, attrs in page.tags if tag == "circle"]
    assert len(circles) == 3 + 2  # the scatter's points and the volcano's
    assert page.count("rect") >= 6 + 2  # six heatmap cells and two track features (plus the white backgrounds)
    assert dashboard["n_panels"] == 4 and dashboard["layout"] == "grid"


def test_dashboard_marks_an_undrawable_panel_and_says_when_it_truncates() -> None:
    huge = {"blob": ["x" * 50] * 100}

    dashboard = d.create_dashboard([{"title": "Raw", "plot_data": huge}, {"title": "Small", "plot_data": {"a": 1}}])

    assert dashboard["html"].count("not drawn: unsupported plot structure") == 2
    assert "... (truncated)" in dashboard["html"]  # the 2000-character cut used to land mid-token without a word
    assert dashboard["html"].count("... (truncated)") == 1  # only the long one


def test_dashboard_refuses_an_unknown_layout_and_honours_panel_width() -> None:
    with pytest.raises(ValueError, match="Unknown dashboard layout 'diagonal'"):
        d.create_dashboard([{"title": "x", "plot_data": _scatter()}], layout="diagonal")

    sized = d.create_dashboard(
        [{"title": "narrow", "plot_data": _scatter(), "width": 300}, {"title": "auto", "plot_data": _scatter()}],
        layout="columns",
    )

    styles = [
        attrs["style"] for tag, attrs in _parse(sized["html"]).tags if tag == "div" and "inline-block" in str(attrs)
    ]
    assert "width:300px" in styles[0].replace(" ", "").replace("\n", "")
    assert "width:50.0%" in styles[1].replace(" ", "").replace("\n", "")
    assert sized["dashboard_data"][0]["size"]["width_px"] == 300
    assert d.create_dashboard([], layout="rows")["n_panels"] == 0


def test_dashboard_escapes_titles_labels_and_data() -> None:
    hostile = "<script>alert(1)</script>"
    scatter = d.create_interactive_scatter([1.0, 2.0], [2.0, 1.0], labels=[hostile, "<img src=x onerror=alert(2)>"])

    dashboard = d.create_dashboard(
        [{"title": hostile, "plot_data": scatter}, {"title": "<b>bold</b>", "plot_data": {"raw": hostile}}],
        title="<i>dash</i>",
    )

    page = _parse(dashboard["html"])
    assert page.count("script") == 0 and page.count("img") == 0 and page.count("b") == 0 and page.count("i") == 0
    assert page.title == "<i>dash</i>"  # the title text, shown as text
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in dashboard["html"]


# ---------------------------------------------------------------------------
# export_to_html
# ---------------------------------------------------------------------------


def test_exported_page_cannot_be_broken_out_of_by_a_label_or_title(tmp_path: Path) -> None:
    """``</script>`` inside a gene name ended the script block and let the rest be parsed as markup."""
    hostile = "</script><img src=x onerror=alert(1)>"
    plot = d.create_interactive_scatter([1.0, 2.0], [1.0, 2.0], labels=[hostile, "ok"])
    path = d.export_to_html(plot, str(tmp_path / "plot.html"), title=f"<u>{hostile}")

    markup = Path(path).read_text(encoding="utf-8")

    page = _parse(markup)
    assert page.count("img") == 0 and page.count("u") == 0
    assert page.count("script") == 1  # the one script of the page itself
    assert markup.count("</script>") == 1
    assert page.title == f"<u>{hostile}"
    embedded = markup.split("const plotData = ", 1)[1].split(";\n", 1)[0]
    assert json.loads(embedded)["plot_data"][0]["label"] == hostile  # the data is intact, only its encoding changed


def test_exported_page_is_utf8_whatever_the_platform_default(tmp_path: Path) -> None:
    """``write_text`` without an encoding used the locale's: under the C locale (ASCII) a Greek title raised."""
    plot = d.create_interactive_scatter([1.0, 2.0], [1.0, 2.0], labels=["α-tubulin", "β-actin"])
    previous = locale.setlocale(locale.LC_CTYPE)
    locale.setlocale(locale.LC_CTYPE, "C")  # a platform whose default text encoding is ASCII
    try:
        path = d.export_to_html(plot, str(tmp_path / "greek.html"), title="Δ expression")
    finally:
        locale.setlocale(locale.LC_CTYPE, previous)

    raw = Path(path).read_bytes()

    assert "Δ expression".encode("utf-8") in raw  # written as UTF-8, as the page's meta charset says
    assert Path(path).is_absolute()


def test_export_draws_heatmap_and_track_data_as_static_svg(tmp_path: Path) -> None:
    """The page's script only knew points with x and y, so heatmaps and tracks exported as an empty box."""
    heat = _parse(Path(d.export_to_html(_heatmap(), str(tmp_path / "heat.html"))).read_text(encoding="utf-8"))
    track = _parse(Path(d.export_to_html(_track(), str(tmp_path / "track.html"))).read_text(encoding="utf-8"))

    assert heat.count("svg") == 1 and heat.count("rect") >= 6
    assert track.count("svg") == 1 and track.count("rect") >= 2
    scatter_page = _parse(Path(d.export_to_html(_scatter(), str(tmp_path / "s.html"))).read_text(encoding="utf-8"))
    assert scatter_page.count("svg") == 0  # a scatter is still drawn by the page's own script (hover tooltips)


# ---------------------------------------------------------------------------
# create_interactive_heatmap
# ---------------------------------------------------------------------------


def _members_are_contiguous(order: List[int], merges: np.ndarray) -> bool:
    n = len(order)
    position = {leaf: index for index, leaf in enumerate(order)}
    members = {i: [i] for i in range(n)}
    for k, (left, right, _, _) in enumerate(merges):
        members[n + k] = members[int(left)] + members[int(right)]
        spots = sorted(position[m] for m in members[n + k])
        if spots[-1] - spots[0] != len(spots) - 1:
            return False
    return True


def test_heatmap_row_order_keeps_every_cluster_contiguous_and_returns_the_dendrogram() -> None:
    """The greedy chain started at row 0 and appended 0.2 (the nearest neighbour of 2.2) at the very end."""
    data = [[2.2], [4.6], [2.9], [0.2], [8.4], [5.6]]

    result = d.create_interactive_heatmap(data, list("abcdef"), ["x"], cluster_cols=False)

    order = result["row_order"]
    assert sorted(order) == list(range(6))
    assert _members_are_contiguous(order, linkage(np.array(data), method="average"))
    rows = result["dendrograms"]["rows"]
    assert rows["order"] == order and rows["clustered"] is True
    merges = np.array(rows["linkage"])
    assert merges.shape == (5, 4) and merges[-1, 3] == 6  # a merge table, ending with a cluster of all six rows
    np.testing.assert_allclose(merges[:, 2], linkage(np.array(data), method="average")[:, 2])
    assert [cell["row"] for cell in result["plot_data"][:1]] == ["d"]  # the cells follow the new order


def test_heatmap_without_clustering_keeps_the_input_order_and_has_no_dendrogram() -> None:
    result = d.create_interactive_heatmap([[3.0, 1.0], [1.0, 2.0]], ["r1", "r2"], ["c1", "c2"], False, False)

    assert result["row_order"] == [0, 1] and result["col_order"] == [0, 1]
    assert result["dendrograms"]["rows"] == {"order": [0, 1], "clustered": False, "linkage": []}
    assert result["color_range"] == (1.0, 3.0)


def test_heatmap_clustering_refuses_non_finite_values() -> None:
    with pytest.raises(ValueError, match="NaN or infinite"):
        d.create_interactive_heatmap([[1.0, float("nan")], [2.0, 3.0]], ["a", "b"], ["x", "y"])
    assert d.create_interactive_heatmap([[1.0, float("nan")], [2.0, 3.0]], ["a", "b"], ["x", "y"], False, False)


# ---------------------------------------------------------------------------
# create_genome_browser_track / create_interactive_volcano
# ---------------------------------------------------------------------------


def test_track_features_missing_required_keys_are_refused_not_invented() -> None:
    for missing in ("start", "end", "name"):
        feature = {"start": 5, "end": 9, "name": "g"}
        del feature[missing]
        with pytest.raises(ValueError, match=f"feature 0 is missing the required key.*{missing}"):
            d.create_genome_browser_track([feature])


def test_track_does_not_guess_a_chromosome_or_a_range() -> None:
    track = d.create_genome_browser_track([{"start": 10, "end": 50, "name": "g"}])

    assert track["track_data"][0]["chromosome"] is None  # not "chr1"
    assert track["coordinate_range"] == (10, 50)
    assert d.create_genome_browser_track([])["coordinate_range"] == (0, 0)  # not an invented 0-1000


def test_volcano_results_without_a_p_value_or_fold_change_are_refused() -> None:
    with pytest.raises(ValueError, match="has no p-value"):
        d.create_interactive_volcano([{"gene": "G1", "log2fc": 3.0}])  # was counted as p = 1.0, "not significant"
    with pytest.raises(ValueError, match="has no fold change"):
        d.create_interactive_volcano([{"gene": "G1", "pvalue": 0.001}])  # was counted as a fold change of 0.0
    ok = d.create_interactive_volcano([{"gene": "G1", "log2_fold_change": 3.0, "adjusted_p": 0.001}])
    assert ok["significance_counts"] == {"up": 1, "down": 0, "not_significant": 0}
