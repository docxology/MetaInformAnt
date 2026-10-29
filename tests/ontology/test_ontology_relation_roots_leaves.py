"""``Ontology.get_roots`` and ``get_leaves`` honour ``relation_type``, and the enrichment dot plot draws its FDR line.

``get_roots(relation_type)`` / ``get_leaves(relation_type)`` computed the ``is_a`` roots and leaves whatever the
argument was, although ``get_parents`` / ``get_children`` / ``get_ancestors`` / ``get_descendants`` all take the same
argument and honour it. ``enrichment_dotplot`` documented ``alpha_threshold`` as the "FDR threshold for dashed-line
annotation" and never drew a line.

The oracles are the parent and child sets of the ontology built in the test (relationship by relationship) and the
line objects of the rendered colour bar. No test doubles.
"""

from __future__ import annotations

import math

import pytest

from metainformant.ontology.core.types import Ontology, Relationship, Term


def _ontology() -> Ontology:
    """B is_a A, C is_a B; D part_of A, E part_of D; F stands alone."""
    ontology = Ontology()
    for term_id in "ABCDEF":
        ontology.add_term(Term(term_id=term_id, name=f"term {term_id}"))
    for source, target, relation in [
        ("B", "A", "is_a"),
        ("C", "B", "is_a"),
        ("D", "A", "part_of"),
        ("E", "D", "part_of"),
    ]:
        ontology.add_relationship(Relationship(source, target, relation))
    return ontology


def test_the_test_ontology_is_what_the_oracles_assume() -> None:
    ontology = _ontology()

    assert ontology.get_parents("B", "is_a") == {"A"} and ontology.get_parents("E", "part_of") == {"D"}
    assert ontology.get_children("A", "part_of") == {"D"} and ontology.get_children("A", "is_a") == {"B"}


def test_roots_follow_the_relation_type() -> None:
    ontology = _ontology()

    assert ontology.get_roots() == {"A", "D", "E", "F"}  # no is_a parent
    assert ontology.get_roots("is_a") == {"A", "D", "E", "F"}
    # the old answer for part_of was the is_a roots; D and E have part_of parents, B and C do not
    assert ontology.get_roots("part_of") == {"A", "B", "C", "F"}
    assert ontology.get_roots("regulates") == set("ABCDEF")  # nothing is related that way


def test_leaves_follow_the_relation_type() -> None:
    ontology = _ontology()

    assert ontology.get_leaves() == {"C", "D", "E", "F"}  # no is_a child
    assert ontology.get_leaves("part_of") == {"B", "C", "E", "F"}
    assert ontology.get_leaves("regulates") == set("ABCDEF")


def test_roots_and_leaves_agree_with_the_parent_and_child_queries_for_every_relation() -> None:
    ontology = _ontology()

    for relation in ("is_a", "part_of", "regulates"):
        roots = {t for t in ontology.terms if not ontology.get_parents(t, relation)}
        leaves = {t for t in ontology.terms if not ontology.get_children(t, relation)}
        assert ontology.get_roots(relation) == roots
        assert ontology.get_leaves(relation) == leaves


def test_the_dotplot_draws_the_fdr_threshold_on_its_colour_bar() -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from metainformant.ontology.visualization.plots import enrichment_dotplot

    results = [
        {"term_id": f"GO:{i}", "term_name": f"term {i}", "adjusted_p": p, "n_overlap": 3 + i, "query_size": 30}
        for i, p in enumerate([1e-8, 1e-4, 0.01, 0.2])
    ]

    figure = enrichment_dotplot(results, alpha_threshold=0.05)
    colour_bar = next(axes for axes in figure.axes if axes.get_ylabel() == "-log10(adj p-value)")

    dashed = [line for line in colour_bar.lines if line.get_linestyle() == "--"]
    assert len(dashed) == 1
    assert list(dashed[0].get_ydata()) == pytest.approx([-math.log10(0.05)] * 2)
    assert any("0.05" in text.get_text() for text in colour_bar.texts)

    # a stricter threshold moves the line
    strict = enrichment_dotplot(results, alpha_threshold=0.001)
    strict_bar = next(axes for axes in strict.axes if axes.get_ylabel() == "-log10(adj p-value)")
    strict_line = next(line for line in strict_bar.lines if line.get_linestyle() == "--")
    assert list(strict_line.get_ydata())[0] == pytest.approx(3.0)


def test_the_dotplot_leaves_the_line_out_when_no_term_reaches_the_threshold_and_rejects_nonsense() -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from metainformant.ontology.visualization.plots import enrichment_dotplot

    results = [{"term_id": "GO:1", "adjusted_p": 0.5, "n_overlap": 3, "query_size": 30}]
    figure = enrichment_dotplot(results, alpha_threshold=0.05)  # -log10(0.5) = 0.3 < 1.3: the bar stops short of it

    colour_bar = next(axes for axes in figure.axes if axes.get_ylabel() == "-log10(adj p-value)")
    assert not [line for line in colour_bar.lines if line.get_linestyle() == "--"]
    for bad in (0.0, -0.1, 1.5):
        with pytest.raises(ValueError, match="alpha_threshold"):
            enrichment_dotplot(results, alpha_threshold=bad)
