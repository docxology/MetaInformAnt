"""Semantics of ``ancestors`` / ``descendants``: the term itself is excluded, and the documentation says so.

The docstrings used to promise "including the term itself" while the code (and every caller and test) excludes it:
``Ontology.get_ancestors`` / ``get_descendants`` never return the start term, ``get_subontology`` and
``validate_go_ontology`` add their roots separately, and the query tests assert ``"C" not in ancestors(onto, "C")``.
"""

from __future__ import annotations

import doctest
import re
import sys
from typing import Callable

import pytest

from metainformant.ontology.core.types import Ontology, Relationship, Term
from metainformant.ontology.query.query import (
    ancestors,
    clear_cache,
    common_ancestors,
    descendants,
    get_subontology,
)


@pytest.fixture(autouse=True)
def fresh_query_cache():
    """The query cache is keyed by ``id(onto)``; start and finish each test with it empty."""
    clear_cache()
    yield
    clear_cache()


def _chain() -> Ontology:
    """A -> B -> C (C is_a B, B is_a A)."""
    onto = Ontology()
    onto.add_term(Term("A", "root"))
    onto.add_term(Term("B", "child", is_a_parents=["A"]))
    onto.add_term(Term("C", "grandchild", is_a_parents=["B"]))
    return onto


class TestTheTermItselfIsExcluded:
    def test_ancestors_exclude_the_term(self) -> None:
        onto = _chain()
        assert ancestors(onto, "C") == {"A", "B"}
        assert ancestors(onto, "B") == {"A"}

    def test_descendants_exclude_the_term(self) -> None:
        onto = _chain()
        assert descendants(onto, "A") == {"B", "C"}
        assert descendants(onto, "B") == {"C"}

    def test_a_root_has_no_ancestors_and_a_leaf_no_descendants(self) -> None:
        onto = _chain()
        assert ancestors(onto, "A") == set()
        assert descendants(onto, "C") == set()

    def test_other_relation_types_follow_the_same_rule(self) -> None:
        onto = _chain()
        onto.add_relationship(Relationship("C", "A", "part_of"))
        assert ancestors(onto, "C", relation_type="part_of") == {"A"}
        assert descendants(onto, "A", relation_type="part_of") == {"C"}

    def test_the_term_can_be_added_back_explicitly(self) -> None:
        onto = _chain()
        assert ancestors(onto, "C") | {"C"} == {"A", "B", "C"}
        assert descendants(onto, "A") | {"A"} == {"A", "B", "C"}

    def test_callers_add_their_roots_themselves(self) -> None:
        # get_subontology unions descendants() with the roots precisely because descendants() leaves them out
        sub = get_subontology(_chain(), ["B"])
        assert set(sub.terms) == {"B", "C"}

    def test_common_ancestors_exclude_the_pair_themselves(self) -> None:
        onto = Ontology()
        onto.add_term(Term("A", "root"))
        onto.add_term(Term("B", "left", is_a_parents=["A"]))
        onto.add_term(Term("C", "right", is_a_parents=["A"]))
        assert common_ancestors(onto, "B", "C") == {"A"}
        # a term that is an ancestor of the other is not reported as a common ancestor
        assert common_ancestors(onto, "B", "A") == set()


def _docstring(func: Callable[..., object]) -> str:
    return func.__doc__ or ""


class TestDocumentationMatchesBehaviour:
    @pytest.mark.parametrize("func", [ancestors, descendants])
    def test_docstring_no_longer_claims_the_term_is_included(self, func: Callable[..., object]) -> None:
        doc = _docstring(func)
        assert not re.search(r"\bincluding the term itself\b", doc, re.IGNORECASE)
        assert re.search(r"\bexcluding the term itself\b", doc, re.IGNORECASE)

    @pytest.mark.parametrize("func", [ancestors, descendants, common_ancestors])
    def test_docstring_examples_run_and_give_the_documented_results(self, func: Callable[..., object]) -> None:
        module = sys.modules[func.__module__]
        tests = doctest.DocTestFinder().find(func, func.__name__, globs=dict(vars(module)))
        runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS)
        report: list[str] = []
        for test in tests:
            runner.run(test, out=report.append)
        assert runner.tries >= 3, "the docstring should carry a runnable example"
        assert runner.failures == 0, "".join(report)
