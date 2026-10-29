r"""Regression tests for OBO ``synonym:`` parsing and writing.

An OBO synonym is a quoted string followed by a scope and a dbxref list: ``synonym: "text" EXACT [GOC:mah]``. The
parser took the text with ``re.search(r'"([^"]*)"', value)`` on the value after its `` ! `` comment had been cut off,
which

* ended the text at the first escaped quote (``"the \"x\" thing"`` became ``the \``),
* dropped a synonym whose text contains `` ! `` altogether (the comment cut removed its closing quote), and
* never decoded escape sequences,

the same flaws ``def:`` had before ``parse_obo_definition``. The writer emitted ``"{synonym}"`` without escaping, so a
synonym containing a quote or backslash could not be read back from a saved OBO file. ``parse_obo_synonym`` now parses
the value the way ``parse_obo_definition`` does, and the writer quotes with ``quote_obo_string``. No test doubles.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.ontology.core import obo as obo_module
from metainformant.ontology.core.obo import parse_obo, quote_obo_string
from metainformant.ontology.core.types import create_ontology, create_term
from metainformant.ontology.query.serialize import load_ontology, save_ontology


@pytest.mark.parametrize(
    ("value", "text", "scope", "xrefs"),
    [
        ('"kid process" EXACT []', "kid process", "EXACT", []),
        ('"the \\"x\\" process" EXACT [GOC:mah]', 'the "x" process', "EXACT", ["GOC:mah"]),
        ('"A ! B" RELATED [GOC:a, PMID:1]', "A ! B", "RELATED", ["GOC:a", "PMID:1"]),
        ('"back\\\\slash, new\\nline" BROAD []', "back\\slash, new\nline", "BROAD", []),
        ('"ends with a backslash\\\\" NARROW []', "ends with a backslash\\", "NARROW", []),
        # a synonym type may sit between the scope and the dbxref list
        ('"typed" EXACT OMO:0003000 [GOC:a]', "typed", "EXACT", ["GOC:a"]),
        # trailing modifiers and comments are not part of the synonym
        ('"x" EXACT [] ! trailing comment', "x", "EXACT", []),
        ('"x" EXACT [GOC:a] {comment="note"}', "x", "EXACT", ["GOC:a"]),
        ('"x" EXACT ! comment with [brackets]', "x", "EXACT", []),
        # no scope, no list
        ('"x" []', "x", None, []),
        ('"just text"', "just text", None, []),
        ('""  EXACT []', "", "EXACT", []),
        # older versions of this package wrote synonyms without escaping inner quotes; they still read back whole
        ('"A "special" process" EXACT []', 'A "special" process', "EXACT", []),
        ('"say "EXACT" twice" EXACT []', 'say "EXACT" twice', "EXACT", []),
        # malformed input is read leniently, not rejected
        ("not quoted ! with a comment", "not quoted", None, []),
        ('"never closed', "never closed", None, []),
    ],
)
def test_parse_obo_synonym(value: str, text: str, scope: str | None, xrefs: list[str]) -> None:
    assert obo_module.parse_obo_synonym(value) == (text, scope, xrefs)


OBO_TEXT = """format-version: 1.2
ontology: demo

[Term]
id: GO:0000001
name: demo process
synonym: "the \\"x\\" thing" EXACT []
synonym: "A ! B" RELATED [GOC:mah]
synonym: "plain" EXACT []
synonym: "path\\\\to\\\\it" NARROW [] ! a comment
synonym: "typed one" BROAD OMO:0003000 [GOC:a]
is_a: GO:0000002 ! other

[Term]
id: GO:0000002
name: other process
synonym: unquoted junk
synonym: "second" EXACT []
"""


def test_parse_obo_reads_synonyms_whole(tmp_path: Path) -> None:
    path = tmp_path / "syn.obo"
    path.write_text(OBO_TEXT)

    onto = parse_obo(path)

    assert onto.terms["GO:0000001"].synonyms == ['the "x" thing', "A ! B", "plain", "path\\to\\it", "typed one"]
    assert onto.terms["GO:0000002"].synonyms == ["second"]  # the unquoted line is not a synonym, as before


SYNONYMS = ['say "hi"', "back\\slash", "multi\nline", "A ! B", "tab\there", "ends with \\", '"leading quote', "[x] {y}"]


def test_synonyms_survive_an_obo_round_trip(tmp_path: Path) -> None:
    term = create_term(id="GO:0000001", name="demo process", synonyms=list(SYNONYMS))
    original = create_ontology(terms={term.id: term})
    path = tmp_path / "saved.obo"

    save_ontology(original, path, format="obo")
    loaded = load_ontology(path, format="obo")

    assert loaded.terms["GO:0000001"].synonyms == SYNONYMS


def test_written_synonym_lines_are_quoted_with_the_obo_escapes(tmp_path: Path) -> None:
    term = create_term(id="GO:0000001", name="demo process", synonyms=['say "hi"', "multi\nline"])
    path = tmp_path / "saved.obo"

    save_ontology(create_ontology(terms={term.id: term}), path, format="obo")

    lines = path.read_text().splitlines()
    quoted = quote_obo_string('say "hi"')
    assert f"synonym: {quoted} EXACT []" in lines
    assert 'synonym: "say \\"hi\\"" EXACT []' in lines
    assert 'synonym: "multi\\nline" EXACT []' in lines
