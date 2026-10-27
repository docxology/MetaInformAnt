"""Regression tests for OBO ``def:`` parsing.

An OBO definition is a quoted string followed by a bracketed dbxref list: ``def: "text" [GOC:mah, PMID:123]``.
The parser used to strip the quote characters from both ends of the whole value, which left ``." [GOC:mah, PMID:123]``
glued to the definition, kept ``\\"`` escapes undecoded, and cut definitions that contain `` ! `` short.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.ontology.core.obo import (
    format_obo_xref_list,
    parse_obo,
    parse_obo_definition,
    quote_obo_string,
)
from metainformant.ontology.core.types import create_ontology, create_term
from metainformant.ontology.query.serialize import load_ontology, save_ontology


@pytest.mark.parametrize(
    ("value", "text", "xrefs"),
    [
        # the example from the bug report: escaped quotes, then a dbxref list
        ('"A process that \\"x\\" does." [GOC:mah, PMID:123]', 'A process that "x" does.', ["GOC:mah", "PMID:123"]),
        # ordinary GO-style definitions
        (
            '"The chemical reactions and pathways." [GOC:go_curators]',
            "The chemical reactions and pathways.",
            ["GOC:go_curators"],
        ),
        (
            '"Catalysis of the reaction: ATP + H2O = ADP + phosphate." [EC:3.6.1.3, RHEA:13065]',
            "Catalysis of the reaction: ATP + H2O = ADP + phosphate.",
            ["EC:3.6.1.3", "RHEA:13065"],
        ),
        # empty or absent dbxref list
        ('"A root term for testing" []', "A root term for testing", []),
        ('"No dbxref list at all"', "No dbxref list at all", []),
        ('  "Padded value"   [ GOC:a ]  ', "Padded value", ["GOC:a"]),
        # a "!" inside the quotes is text, not the start of a comment
        ('"Hmm ! interesting thing." [GOC:x]', "Hmm ! interesting thing.", ["GOC:x"]),
        # what may follow the dbxref list is not part of the definition
        ('"Text." [GOC:x] ! trailing comment', "Text.", ["GOC:x"]),
        ('"Text." [GOC:x] {comment="curator note"}', "Text.", ["GOC:x"]),
        ('"Text." {comment="no dbxrefs"} ! note', "Text.", []),
        # brackets and braces inside the text
        ('"Foo [bar] {baz} qux." [GOC:x]', "Foo [bar] {baz} qux.", ["GOC:x"]),
        # escape sequences: \\ \n \t \W and escaped punctuation
        ('"back\\\\slash, new\\nline, tab\\tx, sp\\Wace" []', "back\\slash, new\nline, tab\tx, sp ace", []),
        ('"ends with a backslash\\\\" [GOC:x]', "ends with a backslash\\", ["GOC:x"]),
        ('"unknown escapes stay: \\d+ and \\q" []', "unknown escapes stay: \\d+ and \\q", []),
        # dbxref details: escaped colon/comma, quoted descriptions, no spaces, stray separators, escaped bracket
        (
            '"Text." [url:http\\://example.org/a\\,b, PMID:1 "a, b"]',
            "Text.",
            ["url:http://example.org/a,b", 'PMID:1 "a, b"'],
        ),
        ('"Text." [GOC:a,GOC:b]', "Text.", ["GOC:a", "GOC:b"]),
        ('"Text." [ , GOC:a, ]', "Text.", ["GOC:a"]),
        ('"Text." [GOC:a\\]b, GOC:c]', "Text.", ["GOC:a]b", "GOC:c"]),
        ('"" []', "", []),
        # older versions of this package wrote definitions without escaping inner quotes; they still read back whole
        ('"A "special" process"', 'A "special" process', []),
        ('"A "special" process" []', 'A "special" process', []),
        # malformed input is read leniently, not rejected
        ("not quoted at all ! with a comment", "not quoted at all", []),
        ('"never closed', "never closed", []),
    ],
)
def test_parse_obo_definition(value: str, text: str, xrefs: list[str]) -> None:
    assert parse_obo_definition(value) == (text, xrefs)


def test_definition_never_keeps_quote_or_dbxref_residue() -> None:
    text, xrefs = parse_obo_definition('"A process that \\"x\\" does." [GOC:mah, PMID:123]')
    assert text == 'A process that "x" does.'
    assert '"' not in text.replace('"x"', "")
    assert "[" not in text and "GOC" not in text and "PMID" not in text
    assert xrefs == ["GOC:mah", "PMID:123"]


OBO_TEXT = """format-version: 1.2
ontology: demo

[Term]
id: GO:0000001
name: demo process
namespace: biological_process
def: "A process that \\"x\\" does." [GOC:mah, PMID:123]

[Term]
id: GO:0000002
name: hydrolase activity
namespace: molecular_function
def: "Catalysis of the reaction: ATP + H2O = ADP + phosphate." [EC:3.6.1.3, RHEA:13065]
is_a: GO:0000001 ! demo process

[Term]
id: GO:0000003
name: bang term
def: "A state ! of mind." [GOC:x]

[Term]
id: GO:0000004
name: bare definition
def: "Bare quoted text"

[Term]
id: GO:0000005
name: padded definition
def: "A root term for testing" []

[Typedef]
id: part_of
name: part of
def: "A typedef \\"quoted\\" definition." [GOC:x]
"""


@pytest.fixture
def parsed(tmp_path: Path):
    path = tmp_path / "defs.obo"
    path.write_text(OBO_TEXT)
    return parse_obo(path)


def test_parse_obo_gives_clean_definitions(parsed) -> None:
    assert parsed.terms["GO:0000001"].definition == 'A process that "x" does.'
    assert parsed.terms["GO:0000002"].definition == "Catalysis of the reaction: ATP + H2O = ADP + phosphate."
    assert parsed.terms["GO:0000003"].definition == "A state ! of mind."
    assert parsed.terms["GO:0000004"].definition == "Bare quoted text"
    assert parsed.terms["GO:0000005"].definition == "A root term for testing"


def test_parse_obo_keeps_the_definition_dbxrefs_apart_from_the_text(parsed) -> None:
    assert parsed.terms["GO:0000001"].metadata["definition_xrefs"] == ["GOC:mah", "PMID:123"]
    assert parsed.terms["GO:0000002"].metadata["definition_xrefs"] == ["EC:3.6.1.3", "RHEA:13065"]
    assert "definition_xrefs" not in parsed.terms["GO:0000004"].metadata
    assert "definition_xrefs" not in parsed.terms["GO:0000005"].metadata
    for term in parsed.terms.values():
        assert term.definition is not None
        assert '" [' not in term.definition and "GOC:" not in term.definition


def test_parse_obo_still_reads_the_rest_of_the_stanza(parsed) -> None:
    assert parsed.terms["GO:0000002"].name == "hydrolase activity"
    assert parsed.terms["GO:0000002"].namespace == "molecular_function"
    assert parsed.parents_of["GO:0000002"] == {"GO:0000001"}


def test_typedef_definition_is_parsed_the_same_way(parsed) -> None:
    typedef = parsed.metadata["typedef:part_of"]
    assert typedef["definition"] == 'A typedef "quoted" definition.'
    assert typedef["definition_xrefs"] == ["GOC:x"]


def test_shipped_mini_fixture_has_a_clean_definition() -> None:
    path = Path(__file__).parent.parent / "data" / "ontology" / "go_mini.obo"
    assert parse_obo(path).terms["GO:0000001"].definition == "A root term for testing"


TRICKY_DEFINITIONS = [
    "Plain definition.",
    'Contains "double quotes" inside.',
    'Ends with a quoted "word"',
    "Back\\slash and a literal backslash-n: \\n",
    "Ends with a backslash\\",
    "Line one\nLine two",
    "Tab\tseparated",
    "Bang ! inside, [brackets] and {braces}",
    "Looks like a dbxref list [GOC:x]",
    "Unicode: naive cafe – “quoted”",
]


@pytest.mark.parametrize("definition", TRICKY_DEFINITIONS)
def test_quote_obo_string_round_trips_through_the_parser(definition: str) -> None:
    quoted = quote_obo_string(definition)
    assert "\n" not in quoted and "\t" not in quoted  # stays on one line
    assert parse_obo_definition(f"{quoted} []") == (definition, [])


@pytest.mark.parametrize("definition", TRICKY_DEFINITIONS)
def test_obo_save_and_load_preserve_the_definition(tmp_path: Path, definition: str) -> None:
    onto = create_ontology(terms={"GO:1": create_term(id="GO:1", name="term", definition=definition)})
    path = tmp_path / "round.obo"
    save_ontology(onto, path, format="obo")
    assert load_ontology(path, format="obo").terms["GO:1"].definition == definition


def test_obo_save_and_load_preserve_definition_dbxrefs(tmp_path: Path) -> None:
    xrefs = ["GOC:mah", "PMID:1,2", 'odd"quote', "a]b", "back\\slash"]
    term = create_term(id="GO:1", name="term", definition="Cited definition.", definition_xrefs=xrefs)
    path = tmp_path / "xrefs.obo"
    save_ontology(create_ontology(terms={"GO:1": term}), path, format="obo")

    content = path.read_text()
    assert 'def: "Cited definition." [GOC:mah, PMID:1\\,2, odd\\"quote, a\\]b, back\\\\slash]' in content
    loaded = load_ontology(path, format="obo").terms["GO:1"]
    assert loaded.definition == "Cited definition."
    assert loaded.metadata["definition_xrefs"] == xrefs


def test_obo_writer_emits_the_dbxref_list_the_format_requires(tmp_path: Path) -> None:
    onto = create_ontology(terms={"GO:1": create_term(id="GO:1", name="term", definition="Needs a list.")})
    path = tmp_path / "plain.obo"
    save_ontology(onto, path, format="obo")
    assert 'def: "Needs a list." []' in path.read_text()


def test_obo_to_obo_copy_keeps_definitions_and_dbxrefs(tmp_path: Path, parsed) -> None:
    path = tmp_path / "copy.obo"
    save_ontology(parsed, path, format="obo")
    copy = load_ontology(path, format="obo")
    for term_id, term in parsed.terms.items():
        assert copy.terms[term_id].definition == term.definition
        assert copy.terms[term_id].metadata.get("definition_xrefs") == term.metadata.get("definition_xrefs")
    # ... and what is preserved is the clean definition and its dbxrefs, not matching residue on both sides
    assert copy.terms["GO:0000001"].definition == 'A process that "x" does.'
    assert copy.terms["GO:0000001"].metadata["definition_xrefs"] == ["GOC:mah", "PMID:123"]


def test_format_obo_xref_list() -> None:
    assert format_obo_xref_list([]) == "[]"
    assert format_obo_xref_list(["GOC:mah", "PMID:123"]) == "[GOC:mah, PMID:123]"
