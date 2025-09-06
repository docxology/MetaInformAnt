from __future__ import annotations

from pathlib import Path
from typing import Iterable, Iterator, List

from .types import Ontology, Term


def _iter_stanzas(lines: Iterable[str]) -> Iterator[List[str]]:
    stanza: List[str] = []
    for raw in lines:
        line = raw.strip("\n")
        if not line:
            # keep blank lines inside stanza
            stanza.append("")
            continue
        if line.startswith("[") and line.endswith("]"):
            # yield previous stanza if any
            if stanza:
                yield stanza
                stanza = []
            stanza.append(line)
        else:
            stanza.append(line)
    if stanza:
        yield stanza


def parse_obo(path: str | Path) -> Ontology:
    """Parse a small subset of OBO format sufficient for GO term hierarchy.

    Supports fields: id, name, namespace, def, alt_id, is_a (parent).
    """
    p = Path(path)
    ontology = Ontology()
    with p.open("rt", encoding="utf-8") as fh:
        for stanza in _iter_stanzas(fh):
            if not stanza:
                continue
            if stanza[0] != "[Term]":
                continue

            term_id: str | None = None
            name: str | None = None
            namespace: str | None = None
            definition: str | None = None
            alt_ids: List[str] = []
            parents: List[str] = []

            for line in stanza[1:]:
                if not line:
                    continue
                if line.startswith("id: "):
                    term_id = line[4:].strip()
                elif line.startswith("name: "):
                    name = line[6:].strip()
                elif line.startswith("namespace: "):
                    namespace = line[11:].strip()
                elif line.startswith("def: "):
                    definition = line[5:].strip()
                elif line.startswith("alt_id: "):
                    alt_ids.append(line[8:].strip())
                elif line.startswith("is_a: "):
                    parents.append(line[6:].split(" ")[0].strip())

            if term_id and name:
                term = Term(
                    term_id=term_id,
                    name=name,
                    namespace=namespace,
                    definition=definition,
                    alt_ids=alt_ids,
                    is_a_parents=parents,
                )
                ontology.add_term(term)
    return ontology
