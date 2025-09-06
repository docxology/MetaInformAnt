from __future__ import annotations

from pathlib import Path
from typing import Optional

from metainformant.core import io as core_io
from metainformant.core.paths import expand_and_resolve

from .obo import parse_obo
from .types import Ontology


def count_go_scripts(go_dir: Path) -> int:
    return sum(1 for p in go_dir.glob("*.py") if p.is_file())


def load_go_obo(path: str | Path) -> Ontology:
    """Load Gene Ontology from an OBO file into an `Ontology`.

    This is a lightweight reader tailored for tests and small workflows.
    """
    return parse_obo(path)


def write_go_summary(onto: Ontology, dest: str | Path | None = None) -> Path:
    """Write a small JSON summary with counts under output/ by default.

    - If `dest` is None, writes to output/ontology/go_summary.json
    """
    if dest is None:
        dest = Path("output/ontology/go_summary.json")
    dest = expand_and_resolve(dest)
    core_io.dump_json({"num_terms": onto.num_terms()}, dest, indent=2)
    return dest
