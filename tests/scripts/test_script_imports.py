"""Every ``metainformant`` import in ``scripts/`` must resolve against the package.

Scripts are not otherwise exercised by the suite, so a renamed or moved module
breaks them silently. Zero-mocks: each import statement is resolved exactly as
Python would import it, without running the scripts.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from metainformant.quality.doc_verification import resolve_import_statement

REPO_ROOT = Path(__file__).resolve().parents[2]

# Scripts whose functionality no longer exists in the package. Each entry must
# still fail to import; once a script resolves again, remove it from this set.
OBSOLETE_SCRIPTS: frozenset[str] = frozenset(
    {
        # rna.engine.discovery returns per-study records; no species-level aggregator exists.
        "scripts/rna/discover_species.py",
        # orchestrate_genome_setup no longer honours dry_run, the skip flags or kmer_size.
        "scripts/rna/setup_genome.py",
    }
)


def _script_import_failures() -> dict[str, list[str]]:
    failures: dict[str, list[str]] = {}
    for path in sorted((REPO_ROOT / "scripts").rglob("*.py")):
        relative = path.relative_to(REPO_ROOT).as_posix()
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"), filename=relative)):
            if not isinstance(node, (ast.Import, ast.ImportFrom)) or getattr(node, "level", 0):
                continue
            module = node.module if isinstance(node, ast.ImportFrom) else node.names[0].name
            if not module or module.split(".")[0] != "metainformant":
                continue
            reason = resolve_import_statement(ast.unparse(node))
            if reason is not None:
                failures.setdefault(relative, []).append(f"line {node.lineno}: {reason}")
    return failures


@pytest.mark.slow
@pytest.mark.timeout(300)  # first import of most of the package under coverage
def test_script_imports_resolve() -> None:
    failures = _script_import_failures()

    unexpected = {script: reasons for script, reasons in failures.items() if script not in OBSOLETE_SCRIPTS}
    assert unexpected == {}, "scripts with unresolvable metainformant imports:\n" + "\n".join(
        f"  {script}: {'; '.join(reasons)}" for script, reasons in sorted(unexpected.items())
    )
    revived = sorted(OBSOLETE_SCRIPTS - failures.keys())
    assert revived == [], f"these scripts import cleanly again; remove them from OBSOLETE_SCRIPTS: {revived}"
