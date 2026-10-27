"""Every ``metainformant`` import quoted in active documentation must resolve.

Zero-mocks: documents are real files and every statement is resolved against the
installed package exactly as Python would import it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.quality.doc_verification import (
    collect_doc_import_statements,
    find_stale_doc_attributes,
    find_stale_doc_imports,
    resolve_import_statement,
)

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_collects_only_complete_active_import_lines(tmp_path: Path) -> None:
    (tmp_path / "guide.md").write_text(
        "Use `from metainformant.core import io` inline.\n"
        "```python\n"
        "from metainformant.core import io\n"
        "import metainformant.core.utils.logging as log\n"
        "from metainformant.dna import ...\n"
        "```\n"
        "```python\n"
        "from metainformant.core.utils import (\n"
        "    config,\n"
        ")\n"
        "from metainformant.core.io import paths  # path helpers\n"
        "```\n",
        encoding="utf-8",
    )
    (tmp_path / "old.md").write_text("> Historical snapshot: kept for provenance.\nfrom metainformant.gone import x\n")
    (tmp_path / "archive").mkdir()
    (tmp_path / "archive" / "report.md").write_text("from metainformant.gone import y\n", encoding="utf-8")

    statements = collect_doc_import_statements([tmp_path])

    assert sorted(statements) == [
        "from metainformant.core import io",
        "from metainformant.core.io import paths",
        "from metainformant.core.utils import config",
        "import metainformant.core.utils.logging as log",
    ]
    assert statements["from metainformant.core import io"] == [f"{tmp_path / 'guide.md'}:3"]
    assert statements["from metainformant.core.utils import config"] == [f"{tmp_path / 'guide.md'}:8"]
    assert statements["from metainformant.core.io import paths"] == [f"{tmp_path / 'guide.md'}:11"]
    assert "from metainformant.gone import x" in collect_doc_import_statements([tmp_path], include_historical=True)


@pytest.mark.parametrize(
    ("statement", "expected"),
    [
        ("from metainformant.core import io", None),
        ("from metainformant.core.io import paths", None),
        ("from metainformant.core.utils import config, logging", None),
        ("import metainformant.core.utils.logging as log", None),
        ("from metainformant.core import not_a_real_name", "ImportError: cannot import name 'not_a_real_name'"),
        ("from metainformant.not_a_package import thing", "ModuleNotFoundError"),
        ("from metainformant.core import io; import os", "not a single absolute import statement"),
    ],
)
def test_statements_resolve_like_python_imports(statement: str, expected: str | None) -> None:
    reason = resolve_import_statement(statement)
    if expected is None:
        assert reason is None
    else:
        assert reason is not None and reason.startswith(expected)


def test_each_location_of_a_stale_import_is_reported(tmp_path: Path) -> None:
    for name in ("a.md", "b.md"):
        (tmp_path / name).write_text("intro\nfrom metainformant.core import missing_helper\n", encoding="utf-8")

    violations = find_stale_doc_imports([tmp_path])

    assert [(Path(v.doc_file).name, v.line_number, v.issue_type) for v in violations] == [
        ("a.md", 2, "StaleImport"),
        ("b.md", 2, "StaleImport"),
    ]


@pytest.mark.slow
@pytest.mark.timeout(300)  # first import of most of the package; ~60 s under coverage
def test_active_documentation_imports_resolve() -> None:
    roots = [REPO_ROOT / name for name in ("docs", "src", "tests", "scripts", "config")]
    roots += sorted(REPO_ROOT.glob("*.md"))
    violations = find_stale_doc_imports([root for root in roots if root.exists()])
    assert not violations, "stale documentation imports:\n" + "\n".join(map(str, violations))


def test_attribute_chains_on_imported_modules_and_classes_are_checked(tmp_path: Path) -> None:
    (tmp_path / "guide.md").write_text(
        "```python\n"
        "from metainformant.core import io\n"
        "from metainformant.core.io import paths\n"
        "io.paths.expand_and_resolve('~')\n"
        "io.no_such_helper.run()\n"
        "io.no_such_helper.other()\n"
        "result = io.load_json('x.json')\n"
        "result.anything_goes\n"
        "```\n"
        "```python-snippet\n"
        "from metainformant.core import io\n"
        "io = object()\n"
        "io.rebound_names_are_not_checked\n"
        "```\n"
        "```python\n"
        "from metainformant.core import io\n"
        "io.missing(  # pseudo-code that does not parse is skipped\n"
        "```\n",
        encoding="utf-8",
    )

    violations = find_stale_doc_attributes([tmp_path])

    assert [(v.code_example, v.issue_type) for v in violations] == [("io.no_such_helper", "StaleAttribute")]
    assert "has no attribute 'no_such_helper'" in violations[0].details


@pytest.mark.slow
@pytest.mark.timeout(300)  # first import of most of the package; ~60 s under coverage
def test_active_documentation_attribute_chains_resolve() -> None:
    roots = [REPO_ROOT / name for name in ("docs", "src", "tests", "scripts", "config")]
    roots += sorted(REPO_ROOT.glob("*.md"))
    violations = find_stale_doc_attributes([root for root in roots if root.exists()])
    assert not violations, "stale documentation attribute chains:\n" + "\n".join(map(str, violations))
