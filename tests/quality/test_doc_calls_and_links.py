"""Quoted calls must fit current signatures, and relative links must resolve.

Zero-mocks: documents are real files, calls are bound against the installed
package's real signatures, and links are resolved against real directory trees.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from metainformant.quality.doc_verification import (
    call_signature_mismatch,
    find_broken_doc_links,
    find_stale_doc_calls,
    markdown_anchors,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
SCAN_ROOTS = ("docs", "src", "tests", "scripts", "config", ".github")


def _call(source: str) -> ast.Call:
    node = ast.parse(source).body[0]
    assert isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
    return node.value


def _documented(a: int, b: int, *, c: int = 1) -> int:
    return a + b + c


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        ("f(1, 2)", None),
        ("f(1, b=2, c=3)", None),
        ("f(1)", "missing a required argument: 'b'"),
        ("f(1, 2, 3)", "too many positional arguments"),
        ("f(1, 2, d=4)", "got an unexpected keyword argument 'd'"),
        ("f(*values)", None),  # unpacking hides the shape
        ("f(1, **options)", None),
    ],
)
def test_call_shapes_are_bound_against_the_signature(source: str, expected: str | None) -> None:
    assert call_signature_mismatch(_documented, _call(source)) == expected


def test_calling_a_non_callable_is_reported() -> None:
    assert call_signature_mismatch(42, _call("value()")) == "'int' object is not callable"


def test_calls_on_imported_names_and_chains_are_checked(tmp_path: Path) -> None:
    (tmp_path / "guide.md").write_text(
        "```python\n"
        "from metainformant.core import io\n"
        "from metainformant.core.io import paths, load_json\n"
        "io.dump_json({'a': 1}, 'out.json', indent=2)\n"
        "load_json('in.json', encoding='utf-8')\n"
        "paths.expand_and_resolve()\n"
        "io.paths.expand_and_resolve('~', strict=True)\n"
        "result = io.load_json('in.json')\n"
        "result.get('key', 1, 2, 3)\n"
        "```\n"
        "```python\n"
        "from metainformant.core.io import load_json\n"
        "load_json = dict\n"
        "load_json(a=1)\n"
        "```\n"
        "```python\n"
        "from metainformant.core import io\n"
        "io.load_json(  # pseudo-code that does not parse is skipped\n"
        "```\n",
        encoding="utf-8",
    )

    violations = find_stale_doc_calls([tmp_path])

    assert [(v.line_number, v.code_example, v.issue_type) for v in violations] == [
        (5, "load_json(...)", "StaleCall"),
        (6, "paths.expand_and_resolve(...)", "StaleCall"),
        (7, "io.paths.expand_and_resolve(...)", "StaleCall"),
    ]
    assert "unexpected keyword argument 'encoding'" in violations[0].details
    assert "missing a required argument: 'path'" in violations[1].details


def test_markdown_anchors_follow_github_heading_slugs() -> None:
    lines = [
        "# The `load_json` helper",
        "## Fine-Mapping & Colocalization",
        "## Usage",
        "## Usage",
        "## _Emphasised_ heading",
        "```bash",
        "# not a heading",
        "```",
        '<a id="Custom-Target"></a>',
    ]

    anchors = markdown_anchors(lines)

    assert {"the-load_json-helper", "fine-mapping--colocalization", "usage", "usage-1"} <= anchors
    assert "emphasised-heading" in anchors
    assert "custom-target" in anchors
    assert "not-a-heading" not in anchors


def test_relative_links_images_and_anchors_are_resolved(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    (repo / "docs" / "guide").mkdir(parents=True)
    (repo / "docs" / "Notes.md").write_text(
        '# Notes\n\n## Fine-Mapping & Colocalization\n\n<a id="custom"></a>\n', encoding="utf-8"
    )
    (repo / "docs" / "guide" / "index.md").write_text(
        "# Guide\n"
        "[ok](../Notes.md) [anchor](../Notes.md#fine-mapping--colocalization) [id](<../Notes.md#custom>)\n"
        "[dir](../guide/) [self](#guide) [web](https://example.org/x) [mail](mailto:someone@example.org)\n"
        "[missing](../missing.md) [case](../notes.md) [stale anchor](../Notes.md#fine-mapping-colocalization)\n"
        "`[quoted](nowhere.md)` and [template]({module}.md)\n"
        "```markdown\n"
        "[fenced](nowhere.md)\n"
        "```\n"
        "![figure](../figure.png) [escape](../../../elsewhere.md)\n"
        "[reference]: ./also-missing.md\n"
        "[^1]: Footnote text is not a link target.\n",
        encoding="utf-8",
    )

    violations = find_broken_doc_links([repo / "docs"], repo_root=repo)

    assert [(v.line_number, v.code_example) for v in violations] == [
        (4, "../missing.md"),
        (4, "../notes.md"),
        (4, "../Notes.md#fine-mapping-colocalization"),
        (9, "../figure.png"),
        (9, "../../../elsewhere.md"),
        (10, "./also-missing.md"),
    ]
    details = {v.code_example: v.details for v in violations}
    assert "on disk: 'Notes.md'" in details["../notes.md"]
    assert "no heading or id" in details["../Notes.md#fine-mapping-colocalization"]
    assert "outside the repository" in details["../../../elsewhere.md"]


def test_links_into_an_uninitialized_submodule_are_not_judged(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    (repo / "projects" / "empty").mkdir(parents=True)
    (repo / "projects" / "present").mkdir()
    (repo / "projects" / "present" / "README.md").write_text("# Present\n", encoding="utf-8")
    (repo / ".gitmodules").write_text(
        '[submodule "projects/empty"]\n\tpath = projects/empty\n'
        '[submodule "projects/present"]\n\tpath = projects/present\n',
        encoding="utf-8",
    )
    (repo / "README.md").write_text(
        "[deferred](projects/empty/README.md) [checked](projects/present/MISSING.md)\n", encoding="utf-8"
    )

    violations = find_broken_doc_links([repo / "README.md"], repo_root=repo)

    assert [v.code_example for v in violations] == ["projects/present/MISSING.md"]


def _scan_roots() -> list[Path]:
    roots = [REPO_ROOT / name for name in SCAN_ROOTS] + sorted(REPO_ROOT.glob("*.md"))
    return [root for root in roots if root.exists()]


def test_active_documentation_links_resolve() -> None:
    violations = find_broken_doc_links(_scan_roots(), repo_root=REPO_ROOT)
    assert not violations, "broken documentation links:\n" + "\n".join(map(str, violations))


@pytest.mark.slow
@pytest.mark.timeout(300)  # first import of most of the package; ~60 s under coverage
def test_active_documentation_calls_fit_current_signatures() -> None:
    violations = find_stale_doc_calls(_scan_roots())
    assert not violations, "stale documentation calls:\n" + "\n".join(map(str, violations))
