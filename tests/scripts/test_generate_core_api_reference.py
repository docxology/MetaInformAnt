"""The generated core API reference stays current and its index links resolve.

Zero-mocks: the generator script runs on the real ``src/metainformant/core``
tree, and its anchors are checked with the same resolver the link gate uses.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

from metainformant.quality.doc_verification import find_broken_doc_links, markdown_anchors

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "package" / "generate_core_api_reference.py"
REFERENCE = REPO_ROOT / "docs" / "core" / "API_REFERENCE.md"

_spec = importlib.util.spec_from_file_location("generate_core_api_reference", SCRIPT)
assert _spec is not None and _spec.loader is not None
generator = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = generator  # dataclasses resolve annotations through sys.modules
_spec.loader.exec_module(generator)


def test_github_anchor_matches_the_link_gate_heading_slug() -> None:
    for heading in ("metainformant.core.io.download_manager", "Fine-Mapping & Colocalization", "Plain heading"):
        assert generator.github_anchor(heading) in markdown_anchors([f"## {heading}"])


def test_rendered_index_links_point_at_rendered_headings(tmp_path: Path) -> None:
    output = tmp_path / "API_REFERENCE.md"
    (tmp_path / "index.md").write_text("# Core\n", encoding="utf-8")

    assert generator.main(["--output", str(output)]) == 0

    assert find_broken_doc_links([output], repo_root=tmp_path) == []


def test_committed_reference_is_current() -> None:
    assert (
        generator.main(["--check", "--output", str(REFERENCE)]) == 0
    ), "docs/core/API_REFERENCE.md is stale; run python scripts/package/generate_core_api_reference.py"
