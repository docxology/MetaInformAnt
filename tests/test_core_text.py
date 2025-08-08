from __future__ import annotations

from metainformant.core import text as core_text


def test_slugify_basic() -> None:
    assert core_text.slugify("Hello, World!") == "hello-world"
    assert core_text.slugify("  A  B  C  ") == "a-b-c"


def test_normalize_whitespace() -> None:
    s = "a\t\t b\n\n c"
    assert core_text.normalize_whitespace(s) == "a b c"


def test_safe_filename() -> None:
    name = "My*Weird:File?.txt"
    safe = core_text.safe_filename(name)
    assert all(ch not in safe for ch in ['*', ':', '?'])
    assert safe.endswith(".txt")


