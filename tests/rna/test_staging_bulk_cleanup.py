"""Cross-volume test pilots never leave directories under the operator's bulk root."""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from tests.rna.conftest import STAGING_BULK_ROOT_VARIABLE, bulk_staging_allocator


def _read_only_tree(root: Path) -> None:
    (root / "locked-quant" / "SRR1").mkdir(parents=True)
    (root / "locked-quant" / "SRR1" / "abundance.tsv").write_text("target_id\ttpm\n")
    os.chmod(root / "locked-quant" / "SRR1" / "abundance.tsv", stat.S_IRUSR)
    os.chmod(root / "locked-quant" / "SRR1", stat.S_IRUSR | stat.S_IXUSR)
    os.chmod(root / "locked-quant", stat.S_IRUSR | stat.S_IXUSR)


def test_allocated_trees_are_removed_even_when_read_only(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    bulk = tmp_path / "bulk"
    bulk.mkdir()
    (bulk / "operator-data").mkdir()
    monkeypatch.setenv(STAGING_BULK_ROOT_VARIABLE, str(bulk))

    with bulk_staging_allocator() as allocate:
        created = allocate("hd-checkpoint")
        created.mkdir()
        _read_only_tree(created)
        never_created = allocate("ssd-publication-pilot")
        assert created.parent == never_created.parent == bulk
        assert created.name != never_created.name

    assert sorted(path.name for path in bulk.iterdir()) == ["operator-data"]


def test_allocated_trees_are_removed_when_the_test_body_fails(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    bulk = tmp_path / "bulk"
    bulk.mkdir()
    monkeypatch.setenv(STAGING_BULK_ROOT_VARIABLE, str(bulk))

    with pytest.raises(RuntimeError, match="interrupted pilot"):
        with bulk_staging_allocator() as allocate:
            allocate("ssd-publication-pilot").mkdir()
            raise RuntimeError("interrupted pilot")

    assert list(bulk.iterdir()) == []


def test_allocation_skips_without_a_declared_bulk_root(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(STAGING_BULK_ROOT_VARIABLE, raising=False)
    with pytest.raises(pytest.skip.Exception):
        with bulk_staging_allocator() as allocate:
            allocate("hd-checkpoint")
