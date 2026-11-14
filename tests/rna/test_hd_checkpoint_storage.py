"""Real quantification/restore controls for the HD-first temporary workspace."""

from __future__ import annotations

import dataclasses
import importlib.util
import json
import shutil
import sys
from collections.abc import Callable
from pathlib import Path
from types import ModuleType
from typing import TYPE_CHECKING

import pytest

from metainformant.rna.engine.acquisition_manifest import sha256_file
from tests.rna.test_ehd_runtime_paths import require_operator_files

pytestmark = pytest.mark.external_tool

if TYPE_CHECKING:
    from tests.rna.conftest import NativeAcquisitionFiles


@pytest.fixture
def checkpoint_module() -> ModuleType:
    path = require_operator_files("checkpoints.py") / "checkpoints.py"
    spec = importlib.util.spec_from_file_location("hd_checkpoint_qualification", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.mark.external_tool
@pytest.mark.parametrize("bulk_restore", [False, True])
def test_real_checkpoint_preserves_restore_contract_and_raw_order(
    tmp_path: Path,
    native_acquisition_factory: Callable[[], NativeAcquisitionFiles],
    checkpoint_module: ModuleType,
    bulk_restore: bool,
    monkeypatch: pytest.MonkeyPatch,
    staging_bulk_dir: Callable[[str], Path],
) -> None:
    # Given a complete native input and an isolated operations/control tree.
    # The bulk copy is allocated (and skipped when the bulk root is unset)
    # before any work; the fixture removes it at teardown and fails if it remains.
    bulk_data = staging_bulk_dir("hd-checkpoint") if bulk_restore else None
    fixture = native_acquisition_factory()
    data, work = fixture.data, fixture.work
    if bulk_data is not None:
        parent = bulk_data.parent
        assert parent.is_dir() and parent.stat().st_dev != tmp_path.stat().st_dev
        data = bulk_data
        shutil.copytree(fixture.data, data)
        work = data / "test_species/work"
        monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data))
        fixture = dataclasses.replace(fixture, data=data, work=work)
    operations = tmp_path / "operations"
    (operations / "tmp").mkdir(parents=True)
    bindings = {
        "test_species/SRR123": checkpoint_module.Binding(
            fixture.config,
            sha256_file(fixture.config),
            sha256_file(fixture.index),
        )
    }
    checkpoint = checkpoint_module.SampleCheckpoint(
        data,
        operations,
        bindings,
        bulk_restore=bulk_restore,
    )
    # When the producer quantifies with Kallisto and the output is locked, restored and checked.
    result = fixture.produce()
    assert result["quantified"] and not result["skipped"], result
    checkpoint("test_species", "SRR123")
    # Then the same restore proof is published and no restore payload stays on SSD.
    proof = json.loads((work / "quant/SRR123/.local_restore_verified.json").read_text())
    assert proof["reference_index_sha256"] == sha256_file(fixture.index)
    assert not list((operations / "tmp").iterdir())
    assert (data / ".checkpoint-verify").is_dir() is bulk_restore
    if bulk_restore:
        assert not list((data / ".checkpoint-verify").iterdir())
        assert (data / ".checkpoint-verify").stat().st_dev != operations.stat().st_dev
    assert (operations / "raw_reclamation.jsonl").is_file()


def test_hd_checkpoint_refuses_symlink_before_sealing(tmp_path: Path, checkpoint_module: ModuleType) -> None:
    # Given an HD workspace redirected outside the selected bulk root.
    data, operations, outside = tmp_path / "data", tmp_path / "operations", tmp_path / "outside"
    data.mkdir()
    outside.mkdir()
    (data / ".checkpoint-verify").symlink_to(outside, target_is_directory=True)
    binding = checkpoint_module.Binding(tmp_path / "absent-config", "config", "index")
    checkpoint = checkpoint_module.SampleCheckpoint(
        data,
        operations,
        {"species/SRR123": binding},
        bulk_restore=True,
    )
    # When checkpointing attempts to admit the altered workspace.
    with pytest.raises(checkpoint_module.CheckpointFailure, match="must not be a symlink"):
        checkpoint("species", "SRR123")
    # Then no external payload or locked output is published.
    assert not list(outside.iterdir())
    assert not (data / "locked-quant").exists()
