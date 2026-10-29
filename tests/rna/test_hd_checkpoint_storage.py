"""Real quantification/restore controls for the HD-first temporary workspace."""

from __future__ import annotations

import importlib.util
import json
import os
import shutil
import sys
from collections.abc import Callable
from pathlib import Path
from types import ModuleType
from typing import TYPE_CHECKING
from uuid import uuid4

import pytest

from metainformant.rna.engine.acquisition_manifest import sha256_file
from metainformant.rna.engine.acquisition_worker import run_manifest
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
) -> None:
    # Given a complete native input and an isolated operations/control tree.
    if bulk_restore and "AMALGKIT_STAGING_TEST_BULK_ROOT" not in os.environ:
        pytest.skip("AMALGKIT_STAGING_TEST_BULK_ROOT is required for cross-volume checkpoint qualification")
    fixture = native_acquisition_factory()
    data, work = fixture.data, fixture.work
    if bulk_restore:
        parent = Path(os.environ["AMALGKIT_STAGING_TEST_BULK_ROOT"])
        assert parent.is_dir() and parent.stat().st_dev != tmp_path.stat().st_dev
        data = parent / ("hd-checkpoint-" + uuid4().hex)
        shutil.copytree(fixture.data, data)
        work = data / "test_species/work"
        monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data))
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
    # When actual Kallisto output is locked, restored and checked before cleanup.
    result = run_manifest(
        manifest_path=fixture.manifest,
        data_root=data,
        config_dir=fixture.config_dir,
        workers=1,
        threads=1,
        fastq_threads=1,
        compression_threads=1,
        validation_slots=1,
        quant_slots=1,
        fasterq_slots=1,
        max_in_flight=1,
        on_quantified=checkpoint,
    )
    # Then the same restore proof is published and no restore payload stays on SSD.
    assert result["counts"]["newly_quantified"] == 1
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
