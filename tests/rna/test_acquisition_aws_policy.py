"""Platform and frozen configuration controls, using real files and values."""

from __future__ import annotations

import hashlib
import json
import tarfile
from pathlib import Path

import pytest

from metainformant.rna.engine.acquisition_aws_policy import (
    WorkerImage,
    validate_worker_image,
)
from metainformant.rna.engine.acquisition_manifest import (
    verify_input_files,
    verify_worker_configs,
)
from metainformant.rna.engine.aws_inputs import _inputs_bundle


def test_default_bootstrap_supports_unlicensed_linux_and_custom_arm_template() -> None:
    validate_worker_image(WorkerImage("available", "x86_64", "Linux/UNIX", False))
    validate_worker_image(WorkerImage("available", "arm64", "Linux/UNIX", False), custom_template=True)


@pytest.mark.parametrize(
    "image",
    [
        WorkerImage("pending", "x86_64", "Linux/UNIX", False),
        WorkerImage("available", "arm64", "Linux/UNIX", False),
        WorkerImage("available", "x86_64", "Windows", False),
        WorkerImage("available", "x86_64", "Linux/UNIX", True),
    ],
)
def test_unsupported_worker_images_refused(image: WorkerImage) -> None:
    with pytest.raises(ValueError):
        validate_worker_image(image)


def test_generic_job_binds_and_stages_worker_configuration(tmp_path: Path) -> None:
    root = tmp_path / "campaign"
    source = root / "inputs/ant_a/work"
    metadata = source / "metadata/metadata_selected.tsv"
    metadata.parent.mkdir(parents=True)
    metadata.write_text("run\tscientific_name\nSRR1\tAnt a\n")
    index = source / "index/Ant_a.idx"
    index.parent.mkdir()
    index.write_bytes(b"checksum fixture")
    configs = tmp_path / "config"
    configs.mkdir()
    config = configs / "amalgkit_ant_a.yaml"
    config.write_text("species_list: [Ant_a]\n")

    def digest(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    species = {
        "species": "ant_a",
        "index_name": index.name,
        "metadata_sha256": digest(metadata),
        "index_sha256": digest(index),
        "config_sha256": digest(config),
    }
    task = {
        "schema": "metainformant.rna.acquisition_task.v1",
        "task_id": "ant_a/SRR1",
        "species": "ant_a",
        "accession": "SRR1",
        "config_name": config.name,
        "batch_index": 1,
    }
    directory = tmp_path / "job"
    bundle, _ = _inputs_bundle(root, species, [task], directory, config_path=config)
    unpack = tmp_path / "unpack"
    unpack.mkdir()
    with tarfile.open(bundle) as archive:
        archive.extractall(unpack, filter="data")
    snapshot = json.loads((unpack / "snapshot.json").read_text())
    assert snapshot["schema"] == "metainformant.rna.acquisition_snapshot.v1"
    verify_input_files(snapshot, unpack)
    verify_worker_configs(snapshot, [task], configs)
    assert (unpack / "config/amalgkit" / config.name).read_bytes() == config.read_bytes()
    config.write_text("species_list: [Different]\n")
    with pytest.raises(ValueError, match="configuration"):
        verify_worker_configs(snapshot, [task], configs)
    with pytest.raises(ValueError, match="configuration"):
        _inputs_bundle(root, species, [task], tmp_path / "bad-job", config_path=config)
