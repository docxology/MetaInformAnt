"""Create and stage portable acquisition envelopes from a frozen inventory."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path

from metainformant.rna.amalgkit import AMALGKIT_RELEASE_TAG, AMALGKIT_SOURCE_REVISION, REQUIRED_AMALGKIT_VERSION
from metainformant.rna.engine.acquisition_manifest import load_snapshot, sha256_file
from metainformant.rna.engine.campaign_status import load_inventory
from metainformant.rna.engine.quant_storage import DirectoryStore


def create_campaign_manifest(campaign_root: Path) -> Path:
    """Publish a generic immutable manifest beside the frozen metadata/index inputs."""
    encoded = (campaign_root / "inventory.json").read_bytes()
    load_inventory(encoded).task_ids()
    inventory = json.loads(encoded)
    tasks = [
        {**task, "config_sha256": species["config_sha256"]}
        for species in inventory["species"]
        for task in species["tasks"]
    ]
    manifest = ("".join(json.dumps(task, sort_keys=True) + "\n" for task in tasks)).encode()
    records = []
    for species in inventory["species"]:
        name = species["index_name"]
        if Path(name).name != name or "\\" in name or name in {"", ".", ".."}:
            raise ValueError("unsafe frozen index name")
        for relative, checksum in (
            ("metadata/metadata_selected.tsv", species["metadata_sha256"]),
            (f"index/{name}", species["index_sha256"]),
        ):
            records.append({"path": f"inputs/{species['species']}/work/{relative}", "sha256": checksum})
    snapshot = {
        "schema": "metainformant.rna.acquisition_snapshot.v1",
        "cloud_launch_policy": "checkpointed",
        "inventory_sha256": hashlib.sha256(encoded).hexdigest(),
        "manifest_sha256": hashlib.sha256(manifest).hexdigest(),
        "task_count": len(tasks),
        "species": [species["species"] for species in inventory["species"]],
        "input_files": records,
        "amalgkit_version": REQUIRED_AMALGKIT_VERSION,
        "amalgkit_release_tag": AMALGKIT_RELEASE_TAG,
        "amalgkit_source_revision": AMALGKIT_SOURCE_REVISION,
    }
    store = DirectoryStore(campaign_root)
    store.put("manifest.jsonl", manifest)
    store.put("snapshot.json", json.dumps(snapshot, sort_keys=True, indent=2).encode() + b"\n")
    return campaign_root / "manifest.jsonl"


def stage_worker_inputs(manifest: Path, data_root: Path, species: frozenset[str]) -> None:
    """Copy only selected species' frozen metadata/index files; never replace conflicts."""
    snapshot, _ = load_snapshot(manifest)
    source_root = manifest.parent.resolve()
    destination_root = data_root.resolve()
    for record in snapshot["input_files"]:
        relative = Path(record["path"])
        parts = relative.parts
        if len(parts) == 3 and parts[:2] == ("config", "amalgkit"):
            source = source_root / relative
            if (
                source.is_symlink()
                or not source.resolve().is_relative_to(source_root)
                or sha256_file(source) != record["sha256"]
            ):
                raise ValueError("frozen configuration input is unavailable or differs from checksum")
            continue
        if (
            relative.is_absolute()
            or ".." in parts
            or len(parts) != 5
            or parts[0] not in {"inputs", "data"}
            or parts[2] != "work"
        ):
            raise ValueError("unsupported staged acquisition input path")
        if parts[1] not in species:
            continue
        source = source_root / relative
        if (
            source.is_symlink()
            or not source.resolve().is_relative_to(source_root)
            or sha256_file(source) != record["sha256"]
        ):
            raise ValueError("frozen input is unowned or differs from its checksum")
        target = destination_root.joinpath(*parts[1:])
        if not target.resolve().is_relative_to(destination_root):
            raise ValueError("local input target escapes data root")
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            if target.is_symlink() or sha256_file(target) != record["sha256"]:
                raise FileExistsError(f"staging would replace different local input: {target}")
            continue
        fd, name = tempfile.mkstemp(prefix=".acquisition-stage-", dir=target.parent)
        temporary = Path(name)
        try:
            with os.fdopen(fd, "wb") as output, source.open("rb") as input_file:
                shutil.copyfileobj(input_file, output)
                output.flush()
                os.fsync(output.fileno())
            if sha256_file(temporary) != record["sha256"]:
                raise ValueError("input changed while staging")
            try:
                os.link(temporary, target)
            except FileExistsError:
                if target.is_symlink() or sha256_file(target) != record["sha256"]:
                    raise FileExistsError(f"concurrent staging conflict: {target}") from None
        finally:
            temporary.unlink(missing_ok=True)
