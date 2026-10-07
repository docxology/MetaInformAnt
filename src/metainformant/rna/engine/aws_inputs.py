"""Immutable source/input archives and shell-safe worker startup rendering."""

from __future__ import annotations
import csv
import hashlib
import json
import shlex
import tarfile
from pathlib import Path
from typing import Any
from metainformant.rna.amalgkit import (
    AMALGKIT_RELEASE_TAG,
    AMALGKIT_SOURCE_REVISION,
    REQUIRED_AMALGKIT_VERSION,
)


def _source_bundle(repo: Path, destination: Path) -> str:
    roots = [
        repo / "src",
        repo / "scripts/rna",
        repo / "projects/hymenoptera_amalgkit/scripts",
        repo / "projects/hymenoptera_amalgkit/config",
        repo / "config/amalgkit",
    ]
    files = [repo / name for name in ("pyproject.toml", "uv.lock", "README.md")]
    for root in roots:
        files.extend(
            p
            for p in root.rglob("*")
            if p.is_file()
            and not p.is_symlink()
            and "__pycache__" not in p.parts
            and not p.name.startswith("._")
            and p.suffix != ".pyc"
        )
    with tarfile.open(destination, "w") as archive:
        for path in sorted(set(files)):
            archive.add(
                path, arcname=path.relative_to(repo).as_posix(), recursive=False
            )
    with tarfile.open(destination) as archive:
        for member in archive:
            if (
                member.islnk()
                or member.issym()
                or Path(member.name).is_absolute()
                or ".." in Path(member.name).parts
            ):
                raise ValueError("unsafe source bundle")
    return hashlib.sha256(destination.read_bytes()).hexdigest()


def _inputs_bundle(
    root: Path, species: dict[str, Any], tasks: list[dict[str, Any]], directory: Path,
    *, config_path: Path | None = None,
) -> tuple[Path, str]:
    directory.mkdir(parents=True, exist_ok=True)
    manifest = directory / "manifest.jsonl"
    manifest.write_text("".join(json.dumps(t, sort_keys=True) + "\n" for t in tasks))
    source = root / "inputs" / species["species"] / "work"
    metadata = source / "metadata" / "metadata_selected.tsv"
    index = source / "index" / species["index_name"]
    if hashlib.sha256(metadata.read_bytes()).hexdigest() != species["metadata_sha256"]:
        raise ValueError("frozen selected metadata changed")
    if hashlib.sha256(index.read_bytes()).hexdigest() != species["index_sha256"]:
        raise ValueError("frozen reference index changed")
    resolved = {t["accession"]: t for t in tasks if t.get("source_evidence_sha256")}
    metadata_override = None
    if resolved:
        with metadata.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields, rows = list(reader.fieldnames or []), list(reader)
        additions = (
            "total_spots",
            "total_bases",
            "size",
            "mean_read_length",
            "source_evidence_sha256",
        )
        fields.extend(name for name in additions if name not in fields)
        seen = set()
        for row in rows:
            task = resolved.get(row["run"])
            if task is not None:
                seen.add(row["run"])
                row.update(
                    total_spots=str(task["total_spots"]),
                    total_bases=str(task["total_bases"]),
                    size=str(task["sra_bytes"]),
                    mean_read_length=str(task["total_bases"] / task["total_spots"]),
                    source_evidence_sha256=task["source_evidence_sha256"],
                )
        if seen != set(resolved):
            raise ValueError("resolved tasks missing from frozen selected metadata")
        metadata_override = directory / "metadata_resolved.tsv"
        with metadata_override.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    records = []
    for path in sorted(source.rglob("*")):
        if not path.is_file():
            continue
        payload_path = (
            metadata_override
            if path == metadata and metadata_override is not None
            else path
        )
        records.append(
            {
                "path": f"data/{species['species']}/work/{path.relative_to(source).as_posix()}",
                "sha256": hashlib.sha256(payload_path.read_bytes()).hexdigest(),
                "size": payload_path.stat().st_size,
            }
        )
    if config_path is not None:
        if config_path.is_symlink():
            raise ValueError("worker configuration must be an owned regular file")
        config_bytes = config_path.read_bytes()
        if hashlib.sha256(config_bytes).hexdigest() != species["config_sha256"]:
            raise ValueError("worker configuration differs from frozen inventory")
        records.append({"path": f"config/amalgkit/{config_path.name}",
                        "sha256": hashlib.sha256(config_bytes).hexdigest(), "size": len(config_bytes)})
    snapshot = {
        "schema": "metainformant.rna.acquisition_snapshot.v1" if config_path is not None else "metainformant.hymenoptera.gcp_snapshot.v1",
        "source_state": "quiescent",
        "cloud_launch_policy": "checkpointed",
        "amalgkit_version": REQUIRED_AMALGKIT_VERSION,
        "amalgkit_release_tag": AMALGKIT_RELEASE_TAG,
        "amalgkit_source_revision": AMALGKIT_SOURCE_REVISION,
        "manifest_sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
        "task_count": len(tasks),
        "species": [species["species"]],
        "input_files": records,
        "raw_reads_included": False,
        "quant_outputs_included": False,
        "source_resolution_tasks": [
            {"task_id": t["task_id"], "evidence_sha256": t["source_evidence_sha256"]}
            for t in resolved.values()
        ],
        "frozen_metadata_sha256": species["metadata_sha256"],
    }
    (directory / "snapshot.json").write_text(
        json.dumps(snapshot, indent=2, sort_keys=True)
    )
    bundle = directory / "inputs.tar"
    with tarfile.open(bundle, "w") as archive:
        archive.add(manifest, arcname="manifest.jsonl")
        archive.add(directory / "snapshot.json", arcname="snapshot.json")
        for record in records:
            if record["path"].startswith("config/"):
                archive.add(config_path, arcname=record["path"], recursive=False)
                continue
            relative = Path(record["path"]).relative_to(
                f"data/{species['species']}/work"
            )
            payload_path = source / relative
            if payload_path == metadata and metadata_override is not None:
                payload_path = metadata_override
            archive.add(payload_path, arcname=record["path"], recursive=False)
    return bundle, hashlib.sha256(bundle.read_bytes()).hexdigest()


def _render_startup(template: Path, replacements: dict[str, Any]) -> str:
    script = template.read_text()
    for key, value in replacements.items():
        script = script.replace(f"@@{key}@@", shlex.quote(str(value)))
    if "@@" in script or len(script.encode()) > 16384:
        raise ValueError("invalid or oversized EC2 user data")
    return script
