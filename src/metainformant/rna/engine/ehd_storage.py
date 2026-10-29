"""Validated EHD write-root admission, independent of unused host storage."""
from __future__ import annotations

import math
import os
import plistlib
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path


def _physical(path: Path, *, directory: bool, required: bool = True) -> Path:
    if (not path.is_absolute() or any(p.is_symlink() for p in (path, *path.parents))
            or path.resolve(strict=required) != path):
        raise ValueError(f"EHD path is absent, relative or redirected: {path}")
    if required and not (path.is_dir() if directory else path.is_file()):
        raise ValueError(f"EHD path has the wrong type: {path}")
    return path


def _contained(path: Path, namespace: Path, *, directory: bool, required: bool = True) -> Path:
    value = _physical(path, directory=directory, required=required)
    if value == namespace or not value.is_relative_to(namespace):
        raise ValueError(f"Managed RNA path escapes EHD Data: {path}")
    existing = next(p for p in (value, *value.parents) if p.exists())
    if existing.stat().st_dev != namespace.stat().st_dev:
        raise ValueError(f"RNA path is on a different filesystem: {path}")
    return value


@dataclass(frozen=True, slots=True)
class EHDWritePolicy:
    """Snapshot the validated roots; refuse redirection or policy drift later."""

    namespace: Path
    operations: Path
    bulk: Path
    config: Path
    log: Path
    database: Path
    root_identities: tuple[tuple[Path, int, int], ...]
    environment: tuple[tuple[str, str | None], ...]

    def revalidate(self, output: Path | None = None) -> None:
        if any(os.environ.get(key) != value for key, value in self.environment):
            raise ValueError("EHD write environment changed after admission")
        if not self.namespace.parent.is_mount():
            raise ValueError("EHD volume is no longer mounted")
        for root, device, inode in self.root_identities:
            _physical(root, directory=True)
            if (root.stat().st_dev, root.stat().st_ino) != (device, inode):
                raise ValueError(f"EHD root identity changed: {root}")
        _contained(self.config, self.namespace, directory=True)
        _contained(self.log, self.namespace, directory=True, required=False)
        _contained(self.database, self.namespace, directory=False)
        for key in _WRITE_ROOTS:
            _contained(Path(os.environ[key]), self.namespace, directory=True)
        if Path(tempfile.gettempdir()).resolve(strict=True) != Path(os.environ["TMPDIR"]):
            raise ValueError("Actual Python temporary directory is outside the EHD binding")
        if output is not None:
            _contained(output, self.namespace, directory=True)
            if not output.is_relative_to(self.bulk):
                raise ValueError("Acquisition output escapes the bound EHD bulk root")

    def receipt(self) -> dict[str, object]:
        return {"policy": "validated_ehd_write_roots_v1", "system_capacity_required": False,
                "namespace": str(self.namespace), "operations": str(self.operations),
                "bulk": str(self.bulk), "external_floor_gib": 256}


_WRITE_ROOTS = {
    "HOME": "home", "TMPDIR": "tmp", "UV_CACHE_DIR": "cache/uv",
    "XDG_CACHE_HOME": "cache", "MPLCONFIGDIR": "cache/matplotlib",
    "NUMBA_CACHE_DIR": "cache/numba", "AMALGKIT_TEMP_ROOT": "tmp",
    "AMALGKIT_CACHE_ROOT": "cache", "AMALGKIT_LOCAL_INDEX_CACHE_DIR": "index-cache",
}


def capture_ehd_write_policy(
    config: Path, log: Path, database: Path | None,
) -> EHDWritePolicy | None:
    """Exempt only a fully validated EHD mode from the legacy host-space gate.

    The production launcher additionally applies the independently qualified
    native write sandbox. This function is a path/capacity contract, not a grant
    to stop a producer, migrate data, or launch an unreviewed profile.
    """
    mode = os.environ.get("AMALGKIT_LOCAL_STORAGE_MODE", "ssd_compute")
    if mode not in {"ssd_compute", "hd_first", "ehd_only"}:
        raise ValueError("Unknown local storage mode")
    if mode != "ehd_only":
        return None
    namespace = _physical(Path(os.environ["AMALGKIT_DATA_NAMESPACE"]), directory=True)
    if namespace.name != "Data" or not namespace.parent.is_mount():
        raise ValueError("EHD namespace must be Data on the mounted volume")
    volume = plistlib.loads(subprocess.check_output(
        ["/usr/sbin/diskutil", "info", "-plist", str(namespace.parent)],
        stdin=subprocess.DEVNULL,
    ))
    if (volume.get("MountPoint") != str(namespace.parent) or volume.get("Locked")
            or volume.get("VolumeUUID") != os.environ["AMALGKIT_DATA_VOLUME_UUID"]):
        raise ValueError("EHD mounted UUID/lock binding differs")
    operations = _contained(Path(os.environ["AMALGKIT_OPERATIONS_ROOT"]), namespace, directory=True)
    bulk = _contained(Path(os.environ["AMALGKIT_DATA_ROOT"]), namespace, directory=True)
    if operations.is_relative_to(bulk) or bulk.is_relative_to(operations):
        raise ValueError("EHD control and bulk roots overlap")
    source = _physical(Path(os.environ["METAINFORMANT_ROOT"]), directory=True)
    source_mount = next((p for p in (source, *source.parents) if p.is_mount()), None)
    if source_mount is None or source_mount in {Path("/"), namespace.parent}:
        raise ValueError("Software must remain on its distinct mounted source volume")
    if not all(Path(p).resolve(strict=True).is_relative_to(source_mount)
               for p in (sys.prefix, sys.base_prefix)):
        raise ValueError("Both Python prefixes must remain on the source volume")
    if not sys.dont_write_bytecode or os.environ.get("PYTHONDONTWRITEBYTECODE") != "1":
        raise ValueError("EHD mode requires source bytecode writes disabled")
    for key, relative in _WRITE_ROOTS.items():
        if os.environ.get(key) != str(operations / relative):
            raise ValueError(f"EHD write binding differs: {key}")
        _contained(Path(os.environ[key]), namespace, directory=True)
    if any(os.environ.get(key) for key in (
        "AMALGKIT_LOCAL_QUANT_SCRATCH_DIR", "AMALGKIT_LOCAL_FASTERQ_SCRATCH_DIR",
        "AMALGKIT_LOCAL_ENA_SCRATCH_DIR",
    )) or os.environ.get("AMALGKIT_LOCAL_FASTQ_STAGE") != "no":
        raise ValueError("EHD mode refuses SSD scratch or FASTQ staging")
    if os.environ.get("AMALGKIT_PROGRESS_DB_JOURNAL_MODE") != "DELETE":
        raise ValueError("EHD progress store requires DELETE/FULL/mmap0")
    if database is None or database != Path(os.environ["AMALGKIT_PROGRESS_DB_PATH"]):
        raise ValueError("EHD progress store must be explicit and match the binding")
    if config != Path(os.environ["AMALGKIT_CONFIG_DIR"]):
        raise ValueError("EHD configuration root differs")
    _contained(config, namespace, directory=True)
    _contained(log, namespace, directory=True, required=False)
    if not (log.is_relative_to(operations) or log.is_relative_to(bulk)):
        raise ValueError("RNA logs must be below bound controls or bulk")
    _contained(database, namespace, directory=False)
    for key in ("AMALGKIT_RUNTIME_SOURCE_MANIFEST", "AMALGKIT_TOOLCHAIN_MANIFEST"):
        _contained(Path(os.environ[key]), namespace, directory=False)
    _contained(Path(os.environ["AMALGKIT_FROZEN_ROOT"]), namespace, directory=True)
    overlay = os.environ.get("AMALGKIT_LOCAL_SIZE_OVERLAY_PATH")
    if overlay:
        _contained(Path(overlay), namespace, directory=False)
    minimum = float(os.environ.get("AMALGKIT_MIN_EXTERNAL_FREE_GB", "256"))
    if not math.isfinite(minimum) or minimum < 256:
        raise ValueError("Validated EHD mode retains at least 256 GiB external headroom")
    keys = {*_WRITE_ROOTS, "AMALGKIT_LOCAL_STORAGE_MODE", "AMALGKIT_DATA_NAMESPACE",
            "AMALGKIT_DATA_VOLUME_UUID", "AMALGKIT_OPERATIONS_ROOT", "AMALGKIT_DATA_ROOT",
            "METAINFORMANT_ROOT", "AMALGKIT_PROGRESS_DB_PATH", "AMALGKIT_CONFIG_DIR",
            "AMALGKIT_FROZEN_ROOT", "AMALGKIT_RUNTIME_SOURCE_MANIFEST",
            "AMALGKIT_TOOLCHAIN_MANIFEST", "AMALGKIT_LOCAL_SIZE_OVERLAY_PATH",
            "AMALGKIT_LOCAL_SIZE_OVERLAY_SHA256", "AMALGKIT_PROGRESS_DB_JOURNAL_MODE",
            "AMALGKIT_LOCAL_QUANT_SCRATCH_DIR", "AMALGKIT_LOCAL_FASTERQ_SCRATCH_DIR",
            "AMALGKIT_LOCAL_ENA_SCRATCH_DIR", "AMALGKIT_LOCAL_FASTQ_STAGE",
            "PYTHONDONTWRITEBYTECODE", "AMALGKIT_MIN_EXTERNAL_FREE_GB"}
    policy = EHDWritePolicy(
        namespace, operations, bulk, config, log, database,
        tuple((p, p.stat().st_dev, p.stat().st_ino) for p in (namespace, operations, bulk)),
        tuple((key, os.environ.get(key)) for key in sorted(keys)),
    )
    policy.revalidate()
    return policy
