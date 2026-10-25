"""Execute one frozen acquisition task using shared local and durable methods."""

from __future__ import annotations

import stat
from pathlib import Path
from typing import Any

from metainformant.rna.engine.acquisition_lanes import Result, Task
from metainformant.rna.engine.acquisition_manifest import sha256_file
from metainformant.rna.engine.durable_quant import ObjectStore, bound_receipt_key
from metainformant.rna.engine.provenance import classify_quantification, read_quant_provenance
from metainformant.rna.engine.streaming_orchestrator import (
    StreamingPipelineOrchestrator,
)


def _regular_owned_file(path: Path, root: Path) -> bool:
    try:
        resolved = path.resolve(strict=True)
        return (
            stat.S_ISREG(path.lstat().st_mode)
            and resolved == path.absolute()
            and resolved.is_relative_to(root.resolve(strict=True))
        )
    except (OSError, RuntimeError):
        return False


def reuse_manifest_quantification(task: Task, data_root: Path, config_dir: Path) -> Result | None:
    species, accession = str(task["species"]), str(task["accession"])
    sample = data_root / species / "work/quant" / accession
    sidecar = sample / ".metainformant_quant_provenance.json"
    if not _regular_owned_file(sidecar, data_root):
        return None
    provenance = read_quant_provenance(sample)
    if provenance is None or provenance.get("species") != species:
        return None
    name = provenance.get("quantification_file")
    if not isinstance(name, str) or Path(name).name != name or name in ("", ".", ".."):
        return None
    payloads = [sample / name, *sample.glob("*.h5")]
    information = [sample / f"{accession}_run_info.json", sample / "run_info.json"]
    present_information = [p for p in information if p.exists() or p.is_symlink()]
    if not present_information:
        return None
    payloads.extend(present_information)
    restored = sample / ".metainformant_restored_inputs"
    if restored.exists() or restored.is_symlink():
        payloads.extend([restored / "manifest.json", restored / "config.yaml"])
        if provenance.get("reference_manifest_sha256"):
            payloads.append(restored / "reference_manifest.json")
    elif not isinstance(provenance.get("config_path"), str) or not _regular_owned_file(
        Path(provenance["config_path"]), config_dir
    ):
        return None
    reference = provenance.get("reference_manifest_path")
    if reference and (not isinstance(reference, str) or not _regular_owned_file(Path(reference), data_root)):
        return None
    if not all(_regular_owned_file(path, data_root) for path in payloads):
        return None
    config = config_dir / str(task["config_name"])
    if not _regular_owned_file(config, config_dir) or sha256_file(config) != task["config_sha256"]:
        return None
    audit = classify_quantification(sample, accession, expected_config_path=config)
    if audit["status"] != "current" or audit["contract_id"] != provenance.get("quant_contract_id"):
        return None
    return {
        "srr": accession,
        "batch": int(task["batch_index"]),
        "downloaded": False,
        "quantified": True,
        "skipped": True,
        "error": None,
    }


def execute_manifest_task(
    task: dict[str, Any],
    orchestrator: StreamingPipelineOrchestrator,
    data_root: Path,
    config_dir: Path,
    durable_store: ObjectStore | None,
    durable_cohort: str,
    quant_threads: int,
    *,
    defer_ncbi_fallback: bool = False,
    acquisition_only: bool = False,
    quantification_only: bool = False,
) -> dict[str, Any]:
    species = str(task["species"])
    accession = str(task["accession"])
    config_path = config_dir / str(task["config_name"])
    if not config_path.is_file():
        raise FileNotFoundError(f"config missing for {species}: {config_path}")
    work_dir = data_root / species / "work"
    fastq_dir = work_dir / "getfastq"
    if durable_store is not None:
        index_hash = task.get("reference_index_sha256")
        if (
            not isinstance(index_hash, str)
            or len(index_hash) != 64
            or any(c not in "0123456789abcdef" for c in index_hash)
        ):
            raise ValueError("durable cloud task lacks a frozen index checksum")
        from metainformant.rna.engine.durable_quant import restore_quantification

        try:
            durable_store.get(bound_receipt_key(durable_cohort, species, accession))
        except FileNotFoundError:
            receipt_present = False
        else:
            receipt_present = True
        if receipt_present:
            restored = restore_quantification(
                durable_store,
                durable_cohort,
                species,
                accession,
                work_dir / "quant" / accession,
                expected_config_sha256=sha256_file(config_path),
                expected_reference_index_sha256=task["reference_index_sha256"],
                verified_config_path=config_path,
            )
            if restored["config_sha256"] != sha256_file(config_path):
                raise ValueError("durable receipt uses a different species configuration")
            orchestrator.db.set_state(species, accession, "quantified")
            return {
                "srr": accession,
                "batch": int(task["batch_index"]),
                "downloaded": False,
                "quantified": True,
                "skipped": True,
                "durable": True,
                "error": None,
            }
    result = orchestrator.process_single_sample(
        accession,
        int(task["batch_index"]),
        fastq_dir,
        config_path,
        species,
        quant_threads,
        task.get("expected_paired"),
        defer_ncbi_fallback=defer_ncbi_fallback,
        acquisition_only=acquisition_only,
        quantification_only=quantification_only,
    )
    if result.get("quantified") and durable_store is not None:
        from metainformant.rna.engine.durable_quant import lock_quantification

        lock_quantification(
            durable_store,
            durable_cohort,
            species,
            accession,
            work_dir / "quant" / accession,
            expected_config_sha256=sha256_file(config_path),
            expected_reference_index_sha256=task["reference_index_sha256"],
        )
        result["durable"] = True
    return result
