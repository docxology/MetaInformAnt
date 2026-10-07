"""Execute one frozen acquisition task using shared local and durable methods."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from metainformant.rna.engine.acquisition_manifest import sha256_file
from metainformant.rna.engine.durable_quant import ObjectStore, bound_receipt_key
from metainformant.rna.engine.streaming_orchestrator import StreamingPipelineOrchestrator


def execute_manifest_task(
    task: dict[str, Any],
    orchestrator: StreamingPipelineOrchestrator,
    data_root: Path,
    config_dir: Path,
    durable_store: ObjectStore | None,
    durable_cohort: str,
    quant_threads: int,
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
