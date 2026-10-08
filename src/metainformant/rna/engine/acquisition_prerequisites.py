"""Resolve native Amalgkit requirements before admitting frozen work."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Literal


@dataclass(frozen=True, slots=True)
class TaskPrerequisite:
    """Admission result without changing a task's scientific semantics."""

    task_id: str
    accession: str
    backend: str
    sequencing_technology: str
    status: Literal["ready", "unresolved"]
    reason: str


class PrerequisiteError(RuntimeError):
    """Selected work cannot run under its existing frozen bindings."""

    def __init__(self, results: Sequence[TaskPrerequisite]) -> None:
        self.results = tuple(results)
        super().__init__(
            "acquisition prerequisites unresolved: "
            + "; ".join(f"{result.task_id}: {result.reason}" for result in results)
        )


def classify_task_prerequisites(
    *, metadata_path: Path, tasks: Sequence[Mapping[str, str | int | bool | None]]
) -> list[TaskPrerequisite]:
    """Use installed native backend rules; reject unsupported frozen bindings.

    The current acquisition envelope binds a Kallisto index checksum. An
    oarfish MMI cannot satisfy that checksum or its provenance contract. Long
    reads therefore require an explicit envelope amendment, even if oarfish is
    installed. This is an admission classification, never a backend override.
    """
    from amalgkit.metadata_utils import load_metadata, parse_bool_flags
    from amalgkit.quant import resolve_oarfish_seq_tech, resolve_quant_backend

    args = SimpleNamespace(metadata=str(metadata_path), quant_backend="auto", oarfish_seq_tech="auto")
    metadata = load_metadata(args)
    sampled = parse_bool_flags(metadata.df["is_sampled"], column_name="is_sampled", default="no")
    batch_runs = metadata.df.loc[sampled, "run"].tolist()
    results: list[TaskPrerequisite] = []
    for task in tasks:
        accession = str(task["accession"])
        task_id = str(task["task_id"])
        backend = ""
        technology = ""
        reason = ""
        batch = int(task["batch_index"])
        if batch < 1 or batch > len(batch_runs) or batch_runs[batch - 1] != accession:
            reason = "frozen batch index does not select the requested accession"
        else:
            try:
                backend = resolve_quant_backend(args, metadata, accession)
                if backend == "oarfish":
                    technology = resolve_oarfish_seq_tech(args, metadata, accession)
                    reason = (
                        "oarfish requires an amended frozen transcript FASTA/MMI binding, "
                        "sequencing-technology provenance and validated native tool bootstrap; "
                        "the existing Kallisto index binding cannot be substituted"
                    )
            except ValueError as exc:
                reason = str(exc)
        results.append(
            TaskPrerequisite(
                task_id,
                accession,
                backend,
                technology,
                "unresolved" if reason else "ready",
                reason,
            )
        )
    return results


def require_worker_prerequisites(results: Sequence[TaskPrerequisite]) -> None:
    """Fail before acquisition and invoke the required native dependency probe."""
    from amalgkit.quant import check_kallisto_dependency

    unresolved = [result for result in results if result.status != "ready"]
    if unresolved:
        raise PrerequisiteError(unresolved)
    if results:
        check_kallisto_dependency()
