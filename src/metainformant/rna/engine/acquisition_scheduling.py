"""Conservative admission scenarios, not measured throughput guarantees."""

from __future__ import annotations

import math
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from typing import NotRequired, TypedDict

from metainformant.rna.engine.acquisition_estimates import AcquisitionEstimateError


def positive_size(value: int | float | str | None, field: str) -> int:
    """Parse an exact positive workload count without accepting unknown/nonfinite values."""
    try:
        number = Decimal(str(value)) if value is not None else Decimal(0)
    except InvalidOperation as exc:
        raise AcquisitionEstimateError(field, "requires a known positive integer") from exc
    if isinstance(value, bool) or not number.is_finite() or number <= 0 or number != number.to_integral_value():
        raise AcquisitionEstimateError(field, "requires a known positive integer")
    return int(number)


@dataclass(frozen=True, slots=True)
class PlanningAssumptions:
    transfer_bytes_per_second: float
    extraction_bases_per_second: float
    quant_bases_per_second: float
    source: str
    extraction_slots: int = 1
    quant_slots: int = 4
    setup_seconds: int = 900
    drain_seconds: int = 300
    task_overhead_seconds: int = 60
    safety_factor: float = 1.5

    def __post_init__(self) -> None:
        for name in (
            "transfer_bytes_per_second",
            "quant_bases_per_second",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not math.isfinite(value) or value <= 0:
                raise AcquisitionEstimateError(name, "requires an explicit finite positive rate")
        if (
            isinstance(self.extraction_bases_per_second, bool)
            or not math.isfinite(self.extraction_bases_per_second)
            or self.extraction_bases_per_second < 0
        ):
            raise AcquisitionEstimateError(
                "extraction_bases_per_second", "requires a finite nonnegative rate; zero disables SRA admission"
            )
        for name in (
            "extraction_slots",
            "quant_slots",
            "setup_seconds",
            "drain_seconds",
            "task_overhead_seconds",
        ):
            value = getattr(self, name)
            if type(value) is not int or value <= 0:
                raise AcquisitionEstimateError(name, "requires a positive integer")
        if (
            isinstance(self.safety_factor, bool)
            or not math.isfinite(self.safety_factor)
            or self.safety_factor < 1
            or not self.source.strip()
        ):
            raise AcquisitionEstimateError("planning", "requires a source and safety factor >= 1")


@dataclass(frozen=True, slots=True)
class Workload:
    task_id: str
    raw_bytes: int
    total_bases: int
    transfer_bytes: int | None = None
    requires_extraction: bool = False

    def __post_init__(self) -> None:
        positive_size(self.raw_bytes, "raw_bytes")
        positive_size(self.total_bases, "total_bases")
        if self.transfer_bytes is not None:
            positive_size(self.transfer_bytes, "transfer_bytes")


class TaskWorkloadInput(TypedDict):
    task_id: str
    fastq_bytes: int | float | str | None
    total_bases: int | float | str | None
    source_evidence_sha256: NotRequired[str]
    sra_bytes: NotRequired[int | float | str | None]


def task_workload(task: TaskWorkloadInput) -> Workload:
    """Share the same source-size interpretation between controller and worker."""
    sra = bool(task.get("source_evidence_sha256"))
    return Workload(
        task["task_id"],
        positive_size(task.get("fastq_bytes"), "fastq_bytes"),
        positive_size(task.get("total_bases"), "total_bases"),
        positive_size(task.get("sra_bytes") if sra else task.get("fastq_bytes"), "transfer_bytes"),
        requires_extraction=sra,
    )


def _stage_seconds(durations: list[float], slots: int) -> float:
    # List-scheduling upper bound does not assume the worker uses the planner's
    # task order. Whole stages remain serialized in the admission scenario.
    if not durations:
        return 0.0
    total = sum(durations)
    return min(total, total / slots + max(durations) * (1 - 1 / slots))


def workload_seconds(
    tasks: list[Workload],
    assumptions: PlanningAssumptions,
    *,
    include_setup: bool = True,
) -> int:
    """Sum stages so no unmeasured download/compute overlap gain is assumed.

    Transfer rate is aggregate per worker. Compute rates are per occupied stage
    slot under the declared profile. Extraction is charged for planned SRA
    tasks only; zero extraction rate refuses those tasks. Unplanned fallback
    is outside the scenario. Rates require workload-matched calibration.
    """
    if not tasks:
        raise AcquisitionEstimateError("tasks", "empty work cannot be admitted")
    transfer = (
        sum(t.transfer_bytes if t.transfer_bytes is not None else t.raw_bytes for t in tasks)
        / assumptions.transfer_bytes_per_second
    )
    extraction_tasks = [t for t in tasks if t.requires_extraction]
    if extraction_tasks and assumptions.extraction_bases_per_second == 0:
        raise AcquisitionEstimateError("extraction", "SRA admission requires calibrated extraction rate")
    extraction = _stage_seconds(
        [t.total_bases / assumptions.extraction_bases_per_second for t in extraction_tasks],
        assumptions.extraction_slots,
    )
    quant = _stage_seconds(
        [t.total_bases / assumptions.quant_bases_per_second for t in tasks],
        assumptions.quant_slots,
    )
    work = assumptions.safety_factor * (transfer + extraction + quant + len(tasks) * assumptions.task_overhead_seconds)
    return math.ceil(work) + (assumptions.setup_seconds + assumptions.drain_seconds if include_setup else 0)


def fits_worker_deadline(
    *,
    now: float,
    deadline: float,
    task_seconds: int,
    reserved_seconds: int,
    drain_seconds: int,
) -> bool:
    """Reject late submissions, including time committed to unfinished tasks."""
    if not math.isfinite(now) or not math.isfinite(deadline):
        raise AcquisitionEstimateError("deadline", "requires finite clock values")
    if (
        type(task_seconds) is not int
        or task_seconds <= 0
        or type(reserved_seconds) is not int
        or reserved_seconds < 0
        or type(drain_seconds) is not int
        or drain_seconds <= 0
    ):
        raise AcquisitionEstimateError("task_seconds", "requires positive duration and nonnegative reservations")
    return now + task_seconds + reserved_seconds + drain_seconds <= deadline
