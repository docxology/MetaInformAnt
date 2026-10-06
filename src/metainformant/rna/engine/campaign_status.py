"""Read-only, inventory-bounded cloud/local campaign reconciliation."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import json
import re
from typing import Final, Literal

STATES: Final = ("pending", "downloading", "downloaded", "quantifying", "quantified", "failed", "quarantined")
CLOUD_COLUMNS: Final = ("locked", *STATES, "worker_unknown", "unassigned")
LOCAL_COLUMNS: Final = (*STATES, "untracked")
COVERAGE_COLUMNS: Final = ("present_locked", "present_unlocked", "partial", "absent", "transfer_gap", "diagnostic_present")


class StatusError(ValueError):
    """An inventory or observation cannot support a coherent report."""


@dataclass(frozen=True, slots=True)
class Task:
    task_id: str
    accession: str


@dataclass(frozen=True, slots=True)
class Species:
    species: str
    tasks: tuple[Task, ...]


@dataclass(frozen=True, slots=True)
class Inventory:
    species: tuple[Species, ...]
    species_count: int
    task_count: int

    def task_ids(self) -> frozenset[str]:
        ids = [t.task_id for s in self.species for t in s.tasks]
        if not ids or len(ids) != self.task_count or len(set(ids)) != len(ids):
            raise StatusError("Empty, duplicate, or inconsistent inventory tasks")
        if len(self.species) != self.species_count or len({s.species for s in self.species}) != self.species_count:
            raise StatusError("Inconsistent species inventory")
        if any(t.task_id != f"{s.species}/{t.accession}" for s in self.species for t in s.tasks):
            raise StatusError("Task identity disagrees with species/accession")
        if any(not re.fullmatch(r"[a-z0-9_]+", s.species) for s in self.species) or any(
            not re.fullmatch(r"[A-Za-z0-9_]+", t.accession) for s in self.species for t in s.tasks
        ):
            raise StatusError("Unsafe species or accession identifier")
        return frozenset(ids)


def load_inventory(data: bytes) -> Inventory:
    """Parse the frozen inventory at the protocol boundary without new dependencies."""
    value = json.loads(data)
    if not isinstance(value, dict) or not isinstance(value.get("species"), list):
        raise StatusError("Inventory must contain a species list")
    species = []
    for row in value["species"]:
        if not isinstance(row, dict) or not isinstance(row.get("species"), str) or not isinstance(row.get("tasks"), list):
            raise StatusError("Malformed inventory species")
        tasks = []
        for task in row["tasks"]:
            if not isinstance(task, dict) or not isinstance(task.get("task_id"), str) or not isinstance(task.get("accession"), str):
                raise StatusError("Malformed inventory task")
            tasks.append(Task(task["task_id"], task["accession"]))
        species.append(Species(row["species"], tuple(tasks)))
    if type(value.get("species_count")) is not int or type(value.get("task_count")) is not int:
        raise StatusError("Inventory totals must be integers")
    result = Inventory(tuple(species), value["species_count"], value["task_count"])
    result.task_ids()
    return result


@dataclass(frozen=True, slots=True)
class Observation:
    task_id: str
    state: str
    source: str


@dataclass(frozen=True, slots=True)
class SampleStatus:
    task_id: str
    species: str
    cloud: str
    local: str
    coverage: Literal["present_locked", "present_unlocked", "partial", "absent"]
    transfer_gap: bool
    diagnostic_present: bool


@dataclass(frozen=True, slots=True)
class SpeciesStatus:
    species: str
    eligible: int
    cloud: dict[str, int]
    local: dict[str, int]
    coverage: dict[str, int]


@dataclass(frozen=True, slots=True)
class StatusReport:
    rows: tuple[SpeciesStatus, ...]
    totals: SpeciesStatus
    samples: tuple[SampleStatus, ...]
    local_outside_inventory: int


def reconcile(
    inventory: Inventory,
    locked: frozenset[str],
    assigned: frozenset[str],
    worker: tuple[Observation, ...],
    local: tuple[Observation, ...],
    present: frozenset[str],
    partial: frozenset[str],
    diagnostic: frozenset[str],
) -> StatusReport:
    """Partition every eligible task exactly once; retain file/DB distinctions."""
    ids = inventory.task_ids()
    if not locked <= ids or not assigned <= ids:
        raise StatusError("Cloud observations contain tasks outside frozen inventory")
    if present & partial:
        raise StatusError("Complete and partial file observations overlap")
    if not (present | partial | diagnostic) <= ids:
        raise StatusError("File observations contain tasks outside frozen inventory")
    workers: dict[str, str] = {}
    locals_: dict[str, str] = {}
    for records, target in ((worker, workers), (local, locals_)):
        for record in records:
            if record.state not in STATES:
                raise StatusError(f"Unknown sample state: {record.state}")
            if record.task_id in target:
                raise StatusError(f"Duplicate sample observation: {record.task_id}")
            target[record.task_id] = record.state
    if not workers.keys() <= assigned:
        raise StatusError("Worker telemetry contains unassigned tasks")
    samples = []
    rows = []
    for species in inventory.species:
        cloud_counts: Counter[str] = Counter()
        local_counts: Counter[str] = Counter()
        coverage_counts: Counter[str] = Counter()
        for task in species.tasks:
            key = task.task_id
            cloud_state = "locked" if key in locked else workers.get(key, "worker_unknown" if key in assigned else "unassigned")
            local_state = locals_.get(key, "untracked")
            coverage: Literal["present_locked", "present_unlocked", "partial", "absent"] = "absent"
            if key in present:
                coverage = "present_locked" if key in locked else "present_unlocked"
            elif key in partial:
                coverage = "partial"
            gap = key in locked and key not in present
            cloud_counts[cloud_state] += 1
            local_counts[local_state] += 1
            coverage_counts[coverage] += 1
            coverage_counts["transfer_gap"] += int(gap)
            coverage_counts["diagnostic_present"] += int(key in diagnostic)
            samples.append(SampleStatus(key, species.species, cloud_state, local_state, coverage, gap, key in diagnostic))
        rows.append(SpeciesStatus(species.species, len(species.tasks), dict(cloud_counts), dict(local_counts), dict(coverage_counts)))
    total = SpeciesStatus("TOTAL", inventory.task_count,
                          {k: sum(r.cloud.get(k, 0) for r in rows) for k in CLOUD_COLUMNS},
                          {k: sum(r.local.get(k, 0) for r in rows) for k in LOCAL_COLUMNS},
                          {k: sum(r.coverage.get(k, 0) for r in rows) for k in COVERAGE_COLUMNS})
    return StatusReport(tuple(rows), total, tuple(samples), len(locals_.keys() - ids))


def markdown_tables(report: StatusReport, *, include_coverage: bool = True) -> str:
    """Render two tables with row and column marginals, without overlapping stages."""
    cloud = ["Species", *CLOUD_COLUMNS, "Total"]
    local = ["Species", *LOCAL_COLUMNS, *(COVERAGE_COLUMNS if include_coverage else ()), "Total"]
    tables = []
    for labels, lane in ((cloud, "cloud"), (local, "local")):
        lines = ["| " + " | ".join(labels) + " |", "| " + " | ".join("---" for _ in labels) + " |"]
        for row in (*report.rows, report.totals):
            counts = row.cloud if lane == "cloud" else row.local
            columns = CLOUD_COLUMNS if lane == "cloud" else LOCAL_COLUMNS
            values = [row.species, *(str(counts.get(k, 0)) for k in columns)]
            if lane == "local" and include_coverage:
                values.extend(str(row.coverage.get(k, 0)) for k in COVERAGE_COLUMNS)
            lines.append("| " + " | ".join((*values, str(row.eligible))) + " |")
        tables.append("\n".join(lines))
    return "\n\n".join(tables) + "\n"
