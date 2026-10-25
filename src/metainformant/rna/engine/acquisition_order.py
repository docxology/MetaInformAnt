from __future__ import annotations

from bisect import bisect_left
from collections import defaultdict, deque
from collections.abc import Iterable, Mapping
from typing import Any, Final

SIZE_BAND_UPPER_BYTES: Final = (
    100 * 1024**2,
    500 * 1024**2,
    1024**3,
    2 * 1024**3,
    5 * 1024**3,
    20 * 1024**3,
)


def plan_acquisition_order(
    pending: Iterable[Mapping[str, Any]],
    checkpoint_ids: Iterable[str],
    *,
    global_size_bands: bool = False,
) -> list[str]:
    """Interleave restore verification with size-ordered, species-fair acquisition.

    Optional global size bands precede species rotation. Every task appears
    once; checkpoint verification remains interleaved. Unknown sizes occupy
    the last band. This changes admission order, not task eligibility.
    """
    by_species: dict[tuple[int, str], deque[str]] = defaultdict(deque)
    seen: set[str] = set()
    checkpoints = deque(sorted(checkpoint_ids))
    for task_id in checkpoints:
        if task_id in seen:
            raise ValueError(f"Duplicate acquisition task: {task_id}")
        seen.add(task_id)
    for task in sorted(
        pending,
        key=lambda entry: (
            int(entry.get("fastq_bytes", 0)) <= 0,
            int(entry.get("fastq_bytes", 0)),
            entry["task_id"],
        ),
    ):
        task_id = str(task["task_id"])
        if task_id in seen:
            raise ValueError(f"Duplicate acquisition task: {task_id}")
        seen.add(task_id)
        size = int(task.get("fastq_bytes", 0))
        band = (
            (bisect_left(SIZE_BAND_UPPER_BYTES, size) if size > 0 else len(SIZE_BAND_UPPER_BYTES) + 1)
            if global_size_bands
            else 0
        )
        by_species[(band, str(task["species"]))].append(task_id)
    ordered = []
    while by_species:
        first_band = min(by_species)[0]
        for species in tuple(key for key in sorted(by_species) if key[0] == first_band):
            if checkpoints:
                ordered.append(checkpoints.popleft())
            ordered.append(by_species[species].popleft())
            if not by_species[species]:
                del by_species[species]
    ordered.extend(checkpoints)
    return ordered
