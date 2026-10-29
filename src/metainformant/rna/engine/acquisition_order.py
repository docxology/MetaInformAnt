from __future__ import annotations

from bisect import bisect_left
from collections import defaultdict, deque
from collections.abc import Iterable, Mapping
from typing import Any, Final

from metainformant.rna.engine.acquisition_size_overlay import effective_fastq_bytes

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
    global_size_order: bool = False,
    size_overrides: Mapping[str, int] | None = None,
) -> list[str]:
    """Interleave restore verification with size-ordered, species-fair acquisition.

    Optional global size bands precede species rotation. Global size order
    instead removes species rotation and sorts all acquisition weights.
    Every task appears once; checkpoint verification remains interleaved.
    Unknown sizes follow positive weights. This changes admission priority,
    not task eligibility or completion order.
    """
    if global_size_bands and global_size_order:
        raise ValueError("Choose global size bands or global size order, not both")
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
            effective_fastq_bytes(entry, size_overrides) <= 0,
            effective_fastq_bytes(entry, size_overrides),
            entry["task_id"],
        ),
    ):
        task_id = str(task["task_id"])
        if task_id in seen:
            raise ValueError(f"Duplicate acquisition task: {task_id}")
        seen.add(task_id)
        size = effective_fastq_bytes(task, size_overrides)
        band = (
            (bisect_left(SIZE_BAND_UPPER_BYTES, size) if size > 0 else len(SIZE_BAND_UPPER_BYTES) + 1)
            if global_size_bands
            else 0
        )
        species = str(task["species"])
        by_species[(band, "" if global_size_order else species)].append(task_id)
    ordered = []
    while by_species:
        first_band = min(by_species)[0]
        for group_key in tuple(key for key in sorted(by_species) if key[0] == first_band):
            if checkpoints:
                ordered.append(checkpoints.popleft())
            ordered.append(by_species[group_key].popleft())
            if not by_species[group_key]:
                del by_species[group_key]
    ordered.extend(checkpoints)
    return ordered
