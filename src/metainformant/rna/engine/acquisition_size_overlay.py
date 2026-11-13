"""Effective compressed FASTQ weights for acquisition ordering, with optional per-task size supplements."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any


def effective_fastq_bytes(task: Mapping[str, Any], sizes: Mapping[str, int] | None = None) -> int:
    """The task's ``fastq_bytes``, or its supplement from ``sizes`` (keyed by ``task_id``).

    A supplement may only fill a zero weight, and only with a positive integer
    below 2**63; anything else raises :class:`ValueError`.
    """

    original = int(task.get("fastq_bytes", 0))
    identity = str(task["task_id"])
    if sizes is None or identity not in sizes:
        return original
    size = sizes[identity]
    if original != 0 or type(size) is not int or not 0 < size < 2**63:
        raise ValueError("Size supplements may only fill a zero weight with positive integer bytes")
    return size
