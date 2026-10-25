from __future__ import annotations

import math
import os
import threading
import time
from collections import deque
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Literal
from urllib.parse import urlsplit


@dataclass
class _Gate:
    limit: int
    active: int = 0
    queue: deque[object] = field(default_factory=deque)


_condition = threading.Condition()
_gates: dict[tuple[str, str, int], _Gate] = {}


def _connection_limit(kind: Literal["file", "metadata"]) -> int:
    variable = "AMALGKIT_ENA_MAX_CONNECTIONS" if kind == "file" else "AMALGKIT_ENA_API_CONNECTIONS"
    limit = int(os.environ.get(variable, "0"))
    if not 0 <= limit <= 256:
        raise ValueError(f"{variable} must be between 0 and 256")
    return limit


def validate_connection_budgets() -> None:
    _connection_limit("file")
    _connection_limit("metadata")


@contextmanager
def ena_request_slot(
    url: str, *, kind: Literal["file", "metadata"], deadline_at: float
) -> Iterator[None]:
    """Share a FIFO request budget across downloader instances in this process.

    Waiting consumes the caller's absolute deadline. A lease covers one native
    request, never its retry backoff; separate origins and metadata have their
    own budgets. Zero keeps the legacy unlimited setting.
    """
    limit = _connection_limit(kind)
    if not math.isfinite(deadline_at):
        raise ValueError("ENA request deadline must be finite")
    if deadline_at <= time.monotonic():
        raise TimeoutError(f"ENA {kind} request deadline exhausted")
    if not limit:
        yield
        return
    parsed = urlsplit(url)
    if parsed.hostname is None or parsed.scheme not in {"http", "https"}:
        raise ValueError("ENA request budget requires an HTTP origin")
    key = kind, parsed.hostname, parsed.port or (443 if parsed.scheme == "https" else 80)
    ticket = object()
    acquired = False
    with _condition:
        gate = _gates.setdefault(key, _Gate(limit))
        if gate.limit != limit:
            raise ValueError("Cannot change an ENA connection budget while requests own it")
        gate.queue.append(ticket)
    try:
        with _condition:
            while gate.active >= gate.limit or gate.queue[0] is not ticket:
                remaining = deadline_at - time.monotonic()
                if remaining <= 0:
                    raise TimeoutError(f"ENA {kind} connection slot deadline exhausted; inputs retained")
                _condition.wait(timeout=remaining)
            if deadline_at <= time.monotonic():
                raise TimeoutError(f"ENA {kind} connection slot deadline exhausted; inputs retained")
            gate.queue.popleft()
            gate.active += 1
            acquired = True
            _condition.notify_all()
        yield
    finally:
        with _condition:
            if acquired:
                gate.active -= 1
            else:
                gate.queue.remove(ticket)
            if not gate.active and not gate.queue:
                del _gates[key]
            _condition.notify_all()
