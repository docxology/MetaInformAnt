"""Registry of the transfer and quantification process groups a writer has started.

Every download, extraction and quantification command runs in its own session
(``start_new_session=True``, see
:func:`metainformant.rna.retrieval.ena_downloader._run_command_in_process_group`),
so a signal sent to the producer's process group never reaches it. This module
records those groups so the producer can stop them before it releases the
campaign writer lock:

* :func:`spawn` starts a command in a new session and registers its group,
  atomically with respect to :func:`begin_shutdown`: once shutdown has begun,
  no new group starts (:class:`ProcessGroupsStopping` is raised instead);
* :func:`signal_all` sends a signal to every registered group;
* :func:`wait_until_empty` waits for the worker threads to reap them;
* :func:`set_inherited_fds` names descriptors (the campaign lock) that every
  new group inherits, so the kernel keeps the writer lock held for as long as
  any of them lives, even after the producer itself was SIGKILLed.

The registry is process-wide; :func:`reset` restores the initial state (a new
run in the same process).
"""

from __future__ import annotations

import os
import signal
import subprocess
import threading
import time
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

_GUARD = threading.Lock()
_CHANGED = threading.Condition(_GUARD)
_GROUPS: set[int] = set()
_INHERITED_FDS: tuple[int, ...] = ()
_STOPPING: str | None = None


class ProcessGroupsStopping(RuntimeError):
    """The writer is stopping; no new transfer or quantification process may start."""


def set_inherited_fds(descriptors: Iterable[int | None]) -> None:
    """Descriptors every subsequently spawned group inherits (``None`` entries are ignored)."""

    global _INHERITED_FDS
    with _GUARD:
        _INHERITED_FDS = tuple(sorted({int(fd) for fd in descriptors if fd is not None}))


def inherited_fds() -> tuple[int, ...]:
    """The descriptors new groups inherit."""

    with _GUARD:
        return _INHERITED_FDS


def spawn(
    command: Sequence[str],
    *,
    env: Mapping[str, str] | None = None,
    stdout: Any = subprocess.PIPE,
    stderr: Any = subprocess.PIPE,
    text: bool = True,
) -> subprocess.Popen[Any]:
    """Start ``command`` as the leader of a new session and register its process group.

    Raises :class:`ProcessGroupsStopping` once :func:`begin_shutdown` has run.
    Call :func:`unregister` after the process has been reaped.
    """

    with _GUARD:
        if _STOPPING is not None:
            raise ProcessGroupsStopping(
                f"writer is stopping ({_STOPPING}); not starting {command[0] if command else ''}"
            )
        process: subprocess.Popen[Any] = subprocess.Popen(
            list(command),
            stdout=stdout,
            stderr=stderr,
            text=text,
            start_new_session=True,
            env=dict(env) if env is not None else None,
            pass_fds=_INHERITED_FDS,
        )
        _GROUPS.add(process.pid)
        _CHANGED.notify_all()
    return process


def unregister(pgid: int) -> None:
    """Forget a group whose leader has been reaped."""

    with _GUARD:
        _GROUPS.discard(pgid)
        _CHANGED.notify_all()


def registered() -> tuple[int, ...]:
    """Process-group ids of the commands started and not yet reaped."""

    with _GUARD:
        return tuple(sorted(_GROUPS))


def begin_shutdown(reason: str) -> bool:
    """Refuse every later :func:`spawn`; return True the first time."""

    global _STOPPING
    with _GUARD:
        if _STOPPING is not None:
            return False
        _STOPPING = reason or "stopping"
        return True


def stopping() -> str | None:
    """The shutdown reason, or ``None`` while new groups may start."""

    with _GUARD:
        return _STOPPING


def signal_all(signum: int = signal.SIGTERM) -> list[int]:
    """Send ``signum`` to every registered group; return the groups that were signalled."""

    signalled: list[int] = []
    for pgid in registered():
        try:
            os.killpg(pgid, signum)
        except ProcessLookupError:
            continue
        except PermissionError:
            continue
        signalled.append(pgid)
    return signalled


def wait_until_empty(timeout: float) -> bool:
    """Wait up to ``timeout`` seconds for every registered group to be reaped; True when none remain."""

    deadline = time.monotonic() + max(0.0, timeout)
    with _GUARD:
        while _GROUPS:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                return False
            _CHANGED.wait(min(remaining, 0.5))
        return True


def reset() -> None:
    """Allow spawning again and forget inherited descriptors (a new writer run in this process)."""

    global _STOPPING, _INHERITED_FDS
    with _GUARD:
        _STOPPING = None
        _INHERITED_FDS = ()
