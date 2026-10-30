"""Conservative native process-state confirmation for recovery observers."""

from __future__ import annotations

import psutil


class ProcessIdentityChangedError(RuntimeError):
    """A numeric PID no longer has its captured native creation time."""


def process_is_confirmed_inert(pid: int, expected_create_time: float | None) -> bool:
    """Confirm current kernel absence or a stable same-identity zombie.

    A live process, an inaccessible process, or an unconfirmed zombie returns
    False. AccessDenied itself is never absence evidence. Call this after an
    observation race using the original captured identity; do not substitute a
    newly observed creation time to accept a reused PID.

    Args:
        pid: Previously observed numeric process identifier.
        expected_create_time: Captured psutil creation time; None permits only
            fresh kernel-confirmed absence, never zombie acceptance.

    Returns:
        True only for a fresh NoSuchProcess or a twice-confirmed zombie with the
        original creation time.

    Raises:
        ValueError: The captured PID is not a positive integer.
        ProcessIdentityChangedError: The PID belongs to a different identity.
    """
    if isinstance(pid, bool) or not isinstance(pid, int) or pid <= 0:
        raise ValueError("A captured positive integer PID is required")
    try:
        current = psutil.Process(pid)
        fresh_created = current.create_time()
        state = current.status()
    except psutil.ZombieProcess:
        return False
    except psutil.NoSuchProcess:
        return True
    except psutil.AccessDenied:
        return False
    if expected_create_time is None:
        return False
    if fresh_created != expected_create_time:
        raise ProcessIdentityChangedError(f"Process identity changed after denied inspection: pid={pid}")
    if state != psutil.STATUS_ZOMBIE:
        return False
    try:
        second = psutil.Process(pid)
        if second.create_time() != expected_create_time:
            raise ProcessIdentityChangedError(f"Zombie identity changed after denied inspection: pid={pid}")
        return second.status() == psutil.STATUS_ZOMBIE
    except psutil.ZombieProcess:
        return False
    except psutil.NoSuchProcess:
        return True
    except psutil.AccessDenied:
        return False
