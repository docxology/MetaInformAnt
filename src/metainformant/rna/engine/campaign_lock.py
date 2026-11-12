"""One exclusive writer per campaign data root, released by the kernel on any death.

A campaign data root (``AMALGKIT_DATA_ROOT``) has exactly one writer at a time:
the streaming producer, a finalization run, or a maintenance command (cohort
migration, quantification import, seal). Two writers on one root corrupt it in
ways nothing detects later: two producers rewrite the same task table while the
other's ``amalgkit quant --batch N`` reads it, and two resets fight over
in-flight rows.

The lock is an ``fcntl.flock(LOCK_EX | LOCK_NB)`` on ``<root>/.campaign.flock``
held for the life of the writer. The kernel drops it when the holder exits for
any reason (SIGKILL, crash, power loss, a USB drop that kills the process), so
there is never a stale lock to clean up by hand. The file's bytes are owner
information for humans (``pid``, ``host``, ``kind``, ``started_at``,
``token``); nothing decides on them except token inheritance (below).

**Children of the holder.** ``python -m metainformant.rna.engine.campaign_lock
run --kind K -- cmd ...`` holds the lock while ``cmd`` runs, passes the locked
descriptor to the child (so the lock survives even if the wrapper itself is
killed while the child still writes) and exports the owner's random token as
``METAINFORMANT_CAMPAIGN_LOCK_TOKEN``. A process whose token equals the
recorded token of the live holder is *inside* the lock: :func:`acquire_writer`
returns an inherited handle and :func:`is_quiescent` reports quiescent.

**In-process nesting.** A second :func:`acquire_writer` of the same ``kind``
in the process that already holds the lock nests (reference counted); a
different ``kind`` is refused, exactly as another process would be.

**Other writers recognized.** Before claiming the root, the lock also refuses
while the legacy shell launcher's directory locks
(``results/.full_campaign.lock``, ``results/.finalization.lock``) name a live
process that is not an ancestor of this one, and while another process holds
the manifest worker's ``<root>/.acquisition.lock`` flock.

Exit codes used by the producer entry points (BSD ``sysexits``):
``0`` drained, ``75`` resumable work remains, ``73`` lock busy, ``69``
environment or preflight failure, ``70`` internal error.
"""

from __future__ import annotations

import argparse
import contextlib
import errno
import fcntl
import json
import os
import secrets
import signal
import socket
import subprocess
import sys
import threading
import time
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Self

LOCK_FILENAME = ".campaign.flock"
LOCK_SCHEMA = "metainformant.rna.campaign_lock.v1"
TOKEN_ENV = "METAINFORMANT_CAMPAIGN_LOCK_TOKEN"
# Directory locks created by the shell launchers (mkdir + pid file).
LEGACY_DIRECTORY_LOCKS = ("results/.full_campaign.lock", "results/.finalization.lock")
# The manifest acquisition worker's flock (acquisition_worker.run_manifest).
FOREIGN_FLOCKS = (".acquisition.lock",)
# A legacy directory lock without a pid file younger than this is a launcher between mkdir and its pid write.
LEGACY_PIDLESS_GRACE_SECONDS = 60.0

EXIT_OK = 0
EXIT_ENVIRONMENT = 69
EXIT_INTERNAL = 70
EXIT_LOCK_BUSY = 73
EXIT_RETRYABLE = 75


class CampaignLockBusy(RuntimeError):
    """Another writer holds (or may hold) the campaign data root."""

    def __init__(self, message: str, *, holder: Mapping[str, Any] | None = None, path: Path | None = None) -> None:
        super().__init__(message)
        self.holder = dict(holder) if holder else None
        self.path = path


@dataclass
class _HeldLock:
    fd: int
    kind: str
    token: str
    owner: dict[str, Any]
    depth: int = 1


_HELD: dict[str, _HeldLock] = {}
_HELD_GUARD = threading.Lock()


@dataclass
class CampaignLock:
    """A held (or inherited) campaign writer lock. Release it, or use it as a context manager."""

    data_root: Path
    path: Path
    kind: str
    token: str
    owner: dict[str, Any]
    inherited: bool
    _key: str | None = field(default=None, repr=False)
    _released: bool = field(default=False, repr=False)

    @property
    def fd(self) -> int | None:
        """The locked descriptor in this process (``None`` for a lock inherited through the token)."""

        if self._key is None:
            return None
        held = _HELD.get(self._key)
        return held.fd if held is not None else None

    def child_environment(self, base: Mapping[str, str] | None = None) -> dict[str, str]:
        """An environment whose processes are inside this lock (carries the holder token)."""

        environment = dict(os.environ if base is None else base)
        environment[TOKEN_ENV] = self.token
        return environment

    def release(self) -> None:
        """Release this handle; the flock is dropped when the outermost handle in the process is released."""

        if self._released:
            return
        self._released = True
        if self._key is None:
            return
        with _HELD_GUARD:
            held = _HELD.get(self._key)
            if held is None:
                return
            held.depth -= 1
            if held.depth > 0:
                return
            del _HELD[self._key]
        try:
            os.close(held.fd)
        except OSError:
            pass

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *_exc: object) -> None:
        self.release()


def lock_path(data_root: str | Path) -> Path:
    """The flock file of a data root."""

    return Path(data_root) / LOCK_FILENAME


def _now() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _key_for(path: Path) -> str:
    return os.path.realpath(path)


def _read_owner_fd(fd: int) -> dict[str, Any] | None:
    try:
        data = os.pread(fd, 65536, 0)
    except OSError:
        return None
    try:
        payload = json.loads(data.decode("utf-8")) if data.strip() else None
    except (UnicodeDecodeError, ValueError):
        return None
    return payload if isinstance(payload, dict) else None


def read_owner(data_root: str | Path) -> dict[str, Any] | None:
    """The owner record last written into the lock file (informational; ``None`` when absent)."""

    try:
        fd = os.open(lock_path(data_root), os.O_RDONLY | os.O_CLOEXEC)
    except OSError:
        return None
    try:
        return _read_owner_fd(fd)
    finally:
        os.close(fd)


def _pid_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _ancestor_pids() -> set[int]:
    """Process ids of this process's ancestors (the shell launcher is one when it started us)."""

    try:
        import psutil

        return {parent.pid for parent in psutil.Process().parents()}
    except Exception:  # noqa: BLE001 - psutil may be unavailable or the table unreadable; fall back to the parent
        return {os.getppid()}


def _legacy_directory_lock_holders(root: Path) -> list[dict[str, Any]]:
    """Live, non-ancestor owners of the shell launchers' directory locks."""

    holders: list[dict[str, Any]] = []
    ancestors: set[int] | None = None
    for relative in LEGACY_DIRECTORY_LOCKS:
        directory = root / relative
        try:
            info = directory.lstat()
        except OSError:
            continue
        pid_text = None
        try:
            pid_text = (directory / "pid").read_text().strip()
        except OSError:
            pass
        if not pid_text:
            if time.time() - info.st_mtime < LEGACY_PIDLESS_GRACE_SECONDS:
                holders.append({"lock": relative, "pid": None, "reason": "launcher lock without a pid file yet"})
            continue
        try:
            pid = int(pid_text)
        except ValueError:
            continue
        if not _pid_alive(pid):
            continue
        if ancestors is None:
            ancestors = _ancestor_pids()
        if pid in ancestors:
            continue
        holders.append({"lock": relative, "pid": pid, "reason": "live launcher process"})
    return holders


def _foreign_flock_holders(root: Path) -> list[dict[str, Any]]:
    """Other lanes' flocks on this root that some process holds right now."""

    holders: list[dict[str, Any]] = []
    for relative in FOREIGN_FLOCKS:
        try:
            fd = os.open(root / relative, os.O_RDONLY | os.O_CLOEXEC)
        except OSError:
            continue
        try:
            try:
                fcntl.flock(fd, fcntl.LOCK_SH | fcntl.LOCK_NB)
            except BlockingIOError:
                holders.append({"lock": relative, "pid": None, "reason": "flock held by another process"})
            else:
                fcntl.flock(fd, fcntl.LOCK_UN)
        finally:
            os.close(fd)
    return holders


def _describe(holder: Mapping[str, Any] | None) -> str:
    if not holder:
        return "unknown holder"
    parts = [f"{key}={holder[key]}" for key in ("kind", "pid", "host", "started_at") if holder.get(key) is not None]
    return ", ".join(parts) or "unknown holder"


def _refuse_other_lanes(root: Path) -> None:
    """Raise :class:`CampaignLockBusy` while a legacy launcher lock or another lane's flock is live."""

    others = _legacy_directory_lock_holders(root) + _foreign_flock_holders(root)
    if others:
        first = others[0]
        raise CampaignLockBusy(
            f"another campaign lane holds {root / first['lock']} ({first['reason']}"
            + (f", pid {first['pid']}" if first.get("pid") else "")
            + ")",
            holder=first,
            path=root / first["lock"],
        )


def acquire(
    data_root: str | Path,
    kind: str,
    *,
    environ: Mapping[str, str] | None = None,
    argv: Sequence[str] | None = None,
) -> CampaignLock:
    """Claim the exclusive writer lock of ``data_root`` for ``kind`` (``producer``, ``import``, ...).

    Never creates the data root: a missing root (an unmounted volume) raises
    :class:`FileNotFoundError`. Raises :class:`CampaignLockBusy` while another
    writer holds the root. Release the returned lock (or use
    :func:`acquire_writer`).
    """

    if not kind or not kind.replace("-", "").replace("_", "").isalnum():
        raise ValueError(f"lock kind must be a short identifier: {kind!r}")
    root = Path(data_root)
    if not root.is_dir():
        raise FileNotFoundError(
            errno.ENOENT, "campaign data root is not a directory (is the volume mounted?)", str(root)
        )
    path = lock_path(root)
    key = _key_for(path)
    environment = os.environ if environ is None else environ

    with _HELD_GUARD:
        held = _HELD.get(key)
        if held is not None:
            if held.kind != kind:
                raise CampaignLockBusy(
                    f"this process already holds {path} as {held.kind!r}; refusing a nested {kind!r} writer",
                    holder=held.owner,
                    path=path,
                )
            held.depth += 1
            return CampaignLock(root, path, kind, held.token, dict(held.owner), inherited=True, _key=key)

    token = environment.get(TOKEN_ENV, "").strip()
    if not token:
        # Refuse before creating anything; re-checked below once the flock is held (a launcher may start between).
        _refuse_other_lanes(root)
    fd = os.open(path, os.O_RDWR | os.O_CREAT | os.O_CLOEXEC, 0o644)
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            owner = _read_owner_fd(fd)
            if token and owner is not None and secrets.compare_digest(str(owner.get("token", "")), token):
                os.close(fd)
                return CampaignLock(root, path, kind, token, owner, inherited=True)
            raise CampaignLockBusy(
                f"campaign data root is locked by another writer ({_describe(owner)}): {path}", holder=owner, path=path
            ) from None
        _refuse_other_lanes(root)
        token = secrets.token_hex(16)
        owner = {
            "schema": LOCK_SCHEMA,
            "kind": kind,
            "pid": os.getpid(),
            "host": socket.gethostname(),
            "started_at": _now(),
            "token": token,
            "argv": list(argv if argv is not None else sys.argv),
        }
        payload = (json.dumps(owner, sort_keys=True) + "\n").encode("utf-8")
        os.ftruncate(fd, 0)
        os.pwrite(fd, payload, 0)
        try:
            os.fsync(fd)
        except OSError:
            pass
    except BaseException:
        os.close(fd)
        raise
    with _HELD_GUARD:
        _HELD[key] = _HeldLock(fd=fd, kind=kind, token=token, owner=owner)
    return CampaignLock(root, path, kind, token, dict(owner), inherited=False, _key=key)


@contextlib.contextmanager
def acquire_writer(
    data_root: str | Path,
    kind: str,
    *,
    environ: Mapping[str, str] | None = None,
) -> Iterator[CampaignLock]:
    """Hold the exclusive writer lock of ``data_root`` for the ``with`` block (see :func:`acquire`)."""

    lock = acquire(data_root, kind, environ=environ)
    try:
        yield lock
    finally:
        lock.release()


def writer_status(data_root: str | Path, *, environ: Mapping[str, str] | None = None) -> dict[str, Any]:
    """Who, if anyone, is writing ``data_root`` right now. Never creates or changes a file.

    ``quiescent`` is true when no other writer holds the root. This process's
    own lock and a lock whose holder token is in ``environ`` (we are its child)
    count as quiescent (``inside_lock`` true).
    """

    root = Path(data_root)
    path = lock_path(root)
    environment = os.environ if environ is None else environ
    status: dict[str, Any] = {
        "data_root": str(root),
        "lock_path": str(path),
        "held": False,
        "inside_lock": False,
        "owner": read_owner(root),
        "other_writers": [],
    }
    key = _key_for(path)
    with _HELD_GUARD:
        own = _HELD.get(key)
    if own is not None:
        status.update(held=True, inside_lock=True, owner=dict(own.owner))
    else:
        try:
            fd = os.open(path, os.O_RDONLY | os.O_CLOEXEC)
        except OSError:
            fd = None
        if fd is not None:
            try:
                try:
                    fcntl.flock(fd, fcntl.LOCK_SH | fcntl.LOCK_NB)
                except BlockingIOError:
                    status["held"] = True
                else:
                    fcntl.flock(fd, fcntl.LOCK_UN)
            finally:
                os.close(fd)
        if status["held"]:
            token = environment.get(TOKEN_ENV, "").strip()
            owner = status["owner"] or {}
            status["inside_lock"] = bool(token) and secrets.compare_digest(str(owner.get("token", "")), token)
    if root.is_dir():
        status["other_writers"] = _legacy_directory_lock_holders(root) + _foreign_flock_holders(root)
    status["quiescent"] = (not status["held"] or status["inside_lock"]) and not status["other_writers"]
    return status


def is_quiescent(data_root: str | Path, *, environ: Mapping[str, str] | None = None) -> bool:
    """True when no other process writes ``data_root`` (see :func:`writer_status`)."""

    return bool(writer_status(data_root, environ=environ)["quiescent"])


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------


def _default_data_root() -> str | None:
    value = os.environ.get("AMALGKIT_DATA_ROOT", "").strip()
    return value or None


def _run_child(lock: CampaignLock, command: Sequence[str]) -> int:
    """Run ``command`` inside ``lock``: the child inherits the locked descriptor and the holder token."""

    descriptor = lock.fd
    pass_fds = (descriptor,) if descriptor is not None else ()
    if descriptor is not None:
        os.set_inheritable(descriptor, True)
    try:
        child = subprocess.Popen(list(command), env=lock.child_environment(), pass_fds=pass_fds)
    finally:
        if descriptor is not None:
            os.set_inheritable(descriptor, False)

    def forward(signum: int, _frame: object) -> None:
        try:
            child.send_signal(signum)
        except ProcessLookupError:
            pass

    previous = {signum: signal.signal(signum, forward) for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP)}
    try:
        returncode = child.wait()
    finally:
        for signum, handler in previous.items():
            signal.signal(signum, handler)
    return 128 - returncode if returncode < 0 else returncode


def main(argv: Sequence[str] | None = None) -> int:
    """``python -m metainformant.rna.engine.campaign_lock {status,run}``."""

    parser = argparse.ArgumentParser(
        prog="python -m metainformant.rna.engine.campaign_lock", description=__doc__.split("\n")[0]
    )
    commands = parser.add_subparsers(dest="command", required=True)
    status_parser = commands.add_parser("status", help="Report the writer of a data root (exit 0 quiescent, 73 busy)")
    status_parser.add_argument(
        "--data-root", default=_default_data_root(), help="Data root (default: AMALGKIT_DATA_ROOT)"
    )
    run_parser = commands.add_parser("run", help="Hold the lock while a command runs (exit 73 when busy)")
    run_parser.add_argument("--data-root", default=_default_data_root(), help="Data root (default: AMALGKIT_DATA_ROOT)")
    run_parser.add_argument("--kind", required=True, help="Writer kind recorded in the lock, e.g. producer")
    run_parser.add_argument("child", nargs=argparse.REMAINDER, help="-- command [args...]")
    args = parser.parse_args(argv)
    if not args.data_root:
        parser.error("--data-root is required (or set AMALGKIT_DATA_ROOT)")

    if args.command == "status":
        status = writer_status(args.data_root)
        print(json.dumps(status, indent=2, sort_keys=True))
        return EXIT_OK if status["quiescent"] else EXIT_LOCK_BUSY

    command = list(args.child)
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        parser.error("run needs a command after --")
    try:
        lock = acquire(args.data_root, args.kind, argv=command)
    except CampaignLockBusy as busy:
        print(f"campaign lock busy: {busy}", file=sys.stderr)
        return EXIT_LOCK_BUSY
    except (OSError, ValueError) as error:
        print(f"cannot take the campaign lock: {error}", file=sys.stderr)
        return EXIT_ENVIRONMENT
    try:
        return _run_child(lock, command)
    except OSError as error:
        print(f"cannot start {command[0]}: {error}", file=sys.stderr)
        return EXIT_ENVIRONMENT
    finally:
        lock.release()


if __name__ == "__main__":
    raise SystemExit(main())
