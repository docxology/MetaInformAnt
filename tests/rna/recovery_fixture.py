"""Deadline-preserving native private fixture runner with exact child cleanup."""

from __future__ import annotations

import json
import os
import subprocess
import time
from pathlib import Path

import psutil


def run_private(command: list[str], root: Path, *, timeout: float = 180) -> subprocess.CompletedProcess[str]:
    if not root.is_absolute() or root.resolve(strict=True) != root:
        raise ValueError("Physical private fixture root required")
    if not any(str(root) in argument for argument in command):
        raise ValueError("Command must name its owned private fixture")
    started = time.time()
    # Reserve every receipt exclusively before a child can execute. Partial
    # reservations are retained if a later reservation refuses.
    with (
        (root / "private-run-start.json").open("x") as start_handle,
        (root / "private-run-finished.json").open("x") as finish_handle,
        (root / "private-timeout-signals.jsonl").open("x") as signal_handle,
    ):
        child = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True,
        )
        first = psutil.Process(child.pid)
        prelaunch = {
            "pid": child.pid,
            "create_time": first.create_time(),
            "argv": command,
            "started_epoch": started,
            "deadline_seconds": timeout,
            "private_root": str(root),
            "receipts_reserved_before_child": True,
        }
        json.dump(prelaunch, start_handle, indent=2)
        start_handle.flush()
        os.fsync(start_handle.fileno())
        try:
            stdout, stderr = child.communicate(timeout=timeout)
        except subprocess.TimeoutExpired as expired:
            targets = []
            for process in [*first.children(recursive=True), first]:
                try:
                    if process.uids().real != os.getuid():
                        raise RuntimeError("Private descendant ownership changed")
                    targets.append(
                        {
                            "pid": process.pid,
                            "create_time": process.create_time(),
                            "argv": process.cmdline(),
                            "uid": process.uids().real,
                        }
                    )
                except psutil.NoSuchProcess:
                    continue
            for target in targets:
                try:
                    process = psutil.Process(target["pid"])
                    if process.create_time() != target["create_time"]:
                        raise RuntimeError("Private descendant identity changed")
                    signal_handle.write(json.dumps({**target, "planned_action": "SIGTERM"}) + "\n")
                    signal_handle.flush()
                    os.fsync(signal_handle.fileno())
                    process.terminate()
                    signal_handle.write(json.dumps({**target, "action": "SIGTERM_RETURNED"}) + "\n")
                except psutil.NoSuchProcess:
                    signal_handle.write(json.dumps({**target, "action": "natively_absent"}) + "\n")
                signal_handle.flush()
                os.fsync(signal_handle.fileno())
            stdout, stderr = child.communicate(timeout=15)
            alive, terminal = [], []
            for target in targets:
                try:
                    process = psutil.Process(target["pid"])
                    if process.create_time() != target["create_time"]:
                        raise RuntimeError("Private descendant PID was reused")
                    state = process.status()
                    terminal.append({**target, "native_state": state})
                    if state != psutil.STATUS_ZOMBIE:
                        alive.append(target)
                except psutil.NoSuchProcess:
                    terminal.append({**target, "native_state": "absent"})
            proof = {
                "timed_out": True,
                "deadline_seconds": timeout,
                "deadline_waived": False,
                "targets": targets,
                "terminal_observations": terminal,
                "remaining_live": alive,
                "returncode": child.returncode,
                "stdout": stdout,
                "stderr": stderr,
            }
            json.dump(proof, finish_handle, indent=2)
            finish_handle.flush()
            os.fsync(finish_handle.fileno())
            if alive:
                raise RuntimeError("Exact private timeout descendants remain; no forced kill")
            raise subprocess.TimeoutExpired(command, timeout, output=stdout, stderr=stderr) from expired
        json.dump(
            {
                "timed_out": False,
                "returncode": child.returncode,
                "elapsed_seconds": time.time() - started,
                "deadline_seconds": timeout,
            },
            finish_handle,
            indent=2,
        )
        finish_handle.flush()
        os.fsync(finish_handle.fileno())
    return subprocess.CompletedProcess(command, child.returncode, stdout, stderr)
