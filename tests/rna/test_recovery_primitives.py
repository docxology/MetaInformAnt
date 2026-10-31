"""Native process-state controls without test doubles."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import cast

import psutil
import pytest

from metainformant.rna.engine.process_observation import (
    ProcessIdentityChangedError,
    process_is_confirmed_inert,
)
from tests.rna.recovery_fixture import run_private


@pytest.mark.skipif(not hasattr(os, "fork"), reason="requires native POSIX fork/zombie observations")
@pytest.mark.parametrize("state", ["live", "zombie", "reaped"])
def test_native_process_inertness_requires_current_kernel_proof(tmp_path: Path, state: str) -> None:
    readfd, writefd = os.pipe()
    pid = os.fork()
    if pid == 0:
        os.close(writefd)
        os.read(readfd, 1)
        os._exit(0)
    os.close(readfd)
    closed = waited = False
    try:
        process = psutil.Process(pid)
        created = process.create_time()
        if state != "live":
            os.write(writefd, b"x")
            os.close(writefd)
            closed = True
            if state == "reaped":
                os.waitpid(pid, 0)
                waited = True
            else:
                deadline = time.monotonic() + 5
                while process.status() != psutil.STATUS_ZOMBIE:
                    assert time.monotonic() < deadline
                    time.sleep(0.01)
        assert process_is_confirmed_inert(pid, created) == (state != "live")
        if state != "reaped":
            assert not process_is_confirmed_inert(pid, None)
    finally:
        if not closed:
            os.write(writefd, b"x")
            os.close(writefd)
        if not waited:
            os.waitpid(pid, 0)


def test_reused_identity_is_not_rebaselined() -> None:
    process = psutil.Process(os.getpid())
    with pytest.raises(ProcessIdentityChangedError, match="identity changed"):
        process_is_confirmed_inert(process.pid, process.create_time() - 1)


@pytest.mark.parametrize("pid", [-1, 0, True, 1.5, "1"])
def test_invalid_process_identity_refuses_before_observation(pid: object) -> None:
    with pytest.raises(ValueError, match="positive integer PID"):
        process_is_confirmed_inert(cast(int, pid), None)


@pytest.mark.skipif(sys.platform != "darwin", reason="requires the real macOS process-info sandbox")
def test_real_denied_process_inspection_cannot_claim_absence(tmp_path: Path) -> None:
    import metainformant.rna.engine.process_observation as module

    source = Path(module.__file__).resolve().parents[3]
    code = r"""import ctypes,sys
sys.path.insert(0,sys.argv[1])
from metainformant.rna.engine.process_observation import process_is_confirmed_inert
import psutil
p=psutil.Process(int(sys.argv[2]));created=p.create_time()
lib=ctypes.CDLL('/usr/lib/libsandbox.dylib',use_errno=True)
lib.sandbox_init.argtypes=[ctypes.c_char_p,ctypes.c_uint64,ctypes.POINTER(ctypes.c_char_p)]
error=ctypes.c_char_p()
profile=b'(version 1)(allow default)(deny process-info*)(deny sysctl-read)'
assert lib.sandbox_init(profile,0,ctypes.byref(error))==0
try:p.cmdline()
except psutil.AccessDenied:pass
else:raise AssertionError('Native process query was not denied')
assert process_is_confirmed_inert(p.pid,created) is False
print('live denied process remains unresolved')
"""
    result = subprocess.run(
        [sys.executable, "-I", "-B", "-c", code, str(source), str(os.getpid())],
        capture_output=True,
        text=True,
        timeout=15,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.strip() == "live denied process remains unresolved"


@pytest.mark.parametrize("code", [0, 7])
def test_owned_fixture_runner_preserves_real_terminal_output(tmp_path: Path, code: int) -> None:
    script = tmp_path / "program.py"
    script.write_text(
        "import sys;print('owned stdout');print('owned stderr',file=sys.stderr);raise SystemExit(" + str(code) + ")"
    )
    result = run_private([sys.executable, "-I", "-B", str(script)], tmp_path, timeout=5)
    assert result.returncode == code
    assert result.stdout.strip() == "owned stdout" and result.stderr.strip() == "owned stderr"


@pytest.mark.skipif(os.name != "posix", reason="fixture ownership uses POSIX UID and session semantics")
def test_owned_fixture_timeout_is_still_failure_and_reaps_descendants(
    tmp_path: Path,
) -> None:
    script, child = tmp_path / "parent.py", tmp_path / "child.py"
    child.write_text("import time;print('ready',flush=True);time.sleep(60)")
    script.write_text(
        "import subprocess,sys,time;subprocess.Popen([sys.executable,'-I','-B',"
        + repr(str(child))
        + "]);print('ready',flush=True);time.sleep(60)"
    )
    with pytest.raises(subprocess.TimeoutExpired):
        run_private([sys.executable, "-I", "-B", str(script)], tmp_path, timeout=5)
    proof = json.loads((tmp_path / "private-run-finished.json").read_text())
    assert proof["timed_out"] and not proof["deadline_waived"]
    assert len(proof["targets"]) >= 2 and not proof["remaining_live"]
    assert len(proof["terminal_observations"]) == len(proof["targets"])
    assert all(
        row["uid"] == os.getuid()
        and row["native_state"] in {"absent", psutil.STATUS_ZOMBIE}
        and any(str(tmp_path) in argument for argument in row["argv"])
        for row in proof["terminal_observations"]
    )
    targets = {(row["pid"], row["create_time"]) for row in proof["targets"]}
    ledger = [json.loads(line) for line in (tmp_path / "private-timeout-signals.jsonl").read_text().splitlines()]
    assert all((row["pid"], row["create_time"]) in targets for row in ledger)
    assert any(row.get("action") == "SIGTERM_RETURNED" for row in ledger)


@pytest.mark.parametrize(
    "receipt",
    [
        "private-run-start.json",
        "private-run-finished.json",
        "private-timeout-signals.jsonl",
    ],
)
def test_existing_private_receipt_refuses_before_any_subprocess(tmp_path: Path, receipt: str) -> None:
    target = tmp_path / receipt
    target.write_bytes(b"retained earlier receipt\n")
    before, stamp = target.read_bytes(), target.stat()
    import tests.rna.recovery_fixture as module

    source = Path(module.__file__).resolve().parents[2]
    code = r"""import sys
from pathlib import Path
sys.path.insert(0,sys.argv[1])
from tests.rna.recovery_fixture import run_private
root=Path(sys.argv[2])
def refuse_exec(event,args):
    if event=="subprocess.Popen":
        raise AssertionError("Child creation preceded receipt reservation")
sys.addaudithook(refuse_exec)
try:
    run_private([sys.executable,"-I","-B",str(root/"unstarted.py")],root,timeout=5)
except FileExistsError:
    print("receipt refused before native subprocess creation")
else:
    raise AssertionError("Existing receipt was accepted")
"""
    result = subprocess.run(
        [sys.executable, "-I", "-B", "-c", code, str(source), str(tmp_path)],
        capture_output=True,
        text=True,
        timeout=15,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.strip() == "receipt refused before native subprocess creation"
    after = target.stat()
    assert target.read_bytes() == before
    assert (
        after.st_dev,
        after.st_ino,
        after.st_mode,
        after.st_uid,
        after.st_gid,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
    ) == (
        stamp.st_dev,
        stamp.st_ino,
        stamp.st_mode,
        stamp.st_uid,
        stamp.st_gid,
        stamp.st_size,
        stamp.st_mtime_ns,
        stamp.st_ctime_ns,
    )
