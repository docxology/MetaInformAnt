"""Real filesystem and launcher contracts for separated RNA runtime roots."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SOURCE_OPS = Path(os.environ.get("RNA_OPERATOR_SOURCE_DIR", ROOT / "output/hymenoptera_local_20261019"))


def require_operator_files(*names: str) -> Path:
    """Require ignored operator assets before their qualification tests.

    Only absent prerequisites skip. Existing malformed assets and failures
    after admission remain ordinary test failures.
    """
    if (SOURCE_OPS.exists() or SOURCE_OPS.is_symlink()) and not SOURCE_OPS.is_dir():
        raise ValueError(f"RNA_OPERATOR_SOURCE_DIR is not a directory: {SOURCE_OPS}")
    missing = [name for name in names if not (SOURCE_OPS / name).exists() and not (SOURCE_OPS / name).is_symlink()]
    if missing:
        pytest.skip(
            "Operator-local qualification assets are absent: "
            + ", ".join(missing)
            + "; provide RNA_OPERATOR_SOURCE_DIR or the ignored campaign output directory"
        )
    for name in names:
        if not (SOURCE_OPS / name).is_file():
            raise ValueError(f"RNA operator asset is not a file: {SOURCE_OPS / name}")
    return SOURCE_OPS


def require_ehd_environment() -> None:
    """Require declared bindings without accepting invalid provided values."""
    missing = [
        name
        for name in ("RNA_EHD_DATA_NAMESPACE", "RNA_EHD_DATA_UUID", "UV_PROJECT_ENVIRONMENT")
        if name not in os.environ
    ]
    if missing:
        pytest.skip("EHD qualification bindings are absent: " + ", ".join(missing))


def profile_environment(tmp_path: Path) -> dict[str, str]:
    require_ehd_environment()
    operations = tmp_path / "runtime with spaces"
    operations.mkdir()
    data = tmp_path / "bulk"
    data.mkdir()
    configs = operations / "configs"
    configs.mkdir()
    frozen = operations / "frozen"
    frozen.mkdir()
    (frozen / "manifest.jsonl").write_text("{}\n")
    source_manifest = operations / "source.json"
    source_manifest.write_text(json.dumps({"files": {}}))
    native = operations / "native.json"
    native.write_text("{}\n")
    db = operations / "progress.db"
    db.touch()
    env = dict(os.environ)
    env.update(
        {
            "AMALGKIT_LOCAL_STORAGE_MODE": "ehd_only",
            "AMALGKIT_OPERATIONS_ROOT": str(operations),
            "METAINFORMANT_ROOT": str(ROOT),
            "AMALGKIT_DATA_NAMESPACE": os.environ["RNA_EHD_DATA_NAMESPACE"],
            "AMALGKIT_DATA_VOLUME_UUID": os.environ["RNA_EHD_DATA_UUID"],
            "AMALGKIT_DATA_ROOT": str(data),
            "AMALGKIT_FROZEN_ROOT": str(tmp_path),
            "AMALGKIT_CONFIG_DIR": str(configs),
            "AMALGKIT_PROGRESS_DB_PATH": str(db),
            "AMALGKIT_PROGRESS_DB_JOURNAL_MODE": "DELETE",
            "AMALGKIT_RUNTIME_SOURCE_MANIFEST": str(source_manifest),
            "AMALGKIT_TOOLCHAIN_MANIFEST": str(native),
            "AMALGKIT_RUNTIME_DIR": os.environ["UV_PROJECT_ENVIRONMENT"],
            "AMALGKIT_TOOLCHAIN_BIN": str(Path(os.environ["UV_PROJECT_ENVIRONMENT"]).parent / "bin"),
        }
    )
    env.pop("AMALGKIT_LOCAL_SIZE_OVERLAY_PATH", None)
    env.pop("AMALGKIT_LOCAL_SIZE_OVERLAY_SHA256", None)
    return env


def inspect(env: dict[str, str]) -> subprocess.CompletedProcess[str]:
    require_operator_files("run.sh")
    return subprocess.run(
        ["/bin/bash", str(SOURCE_OPS / "run.sh"), "print-storage-profile"],
        env=env,
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )


@pytest.mark.external_tool
def test_ehd_launcher_separates_source_from_writable_runtime(tmp_path: Path) -> None:
    env = profile_environment(tmp_path)
    result = inspect(env)
    assert result.returncode == 0, result.stdout + result.stderr
    values = dict(line.split("=", 1) for line in result.stdout.splitlines())
    ops = Path(env["AMALGKIT_OPERATIONS_ROOT"])
    assert values["operations_root"] == str(ops)
    assert values["source_root"] == str(ROOT)
    assert values["index_cache"] == str(ops / "index-cache")
    assert values["tmp_root"] == str(ops / "tmp")
    assert values["fastq_stage"] == "no"
    assert values["quant_scratch"] == values["fasterq_scratch"] == values["ena_scratch"] == ""
    assert values["system_floor_gib"] == "not_applicable_validated_ehd"
    assert not (ops / "tmp").exists()
    assert not (ops / "cache").exists()


@pytest.mark.external_tool
@pytest.mark.parametrize(
    "key,value",
    [
        ("AMALGKIT_DATA_VOLUME_UUID", "00000000-0000-0000-0000-000000000000"),
        ("AMALGKIT_OPERATIONS_ROOT", "relative-runtime"),
        ("AMALGKIT_OPERATIONS_ROOT", "/private/tmp"),
        ("AMALGKIT_DATA_NAMESPACE", "/private/tmp"),
        ("AMALGKIT_PROGRESS_DB_JOURNAL_MODE", "WAL"),
        ("AMALGKIT_PROGRESS_DB_PATH", "/private/tmp/nonexistent-rna-db"),
        ("METAINFORMANT_ROOT", "/private/tmp"),
    ],
)
def test_bad_ehd_roots_refuse_before_creating_controls(tmp_path: Path, key: str, value: str) -> None:
    env = profile_environment(tmp_path)
    original = Path(env["AMALGKIT_OPERATIONS_ROOT"])
    env[key] = value
    result = inspect(env)
    assert result.returncode != 0, result.stdout
    assert not (original / "tmp").exists()
    assert not (original / "active_runtime.json").exists()
    assert not (original / ".run.lock").exists()


@pytest.mark.external_tool
@pytest.mark.parametrize("relative", ["tmp", "home", "cache", "cache/uv", "index-cache"])
def test_ehd_write_root_symlink_refuses(tmp_path: Path, relative: str) -> None:
    env = profile_environment(tmp_path)
    ops = Path(env["AMALGKIT_OPERATIONS_ROOT"])
    outside = tmp_path / "outside"
    outside.mkdir()
    redirect = ops / relative
    redirect.parent.mkdir(parents=True, exist_ok=True)
    redirect.symlink_to(outside, target_is_directory=True)
    result = inspect(env)
    assert result.returncode != 0
    assert not list(outside.iterdir())
    assert not (ops / "active_runtime.json").exists()


def load_paths_module():
    path = require_operator_files("runtime_paths.py") / "runtime_paths.py"
    spec = importlib.util.spec_from_file_location("ehd_runtime_paths", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.mark.external_tool
@pytest.mark.parametrize("bad_name", ["../outside", "/private/tmp/outside"])
def test_source_binding_rejects_traversal(tmp_path: Path, bad_name: str) -> None:
    module = load_paths_module()
    manifest = tmp_path / "source.json"
    manifest.write_text(json.dumps({"files": {bad_name: "0" * 64}}))
    with pytest.raises(ValueError, match="Invalid source binding"):
        module.verify_source_files(ROOT, manifest)


@pytest.mark.external_tool
def test_exact_source_binding_changes_refuse(tmp_path: Path) -> None:
    module = load_paths_module()
    source = tmp_path / "source"
    source.mkdir()
    payload = source / "program.py"
    payload.write_text("print('bound')\n")
    manifest = tmp_path / "source.json"
    manifest.write_text(json.dumps({"files": {"program.py": hashlib.sha256(payload.read_bytes()).hexdigest()}}))
    assert module.verify_source_files(source, manifest) == hashlib.sha256(manifest.read_bytes()).hexdigest()
    payload.write_text("print('changed')\n")
    with pytest.raises(RuntimeError, match="binding changed"):
        module.verify_source_files(source, manifest)


@pytest.mark.external_tool
def test_source_binding_refuses_symlink(tmp_path: Path) -> None:
    module = load_paths_module()
    source = tmp_path / "source"
    source.mkdir()
    target = tmp_path / "payload.py"
    target.write_text("print('outside')\n")
    (source / "program.py").symlink_to(target)
    manifest = tmp_path / "source.json"
    manifest.write_text(json.dumps({"files": {"program.py": hashlib.sha256(target.read_bytes()).hexdigest()}}))
    with pytest.raises(ValueError, match="nonsymlink"):
        module.verify_source_files(source, manifest)


def native_environment(tmp_path: Path, fixture) -> dict[str, str]:
    require_operator_files("continuous.py", "campaign.py", "checkpoints.py")
    env = profile_environment(tmp_path)
    ops = Path(env["AMALGKIT_OPERATIONS_ROOT"])
    env.update(
        {
            "AMALGKIT_DATA_ROOT": str(fixture.data),
            "AMALGKIT_CONFIG_DIR": str(fixture.config_dir),
            "AMALGKIT_LOCAL_QUANT_SLOTS": "8",
            "AMALGKIT_LOCAL_ACQUISITION_ORDER": "smallest_first",
            "AMALGKIT_LOCAL_FRONTIER_SIZE": "48",
            "AMALGKIT_LOCAL_RAW_RESERVATION_BYTES": str(1024**3),
            "AMALGKIT_LOCAL_FASTQ_STAGE": "no",
            "AMALGKIT_LOCAL_QUANT_SCRATCH_DIR": "",
            "AMALGKIT_LOCAL_FASTERQ_SCRATCH_DIR": "",
            "AMALGKIT_LOCAL_ENA_SCRATCH_DIR": "",
            "AMALGKIT_PIPELINE_VALIDATE_FASTQ_RECORDS": "yes",
            "AMALGKIT_PIPELINE_ACQUISITION_SLOTS": "8",
            "AMALGKIT_ENA_MAX_CONNECTIONS": "16",
            "AMALGKIT_ENA_API_CONNECTIONS": "2",
            "AMALGKIT_MIN_EXTERNAL_FREE_GB": "256",
            "AMALGKIT_MIN_SYSTEM_FREE_GB": "0",
            "AMALGKIT_PIPELINE_DEDICATED_RECOVERY": "yes",
            "PYTHONPATH": str(ROOT / "src"),
            "PYTHONDONTWRITEBYTECODE": "1",
        }
    )
    for key, relative in [
        ("HOME", "home"),
        ("TMPDIR", "tmp"),
        ("UV_CACHE_DIR", "cache/uv"),
        ("XDG_CACHE_HOME", "cache"),
        ("MPLCONFIGDIR", "cache/matplotlib"),
        ("NUMBA_CACHE_DIR", "cache/numba"),
        ("AMALGKIT_TEMP_ROOT", "tmp"),
        ("AMALGKIT_CACHE_ROOT", "cache"),
        ("AMALGKIT_LOCAL_INDEX_CACHE_DIR", "index-cache"),
    ]:
        path = ops / relative
        path.mkdir(parents=True, exist_ok=True)
        env[key] = str(path)
    for name in ["manifest.jsonl", "snapshot.json"]:
        shutil.copy2(tmp_path / name, ops / "frozen" / name)
    snapshot = json.loads((tmp_path / "snapshot.json").read_text())
    for row in snapshot["input_files"]:
        dest = ops / "frozen" / row["path"]
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(tmp_path / row["path"], dest)
    source = SOURCE_OPS / "continuous.py"
    Path(env["AMALGKIT_RUNTIME_SOURCE_MANIFEST"]).write_text(
        json.dumps({"files": {str(source.relative_to(ROOT)): hashlib.sha256(source.read_bytes()).hexdigest()}})
    )
    (ops / "preparation_complete.json").write_text("{}\n")
    return env


def continuous(env: dict[str, str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(continuous_command(env), env=env, capture_output=True, text=True, timeout=90, check=False)


def sandbox_command(env: dict[str, str], command: list[str]) -> list[str]:
    if not Path("/usr/bin/sandbox-exec").is_file():
        pytest.skip("The macOS sandbox-exec qualification tool is unavailable")
    private = Path(env["AMALGKIT_OPERATIONS_ROOT"]).parent
    policy = (
        "(version 1)(allow default)(deny file-write*)"
        "(allow file-write* (subpath " + json.dumps(str(private)) + "))"
        '(allow file-write* (literal "/dev/null"))'
    )
    return [
        "/usr/bin/sandbox-exec",
        "-p",
        policy,
        "/bin/bash",
        "-c",
        'export DYLD_LIBRARY_PATH=/opt/homebrew/Cellar/brotli/1.2.0/lib; exec "$@"',
        "--",
        *command,
    ]


def continuous_command(env: dict[str, str]) -> list[str]:
    require_operator_files("continuous.py", "campaign.py", "checkpoints.py")
    return sandbox_command(
        env,
        [
            sys.executable,
            "-B",
            str(SOURCE_OPS / "continuous.py"),
            "--data-root",
            env["AMALGKIT_DATA_ROOT"],
            "--config-dir",
            env["AMALGKIT_CONFIG_DIR"],
            "--workers",
            "12",
            "--quant-slots",
            "8",
        ],
    )


def host_sampler_command(command: list[str]) -> list[str]:
    require_operator_files("host_sample.py", "sample_rate.py")
    sampler = SOURCE_OPS / "sample_rate.py"
    return [
        sys.executable,
        "-B",
        str(SOURCE_OPS / "host_sample.py"),
        "--sampler-sha256",
        hashlib.sha256(sampler.read_bytes()).hexdigest(),
        "--ps-sha256",
        hashlib.sha256(Path("/bin/ps").read_bytes()).hexdigest(),
        *command[3:],
    ]


@pytest.mark.external_tool
@pytest.mark.parametrize("operation", ["write", "mkdir", "rename", "unlink", "chmod", "other_exec", "other_ps_flags"])
def test_read_only_host_collector_rejects_mutation_and_other_commands(tmp_path: Path, operation: str) -> None:
    require_operator_files("host_sample.py")
    source = tmp_path / "original.bin"
    source.write_bytes(b"preserved host-control fixture\n")
    before = (source.stat().st_ino, source.stat().st_mode, source.stat().st_mtime_ns, source.read_bytes())
    code = """
import importlib.util,json,os,subprocess,sys
from pathlib import Path
spec=importlib.util.spec_from_file_location("collector_guard",sys.argv[1])
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
module.install_read_only_guard(sys.argv[2])
path=Path(sys.argv[3]); operation=sys.argv[4]
try:
    if operation=="write":path.write_bytes(b"changed")
    elif operation=="mkdir":path.with_name("new-directory").mkdir()
    elif operation=="rename":path.rename(path.with_name("moved.bin"))
    elif operation=="unlink":path.unlink()
    elif operation=="chmod":path.chmod(0o777)
    elif operation=="other_exec":subprocess.run(["/usr/bin/true"],check=True)
    elif operation=="other_ps_flags":subprocess.run(["/bin/ps","-p",str(os.getpid())],check=True)
    else:raise AssertionError(operation)
except PermissionError as error:
    print(json.dumps({"refused":operation,"reason":str(error)}))
else:raise AssertionError("Host mutation or unbound command was accepted")
"""
    result = subprocess.run(
        [
            sys.executable,
            "-B",
            "-c",
            code,
            str(SOURCE_OPS / "host_sample.py"),
            hashlib.sha256(Path("/bin/ps").read_bytes()).hexdigest(),
            str(source),
            operation,
        ],
        env=dict(os.environ),
        capture_output=True,
        text=True,
        timeout=15,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads(result.stdout)["refused"] == operation
    assert before == (source.stat().st_ino, source.stat().st_mode, source.stat().st_mtime_ns, source.read_bytes())
    assert sorted(p.name for p in tmp_path.iterdir()) == ["original.bin"]


@pytest.mark.external_tool
def test_host_collector_allows_only_pinned_native_ps_despite_path_control(tmp_path: Path) -> None:
    require_operator_files("host_sample.py", "sample_rate.py")
    fake = tmp_path / "ps"
    marker = tmp_path / "unexpected-execution"
    fake.write_text("#!/bin/sh\nprintf unexpected > " + str(marker) + "\n")
    fake.chmod(0o700)
    code = """
import importlib.util,json,os,sys
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
guard=load("collector",sys.argv[1]);sampler=load("sampler",sys.argv[2])
guard.install_read_only_guard(sys.argv[3]);rows=sampler.processes()
assert os.getpid() in rows
print(json.dumps({"native_ps_found_self":True}))
"""
    env = dict(os.environ)
    env["PATH"] = str(tmp_path) + os.pathsep + env["PATH"]
    result = subprocess.run(
        [
            sys.executable,
            "-B",
            "-c",
            code,
            str(SOURCE_OPS / "host_sample.py"),
            str(SOURCE_OPS / "sample_rate.py"),
            hashlib.sha256(Path("/bin/ps").read_bytes()).hexdigest(),
        ],
        env=env,
        capture_output=True,
        text=True,
        timeout=15,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads(result.stdout)["native_ps_found_self"]
    assert not marker.exists()


@pytest.mark.external_tool
@pytest.mark.parametrize("bad_binding", ["sampler", "ps", "output", "optimized"])
def test_host_collector_refuses_bad_eligibility_before_sampling(tmp_path: Path, bad_binding: str) -> None:
    require_operator_files("host_sample.py", "sample_rate.py")
    env = profile_environment(tmp_path)
    ops = Path(env["AMALGKIT_OPERATIONS_ROOT"])
    command = [
        sys.executable,
        "-B",
        str(SOURCE_OPS / "sample_rate.py"),
        "--operations",
        str(ops),
        "--source-root",
        str(ROOT),
        "--manifest",
        env["AMALGKIT_RUNTIME_SOURCE_MANIFEST"],
        "--seconds",
        "1",
    ]
    actual = host_sampler_command(command)
    if bad_binding in {"sampler", "ps"}:
        option = "--sampler-sha256" if bad_binding == "sampler" else "--ps-sha256"
        actual[actual.index(option) + 1] = "0" * 64
    elif bad_binding == "output":
        actual += ["--output", str(ops / "STOP")]
    else:
        actual.insert(1, "-O")
    result = subprocess.run(actual, env=env, capture_output=True, text=True, timeout=15)
    assert result.returncode != 0
    expected = {
        "sampler": "sampler source binding changed",
        "ps": "executable binding changed",
        "output": "unrecognized arguments",
        "optimized": "normal assertions",
    }[bad_binding]
    assert expected in result.stderr, result.stdout + result.stderr
    assert not (ops / "STOP").exists()
    assert not (ops / "active_runtime.json").exists()


@pytest.mark.external_tool
def test_worker_write_sandbox_remains_enforced(tmp_path: Path) -> None:
    env = profile_environment(tmp_path)
    outside = tmp_path.with_name(tmp_path.name + "-outside-worker.bin")
    outside.write_bytes(b"preserve worker boundary\n")
    before = outside.read_bytes(), outside.stat().st_mtime_ns
    code = """
import json,sys
from pathlib import Path
allowed=Path(sys.argv[1]);outside=Path(sys.argv[2]);allowed.write_bytes(b"allowed")
try:outside.write_bytes(b"forbidden")
except PermissionError:print(json.dumps({"allowed":True,"outside_refused":True}))
else:raise AssertionError("Worker write guard lost")
"""
    result = subprocess.run(
        sandbox_command(env, [sys.executable, "-B", "-c", code, str(tmp_path / "allowed.bin"), str(outside)]),
        env=env,
        capture_output=True,
        text=True,
        timeout=15,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads(result.stdout)["outside_refused"]
    assert (tmp_path / "allowed.bin").read_bytes() == b"allowed"
    assert before == (outside.read_bytes(), outside.stat().st_mtime_ns)


@pytest.mark.external_tool
def test_real_ehd_continuous_quant_checkpoint_and_replay(native_acquisition_factory, tmp_path: Path) -> None:
    fixture = native_acquisition_factory()
    env = native_environment(tmp_path, fixture)
    ops = Path(env["AMALGKIT_OPERATIONS_ROOT"])
    first = continuous(env)
    assert first.returncode == 0, first.stdout + first.stderr
    runtime = json.loads((ops / "active_runtime.json").read_text())
    assert runtime["operations_root"] == str(ops) and runtime["source_root"] == str(ROOT)
    assert runtime["storage_mode"] == "ehd_only" and runtime["quant_slots"] == 8
    sample = fixture.work / "quant/SRR123"
    assert (sample / ".local_restore_verified.json").is_file()
    before = {
        p.name: (hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_mtime_ns)
        for p in sample.iterdir()
        if p.is_file() and not p.name.startswith(".local_")
    }
    replay = continuous(env)
    assert replay.returncode == 0, replay.stdout + replay.stderr
    after = {
        p.name: (hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_mtime_ns)
        for p in sample.iterdir()
        if p.is_file() and not p.name.startswith(".local_")
    }
    assert before == after and len(before) >= 3
    assert (ops / "checkpoint_progress.json").is_file()
    assert (fixture.data / "locked-quant").is_dir()
    assert all(
        p.stat().st_dev == Path(env["AMALGKIT_DATA_NAMESPACE"]).stat().st_dev for p in [sample, ops, fixture.data]
    )
    assert not (SOURCE_OPS / "active_runtime.json").exists()


@pytest.mark.external_tool
def test_ehd_existing_stop_is_preserved(native_acquisition_factory, tmp_path: Path) -> None:
    fixture = native_acquisition_factory()
    env = native_environment(tmp_path, fixture)
    ops = Path(env["AMALGKIT_OPERATIONS_ROOT"])
    stop = ops / "STOP"
    stop.write_bytes(b"intentional user pause\n")
    before = stop.stat().st_ino, stop.stat().st_mtime_ns, stop.read_bytes()
    result = continuous(env)
    assert result.returncode == 2, result.stdout + result.stderr
    assert before == (stop.stat().st_ino, stop.stat().st_mtime_ns, stop.read_bytes())
    assert json.loads((ops / "continuous_final.json").read_text())["stopped"]
    assert not (fixture.work / "quant").exists()


@pytest.mark.external_tool
def test_ehd_held_operations_lock_refuses(native_acquisition_factory, tmp_path: Path) -> None:
    import fcntl

    fixture = native_acquisition_factory()
    env = native_environment(tmp_path, fixture)
    ops = Path(env["AMALGKIT_OPERATIONS_ROOT"])
    with (ops / ".run.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        inode = (ops / ".run.lock").stat().st_ino
        result = continuous(env)
        assert result.returncode != 0
        assert "BlockingIOError" in result.stderr
        assert (ops / ".run.lock").stat().st_ino == inode
        assert not (ops / "active_runtime.json").exists()


@pytest.mark.external_tool
@pytest.mark.timeout(180)
def test_ehd_live_sampler_stdout_and_explicit_observer(native_acquisition_factory, tmp_path: Path) -> None:
    require_operator_files("continuous.py", "campaign.py", "checkpoints.py", "host_sample.py", "sample_rate.py")
    observer = ROOT / "projects/hymenoptera_amalgkit/scripts/observe_acquisition.py"
    if not observer.is_file():
        pytest.skip("The operator-local Hymenoptera observer script is unavailable")
    fixture = native_acquisition_factory()
    raw = fixture.work / "getfastq/SRR123/SRR123.fastq.gz"
    payload = raw.read_bytes()
    raw.rename(raw.with_suffix(".preserved"))
    env = native_environment(tmp_path, fixture)
    # A real EHD request must not be blocked by an unrelated host free-space floor.
    env["AMALGKIT_MIN_SYSTEM_FREE_GB"] = "1000000"
    ops = Path(env["AMALGKIT_OPERATIONS_ROOT"])
    requested, release = threading.Event(), threading.Event()

    class Handler(BaseHTTPRequestHandler):
        def do_HEAD(self):
            self.send_response(200)
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()

        def do_GET(self):
            if self.path.startswith("/api"):
                url = f"http://127.0.0.1:{self.server.server_port}/SRR123.fastq.gz"
                md5 = hashlib.md5(payload, usedforsecurity=False).hexdigest()
                body = (
                    "run_accession\tfastq_ftp\tfastq_bytes\tfastq_md5\n" f"SRR123\t{url}\t{len(payload)}\t{md5}\n"
                ).encode()
            else:
                requested.set()
                if not release.wait(30):
                    self.send_error(503, "Private native sample budget expired")
                    return
                body = payload
            self.send_response(200)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, *_):
            return

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever)
    thread.start()
    env["AMALGKIT_ENA_API_BASE_URL"] = f"http://127.0.0.1:{server.server_port}/api"
    process = None
    try:
        with (tmp_path / "live-continuous.log").open("x") as log:
            process = subprocess.Popen(
                continuous_command(env), env=env, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL
            )
            assert requested.wait(25), (tmp_path / "live-continuous.log").read_text()
            command = [
                sys.executable,
                "-B",
                str(SOURCE_OPS / "sample_rate.py"),
                "--operations",
                str(ops),
                "--source-root",
                str(ROOT),
                "--manifest",
                env["AMALGKIT_RUNTIME_SOURCE_MANIFEST"],
                "--seconds",
                "1",
            ]
            refusal = subprocess.run(
                sandbox_command(env, [*command, "--output", str(ops / "STOP")]),
                env=env,
                capture_output=True,
                text=True,
                timeout=15,
            )
            assert refusal.returncode != 0 and "stdout only" in refusal.stderr
            assert not (ops / "STOP").exists()
            # Host process observation is read-only and deliberately outside
            # the worker write sandbox. The worker and artifact writer retain
            # their unchanged kernel-enforced Data-only write restrictions.
            measured = subprocess.run(
                host_sampler_command(command), env=env, capture_output=True, text=True, timeout=15
            )
            assert measured.returncode == 0, measured.stdout + measured.stderr
            rate = json.loads(measured.stdout)
            assert rate["exactly_one_python_producer"] == [process.pid]
            assert rate["seconds"] >= 1 and rate["runtime"]["operations_root"] == str(ops)
            rate_path = tmp_path / "rate.json"
            rate_path.write_text(measured.stdout)
            release.set()
            assert process.wait(timeout=60) == 0, (tmp_path / "live-continuous.log").read_text()
        selection = next((ops / "partitions").glob("continuous_*_primary.json"))
        journal = next(fixture.data.glob("acquisition_runs/*/task_results.jsonl"))
        observe = [
            sys.executable,
            "-B",
            str(ROOT / "projects/hymenoptera_amalgkit/scripts/observe_acquisition.py"),
            "--operations-root",
            str(ops),
            "--data-root",
            str(fixture.data),
            "--manifest",
            str(ops / "frozen/manifest.jsonl"),
            "--task-selection",
            str(selection),
            "--journal",
            str(journal),
            "--transfer-observation",
            str(rate_path),
            "--output-dir",
            str(tmp_path / "observation"),
            "--no-supplement-pdf",
        ]
        observed = subprocess.run(sandbox_command(env, observe), env=env, capture_output=True, text=True, timeout=60)
        assert observed.returncode == 0, observed.stdout + observed.stderr
        observation = json.loads((tmp_path / "observation/observation.json").read_text())
        assert observation["scope"]["total"] == 1
        assert observation["counts"]["new_restore_bound"] == 1
        assert observation["counts"]["reused"] == 0
    finally:
        release.set()
        if process is not None and process.poll() is None:
            process.terminate()  # Exact owned private fixture producer only.
            process.wait(timeout=20)
        server.shutdown()
        server.server_close()
        thread.join(timeout=10)
        assert not thread.is_alive()


def test_broken_operator_source_root_is_an_error_not_a_skip(tmp_path: Path) -> None:
    source = tmp_path / "broken-operator-root"
    target = tmp_path / "missing-operator-root"
    source.symlink_to(target, target_is_directory=True)
    environment = dict(os.environ, RNA_OPERATOR_SOURCE_DIR=str(source))
    environment["PYTHONPATH"] = os.pathsep.join([str(ROOT), str(ROOT / "src")])
    result = subprocess.run(
        [
            sys.executable,
            "-B",
            "-c",
            "from tests.rna.test_ehd_runtime_paths import require_operator_files; require_operator_files('run.sh')",
        ],
        cwd=ROOT,
        env=environment,
        capture_output=True,
        text=True,
        timeout=20,
    )
    (tmp_path / "broken-root-child.json").write_text(
        json.dumps({"exit_code": result.returncode, "stdout": result.stdout, "stderr": result.stderr}, indent=2)
    )
    assert result.returncode == 1
    assert "ValueError: RNA_OPERATOR_SOURCE_DIR is not a directory" in result.stderr
    assert "Skipped:" not in result.stderr
    assert source.is_symlink() and not target.exists()
