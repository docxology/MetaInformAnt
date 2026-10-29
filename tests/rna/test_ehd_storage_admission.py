"""Real EHD admission, path escape and low-capacity controls."""

from __future__ import annotations

import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from tests.rna.test_ehd_runtime_paths import (
    native_environment,
    require_ehd_environment,
    require_operator_files,
    sandbox_command,
)

pytestmark = pytest.mark.external_tool


@pytest.fixture(autouse=True)
def qualified_ehd_inputs() -> None:
    require_operator_files("continuous.py", "campaign.py", "checkpoints.py")
    require_ehd_environment()


def run_case(environment: dict[str, str], code: str) -> subprocess.CompletedProcess[str]:
    prelude = """
import json,os,time,tempfile
from pathlib import Path
from metainformant.rna.engine.ehd_storage import capture_ehd_write_policy
from metainformant.rna.engine.streaming_orchestrator import StreamingPipelineOrchestrator
config=Path(os.environ['AMALGKIT_CONFIG_DIR'])
bulk=Path(os.environ['AMALGKIT_DATA_ROOT'])
db=Path(os.environ['AMALGKIT_PROGRESS_DB_PATH'])
"""
    return subprocess.run(
        sandbox_command(environment, [sys.executable, "-B", "-c", prelude + code]),
        env=environment,
        capture_output=True,
        text=True,
        timeout=30,
    )


def test_valid_EHD_uses_only_Data_capacity_despite_impossible_system_floor(
    native_acquisition_factory,
    tmp_path: Path,
) -> None:
    fixture = native_acquisition_factory()
    env = native_environment(tmp_path, fixture)
    env["AMALGKIT_MIN_SYSTEM_FREE_GB"] = "1000000"
    result = run_case(
        env,
        """
engine=StreamingPipelineOrchestrator(config,bulk/'logs',db_path=db)
try:
    assert engine.system_root is None
    assert engine.min_external_free_gb >= 256
    assert engine._ehd_write_policy.receipt()['system_capacity_required'] is False
    print(json.dumps(engine._ehd_write_policy.receipt()))
finally:engine.db.close()
""",
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads(result.stdout)["policy"] == "validated_ehd_write_roots_v1"


@pytest.mark.parametrize("mode", ["ssd_compute", "hd_first"])
def test_legacy_modes_still_bind_the_real_system_volume(native_acquisition_factory, tmp_path: Path, mode: str) -> None:
    fixture = native_acquisition_factory()
    env = native_environment(tmp_path, fixture)
    env["AMALGKIT_LOCAL_STORAGE_MODE"] = mode
    env["AMALGKIT_MIN_SYSTEM_FREE_GB"] = "32"
    result = run_case(
        env,
        """
engine=StreamingPipelineOrchestrator(config,bulk/'logs',db_path=db)
try:
    assert engine._ehd_write_policy is None
    assert engine.system_root == Path('/System/Volumes/Data')
    assert engine.min_system_free_gb == 32
finally:engine.db.close()
""",
    )
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize(
    "key,value",
    [
        ("HOME", "/private/tmp"),
        ("TMPDIR", "/private/tmp"),
        ("UV_CACHE_DIR", "/private/tmp"),
        ("MPLCONFIGDIR", "/private/tmp"),
        ("AMALGKIT_LOCAL_INDEX_CACHE_DIR", "/private/tmp"),
        ("AMALGKIT_MIN_EXTERNAL_FREE_GB", "0"),
        ("AMALGKIT_MIN_EXTERNAL_FREE_GB", "nan"),
        ("AMALGKIT_LOCAL_QUANT_SCRATCH_DIR", "/private/tmp"),
        ("AMALGKIT_LOCAL_FASTQ_STAGE", "yes"),
        ("AMALGKIT_PROGRESS_DB_JOURNAL_MODE", "WAL"),
        ("AMALGKIT_DATA_VOLUME_UUID", "00000000-0000-0000-0000-000000000000"),
    ],
)
def test_invalid_EHD_declaration_refuses_before_database_or_logs(
    native_acquisition_factory,
    tmp_path: Path,
    key: str,
    value: str,
) -> None:
    fixture = native_acquisition_factory()
    env = native_environment(tmp_path, fixture)
    db = Path(env["AMALGKIT_PROGRESS_DB_PATH"])
    before = db.stat().st_ino, db.stat().st_mtime_ns, db.read_bytes()
    env[key] = value
    result = run_case(
        env,
        """
try:StreamingPipelineOrchestrator(config,bulk/'logs',db_path=db)
except (ValueError,KeyError) as error:print(json.dumps({'refused':str(error)}))
else:raise AssertionError('Invalid write roots admitted')
""",
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "refused" in json.loads(result.stdout)
    assert before == (db.stat().st_ino, db.stat().st_mtime_ns, db.read_bytes())
    assert not (fixture.data / "logs").exists()


def test_actual_low_EHD_headroom_still_blocks_before_transport(native_acquisition_factory, tmp_path: Path) -> None:
    fixture = native_acquisition_factory()
    env = native_environment(tmp_path, fixture)
    env["AMALGKIT_MIN_EXTERNAL_FREE_GB"] = str(math.ceil(shutil.disk_usage(fixture.data).free / 1024**3) + 1)
    env["AMALGKIT_PIPELINE_DOWNLOAD_FIRST_PASS_SECONDS"] = "1"
    env["AMALGKIT_MIN_SYSTEM_FREE_GB"] = "1000000"
    result = run_case(
        env,
        """
from metainformant.rna.retrieval.ena_staging import StagingDeadline
engine=StreamingPipelineOrchestrator(config,bulk/'logs',db_path=db)
output=bulk/'empty-download';output.mkdir()
try:
    try:engine._download_fastq_admitted('SRR999999',output,False,defer_ncbi_fallback=True)
    except StagingDeadline:print(json.dumps({'EHD_low_refused':True}))
    else:raise AssertionError('Low EHD capacity admitted transport')
    assert not list(output.glob('**/*.part'))
finally:engine.db.close()
""",
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads(result.stdout)["EHD_low_refused"]


def test_live_root_replacement_and_environment_drift_refuse(native_acquisition_factory, tmp_path: Path) -> None:
    fixture = native_acquisition_factory()
    env = native_environment(tmp_path, fixture)
    result = run_case(
        env,
        """
policy=capture_ehd_write_policy(config,bulk/'logs',db)
saved=bulk.with_name('preserved-bulk');bulk.rename(saved);bulk.mkdir()
try:policy.revalidate()
except ValueError:pass
else:raise AssertionError('Replacement root silently rebound')
assert saved.is_dir()
""",
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_cached_Python_tempdir_cannot_escape_the_binding(native_acquisition_factory, tmp_path: Path) -> None:
    fixture = native_acquisition_factory()
    env = native_environment(tmp_path, fixture)
    sibling = tmp_path / "old-temp"
    sibling.mkdir()
    env["QUAL_CORRECT_TEMP"] = env["TMPDIR"]
    env["TMPDIR"] = str(sibling)
    result = run_case(
        env,
        """
assert tempfile.gettempdir()==os.environ['TMPDIR']
os.environ['TMPDIR']=os.environ['QUAL_CORRECT_TEMP']
try:capture_ehd_write_policy(config,bulk/'logs',db)
except ValueError as error:assert 'temporary directory' in str(error)
else:raise AssertionError('Cached out-of-binding TMP accepted')
""",
    )
    assert result.returncode == 0, result.stdout + result.stderr
