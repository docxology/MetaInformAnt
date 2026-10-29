from __future__ import annotations
import os
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path
from uuid import uuid4
import pytest

@pytest.mark.parametrize("value", ["0", "9", "invalid"])
@pytest.mark.parametrize("environment", [False, True])
def test_invalid_quant_slots_refuse_before_controls(tmp_path, value, environment):
    root = Path(__file__).resolve().parents[2]
    source = root / "output/hymenoptera_local_20261019"
    private = tmp_path / "operations"
    private.mkdir()
    for name in ("continuous.py", "campaign.py", "checkpoints.py"):
        (private / name).write_bytes((source / name).read_bytes())
    env = dict(os.environ)
    args = [sys.executable, str(private / "continuous.py"), "--data-root", str(tmp_path / "data"), "--config-dir", str(tmp_path / "configs")]
    if environment:
        env["AMALGKIT_LOCAL_QUANT_SLOTS"] = value
    else:
        args.extend(["--quant-slots", value])
    run = subprocess.run(args, env=env, capture_output=True, text=True, timeout=20)
    assert run.returncode == 2, run.stderr
    assert "quant" in run.stderr
    assert not (private / "active_runtime.json").exists()
    assert not (private / ".run.lock").exists()


@pytest.mark.parametrize("slots", [4, 8])
def test_real_cli_scope_reports_selected_quant_profile(native_acquisition_factory, tmp_path, slots):
    fixture = native_acquisition_factory()
    with fixture.metadata.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields, rows = list(reader.fieldnames), list(reader)
    fields.extend(["lib_source", "lib_selection", "exp_title"])
    for row in rows:
        row.update(lib_source="GENOMIC", lib_selection="PCR", exp_title="16S amplicon control")
    with fixture.metadata.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    volume = next(path for path in Path(sys.base_prefix).resolve().parents if path.is_mount())
    assert volume != Path("/") and (volume / "Local").is_dir()
    # The real CLI requires bulk.parent.parent to be a mounted volume. Retain
    # this tiny owned fixture on GitVault; never weaken that production guard.
    bulk = volume / "Local" / ("hymenoptera-cli-profile-" + uuid4().hex)
    fixture.data.rename(bulk)
    private = tmp_path / "operations"
    (private / "frozen").mkdir(parents=True)
    source = Path(__file__).resolve().parents[2] / "output/hymenoptera_local_20261019"
    for name in ("continuous.py", "campaign.py", "checkpoints.py"):
        shutil.copy2(source / name, private / name)
    shutil.copy2(fixture.manifest, private / "frozen/manifest.jsonl")
    shutil.copy2(tmp_path / "snapshot.json", private / "frozen/snapshot.json")
    (private / "preparation_complete.json").write_text("{}\n")
    env = dict(os.environ)
    env.pop("AMALGKIT_RUNTIME_SOURCE_MANIFEST", None)
    env["AMALGKIT_LOCAL_QUANT_SLOTS"] = str(slots)
    result = subprocess.run([sys.executable, str(private / "continuous.py"), "--data-root", str(bulk), "--config-dir", str(fixture.config_dir)], env=env, capture_output=True, text=True, timeout=30)
    assert result.returncode == 1, result.stdout + result.stderr
    runtime = json.loads((private / "active_runtime.json").read_text())
    scope = json.loads((private / "active_scope.json").read_text())
    assert scope["eligible"] == 0 and len(scope["unsupported"]) == 1
    assert "16S" in next(iter(scope["unsupported"].values()))
    assert runtime["quant_slots"] == scope["quant_slots"] == slots
    assert scope["quant_thread_budget"] == 8
    assert scope["quant_threads_per_sample"] == 8 // slots
    assert not list((private / "partitions").glob("*.json"))


def test_real_cli_records_smallest_first_without_changing_profile_or_holds(native_acquisition_factory, tmp_path, monkeypatch):
    monkeypatch.setenv("AMALGKIT_LOCAL_ACQUISITION_ORDER", "smallest_first")
    test_real_cli_scope_reports_selected_quant_profile(native_acquisition_factory, tmp_path, 8)
    private = tmp_path / "operations"
    for name in ("active_runtime.json", "active_scope.json"):
        assert json.loads((private / name).read_text())["acquisition_order"] == "smallest_first"
