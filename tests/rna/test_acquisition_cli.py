"""Real subprocess contracts for generic acquisition entry points."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from metainformant.rna.engine.acquisition_manifest import sha256_file
from metainformant.rna.engine.aws_completion import _render_startup

ROOT = Path(__file__).resolve().parents[2]


def test_plan_and_estimate_cli_use_real_immutable_files(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.jsonl"
    manifest.write_text(
        "".join(
            json.dumps(
                {
                    "schema": "metainformant.rna.acquisition_task.v1",
                    "task_id": f"ant/SRR{i}",
                    "species": "ant",
                    "accession": f"SRR{i}",
                    "config_name": "amalgkit_ant.yaml",
                    "batch_index": i,
                }
            )
            + "\n"
            for i in range(1, 5)
        )
    )
    (tmp_path / "snapshot.json").write_text(
        json.dumps(
            {
                "schema": "metainformant.rna.acquisition_snapshot.v1",
                "manifest_sha256": sha256_file(manifest),
                "task_count": 4,
                "cloud_launch_policy": "checkpointed",
                "input_files": [],
            }
        )
    )
    entry = ROOT / "scripts/rna/acquisition.py"
    command = [
        sys.executable,
        str(entry),
        "plan",
        "--manifest",
        str(manifest),
        "--backend",
        "hybrid",
        "--output-dir",
        str(tmp_path / "plan"),
    ]
    first = subprocess.run(command, check=True, text=True, capture_output=True)
    assert json.loads(first.stdout)["local_pending"] == 2
    subprocess.run(command, check=True, text=True, capture_output=True)
    rates = tmp_path / "rates.json"
    rates.write_text(
        json.dumps(
            {
                lane: {
                    "evidence": {
                        "observed_units": 2,
                        "samples_per_hour_low": 2,
                        "samples_per_hour_high": 4,
                        "source": "controlled arithmetic scenario",
                    },
                    "costs": {"hourly_usd_per_unit": 0 if lane == "local" else 0.55},
                }
                for lane in ("local", "aws")
            }
        )
    )
    output = tmp_path / "estimate.json"
    result = subprocess.run(
        [
            sys.executable,
            str(entry),
            "estimate",
            "--allocation",
            str(tmp_path / "plan/allocation.json"),
            "--rates",
            str(rates),
            "--local-units",
            "2",
            "--aws-units",
            "2",
            "--spent-usd",
            "10",
            "--reserved-usd",
            "3",
            "--ceiling-usd",
            "15",
            "--output",
            str(output),
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    payload = json.loads(result.stdout)
    assert payload["total"]["gross_high_usd"] == pytest.approx(14.1)
    assert payload["total"]["hours_high"] == 1
    assert payload["total"]["fits_conservative_ceiling"] is True
    assert json.loads(output.read_text()) == payload


def test_generic_aws_template_is_complete_shell_safe_and_syntax_valid(
    tmp_path: Path,
) -> None:
    bindings = {
        "BUCKET": "bucket",
        "REGION": "us-east-2",
        "COHORT": "cohort",
        "SPECIES": "ant",
        "SOURCE_KEY": "source.tar",
        "SOURCE_SHA": "a" * 64,
        "INPUT_KEY": "input.tar",
        "INPUT_SHA": "b" * 64,
        "JOB_PREFIX": "jobs/1",
        "LIMIT_SECONDS": 14400,
        "RAW_BYTES": 65536,
        "WORKERS": 16,
        "THREADS": 8,
        "QUANT_SLOTS": 4,
        "FASTQ_SLOTS": 1,
        "MAX_IN_FLIGHT": 12,
        "FASTQ_THREADS": 2,
        "COMPRESSION_THREADS": 2,
        "COMPRESSION_LEVEL": 1,
        "ENA_FILE_WORKERS": 2,
        "VALIDATION_SLOTS": 4,
    }
    text = _render_startup(ROOT / "scripts/rna/aws_acquisition_startup.sh", bindings)
    assert "@@" not in text and "hymenoptera_amalgkit" not in text
    assert "--config-dir /mnt/snapshot/config/amalgkit" in text
    assert "--ena-file-workers 2" in text
    assert "--compression-level 1" in text
    script = tmp_path / "startup.sh"
    script.write_text(text)
    subprocess.run(["bash", "-n", str(script)], check=True, capture_output=True, text=True)


@pytest.mark.parametrize("lane", ["local", "aws"])
@pytest.mark.parametrize(
    "entry", ["scripts/rna/acquisition.py", "projects/hymenoptera_amalgkit/scripts/acquisition.py"]
)
def test_generic_and_project_execution_help_is_callable(lane: str, entry: str) -> None:
    if entry.startswith("projects/") and not (ROOT / "projects/hymenoptera_amalgkit/README.md").is_file():
        pytest.skip("Hymenoptera acquisition adapter is unavailable in this parent-only checkout")
    result = subprocess.run(
        [sys.executable, str(ROOT / entry), lane, "--help"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "--config-dir" in result.stdout
