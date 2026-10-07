"""Local campaign helpers remain usable without the optional AWS runtime."""

from __future__ import annotations

import importlib.util
import json
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

from metainformant.rna.engine.completion_inventory import freeze_inventory
from metainformant.rna.engine.durable_quant import DirectoryStore, receipt_key, restore_quantification


def test_local_status_reads_real_database_without_importing_aws(tmp_path: Path) -> None:
    database = tmp_path / "pipeline_progress.db"
    with sqlite3.connect(database) as connection:
        connection.execute("CREATE TABLE samples(species TEXT,srr_id TEXT,state TEXT)")
        connection.execute("INSERT INTO samples VALUES(?,?,?)", ("ant_a", "SRR1", "pending"))
    before = database.read_bytes()
    code = """import sys
from pathlib import Path
from metainformant.rna.engine.campaign_status_io import read_database
records=read_database(Path(sys.argv[1]))
assert [(r.task_id,r.state) for r in records]==[("ant_a/SRR1","pending")]
assert "boto3" not in sys.modules and "botocore" not in sys.modules
"""
    result = subprocess.run([sys.executable, "-c", code, str(database)], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert database.read_bytes() == before


@pytest.mark.skipif(
    importlib.util.find_spec("boto3") is not None, reason="Requires an environment without the AWS extra"
)
def test_cloud_collection_requires_aws_instead_of_falling_back(tmp_path: Path) -> None:
    from metainformant.rna.engine.campaign_status_io import collect_cloud

    with pytest.raises(ModuleNotFoundError, match="boto3"):
        collect_cloud(tmp_path, "unused-bucket", "unused-profile", "us-east-2")


@pytest.mark.parametrize("value", [None, [], "scalar", 42])
def test_nonobject_quant_receipt_is_refused_before_destination_creation(tmp_path: Path, value: object) -> None:
    store = DirectoryStore(tmp_path / "store")
    store.put(receipt_key("cohort", "test_species", "SRR123"), json.dumps(value).encode())
    destination = tmp_path / "restored"
    with pytest.raises(ValueError, match="receipt must be a JSON object"):
        restore_quantification(store, "cohort", "test_species", "SRR123", destination)
    assert not destination.exists()


@pytest.mark.parametrize("value", [None, [], "scalar", 42])
def test_nonobject_existing_inventory_is_refused_without_discovery(tmp_path: Path, value: object) -> None:
    snapshot = tmp_path / "snapshot"
    snapshot.mkdir()
    inventory = snapshot / "inventory.json"
    inventory.write_text(json.dumps(value))
    before = inventory.read_bytes()
    with pytest.raises(ValueError, match="inventory must be a JSON object"):
        freeze_inventory(tmp_path / "data", tmp_path / "config", snapshot, expected_species_count=None)
    assert inventory.read_bytes() == before
