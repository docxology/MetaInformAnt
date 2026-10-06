"""Read-only SQLite/filesystem probes and AWS telemetry for campaign status."""
from __future__ import annotations

import base64
import json
import os
import shlex
import sqlite3
import time
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Final

import boto3
from botocore.config import Config

from metainformant.rna.engine.campaign_status import Inventory, Observation, StatusError
from metainformant.rna.core.sample_utils import quantification_file_candidates
from metainformant.rna.engine.provenance import QUANT_PROVENANCE_FILENAME

REQUIRED_FILES: Final = ("abundance.tsv", "abundance.h5", "run_info.json", QUANT_PROVENANCE_FILENAME)


@dataclass(frozen=True, slots=True)
class WorkerRow:
    task_id: str
    state: str


@dataclass(frozen=True, slots=True)
class WorkerProbe:
    observed_at: str
    databases: tuple[str, ...]
    rows: tuple[WorkerRow, ...]


@dataclass(frozen=True, slots=True)
class Job:
    job_id: str
    instance_id: str | None
    status: str


@dataclass(frozen=True, slots=True)
class Controller:
    cohort: str
    observed_at: str
    jobs: tuple[Job, ...]


def parse_probe(output: str) -> WorkerProbe:
    """Reject malformed or truncated telemetry rather than substituting zeroes."""
    value = json.loads(output)
    if not isinstance(value, dict) or not isinstance(value.get("observed_at"), str):
        raise StatusError("Malformed worker observation")
    if not isinstance(value.get("databases"), list) or not all(isinstance(p, str) for p in value["databases"]):
        raise StatusError("Malformed worker database inventory")
    if not isinstance(value.get("rows"), list):
        raise StatusError("Missing worker sample observations")
    rows = []
    for row in value["rows"]:
        if not isinstance(row, dict) or not isinstance(row.get("task_id"), str) or not isinstance(row.get("state"), str):
            raise StatusError("Malformed worker sample observation")
        rows.append(WorkerRow(row["task_id"], row["state"]))
    return WorkerProbe(value["observed_at"], tuple(value["databases"]), tuple(rows))


@dataclass(frozen=True, slots=True)
class CloudSnapshot:
    locked: frozenset[str]
    assigned: frozenset[str]
    worker: tuple[Observation, ...]
    diagnostics: tuple[str, ...]
    probes: dict[str, str]
    instances: tuple[str, ...]
    observed_at: str


def read_database(path: Path) -> tuple[Observation, ...]:
    """Open an existing SQLite database read-only; never initialize or reconcile it."""
    with closing(sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True)) as connection:
        rows = connection.execute("SELECT species,srr_id,state FROM samples").fetchall()
    return tuple(Observation(f"{s}/{a}", state, str(path)) for s, a, state in rows)


def file_coverage(inventory: Inventory, root: Path) -> tuple[frozenset[str], frozenset[str]]:
    """Observe expected file presence only; do not claim numerical/hash validation."""
    if not root.is_dir():
        raise StatusError(f"Local root is unavailable: {root}")
    complete = set()
    partial = set()
    for species in inventory.species:
        quant = root / species.species / "work" / "quant"
        if not quant.is_dir():
            continue
        # One directory listing per species avoids 18,200 absent SSD path probes.
        accessions = {task.accession for task in species.tasks}
        with os.scandir(quant) as entries:
            directories = {entry.name: Path(entry.path) for entry in entries if entry.name in accessions}
        for task in species.tasks:
            sample = directories.get(task.accession)
            if sample is None:
                continue
            tables = {p.name for p in quantification_file_candidates(sample, task.accession)}
            infos = {"run_info.json", f"{task.accession}_run_info.json"}
            h5s = {"abundance.h5", f"{task.accession}_abundance.h5"}
            expected = tables | infos | h5s | {QUANT_PROVENANCE_FILENAME}
            with os.scandir(sample) as entries:
                nonempty = {entry.name for entry in entries if entry.name in expected and entry.is_file() and entry.stat().st_size > 0}
            covered = bool(tables & nonempty and infos & nonempty and h5s & nonempty and QUANT_PROVENANCE_FILENAME in nonempty)
            (complete if covered else partial).add(task.task_id)
    return frozenset(complete), frozenset(partial)


def probe_script(task_ids: frozenset[str], root: str = "/mnt/amalgkit") -> str:
    """Build a static read-only worker query; supplied identifiers are encoded data."""
    encoded = base64.b64encode(json.dumps(sorted(task_ids)).encode()).decode()
    return f'''import base64,json,sqlite3
from pathlib import Path
from contextlib import closing
from datetime import datetime,timezone
wanted=set(json.loads(base64.b64decode({encoded!r})))
databases=list(Path({root!r}).glob("**/pipeline_progress.db"))
rows=[]
for path in databases:
    with closing(sqlite3.connect(path.resolve().as_uri()+"?mode=ro",uri=True)) as connection:
        for species,accession,state in connection.execute("SELECT species,srr_id,state FROM samples"):
            key=species+"/"+accession
            if key in wanted:
                rows.append(dict(task_id=key,state=state))
print(json.dumps(dict(observed_at=datetime.now(timezone.utc).isoformat(),databases=[str(p) for p in databases],rows=rows)))
'''


def collect_cloud(campaign_root: Path, bucket: str, profile: str, region: str, *, probe_workers: bool = True) -> CloudSnapshot:
    """List durable receipts and read current worker DBs with bounded SSM commands."""
    raw = json.loads((campaign_root / "aws_controller.json").read_text())
    controller = Controller(raw["cohort"], raw["observed_at"], tuple(Job(j["job_id"], j.get("instance_id"), j["status"]) for j in raw["jobs"]))
    if not controller.cohort or "/" in controller.cohort:
        raise StatusError("Unsafe cohort identifier")
    session = boto3.Session(profile_name=profile, region_name=region)
    config = Config(connect_timeout=10, read_timeout=30, retries={"mode": "standard", "max_attempts": 4})
    ec2 = session.client("ec2", config=config)
    s3 = session.client("s3", config=config)
    ssm = session.client("ssm", config=config)
    instances = tuple(i["InstanceId"] for page in ec2.get_paginator("describe_instances").paginate(
        Filters=[{"Name": "tag:cohort", "Values": [controller.cohort]}, {"Name": "instance-state-name", "Values": ["running"]}])
        for reservation in page["Reservations"] for i in reservation["Instances"])
    known = {j.instance_id for j in controller.jobs}
    if set(instances) - known:
        raise StatusError("Running cohort instances have no controller job binding")
    assignments: dict[str, frozenset[str]] = {}
    for job in controller.jobs:
        if job.instance_id in instances:
            manifest = campaign_root / "jobs" / job.job_id / "manifest.jsonl"
            if manifest.parent.parent.resolve() != (campaign_root / "jobs").resolve():
                raise StatusError("Unsafe job manifest path")
            tasks = [json.loads(line)["task_id"] for line in manifest.read_text().splitlines() if line.strip()]
            if len(set(tasks)) != len(tasks):
                raise StatusError("Duplicate task within worker assignment")
            assignments[job.instance_id] = frozenset(tasks)
    assigned_list = [t for tasks in assignments.values() for t in tasks]
    if len(set(assigned_list)) != len(assigned_list):
        raise StatusError("A task is assigned to multiple live instances")
    diagnostics: list[str] = []
    probes: dict[str, str] = {}

    def inspect(instance: str) -> tuple[str, str]:
        if not probe_workers:
            return instance, ""
        command = "python3 -c " + shlex.quote(probe_script(assignments[instance]))
        try:
            response = ssm.send_command(InstanceIds=[instance], DocumentName="AWS-RunShellScript",
                                        Parameters={"commands": [command], "executionTimeout": ["60"]}, TimeoutSeconds=60)
        except ssm.exceptions.InvalidInstanceId as exc:
            return instance, f"ERROR worker unavailable for SSM: {exc}"
        command_id = response["Command"]["CommandId"]
        deadline = time.monotonic() + 100
        while time.monotonic() < deadline:
            invocations = ssm.list_command_invocations(CommandId=command_id, InstanceId=instance, Details=True)["CommandInvocations"]
            if invocations:
                invocation = invocations[0]
                status = invocation["Status"]
                if status == "Success":
                    output = ssm.get_command_invocation(CommandId=command_id, InstanceId=instance)
                    return instance, output["StandardOutputContent"]
                if status not in {"Pending", "InProgress", "Delayed"}:
                    return instance, f"ERROR SSM command {command_id}: {status}"
            time.sleep(1)
        return instance, f"ERROR SSM command {command_id}: observation timed out"

    records = []
    with ThreadPoolExecutor(max_workers=6) as executor:
        for instance, output in executor.map(inspect, instances):
            probes[instance] = output
            if not output or output.startswith("ERROR"):
                diagnostics.append(f"{instance}: {output or 'worker probe disabled'}")
                continue
            probe = parse_probe(output)
            if not probe.databases:
                diagnostics.append(f"{instance}: no worker database available")
            records.extend(Observation(row.task_id, row.state, instance) for row in probe.rows)
    # Receipts are read after telemetry so newly locked tasks take precedence.
    prefix = f"locked-quant-v1/{controller.cohort}/reference-bound-receipts/"
    locked = set()
    for page in s3.get_paginator("list_objects_v2").paginate(Bucket=bucket, Prefix=prefix):
        for item in page.get("Contents", []):
            relative = item["Key"].removeprefix(prefix)
            if relative.count("/") != 1 or not relative.endswith(".json"):
                raise StatusError(f"Unexpected receipt key: {relative}")
            locked.add(relative[:-5])
    return CloudSnapshot(frozenset(locked), frozenset(assigned_list), tuple(records), tuple(diagnostics), probes, instances, datetime.now(UTC).isoformat())
