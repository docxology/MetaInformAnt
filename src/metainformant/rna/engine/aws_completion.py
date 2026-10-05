"""Budget-bound, resumable EC2 controller for immutable quant receipts."""

from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import math
import shlex
import tarfile
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from metainformant.rna.amalgkit import (
    AMALGKIT_RELEASE_TAG,
    AMALGKIT_SOURCE_REVISION,
    REQUIRED_AMALGKIT_VERSION,
)
from metainformant.rna.engine.durable_quant import S3Store, restore_quantification
from metainformant.rna.engine.source_resolution import (
    SourceTarget,
    load_source_resolutions,
)


def budget_allows(
    spent: float,
    ceiling: float,
    reserved_seconds: int,
    hourly_upper_bound: float,
    reserve: float = 10,
) -> bool:
    """Admit only when the entire hard-bounded job fits, including a storage reserve."""
    values = (spent, ceiling, hourly_upper_bound, reserve)
    if any(not math.isfinite(v) or v < 0 for v in values) or reserved_seconds <= 0:
        raise ValueError("invalid budget inputs")
    return spent + reserved_seconds / 3600 * hourly_upper_bound + reserve <= ceiling


def job_timeout(largest_raw_bytes: int, minimum_seconds: int) -> int:
    """Allow the largest run a bounded 0.5 MiB/s transfer window plus quant time."""
    if largest_raw_bytes <= 0 or minimum_seconds <= 0:
        raise ValueError("task size and minimum deadline must be positive")
    return max(minimum_seconds, math.ceil(largest_raw_bytes / (512 * 1024)) + 3600)


def choose_partition(
    tasks: list[dict[str, Any]],
    locked: set[str],
    *,
    max_bytes: int = 60 * 1024**3,
    max_tasks: int = 120,
) -> list[dict[str, Any]]:
    """Choose one species' byte-bounded missing tasks; never silently admit unknown sizes."""
    selected, size = [], 0
    for task in sorted(tasks, key=lambda t: (int(t["fastq_bytes"]), t["accession"])):
        if task["task_id"] in locked:
            continue
        raw = int(task["fastq_bytes"])
        if raw <= 0:
            continue
        if selected and size + raw > max_bytes:
            continue
        if not selected and raw > max_bytes:
            return [task]
        selected.append(task)
        size += raw
        if len(selected) >= max_tasks:
            break
    return selected


def _write_json(path: Path, payload: Any) -> None:
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True))
    temporary.replace(path)


def verify_locked_campaign(
    inventory: dict[str, Any], store: Any, cohort: str, destination: Path
) -> dict[str, Any]:
    """Require a non-empty complete cohort and validate every restored sample."""
    tasks = [(s, t) for s in inventory["species"] for t in s["tasks"]]
    if (
        not tasks
        or len(tasks) != inventory["task_count"]
        or len({t["task_id"] for _, t in tasks}) != len(tasks)
    ):
        raise ValueError("empty or incomplete completion inventory")
    verified = []
    for species, task in tasks:
        index_hash = species["index_sha256"]
        if (
            not isinstance(index_hash, str)
            or len(index_hash) != 64
            or any(c not in "0123456789abcdef" for c in index_hash)
            or task.get("reference_index_sha256") != index_hash
        ):
            raise ValueError("task reference differs from frozen species index")
        target = destination / species["species"] / "work" / "quant" / task["accession"]
        receipt = restore_quantification(
            store,
            cohort,
            species["species"],
            task["accession"],
            target,
            expected_config_sha256=species["config_sha256"],
            expected_reference_index_sha256=species["index_sha256"],
        )
        if receipt["config_sha256"] != species["config_sha256"]:
            raise ValueError(
                f"completed sample configuration mismatch: {task['task_id']}"
            )
        verified.append(
            {"task_id": task["task_id"], "contract_id": receipt["contract_id"]}
        )
    certificate = {
        "schema": "metainformant.hymenoptera.quant_completion.v1",
        "cohort": cohort,
        "verified_at": datetime.now(UTC).isoformat(),
        "inventory_sha256": hashlib.sha256(
            json.dumps(inventory, sort_keys=True).encode()
        ).hexdigest(),
        "eligible_tasks": len(tasks),
        "verified_tasks": len(verified),
        "all_quant_locked": True,
        "publication_promoted": False,
        "samples": verified,
    }
    _write_json(destination / "quant_completion_certificate.json", certificate)
    return certificate


def _source_bundle(repo: Path, destination: Path) -> str:
    roots = [
        repo / "src",
        repo / "scripts/rna",
        repo / "projects/hymenoptera_amalgkit/scripts",
        repo / "projects/hymenoptera_amalgkit/config",
        repo / "config/amalgkit",
    ]
    files = [repo / name for name in ("pyproject.toml", "uv.lock", "README.md")]
    for root in roots:
        files.extend(
            p
            for p in root.rglob("*")
            if p.is_file()
            and not p.is_symlink()
            and "__pycache__" not in p.parts
            and not p.name.startswith("._")
            and p.suffix != ".pyc"
        )
    with tarfile.open(destination, "w") as archive:
        for path in sorted(set(files)):
            archive.add(
                path, arcname=path.relative_to(repo).as_posix(), recursive=False
            )
    with tarfile.open(destination) as archive:
        for member in archive:
            if (
                member.islnk()
                or member.issym()
                or Path(member.name).is_absolute()
                or ".." in Path(member.name).parts
            ):
                raise ValueError("unsafe source bundle")
    return hashlib.sha256(destination.read_bytes()).hexdigest()


def _inputs_bundle(
    root: Path, species: dict[str, Any], tasks: list[dict[str, Any]], directory: Path
) -> tuple[Path, str]:
    directory.mkdir(parents=True, exist_ok=True)
    manifest = directory / "manifest.jsonl"
    manifest.write_text("".join(json.dumps(t, sort_keys=True) + "\n" for t in tasks))
    source = root / "inputs" / species["species"] / "work"
    metadata = source / "metadata" / "metadata_selected.tsv"
    index = source / "index" / species["index_name"]
    if hashlib.sha256(metadata.read_bytes()).hexdigest() != species["metadata_sha256"]:
        raise ValueError("frozen selected metadata changed")
    if hashlib.sha256(index.read_bytes()).hexdigest() != species["index_sha256"]:
        raise ValueError("frozen reference index changed")
    resolved = {t["accession"]: t for t in tasks if t.get("source_evidence_sha256")}
    metadata_override = None
    if resolved:
        with metadata.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields, rows = list(reader.fieldnames or []), list(reader)
        additions = (
            "total_spots",
            "total_bases",
            "size",
            "mean_read_length",
            "source_evidence_sha256",
        )
        fields.extend(name for name in additions if name not in fields)
        seen = set()
        for row in rows:
            task = resolved.get(row["run"])
            if task is not None:
                seen.add(row["run"])
                row.update(
                    total_spots=str(task["total_spots"]),
                    total_bases=str(task["total_bases"]),
                    size=str(task["sra_bytes"]),
                    mean_read_length=str(task["total_bases"] / task["total_spots"]),
                    source_evidence_sha256=task["source_evidence_sha256"],
                )
        if seen != set(resolved):
            raise ValueError("resolved tasks missing from frozen selected metadata")
        metadata_override = directory / "metadata_resolved.tsv"
        with metadata_override.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    records = []
    for path in sorted(source.rglob("*")):
        if not path.is_file():
            continue
        payload_path = (
            metadata_override
            if path == metadata and metadata_override is not None
            else path
        )
        records.append(
            {
                "path": f"data/{species['species']}/work/{path.relative_to(source).as_posix()}",
                "sha256": hashlib.sha256(payload_path.read_bytes()).hexdigest(),
                "size": payload_path.stat().st_size,
            }
        )
    snapshot = {
        "schema": "metainformant.hymenoptera.gcp_snapshot.v1",
        "source_state": "quiescent",
        "cloud_launch_policy": "checkpointed",
        "amalgkit_version": REQUIRED_AMALGKIT_VERSION,
        "amalgkit_release_tag": AMALGKIT_RELEASE_TAG,
        "amalgkit_source_revision": AMALGKIT_SOURCE_REVISION,
        "manifest_sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
        "task_count": len(tasks),
        "species": [species["species"]],
        "input_files": records,
        "raw_reads_included": False,
        "quant_outputs_included": False,
        "source_resolution_tasks": [
            {"task_id": t["task_id"], "evidence_sha256": t["source_evidence_sha256"]}
            for t in resolved.values()
        ],
        "frozen_metadata_sha256": species["metadata_sha256"],
    }
    (directory / "snapshot.json").write_text(
        json.dumps(snapshot, indent=2, sort_keys=True)
    )
    bundle = directory / "inputs.tar"
    with tarfile.open(bundle, "w") as archive:
        archive.add(manifest, arcname="manifest.jsonl")
        archive.add(directory / "snapshot.json", arcname="snapshot.json")
        for record in records:
            relative = Path(record["path"]).relative_to(
                f"data/{species['species']}/work"
            )
            payload_path = source / relative
            if payload_path == metadata and metadata_override is not None:
                payload_path = metadata_override
            archive.add(payload_path, arcname=record["path"], recursive=False)
    return bundle, hashlib.sha256(bundle.read_bytes()).hexdigest()


def _render_startup(template: Path, replacements: dict[str, Any]) -> str:
    script = template.read_text()
    for key, value in replacements.items():
        script = script.replace(f"@@{key}@@", shlex.quote(str(value)))
    if "@@" in script or len(script.encode()) > 16384:
        raise ValueError("invalid or oversized EC2 user data")
    return script


def _run_controller_locked(args: argparse.Namespace, owned_lock: Any) -> dict[str, Any]:
    """Reconcile AWS receipts and own exactly one bounded instance at a time."""
    import boto3
    from botocore.config import Config

    root, repo = args.campaign_root.resolve(), args.repo.resolve()
    root.mkdir(parents=True, exist_ok=True)
    inventory_bytes = (root / "inventory.json").read_bytes()
    inventory = json.loads(inventory_bytes)
    targets = [
        SourceTarget(t["accession"], sp["species"], int(sp["taxid"]))
        for sp in inventory["species"]
        for t in sp["tasks"]
    ]
    resolutions = load_source_resolutions(
        root / "source_resolutions.json",
        hashlib.sha256(inventory_bytes).hexdigest(),
        targets,
    )
    by_accession = {r.accession: r for r in resolutions}
    for species in inventory["species"]:
        for task in species["tasks"]:
            resolution = by_accession.get(task["accession"])
            if resolution is not None:
                task.update(
                    fastq_bytes=resolution.raw_bound_bytes,
                    total_spots=resolution.total_spots,
                    total_bases=resolution.total_bases,
                    sra_bytes=resolution.sra_bytes,
                    source_evidence_sha256=resolution.evidence_sha256,
                    resource_bound_basis="sra_bytes+3*total_bases+512*total_spots",
                )
    inventory["source_resolutions"] = [r.accession for r in resolutions]
    for species in inventory["species"]:
        config_path = (
            repo
            / "projects/hymenoptera_amalgkit/config/amalgkit"
            / species["config_name"]
        )
        if (
            hashlib.sha256(config_path.read_bytes()).hexdigest()
            != species["config_sha256"]
        ):
            raise ValueError(
                f"frozen species configuration changed: {species['species']}"
            )
    state_path = root / "aws_controller.json"
    session = boto3.Session(profile_name=args.profile, region_name=args.region)
    if args.region != "us-east-2":
        raise ValueError("the completion pricing guard currently requires us-east-2")
    pricing = session.client("pricing", region_name="us-east-1")
    prices = pricing.get_products(
        ServiceCode="AmazonEC2",
        Filters=[
            {"Type": "TERM_MATCH", "Field": key, "Value": value}
            for key, value in {
                "instanceType": args.instance_type,
                "location": "US East (Ohio)",
                "operatingSystem": "Linux",
                "tenancy": "Shared",
                "preInstalledSw": "NA",
                "capacitystatus": "Used",
            }.items()
        ],
        MaxResults=10,
    )
    hourly_prices = [
        float(d["pricePerUnit"]["USD"])
        for product in prices["PriceList"]
        for term in json.loads(product)["terms"]["OnDemand"].values()
        for d in term["priceDimensions"].values()
        if d["unit"] == "Hrs"
    ]
    if len(hourly_prices) != 1 or args.hourly_upper_bound <= hourly_prices[0] + 0.10:
        raise ValueError(
            "hourly cost guard must exceed current compute plus storage/IP allowance"
        )
    ec2 = session.client(
        "ec2",
        config=Config(
            connect_timeout=10,
            read_timeout=30,
            retries={"mode": "standard", "max_attempts": 5},
        ),
    )
    store = S3Store(
        args.bucket, "locked-quant-v1", profile=args.profile, region=args.region
    )
    objects = S3Store(
        args.bucket, "completion-v1", profile=args.profile, region=args.region
    )
    template = repo / "projects/hymenoptera_amalgkit/scripts/cloud/aws_startup.sh"
    source = root / "source.tar"
    source_sha = _source_bundle(repo, source)
    source_key = f"{args.cohort}/sources/{source_sha}.tar"
    objects.put(source_key, source.read_bytes())
    source_bound = "completion-v1/" + source_key
    state = (
        json.loads(state_path.read_text())
        if state_path.exists()
        else {
            "schema": "metainformant.hymenoptera.aws_completion.v1",
            "cohort": args.cohort,
            "budget_ceiling": args.budget,
            "historical_gross": args.historical_gross,
            "hourly_upper_bound": args.hourly_upper_bound,
            "jobs": [],
            "status": "prepared",
        }
    )
    if state["cohort"] != args.cohort or args.budget < state["budget_ceiling"]:
        raise ValueError("controller scope/budget transition is inconsistent")
    state["budget_ceiling"] = args.budget
    try:
        while True:
            now = time.time()
            for job in state["jobs"]:
                if job["status"] == "admitting":
                    response = ec2.run_instances(**job["request"])
                    instance = response["Instances"][0]
                    started = instance["LaunchTime"].timestamp()
                    job.update(
                        instance_id=instance["InstanceId"],
                        started_at=started,
                        deadline=started + job["reserved_seconds"],
                        status="running",
                    )
                    _write_json(state_path, state)
            active = [j for j in state["jobs"] if j["status"] == "running"]
            for job in active:
                try:
                    response = ec2.describe_instances(InstanceIds=[job["instance_id"]])
                except ec2.exceptions.ClientError as exc:
                    if exc.response["Error"]["Code"] != "InvalidInstanceID.NotFound":
                        raise
                    response = {"Reservations": []}
                instances = [
                    i for r in response["Reservations"] for i in r["Instances"]
                ]
                if not instances or instances[0]["State"]["Name"] == "terminated":
                    job.update(status="terminated", finished_at=now)
                elif now >= job["deadline"]:
                    ec2.terminate_instances(InstanceIds=[job["instance_id"]])
                    job.update(status="terminated_by_deadline", finished_at=now)
            spent = state["historical_gross"] + sum(
                (j.get("finished_at", now) - j["started_at"])
                / 3600
                * state["hourly_upper_bound"]
                for j in state["jobs"]
            )
            state["spent_upper_bound"] = spent
            state["observed_at"] = datetime.now(UTC).isoformat()
            locked: set[str] = set()
            prefix = f"locked-quant-v1/{args.cohort}/reference-bound-receipts/"
            for page in store.client.get_paginator("list_objects_v2").paginate(
                Bucket=args.bucket, Prefix=prefix
            ):
                for record in page.get("Contents", []):
                    relative = record["Key"].removeprefix(prefix)
                    pieces = relative.split("/")
                    if len(pieces) == 2 and pieces[1].endswith(".json"):
                        locked.add(f"{pieces[0]}/{pieces[1][:-5]}")
            all_ids = {t["task_id"] for s in inventory["species"] for t in s["tasks"]}
            if not locked.issubset(all_ids):
                raise ValueError(
                    "durable receipt inventory contains tasks outside the frozen cohort"
                )
            state.update(
                locked_count=len(locked),
                eligible_count=len(all_ids),
                missing_count=len(all_ids - locked),
            )
            active = [j for j in state["jobs"] if j["status"] == "running"]
            if not all_ids - locked:
                state["status"] = "verifying_all_outputs"
                _write_json(state_path, state)
                verify_locked_campaign(
                    inventory, store, args.cohort, root / "completed_quant"
                )
                state["status"] = "all_quant_locked"
                _write_json(state_path, state)
                return state
            if active:
                state["status"] = "processing"
                _write_json(state_path, state)
                print(
                    f"locked={len(locked)}/{len(all_ids)} upper-spend=${spent:.2f} instance={active[0]['instance_id']}",
                    flush=True,
                )
                if args.once:
                    return state
                time.sleep(args.poll_seconds)
                continue
            attempts: dict[str, int] = {}
            for job in state["jobs"]:
                for task_id in job["task_ids"]:
                    attempts[task_id] = attempts.get(task_id, 0) + 1
            species_options = sorted(
                inventory["species"],
                key=lambda s: (s["species"] != "nasonia_vitripennis", s["species"]),
            )
            selected_species, partition = None, []
            for species in species_options:
                eligible = [
                    t
                    for t in species["tasks"]
                    if attempts.get(t["task_id"], 0) < args.max_attempts
                ]
                partition = choose_partition(
                    eligible, locked, max_tasks=20 if not state["jobs"] else 120
                )
                if partition:
                    selected_species = species
                    break
            if selected_species is None:
                state["status"] = "unresolved_tasks_require_source_or_size_diagnosis"
                _write_json(state_path, state)
                return state
            largest = max(int(t["fastq_bytes"]) for t in partition)
            limit_seconds = job_timeout(largest, args.job_seconds)
            reserved_seconds = limit_seconds + 900
            if not budget_allows(
                spent, args.budget, reserved_seconds, args.hourly_upper_bound
            ):
                state["status"] = "budget_exhausted"
                _write_json(state_path, state)
                return state
            job_id = f"job-{len(state['jobs']) + 1:05d}"
            directory = root / "jobs" / job_id
            bundle, input_sha = _inputs_bundle(
                root, selected_species, partition, directory
            )
            input_key = f"{args.cohort}/inputs/{input_sha}.tar"
            objects.put(input_key, bundle.read_bytes())
            raw_bound = max(60 * 1024**3, largest)
            disk_gib = max(600, math.ceil(raw_bound * 6 / 1024**3) + 100)
            if disk_gib > 2000:
                raise ValueError(
                    "partition requires an independently reviewed disk profile"
                )
            prefix = f"completion-v1/{args.cohort}/jobs/{job_id}"
            script = _render_startup(
                template,
                {
                    "BUCKET": args.bucket,
                    "REGION": args.region,
                    "COHORT": args.cohort,
                    "SPECIES": selected_species["species"],
                    "SOURCE_KEY": source_bound,
                    "SOURCE_SHA": source_sha,
                    "INPUT_KEY": "completion-v1/" + input_key,
                    "INPUT_SHA": input_sha,
                    "JOB_PREFIX": prefix,
                    "LIMIT_SECONDS": limit_seconds,
                    "RAW_BYTES": raw_bound,
                },
            )
            (directory / "user_data.sh").write_text(script)
            task_ids = [t["task_id"] for t in partition]
            token = hashlib.sha256(
                f"{args.cohort}/{job_id}/{input_sha}/{source_sha}".encode()
            ).hexdigest()
            request = {
                "ImageId": args.ami,
                "InstanceType": args.instance_type,
                "MinCount": 1,
                "MaxCount": 1,
                "ClientToken": token,
                "UserData": script,
                "IamInstanceProfile": {"Name": args.instance_profile},
                "InstanceInitiatedShutdownBehavior": "terminate",
                "MetadataOptions": {"HttpTokens": "required"},
                "BlockDeviceMappings": [
                    {
                        "DeviceName": "/dev/xvda",
                        "Ebs": {
                            "VolumeSize": disk_gib,
                            "VolumeType": "gp3",
                            "Encrypted": True,
                            "DeleteOnTermination": True,
                        },
                    }
                ],
                "TagSpecifications": [
                    {
                        "ResourceType": "instance",
                        "Tags": [
                            {"Key": "Name", "Value": f"hym-lock-{job_id}"},
                            {"Key": "daf-cloud", "Value": "managed"},
                            {"Key": "cohort", "Value": args.cohort},
                            {"Key": "job-id", "Value": job_id},
                        ],
                    }
                ],
            }
            job = {
                "job_id": job_id,
                "task_ids": task_ids,
                "status": "admitting",
                "client_token": token,
                "input_sha256": input_sha,
                "source_sha256": source_sha,
                "disk_gib": disk_gib,
                "request": request,
                "reserved_seconds": reserved_seconds,
            }
            state["jobs"].append(job)
            _write_json(state_path, state)
            response = ec2.run_instances(**request)
            instance = response["Instances"][0]
            started = instance["LaunchTime"].timestamp()
            job.update(
                instance_id=instance["InstanceId"],
                started_at=started,
                deadline=started + reserved_seconds,
                status="running",
            )
            state["status"] = "processing"
            _write_json(state_path, state)
            print(
                f"Launched {job_id} {instance['InstanceId']}: {len(partition)} missing tasks; disk={disk_gib} GiB",
                flush=True,
            )
            if args.once:
                return state
            time.sleep(args.poll_seconds)
    finally:
        owned_lock.close()


def run_controller(args: argparse.Namespace) -> dict[str, Any]:
    """Acquire ownership before constructing or publishing any run artifact."""
    root = args.campaign_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    with (root / "aws_controller.lock").open("a") as owned_lock:
        fcntl.flock(owned_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return _run_controller_locked(args, owned_lock)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-root", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--bucket", required=True)
    parser.add_argument("--cohort", required=True)
    parser.add_argument("--budget", type=float, required=True)
    parser.add_argument("--historical-gross", type=float, required=True)
    parser.add_argument("--hourly-upper-bound", type=float, default=0.55)
    parser.add_argument("--profile", default="dev-agent")
    parser.add_argument("--region", default="us-east-2")
    parser.add_argument("--ami", required=True)
    parser.add_argument("--instance-profile", required=True)
    parser.add_argument("--instance-type", default="c7i.2xlarge")
    parser.add_argument("--job-seconds", type=int, default=14400)
    parser.add_argument("--poll-seconds", type=int, default=60)
    parser.add_argument("--max-attempts", type=int, default=3)
    parser.add_argument("--once", action="store_true")
    args = parser.parse_args(argv)
    if args.job_seconds <= 0 or args.poll_seconds <= 0 or args.max_attempts <= 0:
        parser.error("job, polling and attempt bounds must be positive")
    print(json.dumps(run_controller(args), indent=2))
    return 0
