"""Budget-bound, resumable EC2 controller for immutable quant receipts."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from metainformant.rna.engine.acquisition_allocation import aws_allocation_ids
from metainformant.rna.engine.acquisition_aws_policy import WorkerImage, validate_worker_image
from metainformant.rna.engine.aws_fleet import in_flight_tasks, reserved_future_charge
from metainformant.rna.engine.aws_inputs import (
    _inputs_bundle as _inputs_bundle,
    _render_startup as _render_startup,
    _source_bundle as _source_bundle,
)
from metainformant.rna.engine.aws_resources import (
    AccountedJob,
    CampaignBilling,
    WorkerPrices,
    campaign_charge,
)
from metainformant.rna.engine.campaign_status import load_inventory
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
    selected: list[dict[str, Any]] = []
    size = 0
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
    inventory: dict[str, Any],
    store: Any,
    cohort: str,
    destination: Path,
    *,
    config_dir: Path | None = None,
) -> dict[str, Any]:
    """Require a non-empty complete cohort and validate every restored sample."""
    tasks = [(s, t) for s in inventory["species"] for t in s["tasks"]]
    if not tasks or len(tasks) != inventory["task_count"] or len({t["task_id"] for _, t in tasks}) != len(tasks):
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
        verified_config = None
        if config_dir is not None:
            name = species.get("config_name")
            if not isinstance(name, str) or Path(name).name != name or name in ("", ".", "..") or "\\" in name:
                raise ValueError("unsafe frozen configuration filename")
            verified_config = config_dir / name
        target = destination / species["species"] / "work" / "quant" / task["accession"]
        receipt = restore_quantification(
            store,
            cohort,
            species["species"],
            task["accession"],
            target,
            expected_config_sha256=species["config_sha256"],
            expected_reference_index_sha256=species["index_sha256"],
            verified_config_path=verified_config,
        )
        if receipt["config_sha256"] != species["config_sha256"]:
            raise ValueError(f"completed sample configuration mismatch: {task['task_id']}")
        verified.append({"task_id": task["task_id"], "contract_id": receipt["contract_id"]})
    certificate = {
        "schema": (
            "metainformant.rna.quant_completion.v1"
            if inventory.get("schema") == "metainformant.rna.acquisition_inventory.v1"
            else "metainformant.hymenoptera.quant_completion.v1"
        ),
        "cohort": cohort,
        "verified_at": datetime.now(UTC).isoformat(),
        "inventory_sha256": hashlib.sha256(json.dumps(inventory, sort_keys=True).encode()).hexdigest(),
        "eligible_tasks": len(tasks),
        "verified_tasks": len(verified),
        "all_quant_locked": True,
        "publication_promoted": False,
        "samples": verified,
    }
    _write_json(destination / "quant_completion_certificate.json", certificate)
    return certificate


def _wait_for_fleet(state: dict[str, Any], state_path: Path, args: argparse.Namespace, status: str) -> None:
    """Persist a capacity/task/budget wait without abandoning live workers."""
    state["status"] = status
    _write_json(state_path, state)
    count = sum(job["status"] in {"running", "terminating"} for job in state["jobs"])
    print(
        f"locked={state['locked_count']}/{state['eligible_count']} "
        f"upper-spend=${state['spent_upper_bound']:.2f} workers={count}/{args.max_workers} "
        f"reserved=${state['reserved_future_upper_bound']:.2f} status={status}",
        flush=True,
    )
    if not args.once:
        time.sleep(args.poll_seconds)


def species_order(species: list[dict[str, Any]], priority: str = "", last: Any = ()) -> list[dict[str, Any]]:
    """Order species for partitioning: the priority species first, `last` species at the end, else lexical."""
    deferred = set(last)
    return sorted(species, key=lambda s: (s["species"] in deferred, s["species"] != priority, s["species"]))


def _run_controller_locked(args: argparse.Namespace, owned_lock: Any) -> dict[str, Any]:
    """Reconcile receipts and own a disjoint, fully reserved worker fleet."""
    import boto3
    from botocore.config import Config

    root, repo = args.campaign_root.resolve(), args.repo.resolve()
    root.mkdir(parents=True, exist_ok=True)
    inventory_bytes = (root / "inventory.json").read_bytes()
    inventory_ids = load_inventory(inventory_bytes).task_ids()
    allocation_path = getattr(args, "task_allocation", None)
    aws_ids = (
        aws_allocation_ids(allocation_path, inventory_ids, hashlib.sha256(inventory_bytes).hexdigest())
        if allocation_path
        else inventory_ids
    )
    config_dir = getattr(args, "config_dir", None) or repo / "projects/hymenoptera_amalgkit/config/amalgkit"
    generic_worker = getattr(args, "config_dir", None) is not None
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
        name = species["config_name"]
        if not isinstance(name, str) or Path(name).name != name or "\\" in name or not name.endswith(".yaml"):
            raise ValueError("unsafe frozen configuration name")
        config_path = config_dir / name
        if hashlib.sha256(config_path.read_bytes()).hexdigest() != species["config_sha256"]:
            raise ValueError(f"frozen species configuration changed: {species['species']}")
    state_path = root / "aws_controller.json"
    session = boto3.Session(profile_name=args.profile, region_name=args.region)
    from metainformant.rna.engine.acquisition_pricing import quote_aws_worker

    quote = quote_aws_worker(
        region=args.region,
        instance_type=args.instance_type,
        disk_gib=getattr(args, "min_disk_gib", 600),
        profile=args.profile,
        hourly_floor=args.hourly_upper_bound,
        disk_throughput_mibps=args.disk_throughput_mibps,
    )
    worker_prices = WorkerPrices(
        quote.compute_hourly_usd,
        quote.gp3_gib_month_usd,
        args.hourly_upper_bound,
        gp3_mibps_month=quote.gp3_mibps_month_usd,
    )
    ec2 = session.client(
        "ec2",
        config=Config(
            connect_timeout=10,
            read_timeout=30,
            retries={"mode": "standard", "max_attempts": 5},
        ),
    )
    if generic_worker:
        images = ec2.describe_images(ImageIds=[args.ami])["Images"]
        if len(images) != 1:
            raise ValueError("generic acquisition requires exactly one available worker image")
        image = images[0]
        validate_worker_image(
            WorkerImage(
                image.get("State", ""),
                image.get("Architecture", ""),
                image.get("PlatformDetails", ""),
                bool(image.get("ProductCodes")),
            ),
            custom_template=getattr(args, "startup_template", None) is not None,
        )
    store = S3Store(args.bucket, "locked-quant-v1", profile=args.profile, region=args.region)
    objects = S3Store(args.bucket, "completion-v1", profile=args.profile, region=args.region)
    template = getattr(args, "startup_template", None) or (
        repo / "scripts/rna/aws_acquisition_startup.sh"
        if generic_worker
        else repo / "projects/hymenoptera_amalgkit/scripts/cloud/aws_startup.sh"
    )
    source = root / "source.tar"
    source_sha = _source_bundle(repo, source)
    source_key = f"{args.cohort}/sources/{source_sha}.tar"
    objects.put(source_key, source.read_bytes())
    source_bound = "completion-v1/" + source_key
    state = (
        json.loads(state_path.read_text())
        if state_path.exists()
        else {
            "schema": (
                "metainformant.rna.aws_acquisition.v1"
                if generic_worker
                else "metainformant.hymenoptera.aws_completion.v1"
            ),
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
    allocation_sha = hashlib.sha256(allocation_path.read_bytes()).hexdigest() if allocation_path else None
    if state.get("allocation_sha256") not in {None, allocation_sha} or (
        state.get("allocation_sha256") and not allocation_path
    ):
        raise ValueError("persisted acquisition allocation cannot be changed or omitted on resume")
    if allocation_path:
        state["allocation_sha256"] = allocation_sha
    state["max_workers"] = args.max_workers
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
            active = [j for j in state["jobs"] if j["status"] in {"running", "terminating"}]
            for job in active:
                try:
                    response = ec2.describe_instances(InstanceIds=[job["instance_id"]])
                except ec2.exceptions.ClientError as exc:
                    if exc.response["Error"]["Code"] != "InvalidInstanceID.NotFound":
                        raise
                    response = {"Reservations": []}
                instances = [i for r in response["Reservations"] for i in r["Instances"]]
                if not instances or instances[0]["State"]["Name"] == "terminated":
                    job.update(status="terminated", finished_at=now)
                elif now >= job["deadline"] and job["status"] == "running":
                    ec2.terminate_instances(InstanceIds=[job["instance_id"]])
                    job.update(status="terminating", termination_requested_at=now)
            now = time.time()
            billing_jobs: list[AccountedJob] = [
                {
                    "started_at": job["started_at"],
                    "finished_at": job.get("finished_at", max(now, job["started_at"])),
                    "hourly_upper_bound": job.get("hourly_upper_bound", state["hourly_upper_bound"]),
                }
                for job in state["jobs"]
            ]
            billing: CampaignBilling = {
                "historical_gross": state["historical_gross"],
                "hourly_upper_bound": state["hourly_upper_bound"],
                "jobs": billing_jobs,
            }
            spent = campaign_charge(billing, now)
            state["spent_upper_bound"] = spent
            state["observed_at"] = datetime.now(UTC).isoformat()
            locked: set[str] = set()
            prefix = f"locked-quant-v1/{args.cohort}/reference-bound-receipts/"
            for page in store.client.get_paginator("list_objects_v2").paginate(Bucket=args.bucket, Prefix=prefix):
                for record in page.get("Contents", []):
                    relative = record["Key"].removeprefix(prefix)
                    pieces = relative.split("/")
                    if len(pieces) == 2 and pieces[1].endswith(".json"):
                        locked.add(f"{pieces[0]}/{pieces[1][:-5]}")
            all_ids = {t["task_id"] for s in inventory["species"] for t in s["tasks"]}
            if not locked.issubset(all_ids):
                raise ValueError("durable receipt inventory contains tasks outside the frozen cohort")
            state.update(
                locked_count=len(locked),
                eligible_count=len(all_ids),
                missing_count=len(all_ids - locked),
            )
            active = [j for j in state["jobs"] if j["status"] in {"running", "terminating"}]
            occupied = in_flight_tasks(state["jobs"])
            if allocation_path and occupied - aws_ids - locked:
                raise ValueError("new allocation conflicts with existing live AWS ownership")
            future_charge = reserved_future_charge(state["jobs"], now, state["hourly_upper_bound"])
            state["reserved_future_upper_bound"] = future_charge
            if not all_ids - locked and not active:
                state["status"] = "verifying_all_outputs"
                _write_json(state_path, state)
                verify_locked_campaign(
                    inventory,
                    store,
                    args.cohort,
                    root / "completed_quant",
                    config_dir=config_dir,
                )
                state["status"] = "all_quant_locked"
                _write_json(state_path, state)
                return state
            if not all_ids - locked:
                for job in active:
                    if job["status"] == "running":
                        ec2.terminate_instances(InstanceIds=[job["instance_id"]])
                        job.update(status="terminating", termination_requested_at=now)
                _wait_for_fleet(state, state_path, args, "draining_completed_workers")
                if args.once:
                    return state
                continue
            if len(active) >= args.max_workers:
                _wait_for_fleet(state, state_path, args, "processing")
                if args.once:
                    return state
                continue
            if not aws_ids - locked:
                _wait_for_fleet(state, state_path, args, "waiting_for_local_acquisition")
                if args.once:
                    return state
                continue
            attempts: dict[str, int] = {}
            for job in state["jobs"]:
                for task_id in job["task_ids"]:
                    attempts[task_id] = attempts.get(task_id, 0) + 1
            species_options = species_order(
                inventory["species"],
                priority=getattr(args, "priority_species", "nasonia_vitripennis"),
                last=getattr(args, "last_species", ()),
            )
            selected_species, partition = None, []
            for species in species_options:
                eligible = [
                    t
                    for t in species["tasks"]
                    if t["task_id"] in aws_ids and attempts.get(t["task_id"], 0) < args.max_attempts
                ]
                partition = choose_partition(
                    eligible,
                    locked | occupied,
                    max_bytes=getattr(args, "partition_bytes", 60 * 1024**3),
                    max_tasks=(
                        getattr(args, "partition_tasks", 120)
                        if state["jobs"]
                        else min(20, getattr(args, "partition_tasks", 120))
                    ),
                )
                if partition:
                    selected_species = species
                    break
            if selected_species is None:
                if active:
                    _wait_for_fleet(state, state_path, args, "waiting_for_reserved_tasks")
                    if args.once:
                        return state
                    continue
                state["status"] = "unresolved_tasks_require_source_or_size_diagnosis"
                _write_json(state_path, state)
                return state
            largest = max(int(t["fastq_bytes"]) for t in partition)
            raw_bound = max(getattr(args, "partition_bytes", 60 * 1024**3), largest)
            disk_gib = max(
                getattr(args, "min_disk_gib", 600),
                math.ceil(raw_bound * getattr(args, "disk_expansion_factor", 6) / 1024**3)
                + getattr(args, "disk_reserve_gib", 100),
            )
            if disk_gib > getattr(args, "max_disk_gib", 2000):
                raise ValueError("partition requires an independently reviewed disk profile")
            job_rate = worker_prices.hourly_bound(disk_gib, throughput_mibps=args.disk_throughput_mibps)
            limit_seconds = job_timeout(largest, args.job_seconds)
            reserved_seconds = limit_seconds + 900
            if not budget_allows(spent + future_charge, args.budget, reserved_seconds, job_rate):
                if active:
                    _wait_for_fleet(state, state_path, args, "waiting_for_budget_reservations")
                    if args.once:
                        return state
                    continue
                state["status"] = "budget_exhausted"
                _write_json(state_path, state)
                return state
            job_id = f"job-{len(state['jobs']) + 1:05d}"
            directory = root / "jobs" / job_id
            bundle, input_sha = _inputs_bundle(
                root,
                selected_species,
                partition,
                directory,
                config_path=config_dir / selected_species["config_name"] if generic_worker else None,
            )
            input_key = f"{args.cohort}/inputs/{input_sha}.tar"
            objects.put(input_key, bundle.read_bytes())
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
                    **(
                        {
                            "WORKERS": args.worker_workers,
                            "THREADS": args.worker_threads,
                            "QUANT_SLOTS": args.worker_quant_slots,
                            "FASTQ_SLOTS": args.worker_fastq_slots,
                            "MAX_IN_FLIGHT": args.worker_max_in_flight,
                            "FASTQ_THREADS": args.worker_fastq_threads,
                            "COMPRESSION_THREADS": args.worker_compression_threads,
                            "COMPRESSION_LEVEL": args.worker_compression_level,
                            "VALIDATION_SLOTS": args.worker_validation_slots,
                        }
                        if generic_worker
                        else {}
                    ),
                },
            )
            (directory / "user_data.sh").write_text(script)
            task_ids = [t["task_id"] for t in partition]
            token = hashlib.sha256(f"{args.cohort}/{job_id}/{input_sha}/{source_sha}".encode()).hexdigest()
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
                            "Throughput": args.disk_throughput_mibps,
                            "Iops": 3000,
                            "Encrypted": True,
                            "DeleteOnTermination": True,
                        },
                    }
                ],
                "TagSpecifications": [
                    {
                        "ResourceType": "instance",
                        "Tags": [
                            {"Key": "Name", "Value": f"amalgkit-{args.cohort}-{job_id}"},
                            {"Key": "daf-cloud", "Value": "managed"},
                            {"Key": "cohort", "Value": args.cohort},
                            {"Key": "job-id", "Value": job_id},
                        ],
                    }
                ],
            }
            if args.instance_type.split(".")[0] in {"t2", "t3", "t3a", "t4g"}:
                request["CreditSpecification"] = {"CpuCredits": "standard"}
            job = {
                "job_id": job_id,
                "task_ids": task_ids,
                "status": "admitting",
                "client_token": token,
                "input_sha256": input_sha,
                "source_sha256": source_sha,
                "compression_level": args.worker_compression_level,
                "disk_gib": disk_gib,
                "disk_throughput_mibps": args.disk_throughput_mibps,
                "request": request,
                "reserved_seconds": reserved_seconds,
                "hourly_upper_bound": job_rate,
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
            if len(active) + 1 >= args.max_workers:
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
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--max-attempts", type=int, default=3)
    parser.add_argument(
        "--config-dir", type=Path, help="Use the generic worker with these frozen species configurations"
    )
    parser.add_argument("--startup-template", type=Path)
    parser.add_argument("--task-allocation", type=Path)
    parser.add_argument(
        "--priority-species",
        default="nasonia_vitripennis",
        help="Legacy priority; set an empty value for generic lexical order",
    )
    parser.add_argument(
        "--last-species",
        action="append",
        default=[],
        help="Schedule this species after all others (repeatable); overrides --priority-species",
    )
    parser.add_argument("--partition-bytes", type=int, default=60 * 1024**3)
    parser.add_argument("--partition-tasks", type=int, default=120)
    parser.add_argument("--min-disk-gib", type=int, default=600)
    parser.add_argument("--disk-throughput-mibps", type=int, default=125)
    parser.add_argument("--max-disk-gib", type=int, default=2000)
    parser.add_argument("--disk-expansion-factor", type=float, default=6)
    parser.add_argument("--disk-reserve-gib", type=int, default=100)
    for flag, default in (
        ("workers", 16),
        ("threads", 8),
        ("quant-slots", 4),
        ("fastq-slots", 1),
        ("max-in-flight", 12),
        ("fastq-threads", 2),
        ("compression-threads", 2),
        ("validation-slots", 4),
    ):
        parser.add_argument(f"--worker-{flag}", type=int, default=default)
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--worker-compression-level", type=int, choices=range(1, 10), default=6)
    args = parser.parse_args(argv)
    if not 125 <= args.disk_throughput_mibps <= 750:
        parser.error("disk throughput must be 125–750 MiB/s at baseline 3000 IOPS")
    if args.worker_compression_level != 6 and args.config_dir is None:
        parser.error("nondefault compression levels require the generic worker (--config-dir)")
    if (
        min(args.partition_bytes, args.partition_tasks, args.min_disk_gib, args.max_disk_gib, args.disk_reserve_gib)
        <= 0
        or not math.isfinite(args.disk_expansion_factor)
        or args.disk_expansion_factor < 1
        or args.min_disk_gib > args.max_disk_gib
        or args.max_disk_gib > 16384
    ):
        parser.error("invalid partition/disk bounds; this controller caps configured volumes at 16384 GiB")
    if (
        min(args.job_seconds, args.poll_seconds, args.max_attempts, args.max_workers) <= 0
        or min(
            args.worker_workers,
            args.worker_threads,
            args.worker_quant_slots,
            args.worker_fastq_slots,
            args.worker_max_in_flight,
            args.worker_fastq_threads,
            args.worker_compression_threads,
            args.worker_validation_slots,
        )
        <= 0
    ):
        parser.error("job, polling, attempt and worker bounds must be positive")
    print(json.dumps(run_controller(args), indent=2))
    return 0
