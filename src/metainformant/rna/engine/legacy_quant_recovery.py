"""Recover complete quant members from indexed, truncated S3 tar objects."""

from __future__ import annotations

import concurrent.futures
import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Any

from metainformant.rna.engine.durable_quant import (
    DirectoryStore,
    S3Store,
    lock_quantification,
    restore_quantification,
)


def recover_indexed_archives(
    review_dir: Path,
    output_dir: Path,
    config_dir: Path,
    *,
    bucket: str,
    cohort: str,
    profile: str,
    region: str,
    expected_reference_index_sha256: str,
    workers: int = 6,
) -> dict[str, Any]:
    """Use bounded immutable-object reads; failures remain explicit and resumable."""
    import boto3
    from botocore.config import Config

    output_dir.mkdir(parents=True, exist_ok=True)
    local = DirectoryStore(output_dir / "locked")
    remote = S3Store(bucket, "locked-quant-v1", profile=profile, region=region)
    client = boto3.Session(profile_name=profile, region_name=region).client(
        "s3",
        config=Config(
            connect_timeout=10,
            read_timeout=60,
            retries={"mode": "standard", "max_attempts": 5},
        ),
    )
    sources: dict[str, list[dict[str, Any]]] = {}
    archives = json.loads((review_dir / "archive_integrity_summary.json").read_text())
    for archive in sorted(archives, key=lambda a: a["key"]):
        key = archive["key"]
        instance = key.split("/")[1]
        head = client.head_object(Bucket=bucket, Key=key)
        if head["ContentLength"] != archive["size"]:
            raise ValueError(f"archive changed since inspection: {key}")
        index = json.loads((review_dir / f"{instance}_archive_members.json").read_text())
        by_accession: dict[str, list[dict[str, Any]]] = {}
        for member in index:
            if "/quant/" not in member["name"] or member["size"] <= 0:
                continue
            _, relative = member["name"].split("/quant/", 1)
            pieces = relative.split("/")
            if len(pieces) != 2 or member["offset_data"] + member["size"] > head["ContentLength"]:
                continue
            by_accession.setdefault(pieces[0], []).append({**member, "filename": pieces[1]})
        for accession in archive["quant_ids"]:
            sources.setdefault(accession, []).append(
                {"key": key, "etag": head["ETag"], "members": by_accession[accession]}
            )
    journal = output_dir / "recovery_results.jsonl"
    config = config_dir / "amalgkit_nasonia_vitripennis.yaml"
    expected_config = hashlib.sha256(config.read_bytes()).hexdigest()

    def recover(accession: str) -> dict[str, Any]:
        with tempfile.TemporaryDirectory(prefix=f"retry-{accession}-", dir=output_dir) as temporary:
            sample = Path(temporary) / accession
            try:
                receipt = restore_quantification(
                    local,
                    cohort,
                    "nasonia_vitripennis",
                    accession,
                    sample,
                    expected_config_sha256=expected_config,
                    expected_reference_index_sha256=expected_reference_index_sha256,
                )
            except FileNotFoundError:
                pass
            else:
                lock_quantification(
                    remote,
                    cohort,
                    "nasonia_vitripennis",
                    accession,
                    sample,
                    expected_config_sha256=expected_config,
                    expected_reference_index_sha256=expected_reference_index_sha256,
                )
                return {
                    "accession": accession,
                    "status": "locked",
                    "source": "local_verified_retry",
                    "contract_id": receipt["contract_id"],
                }
        errors = []
        for source in sources[accession]:
            try:
                with tempfile.TemporaryDirectory(prefix=f"{accession}-", dir=output_dir) as temporary:
                    sample = Path(temporary)
                    for member in source["members"]:
                        if member["size"] > 128 * 1024**2:
                            raise ValueError("quant member exceeds recovery bound")
                        start = member["offset_data"]
                        response = client.get_object(
                            Bucket=bucket,
                            Key=source["key"],
                            IfMatch=source["etag"],
                            Range=f"bytes={start}-{start + member['size'] - 1}",
                        )
                        with response["Body"] as body:
                            data = body.read()
                        if len(data) != member["size"]:
                            raise ValueError("incomplete quant member range")
                        name = member["filename"]
                        if Path(name).name != name:
                            raise ValueError("unsafe quant member name")
                        (sample / name).write_bytes(data)
                    receipt = lock_quantification(
                        local,
                        cohort,
                        "nasonia_vitripennis",
                        accession,
                        sample,
                        expected_config_sha256=expected_config,
                        expected_reference_index_sha256=expected_reference_index_sha256,
                    )
                    lock_quantification(
                        remote,
                        cohort,
                        "nasonia_vitripennis",
                        accession,
                        sample,
                        expected_config_sha256=expected_config,
                        expected_reference_index_sha256=expected_reference_index_sha256,
                    )
                    return {
                        "accession": accession,
                        "status": "locked",
                        "source": source["key"],
                        "contract_id": receipt["contract_id"],
                    }
            except (OSError, ValueError, KeyError) as exc:
                errors.append({"source": source["key"], "error": str(exc)})
        return {"accession": accession, "status": "rejected", "errors": errors}

    previous = {}
    if journal.exists():
        for line in journal.read_text().splitlines():
            record = json.loads(line)
            previous[record["accession"]] = record
    results = [
        r
        for accession, r in previous.items()
        if accession in sources
        and r["status"] == "locked"
        and (
            output_dir / "locked" / cohort / "reference-bound-receipts" / "nasonia_vitripennis" / f"{accession}.json"
        ).is_file()
    ]
    remaining = [accession for accession in sorted(sources) if accession not in {r["accession"] for r in results}]
    with (
        concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor,
        journal.open("a") as handle,
    ):
        futures = {executor.submit(recover, accession): accession for accession in remaining}
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            results.append(result)
            handle.write(json.dumps(result, sort_keys=True) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
            if len(results) % 25 == 0:
                print(
                    f"Recovered {len(results)}/{len(sources)} candidates; "
                    f"locked={sum(r['status'] == 'locked' for r in results)}",
                    flush=True,
                )
    summary = {
        "cohort": cohort,
        "candidates": len(sources),
        "locked": sum(r["status"] == "locked" for r in results),
        "rejected": [r for r in results if r["status"] != "locked"],
        "publication_promoted": False,
    }
    (output_dir / "recovery_summary.json").write_text(json.dumps(summary, indent=2))
    return summary
