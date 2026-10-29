from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

from metainformant.rna.engine.campaign_observation import CampaignObservation
from metainformant.visualization.plots.rna_campaign import plot_acquisition_observation


def verify_observation_bundle(output_dir: Path) -> dict:
    output = output_dir.resolve(strict=True)
    manifest = json.loads((output / "bundle_manifest.json").read_text())
    if (
        manifest.get("schema") != "metainformant.rna.operational_bundle.v1"
        or manifest.get("scientific_release_ready") is not False
    ):
        raise ValueError("invalid operational bundle contract")
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, dict) or "observation.json" not in artifacts:
        raise ValueError("bundle has no observation artifact")
    actual = set()
    for path in output.rglob("*"):
        if path.is_symlink():
            raise ValueError("bundle contains a symlink")
        if path.is_file() and path != output / "bundle_manifest.json":
            actual.add(str(path.relative_to(output)))
    if set(artifacts) != actual:
        raise ValueError("bundle artifact inventory mismatch")
    for name, expected in artifacts.items():
        path = output / name
        if not path.resolve().is_relative_to(output):
            raise ValueError("unsafe bundle artifact path")
        raw = path.read_bytes()
        if expected.get("sha256") != hashlib.sha256(raw).hexdigest() or expected.get(
            "bytes"
        ) != len(raw):
            raise ValueError(f"bundle artifact checksum mismatch: {name}")
    observation_raw = (output / "observation.json").read_bytes()
    if hashlib.sha256(observation_raw).hexdigest() != manifest["observation_sha256"]:
        raise ValueError("bundle observation checksum mismatch")
    observation = json.loads(observation_raw)
    if observation.get("scientific_release_ready") is not False:
        raise ValueError("observation cannot grant scientific release readiness")
    for name, source in observation["source_files"].items():
        if source["captured"]:
            record = artifacts.get(f"inputs/{name}")
            if (
                record is None
                or record["sha256"] != source["sha256"]
                or record["bytes"] != source["bytes"]
            ):
                raise ValueError(f"captured source binding mismatch: {name}")
    return manifest


def write_observation_bundle(
    observation: CampaignObservation,
    output_dir: Path,
    *,
    render_pdf: bool = True,
    title: str = "RNA-seq acquisition",
) -> Path:
    output = output_dir.resolve()
    payload = observation.payload
    if any(
        output.is_relative_to(Path(payload[name]).resolve())
        for name in ("data_root", "operations_root")
    ):
        raise ValueError(
            "observation bundles must be isolated from live data and operations roots"
        )
    if payload.get("scientific_release_ready") is not False:
        raise ValueError(
            "working observations must explicitly deny scientific release readiness"
        )
    for name, raw in observation.sources.items():
        if Path(name).is_absolute() or ".." in Path(name).parts:
            raise ValueError("unsafe captured source name")
        if hashlib.sha256(raw).hexdigest() != payload["source_files"][name]["sha256"]:
            raise ValueError("captured source checksum mismatch")
    output.mkdir(parents=True, exist_ok=False)
    for name, raw in observation.sources.items():
        path = output / "inputs" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
    observation_bytes = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
    digest = hashlib.sha256(observation_bytes).hexdigest()
    (output / "observation.json").write_bytes(observation_bytes)
    counts, scope = payload["counts"], payload["scope"]
    with (output / "counts.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "count", "observation_sha256"])
        writer.writerows(
            (key, value, digest) for key, value in {**scope, **counts}.items()
        )
    plot_acquisition_observation(
        payload, output / "acquisition_observation", title=title
    )
    caption = (
        f"Observation {digest}. Frozen scope: {scope['total']:,} runs. The current lower bound combines "
        "the producer-reported startup audit and restore-bound new successes. Panel B distinguishes reused, "
        "new and unconfirmed successes from deferred and failed attempts. These counts do not establish biological "
        "replication or scientific release readiness.\n"
    )
    (output / "figure_caption.md").write_text(caption)
    (output / "figure_alt.txt").write_text(caption)
    observed = datetime.fromtimestamp(payload["observed_at_epoch"], UTC).isoformat()
    transfer = payload["transfer"]
    rate = (
        f"{transfer['MiB_s_lower_bound']:.3f} MiB/s retained-growth lower bound over {transfer['seconds']:.3f} seconds"
        if transfer["available"]
        else f"Unavailable: {transfer['reason']}"
    )
    report = [
        f"# {title}: working operational supplement",
        "",
        "**Working draft. Operational evidence only; scientific release is not authorized by this bundle.**",
        "",
        f"Observation: {observed}",
        "",
        f"Observation SHA256: `{digest}`",
        "",
        f"Selected data root: `{payload['data_root']}`",
        "",
        "## Methods and reproducibility",
        "",
        (
            "The parent MetaInformAnt observation API reads explicit frozen manifest, partition and task journal inputs. "
            "It captures stable control records and newline-complete append-only prefixes without modifying the producer. "
            "Startup counts are the producer's audited report. Unique new successes increment that baseline only when "
            "the current restore marker, frozen reference/configuration bindings, provenance contract and abundance "
            "checksum agree. Reused outputs do not increase the baseline. Missing or inconsistent evidence remains unconfirmed."
        ),
        "",
        (
            "The acquisition method retains bounded primary/recovery ownership, resumable ENA transfers and optional "
            "verified NCBI prefetch. Configured SSD extraction stages authoritative bulk archives under a shared space "
            "reservation. Quantification output alone cannot authorize raw reclamation: the checkpoint must verify "
            "restoration under the same reference/configuration contract."
        ),
        "",
        "## Operational observation",
        "",
        "| Quantity | Value |",
        "| --- | ---: |",
        f"| Frozen run scope | {scope['total']:,} |",
        f"| Current lower bound | {counts['current_lower_bound']:,} |",
        f"| Supported, unconfirmed | {counts['supported_unconfirmed']:,} |",
        f"| Unresolved prerequisites, retained in scope | {scope['unsupported']:,} |",
        f"| New, restore-bound in this phase | {counts['new_restore_bound']:,} |",
        f"| New, unconfirmed in this phase | {counts['new_unconfirmed']:,} |",
        f"| Reused successful attempts | {counts['reused']:,} |",
        f"| Unsuccessful attempts, including deferrals | {counts['failed_attempts']:,} |",
        f"| Deferred download attempts | {counts['deferred_attempts']:,} |",
        f"| Failed attempts, excluding deferrals | {counts['failed_nondeferred_attempts']:,} |",
        "",
        f"Transfer observation: {rate}.",
        "",
        "![Frozen acquisition observation](acquisition_observation.png)",
        "",
        caption.strip(),
        "",
        "## Interpretation and next scientific gates",
        "",
        (
            "Run accessions are computational records, not independent biological replicates. Current quantification "
            "does not establish finalized matrices, harmonized biological units, validated orthology, study-aware "
            "contrasts, phylogenetic inference or a completed cohort. No biological effect or significance claim is made. "
            "Long-read support is a frozen reference/backend/output-contract decision, not a claim that contemporary "
            "quantifiers lack long-read methods. A new backend requires an explicitly amended, validated contract."
        ),
        "",
        "## Scholarly context",
        "",
        (
            "SRA and ENA provide the source archive context ([Leinonen et al., 2011](https://doi.org/10.1093/nar/gkq1019); "
            "[Cummins et al., 2022](https://doi.org/10.1093/nar/gkab1051)). The original Kallisto algorithm "
            "([Bray et al., 2016](https://doi.org/10.1038/nbt.3519)) and the distinct lr-kallisto method "
            "([Loving et al., 2025](https://doi.org/10.1371/journal.pcbi.1013692)) motivate explicit backend/reference "
            "contracts. Study-aware RNA-seq synthesis ([Rau et al., 2014](https://doi.org/10.1186/1471-2105-15-91)) "
            "does not treat downloaded runs as automatically independent replicates. These sources support method "
            "context, not the present observation counts or a biological conclusion."
        ),
        "",
        "## Attribution limits",
        "",
        *[f"- {item}" for item in payload["limitations"]],
        "",
        "## Failure evidence",
        "",
        *[
            f"- {number} attempt(s): {reason}"
            for reason, number in payload["failure_reasons"].items()
        ],
        "",
        "## Verification",
        "",
        (
            "`bundle_manifest.json` binds every captured input and generated artifact. "
            "The source table, caption and figure use this observation hash. The full project manuscript and its "
            "quiescent scientific-release gates are separate from this working supplement. Control records and restore "
            "evidence are copied into inputs; abundance tables remain at their source paths, with content digests "
            "and checked stat signatures recorded in observation.json."
        ),
        "",
    ]
    report_path = output / "working_supplement.md"
    report_path.write_text("\n".join(report))
    if render_pdf:
        result = subprocess.run(
            [
                "pandoc",
                str(report_path),
                "--from=gfm",
                "--pdf-engine=xelatex",
                "--resource-path",
                str(output),
                "-V",
                "geometry:margin=0.75in",
                "-V",
                "fontsize=10pt",
                "-o",
                str(output / "working_supplement.pdf"),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode:
            raise RuntimeError(f"working supplement render failed: {result.stderr}")
        (output / "render.log").write_text(result.stdout + result.stderr)
    artifacts = {
        str(path.relative_to(output)): {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "bytes": path.stat().st_size,
        }
        for path in sorted(output.rglob("*"))
        if path.is_file()
    }
    manifest = {
        "schema": "metainformant.rna.operational_bundle.v1",
        "observation_sha256": digest,
        "scientific_release_ready": False,
        "artifacts": artifacts,
    }
    (output / "bundle_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    verify_observation_bundle(output)
    return output
