from __future__ import annotations

import hashlib
import json
import os
import sys
import threading
from pathlib import Path

import pytest

from metainformant.rna.engine.acquisition_manifest import sha256_file
from metainformant.rna.engine.campaign_observation import (
    ObservationError,
    capture_campaign_observation,
)
from metainformant.rna.engine.provenance import (
    QUANT_PROVENANCE_SCHEMA,
    quantification_contract_id,
)


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) + "\n")


@pytest.fixture
def observation_files(tmp_path: Path) -> dict[str, Path]:
    operations, data = tmp_path / "operations", tmp_path / "data"
    manifest = operations / "frozen/manifest.jsonl"
    manifest.parent.mkdir(parents=True)
    tasks = [
        {
            "schema": "metainformant.rna.acquisition_task.v1",
            "task_id": f"species/SRR{i}",
            "species": "species",
            "accession": f"SRR{i}",
            "batch_index": i,
            "config_name": "test.yaml",
            "reference_index_sha256": "a" * 64,
            "config_sha256": "b" * 64,
        }
        for i in range(1, 5)
    ]
    manifest.write_text("".join(json.dumps(task) + "\n" for task in tasks))
    write_json(
        manifest.with_name("snapshot.json"),
        {
            "schema": "metainformant.rna.acquisition_snapshot.v1",
            "cloud_launch_policy": "checkpointed",
            "manifest_sha256": sha256_file(manifest),
            "task_count": 4,
        },
    )
    write_json(
        operations / "active_runtime.json", {"pid": 123, "started_at_epoch": 100}
    )
    write_json(
        operations / "active_scope.json",
        {
            "scope_total": 4,
            "eligible": 3,
            "verified_existing": 1,
            "unsupported": {"species/SRR4": "requires amended backend"},
        },
    )
    write_json(
        operations / "live_progress.json",
        {"observed_at_epoch": 250, "checkpointed_this_phase": 1},
    )
    selection = operations / "partitions/explicit_primary.json"
    write_json(
        selection,
        {
            "schema": "metainformant.rna.acquisition_partition.v1",
            "task_count": 3,
            "manifest_sha256": sha256_file(manifest),
            "snapshot_sha256": sha256_file(manifest.with_name("snapshot.json")),
            "task_ids": [task["task_id"] for task in tasks[:3]],
        },
    )
    journal = data / "acquisition_runs/explicit_run/task_results.jsonl"
    journal.parent.mkdir(parents=True)
    rows = [
        {
            "task_id": f"species/SRR{i}",
            "species": "species",
            "srr": f"SRR{i}",
            "quantified": i != 3,
            "skipped": i == 1,
            "error": None if i != 3 else "download failed",
        }
        for i in range(1, 4)
    ]
    journal.write_text("".join(json.dumps(row) + "\n" for row in rows))
    sample = data / "species/work/quant/SRR2"
    sample.mkdir(parents=True)
    abundance = sample / "SRR2_abundance.tsv"
    abundance.write_text(
        "target_id\tlength\teff_length\test_counts\ttpm\ngene1\t100\t75\t10\t1000000\n"
    )
    proof = sample / ".local_restore_verified.json"
    provenance = sample / ".metainformant_quant_provenance.json"
    provenance_payload = {
        "schema": QUANT_PROVENANCE_SCHEMA,
        "command": ["amalgkit", "quant", "--batch", "2"],
        "config_sha256": "b" * 64,
        "run_accession": "SRR2",
        "species": "species",
        "quantification_file": abundance.name,
        "quantification_file_sha256": sha256_file(abundance),
    }
    provenance_payload["quant_contract_id"] = quantification_contract_id(
        provenance_payload
    )
    write_json(provenance, provenance_payload)
    write_json(
        proof,
        {
            "contract_id": provenance_payload["quant_contract_id"],
            "reference_index_sha256": "a" * 64,
            "verified_at_epoch": 240,
        },
    )
    write_json(
        operations / "profiling-20261019/live-transfer-observation.json",
        {
            "pid": 123,
            "observed_at_epoch": 250,
            "seconds": 10,
            "retained_transfer_growth_bytes": 10 * 1024**2,
            "tracked_curl_transfers": 2,
        },
    )
    (operations / "raw_reclamation.jsonl").write_text(
        json.dumps({"bytes_freed": 100, "errors": []}) + "\n"
    )
    return {
        "operations": operations,
        "data": data,
        "selection": selection,
        "journal": journal,
        "manifest": manifest,
        "proof": proof,
        "provenance": provenance,
        "abundance": abundance,
    }


def capture(files: dict[str, Path], **kwargs):
    return capture_campaign_observation(
        files["operations"],
        files["data"],
        selection_path=files["selection"],
        journal_path=files["journal"],
        transfer_observation_path=files["operations"]
        / "profiling-20261019/live-transfer-observation.json",
        observed_at=kwargs.pop("observed_at", 300),
        **kwargs,
    )


def test_confirmed_lower_bound_separates_new_reused_and_journal_successes(
    observation_files,
) -> None:
    observation = capture(observation_files)
    counts = observation.payload["counts"]
    assert counts["current_lower_bound"] == 2
    assert counts["journal_successes"] == 2
    assert (
        counts["reused"]
        == counts["new_restore_bound"]
        == counts["failed_attempts"]
        == 1
    )
    assert observation.payload["transfer"]["MiB_s_lower_bound"] == 1
    assert observation.payload["scientific_release_ready"] is False
    assert counts["deferred_attempts"] == 0
    assert counts["failed_nondeferred_attempts"] == 1
    for name, raw in observation.sources.items():
        assert (
            observation.payload["source_files"][name]["sha256"]
            == hashlib.sha256(raw).hexdigest()
        )


def test_download_deferrals_remain_unsuccessful_but_are_reported_separately(observation_files) -> None:
    journal = observation_files["journal"]
    rows = [json.loads(line) for line in journal.read_text().splitlines()]
    rows[-1].update(download_deferred=True, status="deferred_download")
    journal.write_text("".join(json.dumps(row) + "\n" for row in rows))
    counts = capture(observation_files).payload["counts"]
    assert counts["failed_attempts"] == counts["deferred_attempts"] == 1
    assert counts["failed_nondeferred_attempts"] == counts["new_unconfirmed"] == 0
    assert counts["current_lower_bound"] == 2


def test_nonboolean_deferral_flag_is_rejected(observation_files) -> None:
    journal = observation_files["journal"]
    rows = [json.loads(line) for line in journal.read_text().splitlines()]
    rows[-1]["download_deferred"] = "yes"
    journal.write_text("".join(json.dumps(row) + "\n" for row in rows))
    with pytest.raises(ObservationError, match="booleans"):
        capture(observation_files)


@pytest.mark.parametrize(
    "mutation", ["missing", "contract", "reference", "config", "time", "abundance"]
)
def test_unconfirmed_success_never_increments_lower_bound(
    observation_files, mutation: str
) -> None:
    files = observation_files
    if mutation == "missing":
        files["proof"].unlink()
    elif mutation == "abundance":
        files["abundance"].write_bytes(b"corrupted abundance bytes")
    elif mutation == "config":
        value = json.loads(files["provenance"].read_text())
        value["config_sha256"] = "d" * 64
        write_json(files["provenance"], value)
    else:
        value = json.loads(files["proof"].read_text())
        key, replacement = {
            "contract": ("contract_id", "d" * 64),
            "reference": ("reference_index_sha256", "d" * 64),
            "time": ("verified_at_epoch", 99),
        }[mutation]
        value[key] = replacement
        write_json(files["proof"], value)
    payload = capture(files).payload
    assert payload["counts"]["current_lower_bound"] == 1
    assert payload["counts"]["new_unconfirmed"] == 1
    assert payload["unconfirmed_new_tasks"]["species/SRR2"]


def test_abundance_and_sidecar_mutation_cannot_reuse_old_restore_contract(
    observation_files,
) -> None:
    files = observation_files
    files["abundance"].write_bytes(b"replacement abundance bytes")
    provenance = json.loads(files["provenance"].read_text())
    provenance["quantification_file_sha256"] = sha256_file(files["abundance"])
    write_json(files["provenance"], provenance)
    payload = capture(files).payload
    assert payload["counts"]["current_lower_bound"] == 1
    assert payload["counts"]["new_unconfirmed"] == 1


@pytest.mark.parametrize(
    "mutation", ["duplicate", "unknown", "malformed", "identity", "boolean"]
)
def test_invalid_complete_journal_fails_closed(
    observation_files, mutation: str
) -> None:
    path = observation_files["journal"]
    rows = path.read_text().splitlines()
    if mutation == "duplicate":
        rows.append(rows[0])
    elif mutation == "malformed":
        rows.insert(1, "not json")
    else:
        row = json.loads(rows[0])
        if mutation == "unknown":
            row["task_id"] = "species/SRR999"
        elif mutation == "identity":
            row["srr"] = "SRR999"
        else:
            row["quantified"] = "true"
        rows[0] = json.dumps(row)
    path.write_text("\n".join(rows) + "\n")
    with pytest.raises(ValueError):
        capture(observation_files)


def test_partial_tail_is_captured_but_never_counted(observation_files) -> None:
    path = observation_files["journal"]
    with path.open("ab") as handle:
        handle.write(b'{"task_id":')
    observation = capture(observation_files)
    assert observation.payload["journal_incomplete_tail"] is True
    assert observation.payload["counts"]["current_lower_bound"] == 2
    assert observation.sources["task_results.jsonl"].endswith(b'{"task_id":')


@pytest.mark.parametrize("change", ["stale", "pid", "future", "missing"])
def test_transfer_unavailable_is_not_zero(observation_files, change: str) -> None:
    path = (
        observation_files["operations"]
        / "profiling-20261019/live-transfer-observation.json"
    )
    if change == "missing":
        path.unlink()
    else:
        value = json.loads(path.read_text())
        if change == "pid":
            value["pid"] = 999
        else:
            value["observed_at_epoch"] = 299 if change == "stale" else 301
        write_json(path, value)
    result = capture(observation_files, max_rate_age_seconds=0.5).payload["transfer"]
    assert result["available"] is False
    assert "MiB_s_lower_bound" not in result
    assert result["reason"]


def test_transfer_window_must_fit_current_runtime(observation_files) -> None:
    path = (
        observation_files["operations"]
        / "profiling-20261019/live-transfer-observation.json"
    )
    transfer = json.loads(path.read_text())
    transfer["seconds"] = 200
    write_json(path, transfer)
    assert capture(observation_files).payload["transfer"]["available"] is False


def test_capture_interval_and_exact_counted_prefix_are_recorded(
    observation_files,
) -> None:
    with observation_files["journal"].open("ab") as handle:
        handle.write(b'{"partial":')
    observation = capture(observation_files)
    payload = observation.payload
    assert payload["capture_started_at_epoch"] <= payload["capture_finished_at_epoch"]
    assert payload["capture_elapsed_seconds"] >= 0
    source = payload["source_files"]["task_results.jsonl"]
    assert source["excluded_tail_bytes"] == len(b'{"partial":')
    assert (
        source["complete_prefix_bytes"] + source["excluded_tail_bytes"]
        == source["bytes"]
    )
    assert observation.sources["task_results.jsonl"][
        : source["complete_prefix_bytes"]
    ].endswith(b"\n")


def test_scope_and_partition_manifest_mismatches_refuse(observation_files) -> None:
    path = observation_files["selection"]
    value = json.loads(path.read_text())
    value["manifest_sha256"] = "0" * 64
    write_json(path, value)
    with pytest.raises(ValueError, match="manifest hash"):
        capture(observation_files)


def test_proof_symlink_escape_remains_unconfirmed(
    observation_files, tmp_path: Path
) -> None:
    path = observation_files["proof"]
    outside = tmp_path / "unowned_proof.json"
    path.rename(outside)
    path.symlink_to(outside)
    assert capture(observation_files).payload["counts"]["current_lower_bound"] == 1


def test_bundle_replays_from_captured_bytes_without_rereading_live_sources(
    observation_files, tmp_path: Path
) -> None:
    from metainformant.rna.engine.campaign_observation_report import (
        write_observation_bundle,
    )

    observation = capture(observation_files)
    observation_files["journal"].write_text("changed after capture\n")
    first = write_observation_bundle(observation, tmp_path / "first", render_pdf=False)
    second = write_observation_bundle(
        observation, tmp_path / "second", render_pdf=False
    )
    assert (first / "observation.json").read_bytes() == (
        second / "observation.json"
    ).read_bytes()
    assert (first / "acquisition_observation.png").read_bytes() == (
        second / "acquisition_observation.png"
    ).read_bytes()
    assert (first / "acquisition_observation.pdf").read_bytes() == (
        second / "acquisition_observation.pdf"
    ).read_bytes()
    assert (first / "inputs/task_results.jsonl").read_bytes() == observation.sources[
        "task_results.jsonl"
    ]
    manifest = json.loads((first / "bundle_manifest.json").read_text())
    assert manifest["scientific_release_ready"] is False
    for name, record in manifest["artifacts"].items():
        assert sha256_file(first / name) == record["sha256"]
    with pytest.raises(FileExistsError):
        write_observation_bundle(observation, first, render_pdf=False)


def test_bundle_cannot_write_into_active_roots(observation_files) -> None:
    from metainformant.rna.engine.campaign_observation_report import (
        write_observation_bundle,
    )

    observation = capture(observation_files)
    for key in ("operations", "data"):
        destination = observation_files[key] / "forbidden_bundle"
        with pytest.raises(ValueError, match="isolated"):
            write_observation_bundle(observation, destination, render_pdf=False)
        assert not destination.exists()


@pytest.mark.parametrize("change", ["release", "unsafe_source", "source_digest"])
def test_writer_refuses_invalid_observation_before_creating_output(
    observation_files, tmp_path: Path, change: str
) -> None:
    from metainformant.rna.engine.campaign_observation_report import (
        write_observation_bundle,
    )

    observation = capture(observation_files)
    if change == "release":
        observation.payload["scientific_release_ready"] = True
    elif change == "unsafe_source":
        observation.sources["../escape"] = b"invalid"
    else:
        observation.sources["runtime.json"] = b"different bytes"
    output = tmp_path / "refused"
    with pytest.raises(ValueError):
        write_observation_bundle(observation, output, render_pdf=False)
    assert not output.exists()


@pytest.mark.parametrize(
    "change",
    [
        "content",
        "extra",
        "nested_manifest",
        "missing",
        "symlink",
        "release",
        "observation_hash",
    ],
)
def test_bundle_verifier_rejects_tampering(
    observation_files, tmp_path: Path, change: str
) -> None:
    from metainformant.rna.engine.campaign_observation_report import (
        verify_observation_bundle,
        write_observation_bundle,
    )

    output = write_observation_bundle(
        capture(observation_files), tmp_path / "bundle", render_pdf=False
    )
    assert verify_observation_bundle(output)["scientific_release_ready"] is False
    counts = output / "counts.csv"
    if change == "content":
        counts.write_text("corrupt\n")
    elif change == "extra":
        (output / "unlisted.txt").write_text("unbound\n")
    elif change == "nested_manifest":
        (output / "unlisted").mkdir()
        (output / "unlisted/bundle_manifest.json").write_text("{}\n")
    elif change == "missing":
        counts.rename(tmp_path / "preserved.csv")
    elif change == "symlink":
        counts.rename(tmp_path / "preserved.csv")
        counts.symlink_to(tmp_path / "preserved.csv")
    else:
        manifest_path = output / "bundle_manifest.json"
        manifest = json.loads(manifest_path.read_text())
        if change == "release":
            manifest["scientific_release_ready"] = True
        else:
            manifest["observation_sha256"] = "0" * 64
        write_json(manifest_path, manifest)
    with pytest.raises(ValueError):
        verify_observation_bundle(output)


def test_source_change_during_capture_refuses_mixed_observation(
    observation_files,
) -> None:
    path = observation_files["operations"] / "active_runtime.json"
    abundance = observation_files["abundance"]
    read_started, changed = threading.Event(), threading.Event()
    active = True
    reader = threading.get_ident()

    def synchronize_read(event: str, args: tuple) -> None:
        if active and threading.get_ident() == reader and event == "open" and args[0] == str(abundance):
            read_started.set()
            assert changed.wait(5), "Real writer did not change the captured runtime"

    sys.addaudithook(synchronize_read)

    def writer() -> None:
        if not read_started.wait(5):
            return
        temporary = path.with_suffix(".next")
        write_json(temporary, {"pid": 123, "started_at_epoch": 100, "nonce": "changed-during-capture"})
        os.replace(temporary, path)
        changed.set()

    thread = threading.Thread(target=writer)
    thread.start()
    try:
        with pytest.raises(ObservationError, match="source changed"):
            capture(observation_files, stability_attempts=1)
    finally:
        active = False
        thread.join(5)
    assert not thread.is_alive()
    assert read_started.is_set() and changed.is_set()


def test_nonregular_input_is_refused_without_blocking(observation_files) -> None:
    path = observation_files["journal"]
    path.rename(path.with_suffix(".preserved"))
    os.mkfifo(path)
    with pytest.raises(ObservationError, match="regular file"):
        capture(observation_files)


def test_pdf_export_does_not_depend_on_build_clock(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import matplotlib.pyplot as plt

    from metainformant.visualization.config.conventions import save_figure_deterministic

    figure, axis = plt.subplots()
    axis.plot([0, 1], [1, 0])
    try:
        first, second = tmp_path / "first.pdf", tmp_path / "second.pdf"
        monkeypatch.setenv("SOURCE_DATE_EPOCH", "100")
        save_figure_deterministic(figure, first)
        monkeypatch.setenv("SOURCE_DATE_EPOCH", "200")
        save_figure_deterministic(figure, second)
        assert first.read_bytes() == second.read_bytes()
        assert b"/CreationDate" not in first.read_bytes()
    finally:
        plt.close(figure)


@pytest.mark.parametrize(
    "change",
    [
        "pid_boolean",
        "pid_zero",
        "future_start",
        "nonfinite_start",
        "unsupported",
        "scope_count",
        "baseline_count",
        "old_scope",
        "old_partition",
        "unsupported_selection",
        "outside_journal",
        "array_journal",
        "oversized_journal",
        "zero_window",
        "negative_age",
        "invalid_attempts",
    ],
)
def test_invalid_observation_controls_fail_closed(
    observation_files, tmp_path: Path, change: str
) -> None:
    files = observation_files
    kwargs = {}
    if change in {"pid_boolean", "pid_zero", "future_start", "nonfinite_start"}:
        runtime = {"pid": 123, "started_at_epoch": 100}
        key, value = {
            "pid_boolean": ("pid", True),
            "pid_zero": ("pid", 0),
            "future_start": ("started_at_epoch", 301),
            "nonfinite_start": ("started_at_epoch", float("nan")),
        }[change]
        runtime[key] = value
        write_json(files["operations"] / "active_runtime.json", runtime)
    elif change in {"unsupported", "scope_count", "baseline_count"}:
        path = files["operations"] / "active_scope.json"
        scope = json.loads(path.read_text())
        key, value = {
            "unsupported": ("unsupported", {"species/SRR999": "unknown"}),
            "scope_count": ("scope_total", 5),
            "baseline_count": ("verified_existing", 4),
        }[change]
        scope[key] = value
        write_json(path, scope)
    elif change in {"old_scope", "old_partition"}:
        path = (
            files["selection"]
            if change == "old_partition"
            else files["operations"] / "active_scope.json"
        )
        os.utime(path, (50, 50))
    elif change == "unsupported_selection":
        selection = json.loads(files["selection"].read_text())
        selection["task_ids"][-1] = "species/SRR4"
        write_json(files["selection"], selection)
    elif change == "outside_journal":
        outside = tmp_path / "outside.jsonl"
        outside.write_bytes(files["journal"].read_bytes())
        files["journal"] = outside
    elif change == "array_journal":
        files["journal"].write_text("[]\n")
    elif change == "oversized_journal":
        with files["journal"].open("wb") as handle:
            handle.truncate(64 * 1024**2 + 1)
    elif change == "zero_window":
        path = files["operations"] / "profiling-20261019/live-transfer-observation.json"
        value = json.loads(path.read_text())
        value["seconds"] = 0
        write_json(path, value)
    elif change == "negative_age":
        kwargs["max_rate_age_seconds"] = -1
    else:
        kwargs["stability_attempts"] = 0
    with pytest.raises(ObservationError):
        capture(files, **kwargs)


@pytest.mark.parametrize(
    "change", ["empty", "escape_name", "abundance_symlink", "missing_schema"]
)
def test_invalid_abundance_contract_remains_unconfirmed(
    observation_files, tmp_path: Path, change: str
) -> None:
    files = observation_files
    provenance = json.loads(files["provenance"].read_text())
    if change == "empty":
        files["abundance"].write_bytes(b"")
        provenance["quantification_file_sha256"] = sha256_file(files["abundance"])
    elif change == "escape_name":
        provenance["quantification_file"] = "../outside.tsv"
    elif change == "missing_schema":
        del provenance["schema"]
    else:
        preserved = tmp_path / "outside.tsv"
        files["abundance"].rename(preserved)
        files["abundance"].symlink_to(preserved)
    if change != "missing_schema":
        provenance["quant_contract_id"] = quantification_contract_id(provenance)
        proof = json.loads(files["proof"].read_text())
        proof["contract_id"] = provenance["quant_contract_id"]
        write_json(files["proof"], proof)
    write_json(files["provenance"], provenance)
    assert capture(files).payload["counts"]["new_unconfirmed"] == 1


@pytest.mark.external_tool
def test_thin_cli_renders_and_verifies_real_pdf(
    observation_files, tmp_path: Path, capsys
) -> None:
    from metainformant.rna.engine.campaign_observation_cli import main
    from metainformant.rna.engine.campaign_observation_report import (
        verify_observation_bundle,
    )

    files = observation_files
    output = tmp_path / "cli-bundle"
    assert (
        main(
            [
                "--operations-root",
                str(files["operations"]),
                "--data-root",
                str(files["data"]),
                "--task-selection",
                str(files["selection"]),
                "--journal",
                str(files["journal"]),
                "--transfer-observation",
                str(
                    files["operations"]
                    / "profiling-20261019/live-transfer-observation.json"
                ),
                "--output-dir",
                str(output),
            ]
        )
        == 0
    )
    assert json.loads(capsys.readouterr().out)["scientific_release_ready"] is False
    assert (output / "working_supplement.pdf").read_bytes().startswith(b"%PDF")
    assert (
        verify_observation_bundle(output)["artifacts"]["working_supplement.pdf"][
            "bytes"
        ]
        > 1000
    )


@pytest.mark.parametrize("field", ["current_lower_bound", "failed_attempts"])
def test_plot_refuses_inconsistent_counts(
    observation_files, tmp_path: Path, field: str
) -> None:
    from metainformant.visualization.plots.rna_campaign import (
        plot_acquisition_observation,
    )

    observation = capture(observation_files)
    observation.payload["counts"][field] = -1
    with pytest.raises(ValueError):
        plot_acquisition_observation(observation.payload, tmp_path / "invalid")
    assert not list(tmp_path.glob("invalid.*"))
