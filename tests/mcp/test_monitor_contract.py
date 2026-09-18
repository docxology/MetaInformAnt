"""Contract tests for the standalone Amalgkit monitor surface.

The checkout ships a standalone monitor (``amalgkit_monitor.build_status``),
not an MCP transport. These tests pin the adapter-ready contract: explicit
data roots, database/receipt-backed readiness states, resolved evidence
paths, and the permanent withholding of biological inference — even when a
fully finalized campaign fixture is present.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from metainformant.mcp import tool_adapters
from metainformant.mcp.tools import amalgkit_monitor
from metainformant.rna.engine.progress_db import ProgressDB
from metainformant.rna.engine.provenance import (
    DOWNSTREAM_PROVENANCE_FILENAME,
    DOWNSTREAM_PROVENANCE_SCHEMA,
    DOWNSTREAM_STEPS,
    write_downstream_provenance,
)

FULL_STEPS = list(DOWNSTREAM_STEPS)


def _seed_cohort(data_root: Path, species: str, states: list[str]) -> None:
    """Create a real progress database with one sample per requested state."""
    db = ProgressDB(data_root / "pipeline_progress.db")
    db.init_species(species, [f"SRR{index}" for index in range(len(states))])
    for index, state in enumerate(states):
        db.set_state(species, f"SRR{index}", state)
    db.close()


def _write_receipt(data_root: Path, species: str) -> Path:
    """Write a real downstream provenance receipt via the production writer."""
    return write_downstream_provenance(
        data_root / species / "work",
        species=species,
        config_path=data_root / "config.yaml",
        quantified_samples=2,
        steps=FULL_STEPS,
    )


def _snapshot(data_root: Path, **overrides: Any) -> dict[str, object]:
    """Build a monitor snapshot without touching the live process table."""
    return amalgkit_monitor.build_status(data_root=data_root, inspect_processes=False, **overrides)


def _assert_no_inference_payload(snapshot: dict[str, object]) -> None:
    """Assert the snapshot offers no biological-inference state anywhere."""

    forbidden_keys = {
        "expression",
        "matrix",
        "causal",
        "significance",
        "p_value",
        "correlation",
        "enrichment",
        "prediction",
        "conclusion",
        "biological_inference_offered",
    }

    def walk(node: object) -> None:
        if isinstance(node, dict):
            for key, value in node.items():
                assert str(key) not in forbidden_keys
                if "infer" in str(key).lower():
                    assert value == "withheld", f"readiness key {key!r} must be withheld, got {value!r}"
                walk(value)
        elif isinstance(node, list):
            for item in node:
                walk(item)

    walk(snapshot)
    serialized = json.dumps(snapshot)
    assert "biological_inference" in serialized  # the withholding state itself is exposed
    assert '"biological_inference":"offered"' not in serialized.replace(" ", "")


# --- explicit roots -----------------------------------------------------------


def test_build_status_requires_explicit_data_root() -> None:
    """The monitor never guesses a root: data_root is keyword-only and required."""

    with pytest.raises(TypeError):
        amalgkit_monitor.build_status(inspect_processes=False)  # type: ignore[call-arg]
    with pytest.raises(TypeError):
        amalgkit_monitor.build_status(Path("."))  # type: ignore[misc]


def test_evidence_paths_are_resolved(tmp_path: Path) -> None:
    """Evidence records the resolved absolute root and log path, not inputs."""

    root = tmp_path / "campaign"
    root.mkdir()
    snapshot = _snapshot(root)

    evidence = snapshot["evidence"]
    assert isinstance(evidence, dict)
    assert evidence["data_root"] == str(root.resolve())
    assert evidence["log_file"] == str((root / "results" / "full_campaign.log").resolve())
    assert Path(str(evidence["data_root"])).is_absolute()


def test_explicit_log_file_override_is_reflected_in_evidence(tmp_path: Path) -> None:
    """A custom log path is resolved and used for progress parsing."""

    root = tmp_path / "campaign"
    log = root / "logs" / "custom.log"
    log.parent.mkdir(parents=True)
    log.write_text("[7/20] SRR000007 complete\n", encoding="utf-8")
    root.mkdir(exist_ok=True)

    snapshot = _snapshot(root, log_file=log)

    evidence = snapshot["evidence"]
    assert isinstance(evidence, dict)
    assert evidence["log_file"] == str(log.resolve())
    progress = snapshot["progress"]
    assert isinstance(progress, dict)
    assert progress["processed"] == 7
    assert progress["total"] == 20


# --- database/receipt-backed readiness ----------------------------------------


def test_readiness_without_database_reports_unknown_cohort(tmp_path: Path) -> None:
    """No database means unknown readiness and a recorded database path."""

    root = tmp_path / "campaign"
    root.mkdir()

    snapshot = _snapshot(root)
    readiness = snapshot["readiness"]
    assert isinstance(readiness, dict)
    assert readiness["cohort"] == "unknown"
    assert readiness["descriptive_analysis"] == "withheld"
    assert readiness["biological_inference"] == "withheld"
    assert readiness["database"] == str(root / "pipeline_progress.db")
    assert snapshot["status"] == "stopped"
    _assert_no_inference_payload(snapshot)


def test_unresolved_states_keep_cohort_partial(tmp_path: Path) -> None:
    """Pending, in-flight, and failed samples all block cohort readiness."""

    for unresolved_state in ("pending", "downloading", "quantifying", "failed", "quarantined"):
        root = tmp_path / f"campaign-{unresolved_state}"
        root.mkdir()
        _seed_cohort(root, "Apis_mellifera", ["quantified", unresolved_state])

        snapshot = _snapshot(root)
        readiness = snapshot["readiness"]
        assert isinstance(readiness, dict)
        assert readiness["cohort"] == "partial_or_unresolved", unresolved_state


def test_fully_quantified_cohort_is_ready(tmp_path: Path) -> None:
    """A cohort with every sample quantified and no failures is ready."""

    root = tmp_path / "campaign"
    root.mkdir()
    _seed_cohort(root, "Apis_mellifera", ["quantified", "quantified"])
    db = ProgressDB(root / "pipeline_progress.db")
    db.init_species("Bombus_terrestris", ["SRR1"])
    db.set_state("Bombus_terrestris", "SRR1", "quantified")
    db.close()

    snapshot = _snapshot(root)
    readiness = snapshot["readiness"]
    assert isinstance(readiness, dict)
    assert readiness["cohort"] == "ready"
    assert readiness["descriptive_analysis"] == "withheld"  # no receipts yet
    assert readiness["biological_inference"] == "withheld"


def test_receipt_backed_descriptive_analysis_states(tmp_path: Path) -> None:
    """Descriptive analysis follows valid per-species provenance receipts."""

    root = tmp_path / "campaign"
    root.mkdir()
    _seed_cohort(root, "Apis_mellifera", ["quantified", "quantified"])
    _write_receipt(root, "Apis_mellifera")

    snapshot = _snapshot(root)
    readiness = snapshot["readiness"]
    assert isinstance(readiness, dict)
    assert readiness["descriptive_analysis"] == "receipt_present"

    # A second species without a receipt demotes the aggregate state.
    db = ProgressDB(root / "pipeline_progress.db")
    db.init_species("Bombus_terrestris", ["SRR1"])
    db.set_state("Bombus_terrestris", "SRR1", "quantified")
    db.close()
    snapshot = _snapshot(root)
    readiness = snapshot["readiness"]
    assert isinstance(readiness, dict)
    assert readiness["descriptive_analysis"] == "partial_or_stale"


def test_receipt_without_required_steps_is_not_evidence(tmp_path: Path) -> None:
    """A receipt missing required steps or the v2 schema is not counted."""

    root = tmp_path / "campaign"
    root.mkdir()
    _seed_cohort(root, "Apis_mellifera", ["quantified"])

    incomplete = root / "Apis_mellifera" / "work" / DOWNSTREAM_PROVENANCE_FILENAME
    incomplete.parent.mkdir(parents=True)
    incomplete.write_text(
        json.dumps(
            {
                "schema": DOWNSTREAM_PROVENANCE_SCHEMA,
                "steps": ["merge"],
                "amalgkit_version": "1.0.0",
                "amalgkit_release_tag": "v1.0.0",
                "amalgkit_source_revision": "abc123",
            }
        ),
        encoding="utf-8",
    )
    readiness = _snapshot(root)["readiness"]
    assert isinstance(readiness, dict)
    assert readiness["descriptive_analysis"] == "withheld"

    wrong_schema_payload = json.loads(incomplete.read_text(encoding="utf-8"))
    wrong_schema_payload["steps"] = FULL_STEPS
    wrong_schema_payload["schema"] = "metainformant.rna.downstream.v1"
    incomplete.write_text(json.dumps(wrong_schema_payload), encoding="utf-8")
    readiness = _snapshot(root)["readiness"]
    assert isinstance(readiness, dict)
    assert readiness["descriptive_analysis"] == "withheld"


def test_unreadable_database_is_reported_not_raised(tmp_path: Path) -> None:
    """A corrupt database degrades to an explicit unreadable state."""

    root = tmp_path / "campaign"
    root.mkdir()
    (root / "pipeline_progress.db").write_bytes(b"this is not a sqlite database" * 16)

    snapshot = _snapshot(root)
    readiness = snapshot["readiness"]
    assert isinstance(readiness, dict)
    assert readiness["cohort"] == "unreadable_database"
    assert readiness["biological_inference"] == "withheld"
    _assert_no_inference_payload(snapshot)


# --- operational status --------------------------------------------------------


def test_writer_lock_reports_running_without_process_scan(tmp_path: Path) -> None:
    """Writer locks flip status deterministically with process inspection off."""

    root = tmp_path / "campaign"
    results = root / "results"
    results.mkdir(parents=True)

    assert _snapshot(root)["status"] == "stopped"
    (results / ".full_campaign.lock").write_text("", encoding="utf-8")
    snapshot = _snapshot(root)
    assert snapshot["status"] == "running"
    evidence = snapshot["evidence"]
    assert isinstance(evidence, dict)
    assert evidence["writer_lock"] is True


# --- negative contract: biological inference withheld --------------------------


def test_biological_inference_is_withheld_even_with_finalized_matrix(tmp_path: Path) -> None:
    """The maximally complete fixture still gets no biological-inference state.

    Every readiness ingredient is present: fully quantified cohort, valid v2
    provenance receipts, and finalized stage outputs (the merged/finalized
    expression matrix) on disk. The monitor must still expose only
    operational/descriptive readiness and withhold inference.
    """

    root = tmp_path / "campaign"
    root.mkdir()
    for species in ("Apis_mellifera", "Bombus_terrestris"):
        _seed_cohort(root, species, ["quantified", "quantified"])
        _write_receipt(root, species)
        work = root / species / "work"
        (work / "merge" / "merged.tsv").parent.mkdir(parents=True)
        (work / "merge" / "merged.tsv").write_text("gene\tSRR0\tSRR1\ngene1\t1\t2\n", encoding="utf-8")
        finalize = work / "finalize"
        finalize.mkdir(parents=True)
        (finalize / f"{species}_expression.tsv").write_text(
            f"gene\t{species}_SRR0\t{species}_SRR1\ngene1\t10.5\t12.25\n", encoding="utf-8"
        )

    snapshot = _snapshot(root)
    readiness = snapshot["readiness"]
    assert isinstance(readiness, dict)
    # The descriptive readiness is fully earned...
    assert readiness["cohort"] == "ready"
    assert readiness["descriptive_analysis"] == "receipt_present"
    # ...and biological inference is still refused.
    assert readiness["biological_inference"] == "withheld"
    assert readiness["biological_inference"] != "ready"
    _assert_no_inference_payload(snapshot)


def test_registered_monitor_tool_offers_no_inference_surface() -> None:
    """The advertised adapter surface is the read-only monitor, nothing more."""

    assert [tool.name for tool in tool_adapters.TOOLS] == ["amalgkit_monitor"]
    properties = tool_adapters.TOOL.input_schema["properties"]
    assert set(properties) == {"data_root", "log_file", "inspect_processes"}
    assert all("infer" not in str(name).lower() for name in properties)
    assert tool_adapters.TOOL.metadata["read_only"] is True
