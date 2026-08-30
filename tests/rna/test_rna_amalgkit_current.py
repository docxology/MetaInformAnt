"""Current Amalgkit command-contract tests.

These tests exercise the supported Amalgkit 0.16.60 command surface without network
access. They deliberately keep optional cross-species commands separate from
the per-species default workflow.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from metainformant.rna import amalgkit
from metainformant.rna.amalgkit.commands import (
    CURRENT_AMALGKIT_STEPS,
    PER_SPECIES_STEPS,
    SUPPORTED_CLI_OPTIONS,
)
from metainformant.rna.engine.provenance import (
    is_current_downstream,
    write_downstream_provenance,
    write_metadata_provenance,
    write_quant_provenance,
)
from metainformant.rna.steps import STEP_RUNNERS

REPO_ROOT = Path(__file__).resolve().parents[2]
RETIRED_PATHS = (
    REPO_ROOT / "scripts" / "rna" / "run_workflow.py",
    REPO_ROOT / "scripts" / "rna" / "run_workflow_tui.py",
    REPO_ROOT / "src" / "metainformant" / "rna" / "engine" / "orchestration.py",
    REPO_ROOT / "src" / "metainformant" / "rna" / "engine" / "orchestration_multi_species.py",
    REPO_ROOT / "src" / "metainformant" / "rna" / "engine" / "orchestrator.py",
    REPO_ROOT / "src" / "metainformant" / "rna" / "engine" / "monitoring.py",
    REPO_ROOT / "src" / "metainformant" / "rna" / "engine" / "progress_tracker.py",
    REPO_ROOT / "src" / "metainformant" / "rna" / "engine" / "tui_workflow.py",
)
ACTIVE_SURFACE = (
    REPO_ROOT / "scripts" / "rna",
    REPO_ROOT / "src" / "metainformant" / "rna",
    REPO_ROOT / "docs" / "rna",
    REPO_ROOT / "config" / "amalgkit",
    REPO_ROOT / "docs" / "agents" / "rules" / "rna.md",
    REPO_ROOT / "cursorrules" / "rna.cursorrules",
)
RETIRED_REFERENCE = re.compile(
    r"run_workflow(?:_tui)?\.py|run_workflow_for_species|"
    r"metainformant\.rna\.orchestration|orchestration_multi_species|"
    r"metainformant\.rna\.engine\.(?:monitoring|progress_tracker|tui_workflow)|"
    r"(?<!Streaming)\bPipelineOrchestrator\b|\bAK_(?:THREADS|WORK_DIR|LOG_DIR)\b"
)


def _active_text_files() -> list[Path]:
    """Yield tracked-looking source and documentation files in active areas."""

    files: list[Path] = []
    for root in ACTIVE_SURFACE:
        if root.is_file():
            files.append(root)
            continue
        if root.is_dir():
            files.extend(
                path
                for path in root.rglob("*")
                if path.is_file()
                and path.suffix.lower() in {".py", ".md", ".rst", ".yaml", ".yml", ".toml"}
                and "__pycache__" not in path.parts
            )
    return sorted(set(files))


def test_current_step_registry_matches_supported_runner_surface() -> None:
    expected = set(PER_SPECIES_STEPS) | {"integrate", "cstmm", "csfilter"}
    assert set(STEP_RUNNERS) == expected
    assert set(STEP_RUNNERS).issubset(set(CURRENT_AMALGKIT_STEPS))
    assert all(callable(runner) for runner in STEP_RUNNERS.values())


def test_retired_entrypoint_files_are_absent() -> None:
    """Deleted execution paths cannot silently return through a later refactor."""

    assert [path for path in RETIRED_PATHS if path.exists()] == []


def test_legacy_quant_cleanup_selects_only_unprovenanced_plain_outputs(tmp_path: Path) -> None:
    """Legacy archival must not select current, partial, or hidden quant work."""

    from scripts.rna.clean_external_artifacts import _legacy_quant_candidates

    quant_root = tmp_path / "apis_mellifera_all" / "work" / "quant"
    legacy = quant_root / "SRR_legacy"
    legacy.mkdir(parents=True)
    (legacy / "abundance.tsv").write_text("target_id\ttpm\n")

    current = quant_root / "SRR_current"
    current.mkdir()
    (current / "SRR_current_abundance.tsv").write_text("target_id\ttpm\n")
    (current / ".metainformant_quant_provenance.json").write_text("{}")

    partial = quant_root / "SRR_partial"
    partial.mkdir()
    (partial / "abundance.tsv.part").write_text("partial")

    hidden = quant_root / ".SRR_promote"
    hidden.mkdir()
    (hidden / "abundance.tsv").write_text("temporary")

    assert _legacy_quant_candidates(tmp_path, set()) == [legacy]


def test_legacy_quant_cleanup_is_idempotent_and_lock_guarded(tmp_path: Path) -> None:
    """Archival is reversible/idempotent and fails closed around active writers."""

    from scripts.rna.clean_external_artifacts import main as cleanup_main

    quant_root = tmp_path / "apis_mellifera_all" / "work" / "quant"
    legacy = quant_root / "SRR_legacy"
    legacy.mkdir(parents=True)
    (legacy / "abundance.tsv").write_text("target_id\ttpm\n")
    results = tmp_path / "results"
    results.mkdir()
    (results / ".full_campaign.lock").mkdir()

    with pytest.raises(SystemExit, match="active writer lock"):
        cleanup_main(
            [
                "--data-root",
                str(tmp_path),
                "--archive-legacy-quant",
                "--legacy-quant-only",
                "--execute",
                "--manifest",
                str(tmp_path / "blocked.json"),
            ]
        )

    (results / ".full_campaign.lock").rmdir()
    manifest = tmp_path / "archived.json"
    assert (
        cleanup_main(
            [
                "--data-root",
                str(tmp_path),
                "--archive-legacy-quant",
                "--legacy-quant-only",
                "--execute",
                "--manifest",
                str(manifest),
            ]
        )
        == 0
    )
    payload = json.loads(manifest.read_text())
    assert payload["candidate_count"] == 1
    assert not legacy.exists()

    second_manifest = tmp_path / "second.json"
    assert (
        cleanup_main(
            [
                "--data-root",
                str(tmp_path),
                "--archive-legacy-quant",
                "--legacy-quant-only",
                "--execute",
                "--manifest",
                str(second_manifest),
            ]
        )
        == 0
    )
    assert json.loads(second_manifest.read_text())["candidate_count"] == 0


def test_active_surface_contains_no_retired_api_references() -> None:
    """Current code and docs expose only the current producer contract."""

    violations = []
    for path in _active_text_files():
        match = RETIRED_REFERENCE.search(path.read_text(encoding="utf-8", errors="replace"))
        if match:
            violations.append(f"{path.relative_to(REPO_ROOT)}: {match.group(0)}")
    assert violations == []


def test_current_wrappers_are_callable() -> None:
    for name in ("wsfilter", "cstmm", "csfilter", "finalize", "busco", "rerun", "dataset"):
        assert callable(getattr(amalgkit, name))


def test_builder_filters_orchestration_values_by_command() -> None:
    params = {
        "out_dir": Path("work"),
        "threads": 4,
        "input_dir": Path("work/merge"),
        "norm": "log2p1-fpkm",
        "taxon_id": 7460,
        "tissue": ["brain"],
        "batch_effect_alg": "no",
    }

    metadata = amalgkit.build_cli_args(params, subcommand="metadata")
    assert "--out_dir" in metadata
    assert "--threads" in metadata
    assert "--input_dir" not in metadata
    assert "--taxon_id" not in metadata

    finalize = amalgkit.build_cli_args(params, subcommand="finalize")
    assert "--input_dir" in finalize
    assert "--norm" in finalize
    assert "--batch_effect_alg" in finalize
    assert "--taxon_id" not in finalize
    assert "--tissue" not in finalize


def test_cstmm_is_single_process_and_accepts_current_inputs() -> None:
    args = amalgkit.build_cli_args(
        {
            "out_dir": "cross_species/work",
            "threads": 8,
            "orthogroup_table": "orthogroups_transcript.tsv",
            "dir_count": "cross_species/merge_inputs",
            "tmm_backend": "python",
        },
        subcommand="cstmm",
    )
    assert "--threads" not in args
    assert "--orthogroup_table" in args
    assert "--dir_count" in args
    assert "--tmm_backend" in args


def test_current_commands_have_allow_lists() -> None:
    for step in PER_SPECIES_STEPS:
        assert step in SUPPORTED_CLI_OPTIONS
    assert SUPPORTED_CLI_OPTIONS["finalize"] >= {"out_dir", "metadata", "input_dir", "norm"}
    assert SUPPORTED_CLI_OPTIONS["wsfilter"] >= {"out_dir", "metadata", "input_dir", "norm"}
    assert SUPPORTED_CLI_OPTIONS["busco"] >= {"out_dir", "metadata", "species"}
    assert SUPPORTED_CLI_OPTIONS["rerun"] >= {"out_dir", "metadata"}
    assert SUPPORTED_CLI_OPTIONS["dataset"] >= {"out_dir", "name", "rule_set"}
    assert SUPPORTED_CLI_OPTIONS["select"] >= {"select_rules_tsv", "random_seed"}
    assert SUPPORTED_CLI_OPTIONS["getfastq"] >= {
        "sra_download_wait_timeout_seconds",
        "sra_download_transfer_timeout_seconds",
    }


def test_version_parser_accepts_installed_output() -> None:
    valid, message = amalgkit.parse_and_check_version("amalgkit version 0.16.60", "0.16.60")
    assert valid, message
    valid, message = amalgkit.parse_and_check_version("amalgkit version 0.16.60", "0.16.60", exact=True)
    assert valid, message
    invalid, _ = amalgkit.parse_and_check_version("amalgkit version 0.16.99", "0.16.60", exact=True)
    assert not invalid
    invalid, _ = amalgkit.parse_and_check_version("amalgkit development build", "0.16.60", exact=True)
    assert not invalid
    assert amalgkit.REQUIRED_AMALGKIT_VERSION == "0.16.60"
    assert amalgkit.AMALGKIT_RELEASE_TAG == "v0.16.60"
    assert amalgkit.AMALGKIT_SOURCE_REVISION == "c656a52aacdcee6fd3bf7e8031769ca957204ebc"


def test_project_selection_rules_use_current_parser() -> None:
    from amalgkit.select import read_select_config

    rules_path = Path(__file__).parents[2] / "config" / "amalgkit" / "select_rules.tsv"
    parsed = read_select_config(str(rules_path))
    assert set(parsed["parameters"]) >= {
        "min_nspots",
        "max_sample",
        "mark_missing_rank",
        "mark_redundant_biosamples",
    }


def test_command_prefix_is_current() -> None:
    command = amalgkit.build_amalgkit_command(
        "finalize",
        {"out_dir": "work", "metadata": "work/wsfilter/metadata.tsv", "input_dir": "work/wsfilter"},
    )
    assert command[:2] == ["amalgkit", "finalize"]
    assert "--input_dir" in command


def test_downstream_provenance_requires_explicit_steps(tmp_path: Path) -> None:
    """A contract-valid sidecar without stages is not completion evidence."""

    work_dir = tmp_path / "work"
    write_downstream_provenance(
        work_dir,
        species="species_a",
        config_path=tmp_path / "config.yaml",
        quantified_samples=2,
        steps=["merge", "finalize"],
    )
    assert is_current_downstream(work_dir)
    assert is_current_downstream(work_dir, {"finalize"})
    assert not is_current_downstream(work_dir, {"sanity"})

    sidecar = work_dir / ".metainformant_downstream_provenance.json"
    payload = json.loads(sidecar.read_text(encoding="utf-8"))
    payload.pop("steps")
    sidecar.write_text(json.dumps(payload), encoding="utf-8")
    assert not is_current_downstream(work_dir)


def _quant_provenance_sidecar_is_byte_idempotent(tmp_path: Path) -> None:
    """Two writes over unchanged facts produce a byte-identical sidecar.

    Quantification sidecars are re-verified fail-closed on every orchestrator
    resume.  (Companion regression to
    ``test_reference_alias_manifest_is_byte_idempotent``.)
    """
    import time

    config = tmp_path / "config.yaml"
    config.write_bytes(b"species: test")
    sample_dir = tmp_path / "quant" / "SRRX"
    sample_dir.mkdir(parents=True)
    abundance = sample_dir / "abundance.tsv"
    abundance.write_bytes(b"target_id\teff_counts\nENST1\t10\n")
    command = ["amalgkit", "quant", "--sample", "SRRX"]

    write_quant_provenance(
        sample_dir,
        species="test_species",
        run_accession="SRRX",
        config_path=config,
        command=command,
        quantification_file=abundance.name,
    )
    first_bytes = (sample_dir / ".metainformant_quant_provenance.json").read_bytes()
    time.sleep(1.1)
    write_quant_provenance(
        sample_dir,
        species="test_species",
        run_accession="SRRX",
        config_path=config,
        command=command,
        quantification_file=abundance.name,
    )
    assert (sample_dir / ".metainformant_quant_provenance.json").read_bytes() == first_bytes


def test_downstream_provenance_sidecar_is_byte_idempotent(tmp_path: Path) -> None:
    """Two writes over unchanged facts produce a byte-identical sidecar.

    Downstream checkpoints are replayed from disk on resume to skip completed
    stages; any restart-varying byte would make an unchanged checkpoint look
    rewritten.  (Companion regression to
    ``test_reference_alias_manifest_is_byte_idempotent``.)
    """
    import time

    work_dir = tmp_path / "work"
    config = tmp_path / "config.yaml"
    config.write_bytes(b"species: test")

    write_downstream_provenance(
        work_dir,
        species="test_species",
        config_path=config,
        quantified_samples=2,
        steps=["merge", "finalize"],
    )
    first_bytes = (work_dir / ".metainformant_downstream_provenance.json").read_bytes()
    time.sleep(1.1)
    write_downstream_provenance(
        work_dir,
        species="test_species",
        config_path=config,
        quantified_samples=2,
        steps=["merge", "finalize"],
    )
    assert (work_dir / ".metainformant_downstream_provenance.json").read_bytes() == first_bytes


def test_metadata_provenance_sidecar_is_byte_idempotent(tmp_path: Path) -> None:
    """Two writes over unchanged metadata produce a byte-identical sidecar.

    (Companion regression to
    ``test_reference_alias_manifest_is_byte_idempotent``.)
    """
    import time

    work_dir = tmp_path / "work"
    metadata_dir = work_dir / "metadata"
    metadata_dir.mkdir(parents=True)
    (metadata_dir / "metadata.tsv").write_bytes(b"srr\tspecies\nSRR1\ttest\n")
    (metadata_dir / "metadata_selected.tsv").write_bytes(b"srr\tspecies\nSRR1\ttest\n")
    config = tmp_path / "config.yaml"
    config.write_bytes(b"species: test")
    rules = tmp_path / "rules.txt"
    rules.write_bytes(b"keep all")

    write_metadata_provenance(
        work_dir,
        species="test_species",
        config_path=config,
        selection_rules_path=rules,
    )
    first_bytes = (metadata_dir / ".metainformant_metadata_provenance.json").read_bytes()
    time.sleep(1.1)
    write_metadata_provenance(
        work_dir,
        species="test_species",
        config_path=config,
        selection_rules_path=rules,
    )
    assert (metadata_dir / ".metainformant_metadata_provenance.json").read_bytes() == first_bytes
