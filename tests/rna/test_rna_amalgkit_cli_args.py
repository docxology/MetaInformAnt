"""Focused tests for the current Amalgkit CLI argument contract."""

from __future__ import annotations

from pathlib import Path

from metainformant.rna.amalgkit.amalgkit import AmalgkitParams, build_cli_args


def test_builder_normalizes_common_types_without_a_subcommand() -> None:
    args = build_cli_args(
        {
            "threads": 8,
            "dry_run": True,
            "species_list": ["A", "B"],
            "output": Path("output/x"),
            "skip": False,
            "none_val": None,
        }
    )
    assert "--threads" in args and "8" in args
    assert "--dry_run" in args
    assert args.count("--species") == 2
    assert "output/x" in args


def test_builder_preserves_current_metadata_options() -> None:
    args = build_cli_args(
        {"out_dir": "work", "search_string": "RNA-Seq", "resolve_names": True},
        subcommand="metadata",
    )
    assert "--search_string" in args
    assert "--resolve_names" in args
    assert args[args.index("--resolve_names") + 1] == "yes"


def test_builder_removes_orchestration_only_values_from_select() -> None:
    args = build_cli_args(
        {
            "out_dir": "work",
            "metadata": "work/metadata/metadata.tsv",
            "taxon_id": 7460,
            "tissue": ["brain"],
            "threads": 4,
        },
        subcommand="select",
    )
    assert "--metadata" in args
    assert "--threads" in args
    assert "--taxon_id" not in args
    assert "--tissue" not in args


def test_structured_params_emit_explicit_yes_no_values() -> None:
    args = build_cli_args(
        AmalgkitParams("work", threads=2, redo=False, fasterq_size_check=False),
        subcommand="getfastq",
    )
    assert args[args.index("--redo") + 1] == "no"
    assert args[args.index("--fasterq_size_check") + 1] == "no"
