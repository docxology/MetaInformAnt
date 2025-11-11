"""Tests for CLI argument building in metainformant.rna.amalgkit.

Tests parameter normalization, flag ordering, type conversion, and v0.12.20 features.
"""

from pathlib import Path

from metainformant.rna.amalgkit import build_amalgkit_command, build_cli_args


def test_build_cli_args_basic():
    """Test basic CLI argument building with various parameter types."""
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
    # flags normalized and ordered per input iteration
    assert "--threads" in args and "8" in args
    assert "--dry-run" in args
    # repeat flags for lists
    assert args.count("--species-list") == 2
    # path stringified
    assert "output/x" in args


def test_build_cli_args_resolve_names_string():
    """Test resolve_names parameter as string (v0.12.20 feature)."""
    args = build_cli_args({"resolve_names": "yes", "search_string": "test"}, for_cli=True)
    assert "--resolve_names" in args
    assert "yes" in args
    # Should be flag + value pair
    resolve_idx = args.index("--resolve_names")
    assert args[resolve_idx + 1] == "yes"


def test_build_cli_args_resolve_names_boolean():
    """Test resolve_names parameter as boolean (v0.12.20 feature)."""
    args = build_cli_args({"resolve_names": True, "search_string": "test"}, for_cli=True)
    assert "--resolve_names" in args
    assert "yes" in args
    # Boolean True should convert to "yes"
    resolve_idx = args.index("--resolve_names")
    assert args[resolve_idx + 1] == "yes"


def test_build_cli_args_mark_missing_rank():
    """Test mark_missing_rank parameter (v0.12.20 feature)."""
    args = build_cli_args({"mark_missing_rank": "species", "out_dir": "test"}, for_cli=True)
    assert "--mark_missing_rank" in args
    assert "species" in args
    # Should be flag + value pair
    mark_idx = args.index("--mark_missing_rank")
    assert args[mark_idx + 1] == "species"


def test_build_amalgkit_command_with_v0_12_20_features():
    """Test command building with v0.12.20 features."""
    # Test metadata with resolve_names
    cmd1 = build_amalgkit_command("metadata", {"resolve_names": "yes", "search_string": "test"})
    assert "--resolve_names" in cmd1
    assert "yes" in cmd1
    
    # Test select with mark_missing_rank
    cmd2 = build_amalgkit_command("select", {"mark_missing_rank": "genus", "out_dir": "test"})
    assert "--mark_missing_rank" in cmd2
    assert "genus" in cmd2
    
    # Test select with all new parameters
    cmd3 = build_amalgkit_command(
        "select",
        {
            "mark_missing_rank": "species",
            "min_nspots": 10000000,
            "max_sample": 50,
            "mark_redundant_biosamples": "yes",
        },
    )
    assert "--mark_missing_rank" in cmd3
    assert "--min_nspots" in cmd3
    assert "--max_sample" in cmd3
    assert "--mark_redundant_biosamples" in cmd3
