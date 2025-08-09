from pathlib import Path

from metainformant.rna.amalgkit import build_cli_args


def test_build_cli_args_basic():
    args = build_cli_args({
        "threads": 8,
        "dry_run": True,
        "species_list": ["A", "B"],
        "output": Path("output/x"),
        "skip": False,
        "none_val": None,
    })
    # flags normalized and ordered per input iteration
    assert "--threads" in args and "8" in args
    assert "--dry-run" in args
    # repeat flags for lists
    assert args.count("--species-list") == 2
    # path stringified
    assert "output/x" in args


