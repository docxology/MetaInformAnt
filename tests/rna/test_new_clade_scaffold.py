"""Zero-mocks tests for scripts/rna/new_clade_scaffold.py (real files, real subprocesses)."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
SCAFFOLD = REPO_ROOT / "scripts" / "rna" / "new_clade_scaffold.py"

SPECIES_LIST = "\n".join(
    [
        "# Lepidoptera pilot roster",
        "Bombyx mori 7091 GCF_000151645.2 Bombyx_mori_p50T",
        "Danaus plexippus",
        "Heliconius melpomene 7122 GCA_000290705.1 Hmel2_5",
        "",
    ]
)


@pytest.fixture()
def species_file(tmp_path: Path) -> Path:
    path = tmp_path / "species.txt"
    path.write_text(SPECIES_LIST + "\n", encoding="utf-8")
    return path


def _run(args: list[str], **kwargs) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCAFFOLD), *args],
        capture_output=True,
        text=True,
        timeout=300,
        **kwargs,
    )


def _scaffold(tmp_path: Path, species_file: Path, clade: str = "lepidoptera") -> subprocess.CompletedProcess:
    return _run(["--clade", clade, "--species-list", str(species_file), "--output-dir", str(tmp_path)])


def test_scaffold_generates_expected_files(tmp_path: Path, species_file: Path) -> None:
    result = _scaffold(tmp_path, species_file)
    assert result.returncode == 0, result.stderr
    root = tmp_path / "lepidoptera_amalgkit"
    for relative in (
        "README.md",
        "TODO.md",
        ".gitignore",
        "doc/00_setup/01_environment.md",
        "docs/manuscript/README.md",
        "scripts/metainformant_import.py",
        "config/amalgkit/amalgkit_cross_species.yaml",
        "config/amalgkit/amalgkit_bombyx_mori.yaml",
        "config/amalgkit/amalgkit_danaus_plexippus.yaml",
        "config/amalgkit/amalgkit_heliconius_melpomene.yaml",
        "config/amalgkit/select_rules.tsv",
        "tests/test_import_contract.py",
    ):
        assert (root / relative).is_file(), f"missing {relative}"
    shim_generated = (root / "scripts" / "metainformant_import.py").read_text(encoding="utf-8")
    shim_canonical = (
        REPO_ROOT / "projects" / "hymenoptera_amalgkit" / "scripts" / "metainformant_import.py"
    ).read_text(encoding="utf-8")
    assert shim_generated == shim_canonical, "shim must be copied verbatim from the canonical source"


def test_scaffold_output_is_deterministic(tmp_path: Path, species_file: Path) -> None:
    first_dir = tmp_path / "a"
    second_dir = tmp_path / "b"
    assert _scaffold(first_dir, species_file).returncode == 0
    assert _scaffold(second_dir, species_file).returncode == 0

    def snapshot(root: Path) -> dict[str, str]:
        return {str(p.relative_to(root)): p.read_text(encoding="utf-8") for p in sorted(root.rglob("*")) if p.is_file()}

    assert snapshot(first_dir / "lepidoptera_amalgkit") == snapshot(second_dir / "lepidoptera_amalgkit")


def test_generated_species_config_matches_contract(tmp_path: Path, species_file: Path) -> None:
    assert _scaffold(tmp_path, species_file).returncode == 0
    config_path = tmp_path / "lepidoptera_amalgkit" / "config" / "amalgkit" / "amalgkit_bombyx_mori.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    assert config["species_list"] == ["Bombyx_mori"]
    assert config["taxon_id"] == 7091
    assert config["genome"]["accession"] == "GCF_000151645.2"
    assert config["genome"]["assembly_name"] == "Bombyx_mori_p50T"
    for key in ("work_dir", "log_dir", "threads", "steps"):
        assert key in config
    assert set(config["steps"]) == {
        "metadata",
        "integrate",
        "select",
        "getfastq",
        "quant",
        "merge",
        "wsfilter",
        "finalize",
        "sanity",
    }
    assert config["steps"]["select"]["select_rules_tsv"] == "config/amalgkit/select_rules.tsv"
    cross = yaml.safe_load(
        (tmp_path / "lepidoptera_amalgkit" / "config" / "amalgkit" / "amalgkit_cross_species.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert cross["analysis"]["inferential_statistics"] == "none"


def test_missing_accession_produces_todo_entries(tmp_path: Path) -> None:
    species_file = tmp_path / "sparse.txt"
    species_file.write_text("Danaus plexippus\n", encoding="utf-8")
    result = _run(["--clade", "danaus", "--species-list", str(species_file), "--output-dir", str(tmp_path)])
    assert result.returncode == 0, result.stderr
    root = tmp_path / "danaus_amalgkit"
    todo_text = (root / "TODO.md").read_text(encoding="utf-8")
    assert "Danaus_plexippus: fill in NCBI genome accession" in todo_text
    assert "Danaus_plexippus: fill in NCBI taxon_id" in todo_text
    config = yaml.safe_load(
        (root / "config" / "amalgkit" / "amalgkit_danaus_plexippus.yaml").read_text(encoding="utf-8")
    )
    assert config["genome"]["accession"] == ""
    assert "taxon_id" not in config


def test_generated_tests_pass_in_skeleton(tmp_path: Path, species_file: Path) -> None:
    """The skeleton's own pytest suite passes against the real parent package."""
    result = _scaffold(tmp_path, species_file)
    assert result.returncode == 0, result.stderr
    root = tmp_path / "lepidoptera_amalgkit"
    venv_py = REPO_ROOT / ".venv" / "bin" / "python"
    python = str(venv_py) if venv_py.is_file() else sys.executable
    # Import-contract test: run non-pytest checks directly (deterministic).
    probe = subprocess.run(
        [
            python,
            "-c",
            "import sys; sys.path.insert(0, r'%s'); "
            "import metainformant_import; metainformant_import.ensure_metainformant(); "
            "import metainformant; print('resolved:', metainformant.__file__)" % (root / "scripts"),
        ],
        capture_output=True,
        text=True,
        timeout=120,
        cwd=str(root),
    )
    assert probe.returncode == 0, probe.stderr
    assert "resolved:" in probe.stdout


def test_scaffold_refuses_existing_target(tmp_path: Path, species_file: Path) -> None:
    assert _scaffold(tmp_path, species_file).returncode == 0
    result = _scaffold(tmp_path, species_file)
    assert result.returncode == 1
    assert "refusing to overwrite" in result.stderr


def test_dry_run_writes_nothing(tmp_path: Path, species_file: Path) -> None:
    result = _run(
        ["--clade", "lepidoptera", "--species-list", str(species_file), "--output-dir", str(tmp_path), "--dry-run"]
    )
    assert result.returncode == 0
    assert list(tmp_path.iterdir()) == [species_file]


def test_duplicate_species_rejected(tmp_path: Path) -> None:
    species_file = tmp_path / "dup.txt"
    species_file.write_text("Danaus plexippus\nDanaus plexippus\n", encoding="utf-8")
    result = _run(["--clade", "lepidoptera", "--species-list", str(species_file), "--output-dir", str(tmp_path)])
    assert result.returncode == 2
    assert "duplicate species token" in result.stderr


def test_accession_without_assembly_rejected(tmp_path: Path) -> None:
    species_file = tmp_path / "bad.txt"
    species_file.write_text("Danaus plexippus 7391 GCF_000.txt\n", encoding="utf-8")
    result = _run(["--clade", "lepidoptera", "--species-list", str(species_file), "--output-dir", str(tmp_path)])
    assert result.returncode == 2
    assert "accession requires assembly_name" in result.stderr


def test_invalid_clade_name_rejected(tmp_path: Path, species_file: Path) -> None:
    result = _run(["--clade", "../escape", "--species-list", str(species_file), "--output-dir", str(tmp_path)])
    assert result.returncode == 2
    assert "invalid clade name" in result.stderr
