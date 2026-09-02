#!/usr/bin/env python3
"""Scaffold a complete sibling clade project from the hymenoptera template pattern.

Given a clade name and a species list, generate a standalone
``<clade>_amalgkit`` project skeleton (config species YAMLs, import-contract
shim, minimal docs, TODO) that follows ``docs/rna/CLADE_REPLICATION.md`` and
validates with the same gates where possible. Output is deterministic: the
same inputs always produce byte-identical files (sorted species, no
timestamps, no machine-specific paths).

Species-list format (one per line, ``#`` comments allowed)::

    Genus species [taxon_id [accession assembly_name]]

Examples::

    Bombyx mori 7091 GCF_000151645.2 Bombyx_mori_p50T
    Danaus plexippus

``taxon_id``/``accession``/``assembly_name`` are optional; when an accession
is missing the emitted YAML carries an empty ``genome.accession`` plus a
``TODO`` entry so the resulting project fails fast at its own validation
gates instead of silently downloading nothing.

Usage::

    python scripts/rna/new_clade_scaffold.py \
        --clade lepidoptera --species-list species.txt \
        --output-dir /tmp/clades
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
HYMENOPTERA_SCRIPTS = REPO_ROOT / "projects" / "hymenoptera_amalgkit" / "scripts"
PARENT_CONFIG_DIR = REPO_ROOT / "config" / "amalgkit"

IMPORT_SHIM_NAME = "metainformant_import.py"
SHARED_CONFIG_ASSETS = ("select_rules.tsv", "tissue_mapping.yaml", "tissue_patches.yaml")
NAME_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_-]*$")


@dataclass(frozen=True)
class SpeciesEntry:
    """One species row parsed from the species-list file."""

    genus: str
    epithet: str
    taxon_id: str
    accession: str
    assembly_name: str

    @property
    def token(self) -> str:
        """Genus_species token used in configs and manifests."""
        return f"{self.genus}_{self.epithet}"

    @property
    def display(self) -> str:
        """Display name with the space separator."""
        return f"{self.genus} {self.epithet}"

    @property
    def config_name(self) -> str:
        """amalgkit_<genus_species>.yaml config filename."""
        return f"amalgkit_{self.token.lower()}.yaml"


def parse_species_list(path: Path) -> list[SpeciesEntry]:
    """Parse the species-list file into entries; raise on malformed input."""
    entries: list[SpeciesEntry] = []
    seen: set[str] = set()
    for lineno, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        fields = line.split()
        if len(fields) < 2 or len(fields) > 5:
            raise ValueError(
                f"{path.name}:{lineno}: expected 'Genus species [taxon_id [accession assembly_name]]', got: {line!r}"
            )
        genus, epithet = fields[0], fields[1]
        taxon_id = fields[2] if len(fields) >= 3 else ""
        accession = fields[3] if len(fields) >= 4 else ""
        assembly_name = fields[4] if len(fields) >= 5 else ""
        if taxon_id and not taxon_id.isdigit():
            raise ValueError(f"{path.name}:{lineno}: taxon_id must be numeric: {line!r}")
        if accession and not assembly_name:
            raise ValueError(f"{path.name}:{lineno}: accession requires assembly_name: {line!r}")
        entry = SpeciesEntry(genus, epithet, taxon_id, accession, assembly_name)
        if entry.token in seen:
            raise ValueError(f"{path.name}:{lineno}: duplicate species token: {entry.token}")
        seen.add(entry.token)
        entries.append(entry)
    if not entries:
        raise ValueError(f"{path.name}: no species entries found")
    return entries


def validate_clade_name(clade: str) -> str:
    """Validate the clade token used as project-name prefix."""
    if not NAME_RE.match(clade):
        raise ValueError(f"invalid clade name {clade!r}: use letters, digits, hyphen or underscore (e.g. lepidoptera)")
    return clade


def _species_yaml(entry: SpeciesEntry, clade: str) -> str:
    """Render one species config YAML following the apis_mellifera pattern."""
    token = entry.token
    base = f"output/amalgkit/{token.lower()}"
    genome_lines = [
        "genome:",
        *([f'  accession: "{entry.accession}"'] if not entry.accession else [f"  accession: {entry.accession}"]),
        *(
            [f'  assembly_name: "{entry.assembly_name}"']
            if not entry.assembly_name
            else [f"  assembly_name: {entry.assembly_name}"]
        ),
        f"  dest_dir: {base}/genome",
        "  include:",
    ]
    for include in ("genome", "gff3", "rna", "cds", "protein"):
        genome_lines.append(f"    - {include}")
    if not entry.accession:
        genome_lines.append("  # TODO(new-clade-scaffold): fill in the NCBI accession and assembly_name")
    header_lines = [
        "# METAINFORMANT Amalgkit Configuration",
        f"# Species: {entry.display}",
        f"# Clade: {clade}",
    ]
    if entry.taxon_id:
        header_lines.append(f"# NCBI Taxonomy ID: {entry.taxon_id}")
    else:
        header_lines.append("# NCBI Taxonomy ID: TODO(new-clade-scaffold)")
    header_lines.append(
        f"# Assembly: {entry.accession} ({entry.assembly_name})"
        if entry.accession
        else "# Assembly: TODO(new-clade-scaffold)"
    )
    body = [
        "",
        f"work_dir: {base}/work",
        f"log_dir: {base}/logs",
        "threads: 16",
        "",
        "auto_install_amalgkit: true",
        "",
        "filters:",
        "  require_tissue: false",
        "",
        "species_list:",
        f"  - {token}",
        "",
    ]
    if entry.taxon_id:
        body += [f"taxon_id: {entry.taxon_id}", ""]
    body += genome_lines + ["", "steps:"]

    steps: list[tuple[str, list[tuple[str, str]]]] = [
        (
            "metadata",
            [
                (
                    "search_string",
                    f"'\"{entry.display}\"[Organism] AND RNA-Seq[Strategy] AND Illumina[Platform]'",
                ),
                ("redo", "no"),
            ],
        ),
        ("integrate", [("fastq_dir", f"{base}/fastq/getfastq")]),
        ("select", [("select_rules_tsv", "config/amalgkit/select_rules.tsv")]),
        ("getfastq", [("remove_sra", "yes"), ("threads", "16"), ("redo", "no")]),
        ("quant", [("threads", "16"), ("redo", "no"), ("clean_fastq", "no")]),
        ("merge", [("out_dir", f"{base}/work")]),
        (
            "wsfilter",
            [
                ("out_dir", f"{base}/work"),
                ("input_dir", f"{base}/work/merge"),
                ("norm", "log2p1-fpkm"),
                ("mapping_rate", "0.2"),
            ],
        ),
        (
            "finalize",
            [
                ("out_dir", f"{base}/work"),
                ("input_dir", f"{base}/work/wsfilter"),
                ("norm", "log2p1-fpkm"),
                ("batch_effect_alg", "no"),
            ],
        ),
        ("sanity", []),
    ]
    for step, params in steps:
        body.append(f"  {step}:")
        for key, value in params:
            body.append(f"    {key}: {value}")
        if not params:
            body.append("    {}")
    return "\n".join(header_lines) + "\n" + "\n".join(body) + "\n"


def _cross_species_yaml() -> str:
    """Render the cross-species config (descriptive-only, inferential_stats: none)."""
    return """# Cross-species expression-fingerprint analysis configuration
# (clade scaffold; descriptive-only until the evidence manifest freezes)

work_dir: cross_species/work
log_dir: cross_species/logs
threads: 4

cohort:
  source: config/amalgkit/amalgkit_<species>.yaml
  minimum_samples: 2
  required_matrix_suffix: _expression.tsv

analysis:
  method: expression_fingerprint
  transform: log1p
  pooled_percentile_range: [1, 99]
  fingerprint_bins: 100
  minimum_finite_positive_genes: 100
  distance: 1 - spearman_rho
  clustering: average_linkage
  inferential_statistics: none

cstmm:
  enabled: false
  orthogroup_table: ""
  dir_count: cross_species/merge_inputs

csfilter:
  enabled: false
  input_dir: cross_species/work/cstmm
  orthogroup_table: ""
"""


def _readme(clade: str, entries: list[SpeciesEntry]) -> str:
    """Render the project README."""
    project = f"{clade}_amalgkit"
    lines = [
        f"# {project}",
        "",
        f"Standalone {clade} RNA-seq campaign project built on the MetaInformAnt",
        "platform. Generated by `scripts/rna/new_clade_scaffold.py` from the",
        "hymenoptera_amalgkit template pattern; see",
        "`docs/rna/CLADE_REPLICATION.md` in the parent checkout for the full",
        "replication contract.",
        "",
        "## Species roster",
        "",
        "| Token | Taxon ID | Accession | Assembly |",
        "| ----- | -------- | --------- | -------- |",
    ]
    for entry in sorted(entries, key=lambda e: e.token):
        lines.append(
            f"| {entry.token} | {entry.taxon_id or 'TODO'} | {entry.accession or 'TODO'}"
            f" | {entry.assembly_name or 'TODO'} |"
        )
    lines += [
        "",
        "## Import contract (non-negotiable)",
        "",
        "```python",
        "from metainformant_import import ensure_metainformant",
        "",
        "ensure_metainformant()",
        "```",
        "",
        "## Gates before first campaign",
        "",
        "```bash",
        "uv run python scripts/validate_configs.py --config-dir config/amalgkit",
        "uv run python scripts/validate_project_docs.py",
        "```",
        "",
    ]
    return "\n".join(lines)


def _todo(clade: str, entries: list[SpeciesEntry]) -> str:
    """Render the project TODO with the concrete gaps left by the scaffold."""
    project = f"{clade}_amalgkit"
    gaps: list[str] = []
    for entry in sorted(entries, key=lambda e: e.token):
        if not entry.accession:
            gaps.append(f"- [ ] {entry.token}: fill in NCBI genome accession + assembly_name in `{entry.config_name}`")
        if not entry.taxon_id:
            gaps.append(f"- [ ] {entry.token}: fill in NCBI taxon_id in `{entry.config_name}`")
    lines = [
        f"# {project} TODO",
        "",
        "## Scaffold gaps (must close before first campaign)",
        "",
    ]
    lines += gaps if gaps else ["- [x] No scaffold gaps: every species carries taxon_id, accession and assembly_name."]
    lines += [
        "",
        "## Standard campaign checklist",
        "",
        "- [ ] Export `AMALGKIT_DATA_ROOT` to an empty external data root.",
        "- [ ] Run `scripts/validate_configs.py --config-dir config/amalgkit` clean.",
        "- [ ] Run the parent docs gate `scripts/verify_documentation_code.py` where applicable.",
        "- [ ] Initialize the species manifest and verify tissue coverage.",
        "- [ ] Run `tests/` clean (one pytest directory per invocation).",
        "",
    ]
    return "\n".join(lines)


def _setup_doc(clade: str) -> str:
    """Render the minimal setup doc (no machine-specific absolute paths)."""
    return f"""# {clade}_amalgkit environment setup

1. Parent checkout: `uv sync` in the MetaInformAnt repository, then verify
   `import metainformant` succeeds from any working directory.
2. Data root: create an empty external data directory (never inside either
   repository) and export it as `AMALGKIT_DATA_ROOT`.
3. Import contract: every script/test imports the parent package only through
   `scripts/{IMPORT_SHIM_NAME}` (`ensure_metainformant()`); hand-rolled
   `sys.path.insert` hacks are contract violations.
4. Statistics boundary: cross-species outputs are descriptive-only until the
   evidence manifest freezes; inferential tests stay gated and labeled.
"""


def _manuscript_readme(clade: str) -> str:
    """Render the manuscript placeholder pointing at the campaign runbook."""
    return f"""# {clade}_amalgkit manuscript

Placeholder scaffold. The hymenoptera campaign's full manuscript contract
(`docs/manuscript/` section set, evidence-boundary phrases, claim ledger) is
project-specific; port it deliberately when this clade reaches its own
freeze. Until then, this directory must not claim executed analyses.
"""


def _gitignore() -> str:
    """Render the local-only output hygiene rules."""
    return """output/
logs/
__pycache__/
*.pyc
.pytest_cache/
"""


def _import_shim() -> str:
    """Read the canonical import shim verbatim from the hymenoptera project."""
    shim = HYMENOPTERA_SCRIPTS / IMPORT_SHIM_NAME
    if not shim.is_file():
        raise FileNotFoundError(
            "canonical import shim not found at " f"{shim}; the scaffold must copy it verbatim (single source of truth)"
        )
    return shim.read_text(encoding="utf-8")


def _write_file(path: Path, content: str) -> Path:
    """Write one file deterministically (parents created, newline-terminated)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if not content.endswith("\n"):
        content += "\n"
    path.write_text(content, encoding="utf-8")
    return path


def scaffold(clade: str, species_list: Path, output_dir: Path) -> list[Path]:
    """Generate the clade project skeleton; return the list of files written."""
    clade = validate_clade_name(clade)
    entries = parse_species_list(species_list)
    project_root = output_dir / f"{clade}_amalgkit"
    if project_root.exists():
        raise FileExistsError(f"refusing to overwrite existing target: {project_root}")

    written: list[Path] = []
    config_dir = project_root / "config" / "amalgkit"
    scripts_dir = project_root / "scripts"
    tests_dir = project_root / "tests"
    doc_dir = project_root / "doc" / "00_setup"
    docs_dir = project_root / "docs" / "manuscript"

    written.append(_write_file(config_dir / "amalgkit_cross_species.yaml", _cross_species_yaml()))
    for asset in SHARED_CONFIG_ASSETS:
        source = PARENT_CONFIG_DIR / asset
        if source.is_file():
            shutil.copyfile(source, config_dir / asset)
            written.append(config_dir / asset)

    for entry in sorted(entries, key=lambda e: e.token):
        written.append(_write_file(config_dir / entry.config_name, _species_yaml(entry, clade)))

    written.append(_write_file(scripts_dir / IMPORT_SHIM_NAME, _import_shim()))
    validate_configs_source = REPO_ROOT / "scripts" / "rna" / "validate_configs.py"
    if validate_configs_source.is_file():
        written.append(
            _write_file(scripts_dir / "validate_configs.py", validate_configs_source.read_text(encoding="utf-8"))
        )

    written.append(_write_file(project_root / "README.md", _readme(clade, entries)))
    written.append(_write_file(project_root / "TODO.md", _todo(clade, entries)))
    written.append(_write_file(doc_dir / "01_environment.md", _setup_doc(clade)))
    written.append(_write_file(docs_dir / "README.md", _manuscript_readme(clade)))
    written.append(_write_file(project_root / ".gitignore", _gitignore()))
    written.append(_write_file(tests_dir / "test_import_contract.py", TEST_IMPORT_CONTRACT))

    structure_errors = verify_structure(project_root, clade, entries)
    if structure_errors:
        raise RuntimeError("scaffold self-verification failed: " + "; ".join(structure_errors))
    return sorted(written)


def verify_structure(project_root: Path, clade: str, entries: list[SpeciesEntry]) -> list[str]:
    """Re-derive the scaffold contract; return violations (empty list = pass)."""
    errors: list[str] = []
    config_dir = project_root / "config" / "amalgkit"
    for entry in sorted(entries, key=lambda e: e.token):
        path = config_dir / entry.config_name
        if not path.is_file():
            errors.append(f"missing species config: {entry.config_name}")
            continue
        config = yaml.safe_load(path.read_text(encoding="utf-8"))
        for key in ("work_dir", "species_list", "genome", "steps", "log_dir", "threads"):
            if key not in config:
                errors.append(f"{entry.config_name}: missing required top-level key: {key}")
        steps = config.get("steps") or {}
        if not isinstance(steps, dict) or not steps:
            errors.append(f"{entry.config_name}: steps section must be a non-empty mapping")
        genome = config.get("genome") or {}
        if not isinstance(genome, dict) or "accession" not in genome:
            errors.append(f"{entry.config_name}: genome.accession key missing")
        # NOTE: an empty accession with a TODO entry is an intentional scaffold
        # gap, not a structure violation; TODO.md records it.
    shim = project_root / "scripts" / IMPORT_SHIM_NAME
    if not shim.is_file() or "ensure_metainformant" not in shim.read_text(encoding="utf-8"):
        errors.append(f"missing import-contract shim: scripts/{IMPORT_SHIM_NAME}")
    for relative in (
        "README.md",
        "TODO.md",
        "doc/00_setup/01_environment.md",
        "docs/manuscript/README.md",
    ):
        if not (project_root / relative).is_file():
            errors.append(f"missing doc scaffold: {relative}")
    for text_path in (project_root / "doc").rglob("*.md"):
        text = text_path.read_text(encoding="utf-8")
        if "/Volumes/" in text or "/home/" in text:
            errors.append(f"machine-specific absolute path in {text_path.relative_to(project_root)}")
    return errors


TEST_IMPORT_CONTRACT = '''"""Real-implementation checks for the generated clade project skeleton."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SHIM = PROJECT_ROOT / "scripts" / "metainformant_import.py"


def test_shim_exists_and_defines_contract() -> None:
    assert SHIM.is_file()
    text = SHIM.read_text(encoding="utf-8")
    assert "def ensure_metainformant" in text


def test_shim_resolves_parent_package() -> None:
    bootstrap = (
        "import metainformant_import; metainformant_import.ensure_metainformant(); "
        "import metainformant; print(metainformant.__file__)"
    )
    result = subprocess.run(
        [sys.executable, "-c", bootstrap],
        capture_output=True,
        text=True,
        timeout=120,
        cwd=str(SHIM.parent),
    )
    assert result.returncode == 0, result.stderr
    assert "metainformant" in result.stdout


def test_species_configs_are_wellformed() -> None:
    import yaml

    configs = sorted((PROJECT_ROOT / "config" / "amalgkit").glob("amalgkit_*.yaml"))
    species = [
        p
        for p in configs
        if p.name
        not in {
            "amalgkit_template.yaml",
            "amalgkit_test.yaml",
            "amalgkit_cross_species.yaml",
        }
    ]
    assert species, "scaffold produced no species configs"
    for path in species:
        config = yaml.safe_load(path.read_text(encoding="utf-8"))
        for key in ("work_dir", "species_list", "genome", "steps"):
            assert key in config, path.name
        assert config["species_list"], path.name
'''


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--clade", required=True, help="clade token, e.g. lepidoptera -> lepidoptera_amalgkit")
    parser.add_argument("--species-list", required=True, type=Path, help="species list file (see module docstring)")
    parser.add_argument("--output-dir", required=True, type=Path, help="directory to create <clade>_amalgkit under")
    parser.add_argument("--dry-run", action="store_true", help="validate inputs without writing")
    args = parser.parse_args(argv)

    try:
        clade = validate_clade_name(args.clade)
        entries = parse_species_list(args.species_list)
    except (ValueError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    print(f"clade: {clade} -> {clade}_amalgkit")
    print(f"species ({len(entries)}):")
    for entry in sorted(entries, key=lambda e: e.token):
        gaps = [f for f in ("taxon_id", "accession", "assembly_name") if not getattr(entry, f)]
        status = f" (TODO: {', '.join(gaps)})" if gaps else ""
        print(f"  - {entry.token}{status}")
    if args.dry_run:
        print("dry-run: no files written")
        return 0

    try:
        written = scaffold(clade, args.species_list, args.output_dir)
    except (FileExistsError, FileNotFoundError, RuntimeError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"wrote {len(written)} files under {args.output_dir / (clade + '_amalgkit')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
