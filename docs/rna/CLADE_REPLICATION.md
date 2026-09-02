# Clade Replication Contract — Reusing the MetaInformAnt RNA Stack for a New Clade

This is the owner-facing recipe for replicating the hymenoptera_amalgkit
campaign for another clade (e.g. `lepidoptera_amalgkit`, `diptera_amalgkit`)
on top of the same MetaInformAnt platform. Follow it in order; each step names
the exact artifacts to copy and the contract that must hold.

## 0. Prerequisites

- Parent checkout of MetaInformAnt with an importable package: run `uv sync`,
  then verify that `import metainformant` succeeds from
  any working directory (the project installs editable via
  `[tool.setuptools.packages.find] where = ["src"]`).
- Empty external data root for the new clade, exported as
  `AMALGKIT_DATA_ROOT` (never inside either repository).

## 1. Import contract (non-negotiable)

A clade project's scripts and tests import the parent package **only** through
the shared bootstrap:

```python
from metainformant_import import ensure_metainformant

ensure_metainformant()
from metainformant.rna.analysis.cross_species import ...  # then normal imports
```

Resolution order of `ensure_metainformant()`:

1. `metainformant` already importable (editable install) — no-op.
2. `METAINFORMANT_SRC` env var pointing at the parent `src/` (or checkout root).
3. Parent-checkout layout discovery from the project location
   (`<parent>/projects/<clade_project>` → `<parent>/src`).

Hand-rolled `sys.path.insert` hacks are a contract violation. The canonical
shim lives at `projects/hymenoptera_amalgkit/scripts/metainformant_import.py`;
copy it verbatim into `<new_clade>/scripts/`. Scripts in nested subpackages
(e.g. `scripts/cloud/`) add the project's `scripts/` directory to `sys.path`
once for sibling imports, exactly as the existing cloud scripts do for
`amalgkit_paths`.

## 2. Repository template

Create the new clade as a standalone git repository with this shape (copy the
hymenoptera_amalgkit layout, not its data):

```
<new_clade>/
  config/amalgkit/          # one YAML per species + cross-species YAML
  scripts/                  # thin orchestrators incl. metainformant_import.py
  tests/                    # real-implementation pytest suite
  doc/  docs/               # project docs (docs gate validates these)
  output/                   # working outputs (local-only)
```

Optionally mirror it into the parent checkout under `projects/` for discovery,
as `hymenoptera_amalgkit` is mirrored.

## 3. Species configuration pattern

One YAML per species under `config/amalgkit/`, named
`amalgkit_<genus_species>.yaml`, plus one `amalgkit_cross_species.yaml`. Each
species file declares: `work_dir`, `log_dir`, `threads`, `species_list` (the
`Genus_species` token), `taxon_id`, and a `genome:` block with the NCBI
`accession`, `assembly_name`, and `include:` file set (`genome`, `gff3`, `rna`,
`cds`). Reference aliases (`reference_aliases:`) map public subspecies-label
variants onto the species-level reference contract. Copy
`amalgkit_apis_mellifera.yaml` as the annotated exemplar and replace the
taxonomy/assembly fields.

The authoritative roster lives in a `species_manifest.tsv` consumed by the
pipeline; its exact schema is documented in the parent repo's
`docs/rna/` TSE design integration contract.

## 4. Analysis modules you inherit (no new code required)

- Orthology bridge: `scripts/build_ortholog_bridge.py` chaining
  OrthoDB → NCBI Protein → NCBI RNA
  (`metainformant.rna.analysis.ortholog_mapping`).
- Cross-species descriptive analysis: `run_cross_species_analysis.py`
  (`metainformant.rna.analysis.cross_species`) — fingerprint stability /
  divergence matrices.
- Tissue specificity (tau), orthology classes (1:1 vs multicopy),
  duplication–specialization coupling, and phylogeny-aware specificity
  gain/loss: `metainformant.rna.analysis` tissue-specificity evolution module
  (Mantica 2024 adaptation; see `docs/rna/METHODS_LITERATURE_ALIGNMENT.md`).
- Atlas-style figures/stat regeneration: idempotent, descriptive-only.
- Evidence manifest + downstream provenance: `generate_evidence_manifest.py`,
  `record_downstream_provenance.py`
  (`metainformant.rna.engine.provenance`).

All of these are clade-agnostic given the config inputs of step 3; none
embed Hymenoptera-specific logic.

## 5. Hard boundaries (unchanged)

- Descriptive-only statistics in cross-species outputs until the evidence
  manifest freezes; inferential tests (Wilcoxon/chi²-style) are gated and
  labeled for post-freeze use only.
- Never write into the data root from analysis code; the data root is an
  input volume only.
- Real-implementation tests: real files, real subprocesses; no test doubles.

## 6. Gates before first campaign

```bash
cd <new_clade>
uv run python scripts/validate_project_docs.py          # docs gate
../../.venv/bin/python -m pytest tests/ -q              # per-project suite (one dir per pytest invocation)
```

In the parent checkout:

```bash
uv run python scripts/verify_documentation_code.py
```

## 7. Scaffold generator (single source of truth)

`scripts/rna/new_clade_scaffold.py` in the parent checkout automates steps 2-3
and part of step 1 from the actual template pattern. It is the authoritative
description of the skeleton layout: this section records its observed
behavior; the script owns the shape.

```bash
uv run python scripts/rna/new_clade_scaffold.py \
    --clade lepidoptera --species-list species.txt --output-dir /tmp/clades
```

Species-list format (one per line, `#` comments):

```
Genus species [taxon_id [accession assembly_name]]
```

Example (lepidoptera pilot used to prove the generator, 2026-09-01):

```
Bombyx mori 7091 GCF_000151645.2 Bombyx_mori_p50T
Danaus plexippus
Heliconius melpomene 7122 GCA_000290705.1 Hmel2_5
Plutella xylostella 515471 GCA_000330365.1 DbM_p505
```

Generated skeleton (16 files for the 4-species roster above):

```
<clade>_amalgkit/
  README.md                     # roster table + import contract + gates
  TODO.md                       # per-species scaffold gaps + campaign checklist
  .gitignore                    # output/, logs/, caches
  config/amalgkit/
    amalgkit_<genus_species>.yaml   # one per species, apis_mellifera pattern
    amalgkit_cross_species.yaml     # descriptive-only (inferential_statistics: none)
    select_rules.tsv, tissue_mapping.yaml, tissue_patches.yaml  # copied from parent config/amalgkit/
  scripts/
    metainformant_import.py     # copied VERBATIM from the hymenoptera project
    validate_configs.py         # copied from parent scripts/rna/
  tests/test_import_contract.py # real-subprocess contract tests for the skeleton
  doc/00_setup/01_environment.md
  docs/manuscript/README.md
```

Guarantees:

- **Deterministic**: identical inputs produce byte-identical trees (sorted
  species, no timestamps, no machine paths).
- **Fail-fast gaps**: a species without `taxon_id`/`accession` still emits a
  config, plus TODO entries in `TODO.md`; the skeleton's own
  `validate_configs.py` run then flags the empty accession instead of
  silently downloading nothing.
- **Import contract**: the shim is copied verbatim; the script refuses to run
  if the canonical source is missing.
- **Safety**: refuses to overwrite an existing target; `--dry-run` validates
  input without writing; the clade token is validated (no path traversal).

### Proof run (2026-09-01, lepidoptera)

Recorded from an actual scaffold into a temp directory (never `projects/`):

| Gate | Result |
| --- | --- |
| Skeleton self-verification (structure contract) | PASS |
| `scripts/validate_configs.py` vs Amalgkit 0.16.60 | 3/4 species configs PASS; the intentionally TODO-gapped config (`Danaus_plexippus`, no accession) FAILS with `genome section must define accession or files` — the designed fail-fast |
| Skeleton's own `tests/test_import_contract.py` | 3 passed |
| `ensure_metainformant()` from the skeleton | resolves the parent package (editable-install / parent-layout discovery) |
| Child `validate_project_docs.py` against the skeleton | NOT directly runnable (see gap below) |

Known template gaps surfaced by the proof run (candidates for later rounds,
not blockers):

1. `validate_project_docs.py` hard-codes `PROJECT_ROOT` from its own file
   location and encodes hymenoptera-specific counts (27 species configs) and
   manuscript phrases; a new clade cannot reuse it until it accepts a
   project-root argument and a per-clade contract file. The skeleton
   therefore ships `validate_configs.py` + its own pytest suite as its
   executable gates, and the docs gate must be re-derived per clade.
2. `scripts/validate_configs.py` resolves `select.select_rules_tsv` relative
   to `config_path.parents[2]` (repo root); the skeleton keeps the same
   relative path so a copy works unchanged inside the new repo.

## 8. Post-freeze runbook pointer

Once the evidence manifest freezes for the new clade, follow the same
release/inferential runbook used by the hymenoptera campaign (see
`docs/rna/` campaign and validation docs). Until then, report status only via
`scripts/report_campaign_status.py --json` — never prose snapshots of
counts.
