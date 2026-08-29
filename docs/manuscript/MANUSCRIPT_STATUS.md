# Manuscript status — MetaInformAnt

**Repo type:** Tool/library ecosystem. `metainformant` is a multi-domain
ant-science analysis platform (RNA-seq, genomics, eQTL, GWAS, ecology,
epigenome, phylogenetics, multi-omics integration — see `docs/` topic
folders and `docs/CAPABILITY_MATRIX.md`), driven by `uv` and organized as
`src/metainformant/` plus per-study subprojects under `projects/`.

**Evidence checked:** root listing (`README.md`, `pyproject.toml`,
`QUICKSTART.md`, `SPEC.md`, `src/`, `tests/`, `scripts/`, `docs/` with 60+
topic files, `projects/` gitlinks/subprojects, `output/`), `docs/index.md`,
`docs/AGENTS.md`. No `manuscript/` or `docs/manuscript/` directory existed
before this file.

**Why no publication-target manuscript applies today:** the repository is an
analysis platform whose deliverables are code, pipelines, configs, and
generated analyses; it hosts multiple independent subprojects
(`projects/apis_gwas`, `projects/hymenoptera_amalgkit`,
`projects/drosophila_scrna_2026`) rather than a single narrative research
output. Per-study manuscripts, if any, belong to those subprojects.

**What would trigger creating one:** a methods/tool paper describing the
metainformant platform itself (architecture, guarded-config real-data
policy, pipeline orchestration), or promotion of one subproject's findings
into a publication-track manuscript. At that point, add a full `manuscript/`
tree at the repo top level (config.yaml, section files 00–99,
references.bib) following the docxology/template standard.
