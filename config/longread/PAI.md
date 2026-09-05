# PAI - config/longread

## Context & Intent

`config/longread/` configures long-read sequencing pipelines (PacBio/ONT) for assembly, error correction, methylation, SV calling, and phasing. Per `AGENTS.md`: "Long-read sequencing (PacBio/ONT) pipeline configurations for assembly and error correction." It holds three real configs: `longread_template.yaml` (all options), `longread_ont_r10.yaml` (ONT R10.4.1 defaults), and `longread_pacbio_hifi.yaml` (PacBio HiFi/Revio defaults).

## Virtual Hierarchy

- Parent: `config/` (repo-wide YAML config layer).
- Per `README.md`: environment overrides use the `LR_` prefix (e.g. `LR_THREADS=16`, `LR_MIN_LENGTH=5000`); extended module docs in `docs/longread/README.md`.
- Downstream consumers: `src/metainformant/longread/` module code; outputs go to `output/longread/...` as declared in each YAML.

## Maintenance Notes

From `AGENTS.md` (binding rules):

- "Validate with schema before committing new configs."
- "Follow REAL IMPLEMENTATION policy — tests use real config files."
- "Use `uv` for dependency management."
- "Environment overrides use the long-read namespace documented by the module."

## AI Workflows

- **New config**: start from `longread_template.yaml`; its top sections are `input` (format fast5/pod5/bam/fastq, `platform`, `chemistry`), `qc`, `assembly`, `methylation`, `sv_calling`, `phasing`, `visualization`, `reporting`, plus `work_dir`, `log_dir`, `output_dir`, `threads`.
- **Platform-specific defaults**: derive ONT work from `longread_ont_r10.yaml` (e.g. `min_overlap: 3000`, `polish_iterations: 1`, CpG-only 5mC) and PacBio from `longread_pacbio_hifi.yaml` (e.g. `min_length: 5000`, `min_quality: 20.0`, `polish_iterations: 0`).
- **Verification**: load the YAML with `metainformant.core.utils.config` helpers under `uv run python`; run targeted longread tests, never repo-wide (shared-instruction USB constraint).
