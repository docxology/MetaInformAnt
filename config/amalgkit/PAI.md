# Personal AI Infrastructure (PAI) - amalgkit

## 🧠 Context & Intent

- **Path**: `config/amalgkit/`
- **Purpose**: YAML configurations for amalgkit RNA-seq transcript quantification workflows
- **Domain**: config → bioinformatics → RNA-seq

## 🏗️ Virtual Hierarchy

- **Type**: Configuration Directory
- **Parent**: `config`
- **Consumers**: `scripts/rna/`, `src/metainformant/rna/`

## 📊 Production Status

| Config | Samples | Status |
|--------|---------|--------|
| `amalgkit_acromyrmex_echinatior.yaml` | 48 | 🔄 In Progress |
| `amalgkit_apis_mellifera_all.yaml` | ~7,342 | ⏳ Queued |
| All 23 species | Varies | ⏳ Queued (sequential via `run_all_species.sh`) |

## 📝 Maintenance Notes

- **Dependencies**: `amalgkit>=0.12.20`, `kallisto`, `fasterq-dump`
- **Disk Strategy**: Per-sample stream-and-clean (concurrent within chunks of 6)
- **Critical Settings**: `redo: no` for production runs (idempotent), `aws: yes`, `ncbi: no`
- **Shared Resources**: Genome/index in `output/amalgkit/shared/`
- **Threads**: 16 total, dynamically divided across concurrent samples
- **Metadata Filtering**: Always use `AND "RNA-Seq"[Strategy] AND "Illumina"[Platform]`

## 🔄 AI Workflows

- **New Species**: Copy `amalgkit_template.yaml`, adjust paths/taxon
- **Recovery**: Re-run `bash scripts/rna/run_all_species.sh` — it auto-resumes via `redo: no`
- **Monitoring**: `.venv/bin/python scripts/package/generate_custom_summary.py`
- **Tissue Patches**: Add new bioproject annotations to `tissue_patches.yaml`
- **Documentation**: Update this file and `README.md` when adding configs
