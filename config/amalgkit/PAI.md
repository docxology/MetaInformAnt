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
| `amalgkit_pbarbatus.yaml` | ~110 | ✅ Complete (95 valid) |
| `amalgkit_apis_mellifera_all.yaml` | ~7,270 | 🔄 In Progress |

## 📝 Maintenance Notes

- **Dependencies**: `amalgkit>=0.12.20`, `kallisto`, `fastp`
- **Disk Strategy**: Stream-and-clean (minimal persistent footprint)
- **Critical Settings**: `redo: no` for production runs (idempotent)
- **Shared Resources**: Genome/index in `output/amalgkit/shared/`
- **Metadata Filtering**: Always use `AND "RNA-Seq"[Strategy] AND "Illumina"[Platform]`

## 🔄 AI Workflows

- **New Species**: Copy `amalgkit_template.yaml`, adjust paths/taxon
- **Recovery**: Use `scripts/rna/recover_missing_parallel.py` for failed samples
- **Tissue Patches**: Add new bioproject annotations to `tissue_patches.yaml`
- **Documentation**: Update this file and `README.md` when adding configs
