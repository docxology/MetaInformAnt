# Amalgkit Configuration

## Overview

YAML configurations for the **amalgkit** RNA-seq data integration pipeline. These configs control the full workflow: metadata retrieval → FASTQ download → transcript quantification → expression matrix generation → quality curation.

## 📦 Species Inventory (23 Species)

Configuration files for **22 ant species** plus *Apis mellifera*. All use NCBI RefSeq accessions with verified annotations.

| Species | Code | Accession | Config File | Status |
|---|---|---|---|---|
| *Acromyrmex echinatior* | `acromyrmex_echinatior` | GCF_000204515.1 | `amalgkit_acromyrmex_echinatior.yaml` | ✅ Verified |
| *Anoplolepis gracilipes* | `anoplolepis_gracilipes` | GCF_047496725.1 | `amalgkit_anoplolepis_gracilipes.yaml` | ✅ Verified |
| *Apis mellifera* | `apis_mellifera` | GCF_003254395.2 | `amalgkit_apis_mellifera_all.yaml` | ✅ Verified |
| *Atta cephalotes* | `atta_cephalotes` | GCF_000143395.1 | `amalgkit_atta_cephalotes.yaml` | ✅ Verified |
| *Camponotus floridanus* | `camponotus_floridanus` | GCF_003227725.1 | `amalgkit_camponotus_floridanus.yaml` | ✅ Verified |
| *Cardiocondyla obscurior* | `cardiocondyla_obscurior` | GCF_019399895.1 | `amalgkit_cardiocondyla_obscurior.yaml` | ✅ Verified |
| *Dinoponera quadriceps* | `dinoponera_quadriceps` | GCF_001313825.1 | `amalgkit_dinoponera_quadriceps.yaml` | ✅ Verified |
| *Formica exsecta* | `formica_exsecta` | GCF_003651465.1 | `amalgkit_formica_exsecta.yaml` | ✅ Verified |
| *Harpegnathos saltator* | `harpegnathos_saltator` | GCF_003227715.2 | `amalgkit_harpegnathos_saltator.yaml` | ✅ Verified |
| *Linepithema humile* | `linepithema_humile` | GCF_040581485.1 | `amalgkit_linepithema_humile.yaml` | ✅ Verified |
| *Monomorium pharaonis* | `monomorium_pharaonis` | GCF_013373865.1 | `amalgkit_monomorium_pharaonis.yaml` | ✅ Verified |
| *Nylanderia fulva* | `nylanderia_fulva` | GCF_005281655.2 | `amalgkit_nylanderia_fulva.yaml` | ✅ Verified |
| *Odontomachus brunneus* | `odontomachus_brunneus` | GCF_010583005.1 | `amalgkit_odontomachus_brunneus.yaml` | ✅ Verified |
| *Ooceraea biroi* | `ooceraea_biroi` | GCF_003672135.1 | `amalgkit_ooceraea_biroi.yaml` | ✅ Verified |
| *Pogonomyrmex barbatus* | `pbarbatus` | GCF_000187915.1 | `amalgkit_pbarbatus.yaml` | ✅ Verified |
| *Solenopsis invicta* | `solenopsis_invicta` | GCF_016802725.1 | `amalgkit_solenopsis_invicta.yaml` | ✅ Verified |
| *Temnothorax americanus* | `temnothorax_americanus` | GCF_048541705.1 | `amalgkit_temnothorax_americanus.yaml` | ✅ Verified |
| *Temnothorax curvispinosus* | `temnothorax_curvispinosus` | GCF_003070985.1 | `amalgkit_temnothorax_curvispinosus.yaml` | ✅ Verified |
| *Temnothorax longispinosus* | `temnothorax_longispinosus` | GCF_030848805.1 | `amalgkit_temnothorax_longispinosus.yaml` | ✅ Verified |
| *Temnothorax nylanderi* | `temnothorax_nylanderi` | GCF_030848795.1 | `amalgkit_temnothorax_nylanderi.yaml` | ✅ Verified |
| *Vollenhovia emeryi* | `vollenhovia_emeryi` | GCF_000949405.1 | `amalgkit_vollenhovia_emeryi.yaml` | ✅ Verified |
| *Wasmannia auropunctata* | `wasmannia_auropunctata` | GCF_000956235.1 | `amalgkit_wasmannia_auropunctata.yaml` | ✅ Verified |

## 🧪 Template & Test

| File | Purpose |
|---|---|
| `amalgkit_template.yaml` | Full reference template |
| `amalgkit_test.yaml` | Minimal test configuration |
| `amalgkit_cross_species.yaml` | Cross-species TMM normalization |

## 📊 Workflow Overview

```mermaid
graph LR
    A[metadata] --> B[select]
    B --> C["getfastq + quant + cleanup<br/>(per-sample, concurrent)"]
    C --> D[merge]
    D --> E[curate]
```

Samples are processed concurrently within chunks of 16. Each sample flows through `getfastq → quant → cleanup` independently, maximizing CPU/network utilization.

## 🚀 Usage

## Configuration Validation

To ensure all species configurations maintain consistency with the template and schema, use the validation script:

```bash
python3 scripts/rna/validate_configs.py
```

This script checks for:

- Required top-level keys
- Valid genome configuration structure
- Explicit `pfd` (parallel-fastq-dump) settings
- Deprecated or problematic flags

## Usage

```bash
# Run all species sequentially (recommended)
nohup bash scripts/rna/run_all_species.sh > output/amalgkit/run_all_species_incremental.log 2>&1 &

# Run a single species
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_{species}.yaml --stream --chunk-size 16

# Check progress
.venv/bin/python scripts/package/generate_custom_summary.py
```

## ⚙️ Key Configuration Options

```yaml
# Basic settings
work_dir: output/amalgkit/{species}/work
threads: 16

# Species
species_list:
  - Pogonomyrmex_barbatus
taxon_id: 144034

# Step-specific settings
steps:
  metadata:
    search_string: '"Species"[Organism] AND "RNA-Seq"[Strategy] AND "Illumina"[Platform]'
  getfastq:
    redo: no          # Skip already-downloaded samples
    aws: yes          # Prefer AWS downloads
    ncbi: no          # Avoid sralite issues
  quant:
    redo: no          # Skip already-quantified samples
    keep_fastq: no    # Delete FASTQs after quant (saves disk)
    index_dir: ...    # Reuse existing kallisto index
```

## 💾 Disk Management

The workflow uses a **per-sample stream-and-clean** pattern with concurrent processing:

1. Download sample FASTQs (~2-15 GB each) — up to 6 samples concurrently
2. Quantify with kallisto (~30 sec per sample)
3. Delete FASTQs immediately after successful quantification
4. Final abundance file: ~2 MB per sample

This allows processing hundreds of samples with only ~80 GB free disk space.

## 🔗 Related Resources

- [Amalgkit Documentation](https://github.com/kfuku52/amalgkit)
- [Amalgkit FAQ](./amalgkit_faq.md) - Common errors and solutions
- [RNA Scripts](../../../scripts/rna/)
