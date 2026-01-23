# Agent Directives: config/amalgkit

## Role
Active amalgkit RNA-seq workflow configurations for production use.

## Contents
Species-specific workflow configurations:
- `amalgkit_template.yaml` - Full template with all options documented
- `amalgkit_test.yaml` - Minimal test configuration
- `amalgkit_pbarbatus_*.yaml` - Pogonomyrmex barbatus configurations
- `amalgkit_pogonomyrmex_barbatus.yaml` - Full species config

## Configuration Structure
```yaml
work_dir: output/amalgkit/{species}
threads: 16
species:
  - scientific_name: "Species name"
    taxid: 12345
steps:
  - metadata
  - getfastq
  - quant
  - merge
```

## Environment Overrides
Use `AK_` prefix:
- `AK_THREADS=8`
- `AK_WORK_DIR=/path/to/output`

## Adding New Species
1. Copy `amalgkit_template.yaml`
2. Fill in species-specific values (taxid, scientific name)
3. Adjust thread/memory based on dataset size
4. Test with small sample subset first
