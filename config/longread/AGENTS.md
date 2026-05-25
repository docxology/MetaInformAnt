# Agent Directives: config/longread

## Role

Long-read sequencing (PacBio/ONT) pipeline configurations for assembly and error correction.

## Contents

| File | Description |
| :--- | :--- |
| `longread_template.yaml` | Reference template with all options |
| `longread_ont_r10.yaml` | Oxford Nanopore R10 chemistry defaults |
| `longread_pacbio_hifi.yaml` | PacBio HiFi mode defaults |

## Configuration Structure

```yaml
# Long-read pipeline configuration
platform: ont  # or pacbio
chemistry: r10  # r9 | r10 | hifi | clr
assembler: flye
error_correction: medaka
threads: 16
```

## Rules

- Validate with schema before committing new configs
- Follow REAL IMPLEMENTATION policy — tests use real config files
- Use `uv` for dependency management
- Environment overrides use `AK_` prefix
