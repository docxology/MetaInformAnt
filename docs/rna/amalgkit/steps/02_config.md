# amalgkit config: Configuration File Generation

## Purpose

Generates a series of configuration files (`.config`) for metadata selection and quality filtering. These config files define tissue mappings, quality thresholds, and sample selection criteria for downstream analysis.

## Overview

The `config` step:
- Creates `.config` files for tissue/sample group mappings
- Provides templates for custom filtering criteria
- Supports pre-configured sets (vertebrate, plantae, test)
- Enables species-specific customization
- Generates base templates for manual editing

## Usage

### Basic Usage

```bash
amalgkit config \
  --out_dir output/amalgkit/work \
  --config base \
  --overwrite no
```

### Python API

```python
from metainformant.rna import amalgkit

result = amalgkit.config(
    out_dir="output/amalgkit/work",
    config="base",
    overwrite=False
)
```

### Configuration File

```yaml
steps:
  config:
    out_dir: output/amalgkit/amellifera/work
    config: base
    overwrite: no
```

## Parameters

### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--out_dir` | PATH | `./` | Directory where config files are generated. |
| `--config` | STR | `base` | Config dataset to export. Options: `base`, `base_all`, `test`, `vertebrate`, `plantae` |
| `--overwrite` | yes/no | `no` | Allow overwriting existing config files. |

## Config Dataset Options

### `base` - Minimal Template Set

**Purpose**: Minimal set for creating custom configs  
**Files Generated**: Essential `.config` files only  
**Use Case**: Starting point for custom species configuration

```bash
amalgkit config --config base --out_dir output/work
```

**Generated Files**:
```
out_dir/config/
├── tissue.config              # Tissue name mappings
├── sample_group.config        # Sample group definitions
└── quality.config             # Quality thresholds
```

### `base_all` - Complete Empty Template Set

**Purpose**: Complete set of near-empty `.config` files  
**Files Generated**: All possible config files  
**Use Case**: Comprehensive customization, advanced users

```bash
amalgkit config --config base_all --out_dir output/work
```

**Generated Files**: ~20+ config files with all possible parameters

### `test` - Testing Dataset

**Purpose**: Short animal dataset for testing amalgkit  
**Files Generated**: Pre-populated configs for test species  
**Use Case**: Workflow validation, testing, tutorials

```bash
amalgkit config --config test --out_dir output/work
```

**Species Included**: Small curated set (Drosophila, C. elegans, etc.)

### `vertebrate` - Vertebrate Preset

**Purpose**: Pre-configured for vertebrate animal data  
**Files Generated**: Vertebrate-specific tissue mappings  
**Use Case**: Mammal, bird, fish, reptile, amphibian studies

```bash
amalgkit config --config vertebrate --out_dir output/work
```

**Tissue Mappings**: Brain, liver, heart, kidney, muscle, etc.

### `plantae` - Plant Preset

**Purpose**: Pre-configured for plant data  
**Files Generated**: Plant-specific tissue mappings  
**Use Case**: Angiosperm, gymnosperm, fern studies

```bash
amalgkit config --config plantae --out_dir output/work
```

**Tissue Mappings**: Leaf, root, flower, seed, stem, etc.

## Input Requirements

### Prerequisites

- **Output Directory**: Must exist or be creatable
- **Write Permissions**: For config file generation

### No External Dependencies

The `config` step is self-contained and requires only:
- Python (via amalgkit)
- File system write access

## Output Files

### Directory Structure

```
out_dir/config/
└── <config_name>/
    ├── tissue.config
    ├── sample_group.config
    ├── quality.config
    ├── platform.config
    ├── layout.config
    └── [additional .config files]
```

### Config File Format

Each `.config` file is tab-delimited with mappings:

**Example: `tissue.config`**
```
original_name	standard_name	include
brain	brain	yes
whole brain	brain	yes
cerebral cortex	brain	yes
liver tissue	liver	yes
muscle	muscle	yes
unknown tissue	NA	no
```

**Example: `quality.config`**
```
parameter	threshold	operator
min_spots	5000000	>=
min_bases	1000000000	>=
max_spots	999999999	<=
avgLength	50	>=
```

## Workflow Integration

### Position in Pipeline

```mermaid
flowchart LR
    A[metadata] --> B[config]
    B --> C[select]
    C --> D[getfastq]
```

**config** runs **after metadata** and **before select**.

### Downstream Usage

| Step | Usage | Description |
|------|-------|-------------|
| `select` | Reads `.config` files | Applies tissue mappings and quality filters |
| `curate` | Uses sample groups | Groups samples for batch correction |

## Common Use Cases

### 1. Generate Base Templates for Customization

```bash
# Create minimal config files
amalgkit config \
  --out_dir output/amalgkit/amellifera/work \
  --config base

# Edit generated files to customize tissue mappings
nano output/amalgkit/amellifera/work/config/base/tissue.config
```

**Workflow**: Generate → Customize → Use in select

### 2. Use Pre-configured Vertebrate Set

```bash
# For mouse, rat, human, zebrafish, etc.
amalgkit config \
  --out_dir output/amalgkit/vertebrate/work \
  --config vertebrate
```

**Result**: Ready-to-use configs with common vertebrate tissues

### 3. Use Plant-Specific Configs

```bash
# For Arabidopsis, rice, maize, etc.
amalgkit config \
  --out_dir output/amalgkit/plantae/work \
  --config plantae
```

**Result**: Plant-specific tissue terminology and mappings

### 4. Testing/Development Mode

```bash
# Quick test with pre-populated data
amalgkit config \
  --out_dir output/test/work \
  --config test
```

**Result**: Small dataset for rapid workflow testing

## Customizing Config Files

### Editing Tissue Mappings

```bash
# Generate base configs
amalgkit config --config base --out_dir output/work

# Edit tissue.config
nano output/work/config/base/tissue.config
```

**Add Custom Mappings**:
```
# Map various brain tissue descriptions to "brain"
brain	brain	yes
whole brain	brain	yes
brain tissue	brain	yes
cerebral cortex	brain	yes
hippocampus	brain	yes
prefrontal cortex	brain	yes

# Exclude ambiguous tissues
mixed tissue	NA	no
whole organism	NA	no
```

### Setting Quality Thresholds

Edit `quality.config`:
```
parameter	threshold	operator
min_spots	10000000	>=    # At least 10M reads
min_bases	1000000000	>=  # At least 1Gb
avgLength	75	>=            # At least 75bp reads
mapping_rate	0.5	>=        # At least 50% mapping
```

### Defining Sample Groups

Edit `sample_group.config`:
```
tissue	sample_group	include
brain	brain	yes
liver	liver	yes
heart	heart	yes
kidney	kidney	yes
muscle	muscle	yes
```

## Performance Considerations

### Runtime

- **All modes**: <1 second
- **File I/O only**: No computation

### Storage

- **base**: ~10KB total
- **base_all**: ~100KB total  
- **vertebrate/plantae**: ~50KB total

## Troubleshooting

### Issue: Config files already exist

```
Error: Config files exist and --overwrite is 'no'
```

**Solutions**:
1. Use `--overwrite yes` to replace:
   ```bash
   amalgkit config --config base --overwrite yes
   ```

2. Delete existing configs manually:
   ```bash
   rm -rf output/work/config/base
   amalgkit config --config base
   ```

3. Use different output directory:
   ```bash
   amalgkit config --config base --out_dir output/work_v2
   ```

### Issue: Select step can't find config files

```
Error: Config directory not found
```

**Solutions**:
1. Verify config files exist:
   ```bash
   ls output/work/config/base/
   ```

2. Specify config_dir explicitly in select:
   ```bash
   amalgkit select --config_dir output/work/config/base
   ```

3. Check directory structure matches expected format

### Issue: Custom tissue mappings not working

**Diagnosis**:
```bash
# Check tissue.config format
cat output/work/config/base/tissue.config
```

**Solutions**:
1. Ensure tab-delimited format (not spaces)
2. Verify header row exists: `original_name	standard_name	include`
3. Check for typos in include column (must be `yes` or `no`)
4. Save file with Unix line endings (LF, not CRLF)

## Best Practices

### 1. Start with Pre-configured Sets

```bash
# For vertebrates
--config vertebrate

# For plants
--config plantae

# Then customize as needed
```

### 2. Version Control Config Files

```bash
# Keep configs in version control
git add config/
git commit -m "Add custom tissue mappings for Apis mellifera"
```

### 3. Document Custom Mappings

Add comments to config files:
```
# Bee-specific tissue mappings
# Author: Research Team
# Date: 2025-10-29
brain	brain	yes
antenna	antenna	yes  # Bee-specific sensory organ
mandibular_gland	gland	yes  # Social bee gland
```

### 4. Test Configs Before Large-Scale Analysis

```bash
# Generate test configs
amalgkit config --config base --out_dir output/test/work

# Test with select step on small subset
amalgkit select --config_dir output/test/work/config/base --max_sample 5
```

## Real-World Examples

### Example 1: Honey Bee Brain Study

```bash
# Generate base configs
amalgkit config \
  --out_dir output/amalgkit/amellifera/work \
  --config base

# Customize for bee tissues
cat > output/amalgkit/amellifera/work/config/base/tissue.config << 'EOF'
original_name	standard_name	include
brain	brain	yes
whole brain	brain	yes
head	brain	yes
antenna	antenna	yes
mandibular gland	gland	yes
fat body	fat_body	yes
ovary	ovary	yes
EOF
```

### Example 2: Cross-Species Vertebrate Analysis

```bash
# Use vertebrate preset for consistency
for species in human mouse zebrafish; do
  amalgkit config \
    --out_dir output/amalgkit/${species}/work \
    --config vertebrate
done
```

**Result**: Consistent tissue mappings across all three species

### Example 3: Plant Developmental Study

```bash
# Plant-specific configs
amalgkit config \
  --out_dir output/amalgkit/arabidopsis/work \
  --config plantae

# Customize for developmental stages
cat >> output/amalgkit/arabidopsis/work/config/plantae/tissue.config << 'EOF'
seedling_leaf	leaf	yes
mature_leaf	leaf	yes
senescent_leaf	leaf	yes
flower_bud	flower	yes
open_flower	flower	yes
silique	fruit	yes
EOF
```

## Integration with METAINFORMANT Workflow

### Automatic Config Generation

```python
from metainformant.rna.workflow import execute_workflow, load_workflow_config

cfg = load_workflow_config("config/amalgkit_amellifera.yaml")
execute_workflow(cfg)  # config step runs automatically after metadata
```

### Config File Location

The workflow expects configs at:
```
work_dir/config/<config_name>/
```

Default `config_name` is `base`.

## References

- **Tissue Ontology**: http://www.obofoundry.org/ontology/uberon.html
- **METAINFORMANT Workflow**: `docs/rna/workflow.md`

## See Also

- **Previous Step**: [`metadata.md`](metadata.md) - Retrieving SRA metadata
- **Next Step**: [`select.md`](select.md) - Selecting samples with configs
- **Workflow Overview**: [`../amalgkit.md`](../amalgkit.md)
- **Testing**: `tests/test_rna_amalgkit_steps.py::test_config_basic_execution`

---

**Last Updated**: October 29, 2025  
**AMALGKIT Version**: 0.12.19  
**Status**: ✅ Production-ready, comprehensively tested


