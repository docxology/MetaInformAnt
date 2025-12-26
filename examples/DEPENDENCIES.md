# Examples Dependencies Matrix

This document provides a comprehensive overview of dependencies required by METAINFORMANT examples, helping users understand what needs to be installed to run specific examples.

## Dependency Overview

### Core Dependencies (Required by All)
- **Python**: 3.11+ (3.12 recommended)
- **METAINFORMANT**: Core package (automatically installed)
- **Standard Library**: pathlib, json, typing (built-in)

### Optional Dependencies by Category

#### Scientific Computing
- **numpy**: Numerical computing (required by: DNA, GWAS, ML, Math, Information, Networks)
- **scipy**: Scientific computing (required by: Math, Information, Ecology)
- **pandas**: Data manipulation (required by: GWAS, RNA, Multi-omics, Life Events)

#### Machine Learning
- **scikit-learn**: ML algorithms (required by: ML, Multi-omics)
- **matplotlib**: Plotting (required by: Visualization, GWAS, Math, Life Events)

#### Bioinformatics
- **biopython**: Sequence analysis (optional for: DNA, Protein)
- **scanpy**: Single-cell analysis (required by: Single-cell)
- **anndata**: Single-cell data structures (required by: Single-cell)

#### Network Analysis
- **networkx**: Graph algorithms (required by: Networks, Ontology)

#### External Tools
- **amalgkit**: RNA-seq pipeline (required by: RNA examples)

## Domain Dependencies Matrix

| Domain | Required Dependencies | Optional Dependencies | External Tools |
|--------|----------------------|----------------------|----------------|
| **Core** | metainformant.core | pyyaml | None |
| **DNA** | metainformant.dna, numpy | biopython | None |
| **RNA** | metainformant.rna | pandas, numpy | amalgkit |
| **GWAS** | metainformant.gwas, numpy, pandas, matplotlib | scikit-learn | None |
| **Protein** | metainformant.protein | biopython | None |
| **Epigenome** | metainformant.epigenome | pandas, numpy | None |
| **Ontology** | metainformant.ontology | networkx | None |
| **Phenotype** | metainformant.phenotype | None | None |
| **Ecology** | metainformant.ecology, numpy | scipy | None |
| **Math** | metainformant.math, numpy | scipy, matplotlib | None |
| **Information** | metainformant.information, numpy | scipy | None |
| **Life Events** | metainformant.life_events | pandas, numpy, matplotlib | None |
| **Multi-Omics** | metainformant.multiomics | numpy, pandas, scikit-learn | None |
| **Single-Cell** | metainformant.singlecell | scanpy, anndata, numpy, pandas | None |
| **Quality** | metainformant.quality | None | None |
| **Networks** | metainformant.networks | networkx, numpy | None |
| **ML** | metainformant.ml, numpy, scikit-learn | pandas, matplotlib | None |
| **Simulation** | metainformant.simulation | numpy | None |
| **Visualization** | metainformant.visualization, matplotlib | seaborn, plotly | None |

## Example-Specific Dependencies

### Core Examples

#### `examples/core/example_config.py`
- **Required**: metainformant.core
- **Optional**: pyyaml (for YAML config files)
- **Install**: `pip install pyyaml`

#### `examples/core/example_io.py`
- **Required**: metainformant.core
- **Optional**: None
- **Notes**: Tests various file I/O formats

#### `examples/core/example_logging.py`
- **Required**: metainformant.core
- **Optional**: None
- **Notes**: Demonstrates logging patterns

#### `examples/core/example_paths.py`
- **Required**: metainformant.core
- **Optional**: None
- **Notes**: Path validation and manipulation

#### `examples/core/example_workflow.py`
- **Required**: metainformant.core
- **Optional**: pyyaml
- **Install**: `pip install pyyaml`

### DNA Examples

#### `examples/dna/example_sequences.py`
- **Required**: metainformant.dna, numpy
- **Optional**: biopython (for advanced sequence operations)
- **Install**: `pip install numpy biopython`

#### `examples/dna/example_alignment.py`
- **Required**: metainformant.dna, numpy
- **Optional**: None
- **Notes**: Sequence alignment algorithms

#### `examples/dna/example_phylogeny.py`
- **Required**: metainformant.dna, numpy
- **Optional**: None
- **Notes**: Phylogenetic tree construction

#### `examples/dna/example_population.py`
- **Required**: metainformant.dna, numpy
- **Optional**: None
- **Notes**: Population genetics statistics

### RNA Examples

#### `examples/rna/example_amalgkit.py`
- **Required**: metainformant.rna
- **Optional**: pandas, numpy
- **External**: amalgkit CLI tool
- **Install**: Follow RNA documentation for amalgkit installation
- **Notes**: Requires amalgkit to be installed and in PATH

#### `examples/rna/example_quantification.py`
- **Required**: metainformant.rna
- **Optional**: pandas, numpy
- **Install**: `pip install pandas numpy`

### GWAS Examples

#### `examples/gwas/example_association.py`
- **Required**: metainformant.gwas, numpy, pandas
- **Optional**: scikit-learn
- **Install**: `pip install numpy pandas scikit-learn matplotlib`

#### `examples/gwas/example_visualization.py`
- **Required**: metainformant.gwas, numpy, pandas, matplotlib
- **Optional**: seaborn
- **Install**: `pip install numpy pandas matplotlib seaborn`

### ML Examples

#### `examples/ml/example_pipeline.py`
- **Required**: metainformant.ml, numpy, scikit-learn
- **Optional**: pandas, matplotlib
- **Install**: `pip install numpy scikit-learn pandas matplotlib`

### Integration Examples

#### `examples/integration/example_multiomics.py`
- **Required**: metainformant, numpy
- **Optional**: pandas, scikit-learn
- **Install**: `pip install numpy pandas scikit-learn`

#### `examples/integration/example_dna_rna.py`
- **Required**: metainformant.dna, metainformant.rna, numpy
- **Optional**: pandas
- **Install**: `pip install numpy pandas`

#### `examples/integration/example_complete_workflow.py`
- **Required**: metainformant (multiple modules), numpy, pandas
- **Optional**: matplotlib, scikit-learn
- **Install**: `pip install numpy pandas matplotlib scikit-learn`

## Installation Commands

### Minimal Installation (Core Examples Only)
```bash
# Install METAINFORMANT
pip install -e .

# Test core examples
python scripts/test_examples.py --domain core
```

### Full Scientific Stack
```bash
# Install METAINFORMANT with scientific dependencies
pip install -e ".[scientific]"

# Additional bioinformatics packages
pip install biopython scanpy anndata networkx

# Visualization
pip install matplotlib seaborn plotly
```

### Domain-Specific Installations

#### DNA Analysis
```bash
pip install numpy biopython
```

#### RNA Analysis
```bash
pip install pandas numpy
# Install amalgkit separately (see RNA docs)
```

#### GWAS Analysis
```bash
pip install numpy pandas matplotlib scikit-learn
```

#### ML Analysis
```bash
pip install numpy scikit-learn pandas matplotlib
```

#### Single-Cell Analysis
```bash
pip install scanpy anndata numpy pandas matplotlib
```

#### Network Analysis
```bash
pip install networkx numpy
```

### Conda Installation (Alternative)

```bash
# Create environment
conda create -n metainformant python=3.12
conda activate metainformant

# Install conda packages
conda install numpy scipy pandas matplotlib scikit-learn

# Install pip packages
pip install biopython scanpy anndata networkx

# Install METAINFORMANT
pip install -e .
```

## Dependency Checking

### Automated Dependency Validation
```bash
# Check all dependencies
python scripts/check_example_dependencies.py

# Check specific domain
python scripts/check_example_dependencies.py --domain dna

# Check specific example
python scripts/check_example_dependencies.py --example dna/example_sequences.py
```

### Manual Dependency Checking

#### Check Python Packages
```python
# In Python REPL
import metainformant
print(f"METAINFORMANT: {metainformant.__version__}")

try:
    import numpy
    print("numpy: OK")
except ImportError:
    print("numpy: MISSING")

try:
    import pandas
    print("pandas: OK")
except ImportError:
    print("pandas: MISSING")
```

#### Check External Tools
```bash
# Check amalgkit
which amalgkit
amalgkit --version

# Check other bioinformatics tools
which muscle  # For alignments
which seqkit  # For sequence processing
```

## Troubleshooting Dependencies

### Common Issues

#### Import Errors
**Problem**: `ModuleNotFoundError` for scientific packages
**Solution**:
```bash
pip install --upgrade pip
pip install numpy scipy pandas matplotlib scikit-learn
```

#### Version Conflicts
**Problem**: Package version incompatibilities
**Solution**:
```bash
# Create fresh environment
python -m venv fresh_env
source fresh_env/bin/activate  # Linux/Mac
# fresh_env\Scripts\activate   # Windows

# Install in correct order
pip install numpy
pip install scipy
pip install pandas
pip install -e .
```

#### Missing System Libraries
**Problem**: Compilation errors during installation
**Solution** (Ubuntu/Debian):
```bash
sudo apt-get update
sudo apt-get install build-essential python3-dev
```

**Solution** (macOS):
```bash
# Install Xcode command line tools
xcode-select --install

# Or with Homebrew
brew install gcc
```

#### Conda Package Conflicts
**Problem**: Conda solver issues
**Solution**:
```bash
# Use mamba (faster conda alternative)
conda install mamba
mamba install numpy pandas matplotlib scikit-learn

# Or install in smaller groups
conda install numpy
conda install pandas
conda install matplotlib
```

### External Tool Issues

#### Amalgkit Installation
**Problem**: amalgkit not found in PATH
**Solutions**:
1. **Check installation**:
   ```bash
   which amalgkit
   amalgkit --help
   ```

2. **Add to PATH** (if installed in custom location):
   ```bash
   export PATH=$PATH:/path/to/amalgkit
   # Add to ~/.bashrc for permanence
   ```

3. **Reinstall amalgkit** following RNA documentation

#### Bioinformatics Tools
**Problem**: muscle, seqkit, or other tools missing
**Solutions**:
```bash
# Ubuntu/Debian
sudo apt-get install muscle seqkit

# macOS with Homebrew
brew install muscle seqkit

# Or download from official sources
```

## Performance Impact

### Dependency Size Impact

| Dependency Group | Size Impact | Performance Impact |
|------------------|-------------|-------------------|
| Core only | Minimal | Fast startup |
| Scientific stack | ~200MB | Moderate |
| Full bioinformatics | ~500MB+ | Slower startup |
| ML libraries | ~300MB | GPU acceleration possible |
| Single-cell libs | ~400MB | Memory intensive |

### Startup Time Comparison

| Configuration | Import Time | Notes |
|---------------|-------------|--------|
| Core only | < 1 second | Minimal dependencies |
| Scientific | 2-3 seconds | NumPy, SciPy loading |
| Full ML | 3-5 seconds | Scikit-learn loading |
| Single-cell | 5-8 seconds | Scanpy, AnnData loading |

## Testing Dependencies

### CI/CD Considerations
- **Core examples**: Always tested in CI
- **Scientific examples**: Tested with optional dependencies
- **External tools**: Skipped if not available
- **Performance**: Baseline established, regressions caught

### Local Development
```bash
# Test with minimal dependencies
python scripts/test_examples.py --domain core

# Test with scientific dependencies
python scripts/test_examples.py --domain dna --domain math

# Skip examples with missing dependencies
python scripts/test_examples.py --continue-on-error
```

## Future Dependencies

### Planned Additions
- **PyTorch/TensorFlow**: Deep learning examples
- **Dask**: Parallel computing examples
- **Polars**: High-performance DataFrames
- **BioPython 2.0**: Enhanced sequence analysis
- **Scanpy 2.0**: Improved single-cell analysis

### Deprecation Warnings
- **Python 3.8-3.10**: Will be deprecated in future versions
- **Old NumPy versions**: < 1.20 may cause issues
- **Legacy BioPython**: < 1.80 may have compatibility issues

## Contributing Dependencies

### Adding New Dependencies
1. **Update `examples/dependencies.json`**
2. **Add installation instructions**
3. **Update CI/CD configurations**
4. **Test with new dependencies**
5. **Update documentation**

### Dependency Guidelines
- **Prefer pure Python** when possible
- **Document version requirements**
- **Consider cross-platform compatibility**
- **Test optional dependency handling**
- **Update dependency matrix**

## Quick Reference

### Most Common Dependencies
```bash
# Essential for most examples
pip install numpy pandas matplotlib

# For ML examples
pip install scikit-learn

# For bioinformatics
pip install biopython

# For single-cell
pip install scanpy anndata

# For networks
pip install networkx
```

### Dependency Checking Commands
```bash
# Quick dependency check
python -c "import metainformant, numpy, pandas; print('Core deps OK')"

# Full dependency validation
python scripts/check_example_dependencies.py

# Fix missing dependencies
python scripts/check_example_dependencies.py --fix
```
