# Examples Troubleshooting Guide

This guide helps resolve common issues when running METAINFORMANT examples.

## Common Issues

### Import Errors

**Problem**: `ModuleNotFoundError` or `ImportError` when running examples.

**Solutions**:
1. **Install METAINFORMANT**:
   ```bash
   uv pip install -e .
   ```

2. **Check Python path**:
   ```bash
   python -c "import metainformant; print(metainformant.__file__)"
   ```

3. **Missing dependencies**: Run dependency check:
   ```bash
   python scripts/check_example_dependencies.py --example <domain>/example_name.py
   ```

### Permission Errors

**Problem**: `PermissionError` when writing output files.

**Solutions**:
1. **Check output directory permissions**:
   ```bash
   ls -la output/
   chmod -R 755 output/
   ```

2. **Run with proper user permissions** (avoid `sudo` with Python).

### Missing Data Files

**Problem**: Examples fail because data files don't exist.

**Solutions**:
1. **Create sample data**:
   ```python
   # In examples/<domain>/ directory
   python -c "
   import json
   # Create sample data file
   data = {'sample': 'data'}
   with open('sample_data.json', 'w') as f:
       json.dump(data, f)
   "
   ```

2. **Use synthetic data generation** (built into most examples).

### Matplotlib Backend Issues

**Problem**: GUI-related errors when running headless.

**Solutions**:
1. **Set headless backend** (automatically done by test runner):
   ```bash
   export MPLBACKEND=Agg
   ```

2. **Check display environment**:
   ```bash
   echo $DISPLAY  # Should be empty on headless systems
   ```

### Memory Issues

**Problem**: Examples run out of memory on large datasets.

**Solutions**:
1. **Reduce dataset size** in example code
2. **Use streaming processing** for large files
3. **Add memory monitoring**:
   ```python
   import psutil
   print(f"Memory usage: {psutil.virtual_memory().percent}%")
   ```

### Network Timeout Issues

**Problem**: Network-dependent examples timeout.

**Solutions**:
1. **Increase timeout** in example code
2. **Check network connectivity**:
   ```bash
   curl -I https://example.com
   ```
3. **Use local data** instead of remote downloads

## Domain-Specific Issues

### DNA Analysis Examples

**Sequence validation errors**:
- Ensure sequences contain only valid nucleotides (ATCG)
- Check sequence length requirements
- Verify complementary sequences are equal length

**Alignment failures**:
- Check scoring parameters (match, mismatch, gap)
- Ensure sequences are not too dissimilar
- Verify algorithm parameters

### RNA Analysis Examples

**Amalgkit not found**:
```bash
which amalgkit  # Should return path
# If missing, install amalgkit or skip example
```

**Workflow configuration errors**:
- Validate YAML configuration syntax
- Check file paths exist
- Verify species names are supported

### GWAS Examples

**Statistical computation errors**:
- Check genotype matrix format (0/1/2 for diploid)
- Verify phenotype data is numeric
- Ensure sufficient sample size

**Visualization errors**:
- Check matplotlib installation
- Verify data dimensions
- Ensure output directory is writable

### ML Examples

**Scikit-learn import errors**:
```bash
python -c "import sklearn; print(sklearn.__version__)"
# Install if missing: uv pip install scikit-learn
```

**Model training failures**:
- Check data dimensions
- Verify feature matrix is numeric
- Ensure sufficient training samples

## Debugging Techniques

### Enable Verbose Logging

```bash
# Run example with verbose output
python examples/dna/example_sequences.py 2>&1 | tee debug.log
```

### Check Environment

```bash
# Python environment info
python -c "
import sys, platform
print(f'Python: {sys.version}')
print(f'Platform: {platform.platform()}')
import metainformant
print(f'METAINFORMANT: {metainformant.__version__}')
"
```

### Isolate Problems

1. **Minimal reproduction**:
   ```python
   # Create minimal test case
   from metainformant.core import io
   print("Core import successful")

   # Test specific function
   try:
       result = some_function()
       print(f"Function call successful: {result}")
   except Exception as e:
       print(f"Function failed: {e}")
   ```

2. **Step-by-step execution**:
   ```python
   print("Step 1: Importing modules")
   # ... imports ...

   print("Step 2: Loading data")
   # ... data loading ...

   print("Step 3: Processing data")
   # ... processing ...
   ```

### Test Runner Diagnostics

```bash
# Run with detailed error reporting
python scripts/test_examples.py --verbose --domain dna

# Check specific example
python scripts/test_examples.py --example dna/example_sequences.py

# Generate HTML report for visualization
python scripts/test_examples.py --html
```

## Performance Issues

### Slow Examples

**Solutions**:
1. **Profile execution**:
   ```bash
   python -m cProfile examples/dna/example_alignment.py
   ```

2. **Use smaller datasets** in examples
3. **Enable parallel processing** where available

### Memory Leaks

**Detection**:
```python
import tracemalloc
tracemalloc.start()
# Run example code
snapshot = tracemalloc.take_snapshot()
for stat in snapshot.statistics('lineno')[:10]:
    print(stat)
```

## External Dependencies

### Conda Environment

If using conda, ensure proper environment:

```bash
# Create environment
conda create -n metainformant python=3.11
conda activate metainformant

# Install dependencies
conda install numpy scipy matplotlib scikit-learn
uv pip install -e .
```

### Docker Issues

When running in Docker:

```dockerfile
# Ensure matplotlib backend is set
ENV MPLBACKEND=Agg

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*
```

## Getting Help

### Check Documentation

1. **Example README**: `examples/<domain>/README.md`
2. **Main documentation**: `docs/`
3. **API documentation**: Module docstrings

### Report Issues

When reporting problems, include:

1. **Full error traceback**
2. **Environment information**:
   ```bash
   python -c "import sys; print(sys.version)"
   python -c "import metainformant; print(metainformant.__version__)"
   ```
3. **Command used** to run the example
4. **Expected vs actual behavior**

### Community Support

- Check existing GitHub issues
- Create new issue with complete information
- Include example code and error output
