# Examples Best Practices Guide

This guide outlines best practices for writing, maintaining, and using METAINFORMANT examples.

## Writing Examples

### Structure and Organization

#### File Structure
```
examples/
├── domain/
│   ├── README.md              # Domain overview and examples list
│   ├── example_name.py        # Individual example
│   └── example_other.py       # Another example
```

#### Example Template
```python
#!/usr/bin/env python3
"""Brief description of what this example demonstrates.

This example shows:
- Key concept 1
- Key concept 2
- Key concept 3

Usage:
    python examples/domain/example_name.py

Expected output:
    output/examples/domain/results.json
"""

from __future__ import annotations

from pathlib import Path
from metainformant.core import io

def main():
    """Main function."""
    print("Example Name")
    print("=" * 40)

    # Create output directory
    output_dir = Path("output/examples/domain")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "results.json"

    try:
        # Example code here
        results = {"key": "value"}

        # Save results
        io.dump_json(results, output_file)
        print(f"Results saved to: {output_file}")

    except Exception as e:
        print(f"Error: {e}")
        return 1

    return 0

if __name__ == "__main__":
    exit(main())
```

### Code Quality Standards

#### Imports
```python
# Good: Specific imports
from metainformant.dna import sequences
from metainformant.core import io, paths

# Avoid: Wildcard imports
# from metainformant.dna import *
```

#### Error Handling
```python
# Good: Specific exception handling
try:
    result = risky_operation()
except ValueError as e:
    print(f"Invalid value: {e}")
    return 1
except Exception as e:
    print(f"Unexpected error: {e}")
    return 1

# Avoid: Bare except clauses
# except:
#     pass
```

#### Documentation
```python
def process_data(data: list) -> dict:
    """Process biological data and return statistics.

    Args:
        data: List of measurements

    Returns:
        Dictionary with statistics

    Raises:
        ValueError: If data is empty

    Examples:
        >>> process_data([1, 2, 3])
        {'mean': 2.0, 'std': 1.0}
    """
    if not data:
        raise ValueError("Data cannot be empty")
    # ... implementation
```

### Data Management

#### Sample Data
```python
# Good: Small, representative datasets
SAMPLE_SEQUENCES = {
    "gene1": "ATCGATCGATCG",
    "gene2": "GCTAGCTAGCTA"
}

# Avoid: Large datasets in code
# GENOME_DATA = "..."  # 3GB string
```

#### Output Structure
```python
# Consistent output structure
results = {
    "example": "example_name",
    "domain": "dna",
    "description": "Sequence analysis",
    "timestamp": "2024-01-01T12:00:00Z",
    "results": {
        "analysis_type": "gc_content",
        "data": {...}
    }
}
```

### Performance Considerations

#### Efficient Code
```python
# Good: Vectorized operations
import numpy as np
gc_contents = np.mean([seq.count('G') + seq.count('C') for seq in sequences])

# Avoid: Inefficient loops
# gc_contents = []
# for seq in sequences:
#     gc = (seq.count('G') + seq.count('C')) / len(seq)
#     gc_contents.append(gc)
```

#### Memory Management
```python
# Good: Process large files in chunks
with open(large_file, 'r') as f:
    for line in f:
        process_line(line)  # Process one line at a time

# Avoid: Loading entire files into memory
# with open(large_file, 'r') as f:
#     data = f.read()  # May cause memory issues
```

## Testing Examples

### Automated Testing
```bash
# Test all examples
python scripts/test_examples.py

# Test specific domain
python scripts/test_examples.py --domain dna

# Test with HTML report
python scripts/test_examples.py --html
```

### Manual Testing
```bash
# Run individual example
python examples/dna/example_sequences.py

# Check output
ls -la output/examples/dna/

# Validate JSON output
python -c "import json; print(json.load(open('output/examples/dna/results.json')))"
```

### Performance Testing
```bash
# Establish baseline
python scripts/benchmark_examples.py --baseline

# Compare performance
python scripts/benchmark_examples.py --compare --threshold 0.1
```

## Maintenance Practices

### Regular Updates
- **Weekly testing**: Run full example suite
- **Dependency checks**: Verify dependencies are available
- **Performance monitoring**: Check for regressions

### Updating Examples
```bash
# 1. Test current functionality
python scripts/test_examples.py --example domain/example_name.py

# 2. Make changes
# Edit example file...

# 3. Test again
python scripts/test_examples.py --example domain/example_name.py

# 4. Update documentation if needed
# Edit README.md...
```

### Deprecation Handling
```python
# Mark deprecated examples
"""
DEPRECATED: This example uses old API.
Use examples/new_domain/example_name.py instead.
"""

# Update imports for new API versions
# from metainformant.old_module import old_function  # Old
from metainformant.new_module import new_function    # New
```

## Using Examples

### Learning Path
1. **Start with core examples**: Learn fundamental concepts
2. **Progress by domain**: Explore specific analysis types
3. **Try integration examples**: Combine multiple domains
4. **Create custom examples**: Adapt for your use cases

### Integration with Research
```python
# Example: Research workflow
def research_workflow():
    # 1. Load and validate data
    data = load_research_data()

    # 2. Quality control (from examples/quality/)
    qc_results = quality_control(data)

    # 3. Analysis (from examples/dna/, examples/gwas/, etc.)
    analysis_results = perform_analysis(data)

    # 4. Visualization (from examples/visualization/)
    create_plots(analysis_results)

    return analysis_results
```

### Production Adaptation
```python
# Example adaptation for production
def production_analysis(data_path: str, config: dict) -> dict:
    """Production-ready analysis function."""

    # Load configuration (from examples/core/)
    full_config = load_config_with_defaults(config)

    # Validate inputs (from examples/core/)
    validate_input_data(data_path)

    # Perform analysis (adapted from domain examples)
    results = run_analysis(data_path, full_config)

    # Save results (consistent output structure)
    save_results(results, config.get('output_path'))

    return results
```

## Documentation Standards

### README Files
Each domain should have a README.md with:

```markdown
# Domain Name Analysis Examples

Brief description of the domain and available examples.

## Examples

### example_name.py
Description of what this example demonstrates.

```bash
python examples/domain/example_name.py
```

Expected output: `output/examples/domain/results.json`

### example_other.py
Another example description...

## Dependencies

- metainformant.domain
- numpy
- matplotlib (optional)

## Related Documentation

- [Main documentation](../docs/domain/)
- [API reference](../../docs/api/)
```

### Inline Documentation
```python
"""
DNA Sequence Analysis Example

This example demonstrates:
- FASTA file reading
- Sequence validation
- Basic sequence statistics
- GC content calculation
- Reverse complement generation

Usage:
    python examples/dna/example_sequences.py

The example uses sample sequences but can be adapted for real data
by modifying the SAMPLE_SEQUENCES dictionary or loading from files.

Expected output:
    output/examples/dna/sequences_results.json

Output format:
{
    "example": "sequences",
    "domain": "dna",
    "description": "DNA sequence analysis example",
    "results": {
        "seq1": {
            "sequence": "ATCG...",
            "length": 100,
            "gc_content": 0.45,
            "reverse_complement": "CGAT..."
        }
    }
}
"""
```

## Contributing Examples

### Adding New Examples
```bash
# 1. Use the example generator
python scripts/generate_example.py dna sequences \
    --description "DNA sequence analysis demonstration" \
    --features "FASTA reading,GC content,reverse complement"

# 2. Implement the example
# Edit the generated file...

# 3. Test the example
python scripts/test_examples.py --example dna/example_sequences.py

# 4. Update documentation
# Edit examples/dna/README.md...

# 5. Add to dependencies if needed
# Edit examples/dependencies.json...
```

### Quality Checklist
- [ ] **Runs without errors**: `python scripts/test_examples.py --example path`
- [ ] **Follows structure**: Proper imports, error handling, output format
- [ ] **Well documented**: Docstrings, comments, README updates
- [ ] **Dependencies listed**: In `examples/dependencies.json`
- [ ] **Performance reasonable**: Completes in < 30 seconds
- [ ] **Output validated**: JSON structure correct, data reasonable

### Review Process
1. **Automated checks**: Dependencies, syntax, imports
2. **Manual testing**: Run example, check output
3. **Code review**: Style, documentation, best practices
4. **Documentation review**: README updates, examples list
5. **Integration testing**: Full example suite passes

## Troubleshooting

### Common Issues
- **Import errors**: Check METAINFORMANT installation
- **Permission errors**: Ensure output directories are writable
- **Memory issues**: Reduce dataset sizes for examples
- **Network timeouts**: Use local data or increase timeouts

### Getting Help
- Check `examples/TROUBLESHOOTING.md`
- Run dependency checks: `python scripts/check_example_dependencies.py`
- Test with verbose output: `python scripts/test_examples.py --verbose`

## Performance Guidelines

### Target Metrics
- **Execution time**: < 30 seconds for basic examples
- **Memory usage**: < 500MB RAM
- **Disk space**: < 100MB output files

### Optimization Techniques
```python
# Use efficient data structures
from collections import defaultdict
results = defaultdict(list)

# Vectorize computations where possible
import numpy as np
values = np.array(data)

# Use generators for large datasets
def process_large_file(file_path):
    with open(file_path) as f:
        for line in f:
            yield process_line(line)
```

### Profiling Examples
```bash
# CPU profiling
python -m cProfile examples/dna/example_alignment.py

# Memory profiling
python -m memory_profiler examples/dna/example_alignment.py

# Line-by-line timing
python -m line_profiler examples/dna/example_alignment.py
```

## Future Improvements

### Planned Enhancements
- **Interactive examples**: Jupyter notebook versions
- **Web interface**: Online example runner
- **Example gallery**: Visual showcase of outputs
- **Performance dashboard**: Real-time metrics
- **Automated updates**: Keep examples current with API changes

### Community Contributions
- **Domain experts**: Add domain-specific examples
- **Performance optimizations**: Improve slow examples
- **New analysis types**: Cover emerging methods
- **Educational content**: Tutorials and walkthroughs
