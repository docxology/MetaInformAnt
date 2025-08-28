# FASTQ Quality Analysis

The FASTQ quality analysis module provides comprehensive quality assessment for high-throughput sequencing data, offering FastQC-like functionality with enhanced performance and integration.

## Main Analysis Function

### analyze_fastq_quality()

Perform comprehensive quality analysis on FASTQ files:

```python
from metainformant.quality.fastq import analyze_fastq_quality

# Basic analysis
results = analyze_fastq_quality("sample.fastq.gz")

# With custom parameters
results = analyze_fastq_quality(
    "sample.fastq.gz",
    subsample_size=100000,      # Subsample for large files
    quality_encoding='phred33', # Quality score encoding  
    adapter_sequences=None,     # Custom adapter sequences
    output_dir="output/qc"      # Output directory
)

# Print summary
print(f"Analysis complete for {results['filename']}")
print(f"Total sequences: {results['basic_stats']['total_sequences']:,}")
print(f"Mean quality: {results['basic_stats']['mean_quality']:.2f}")
```

**Parameters:**
- `fastq_path`: Path to FASTQ file (supports .gz compression)
- `subsample_size`: Number of sequences to analyze (None = all)
- `quality_encoding`: Quality score encoding ('phred33' or 'phred64')
- `adapter_sequences`: List of custom adapter sequences to detect
- `output_dir`: Directory for output files

**Returns:** Dictionary with comprehensive quality metrics

## Analysis Modules

### Basic Statistics

```python
from metainformant.quality.fastq import basic_statistics

# Get basic statistics
stats = basic_statistics("sample.fastq.gz")
print(stats)
```

**Output:**
```python
{
    'filename': 'sample.fastq.gz',
    'file_type': 'Conventional base calls',
    'encoding': 'Sanger / Illumina 1.9',
    'total_sequences': 1500000,
    'filtered_sequences': 0,
    'sequence_length': '151',  # or range like "35-151"
    'gc_content': 42.5,
    'n_content': 0.1,
    'mean_quality': 35.2
}
```

### Per-Base Quality Scores

```python
from metainformant.quality.fastq import per_base_quality

# Analyze quality scores by position
quality_data = per_base_quality("sample.fastq.gz")
print(f"Quality data shape: {quality_data.shape}")  # (sequence_length, 6)

# Columns: position, mean, median, q1, q3, max, min
```

### Per-Sequence Quality Scores

```python
from metainformant.quality.fastq import per_sequence_quality

# Quality score distribution across sequences
seq_quality = per_sequence_quality("sample.fastq.gz")
print(f"Quality bins: {len(seq_quality)}")

# Returns dictionary: {quality_score: count}
```

### Sequence Length Distribution

```python
from metainformant.quality.fastq import sequence_length_distribution

# Distribution of sequence lengths
length_dist = sequence_length_distribution("sample.fastq.gz")
print(f"Length range: {min(length_dist.keys())} - {max(length_dist.keys())}")

# Returns dictionary: {length: count}
```

### GC Content Distribution

```python
from metainformant.quality.fastq import gc_content_distribution

# GC content distribution across sequences
gc_dist = gc_content_distribution("sample.fastq.gz")
print(f"GC content range: 0-100%")

# Returns dictionary: {gc_percentage: count}
```

### Adapter Content Analysis

```python
from metainformant.quality.fastq import adapter_content

# Detect adapter contamination
adapter_data = adapter_content(
    "sample.fastq.gz",
    adapter_sequences=['AGATCGGAAGAG']  # Illumina adapter
)

print(f"Max adapter content: {adapter_data['max_adapter_content']:.2f}%")
print(f"Adapter positions: {len(adapter_data['position_data'])}")

# Position-wise adapter content
for position, adapters in adapter_data['position_data'].items():
    print(f"Position {position}: {adapters}")
```

**Built-in Adapters:**
- Illumina Universal Adapter: `AGATCGGAAGAG`
- Illumina Small RNA Adapter: `TGGAATTCTCGG`
- Nextera Transposase Sequence: `CTGTCTCTTATA`

### Overrepresented Sequences

```python
from metainformant.quality.fastq import overrepresented_sequences

# Find highly abundant sequences
overrep = overrepresented_sequences(
    "sample.fastq.gz",
    threshold=0.1,      # Sequences appearing in >0.1% of reads
    max_sequences=20    # Report top 20 sequences
)

for seq_info in overrep:
    print(f"Sequence: {seq_info['sequence'][:30]}...")
    print(f"Count: {seq_info['count']}")
    print(f"Percentage: {seq_info['percentage']:.3f}%")
    print(f"Possible source: {seq_info['possible_source']}")
    print()
```

### Duplication Levels

```python
from metainformant.quality.fastq import duplication_levels

# Assess sequence duplication
dup_data = duplication_levels("sample.fastq.gz")

print(f"Unique sequences: {dup_data['unique_sequences']}")
print(f"Duplicate sequences: {dup_data['duplicate_sequences']}")
print(f"Duplication percentage: {dup_data['duplication_percentage']:.2f}%")

# Duplication level distribution
for level, percentage in dup_data['duplication_levels'].items():
    print(f"{level}: {percentage:.2f}%")
```

### N Content per Position

```python
from metainformant.quality.fastq import n_content_per_position

# Analyze N content across sequence positions
n_content = n_content_per_position("sample.fastq.gz")

print(f"Positions analyzed: {len(n_content)}")
print(f"Max N content: {max(n_content.values()):.2f}%")

# Position-wise N content
for position, n_pct in n_content.items():
    if n_pct > 5:  # Report positions with >5% N
        print(f"Position {position}: {n_pct:.2f}% N")
```

### Quality Score Distribution

```python
from metainformant.quality.fastq import quality_score_distribution

# Overall quality score distribution
qual_dist = quality_score_distribution("sample.fastq.gz")

print("Quality Score Distribution:")
for score, count in sorted(qual_dist.items()):
    print(f"Q{score}: {count:,} ({count/sum(qual_dist.values())*100:.1f}%)")
```

## Batch Processing

### Multiple File Analysis

```python
from pathlib import Path

def analyze_multiple_fastq(input_dir, pattern="*.fastq.gz", output_dir="output/qc"):
    """Analyze multiple FASTQ files."""
    
    fastq_files = list(Path(input_dir).glob(pattern))
    results = {}
    
    for fastq_file in fastq_files:
        print(f"Analyzing {fastq_file.name}...")
        
        # Analyze each file
        result = analyze_fastq_quality(
            str(fastq_file),
            output_dir=f"{output_dir}/{fastq_file.stem}"
        )
        
        results[fastq_file.name] = result
    
    return results

# Batch analyze all FASTQ files
batch_results = analyze_multiple_fastq("data/raw_reads/")
```

### Paired-End Analysis

```python
def analyze_paired_end(r1_file, r2_file, output_dir="output/paired_qc"):
    """Analyze paired-end FASTQ files."""
    
    # Analyze R1 and R2 separately
    r1_results = analyze_fastq_quality(r1_file, output_dir=f"{output_dir}/R1")
    r2_results = analyze_fastq_quality(r2_file, output_dir=f"{output_dir}/R2")
    
    # Compare basic statistics
    print("Paired-End Comparison:")
    print(f"R1 sequences: {r1_results['basic_stats']['total_sequences']:,}")
    print(f"R2 sequences: {r2_results['basic_stats']['total_sequences']:,}")
    print(f"R1 mean quality: {r1_results['basic_stats']['mean_quality']:.2f}")
    print(f"R2 mean quality: {r2_results['basic_stats']['mean_quality']:.2f}")
    
    return {'R1': r1_results, 'R2': r2_results}

# Analyze paired files
paired_results = analyze_paired_end(
    "sample_R1.fastq.gz",
    "sample_R2.fastq.gz"
)
```

## Quality Interpretation

### Quality Score Guidelines

```python
def interpret_quality_scores(results):
    """Provide interpretation of quality metrics."""
    
    mean_quality = results['basic_stats']['mean_quality']
    
    if mean_quality >= 35:
        quality_status = "Excellent"
    elif mean_quality >= 30:
        quality_status = "Good"
    elif mean_quality >= 25:
        quality_status = "Acceptable"
    else:
        quality_status = "Poor"
    
    print(f"Overall Quality: {quality_status} (Q{mean_quality:.1f})")
    
    # Check for specific issues
    issues = []
    
    if results['adapter_content']['max_adapter_content'] > 20:
        issues.append("High adapter content detected")
    
    if results['basic_stats']['n_content'] > 5:
        issues.append("High N content")
    
    if results['duplication_levels']['duplication_percentage'] > 50:
        issues.append("High sequence duplication")
    
    # Check per-base quality
    per_base = results.get('per_base_quality', {})
    if per_base and any(pos['mean'] < 25 for pos in per_base.values()):
        issues.append("Low quality at some positions")
    
    if issues:
        print("Potential Issues:")
        for issue in issues:
            print(f"  - {issue}")
    else:
        print("No major quality issues detected")

# Interpret results
interpret_quality_scores(results)
```

### GC Content Assessment

```python
def assess_gc_content(results, expected_gc=None):
    """Assess GC content against expected values."""
    
    observed_gc = results['basic_stats']['gc_content']
    
    if expected_gc is None:
        # Common organism GC contents
        organisms = {
            'human': (40, 42),
            'mouse': (40, 42),
            'yeast': (38, 40),
            'e_coli': (50, 52),
            'arabidopsis': (35, 37)
        }
        
        print(f"Observed GC content: {observed_gc:.1f}%")
        print("Expected ranges for common organisms:")
        for org, (low, high) in organisms.items():
            status = "✓" if low <= observed_gc <= high else "✗"
            print(f"  {status} {org}: {low}-{high}%")
    else:
        deviation = abs(observed_gc - expected_gc)
        if deviation <= 2:
            print(f"GC content within expected range: {observed_gc:.1f}% (expected {expected_gc}%)")
        else:
            print(f"Warning: GC content deviation: {observed_gc:.1f}% vs expected {expected_gc}%")

# Assess GC content
assess_gc_content(results, expected_gc=41)  # Human genome
```

## Visualization

### Quality Plots

```python
import matplotlib.pyplot as plt
import numpy as np

def plot_per_base_quality(results):
    """Plot per-base quality scores."""
    
    per_base = results.get('per_base_quality', {})
    if not per_base:
        print("No per-base quality data available")
        return
    
    positions = []
    means = []
    q1s = []
    q3s = []
    
    for pos, data in sorted(per_base.items(), key=lambda x: int(x[0])):
        positions.append(int(pos))
        means.append(data['mean'])
        q1s.append(data['q1'])
        q3s.append(data['q3'])
    
    plt.figure(figsize=(12, 6))
    
    # Plot quartile range
    plt.fill_between(positions, q1s, q3s, alpha=0.3, color='lightblue', 
                     label='Interquartile Range')
    
    # Plot mean quality
    plt.plot(positions, means, 'b-', linewidth=2, label='Mean Quality')
    
    # Quality thresholds
    plt.axhline(y=30, color='g', linestyle='--', alpha=0.7, label='Q30 (99.9% accuracy)')
    plt.axhline(y=20, color='orange', linestyle='--', alpha=0.7, label='Q20 (99% accuracy)')
    plt.axhline(y=10, color='r', linestyle='--', alpha=0.7, label='Q10 (90% accuracy)')
    
    plt.xlabel('Position in Read')
    plt.ylabel('Quality Score')
    plt.title('Per Base Quality Scores')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()

def plot_gc_content(results):
    """Plot GC content distribution."""
    
    gc_dist = results.get('gc_content_distribution', {})
    if not gc_dist:
        print("No GC content distribution data available")
        return
    
    gc_values = list(gc_dist.keys())
    counts = list(gc_dist.values())
    
    plt.figure(figsize=(10, 6))
    plt.hist(gc_values, bins=50, weights=counts, alpha=0.7, color='skyblue', edgecolor='black')
    
    # Mark mean GC content
    mean_gc = results['basic_stats']['gc_content']
    plt.axvline(mean_gc, color='red', linestyle='--', linewidth=2, 
                label=f'Mean GC: {mean_gc:.1f}%')
    
    plt.xlabel('GC Content (%)')
    plt.ylabel('Number of Sequences')
    plt.title('GC Content Distribution')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()

# Generate plots
plot_per_base_quality(results)
plot_gc_content(results)
```

## Advanced Features

### Custom Adapter Detection

```python
def detect_custom_adapters(fastq_path, adapter_file):
    """Detect custom adapters from a file."""
    
    # Load custom adapters
    custom_adapters = []
    with open(adapter_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                custom_adapters.append(line.strip())
    
    # Run adapter analysis
    results = adapter_content(fastq_path, adapter_sequences=custom_adapters)
    return results

# Example adapter file format:
# AGATCGGAAGAGCACACGTCTGAACTCCAGTCA    # Illumina TruSeq
# AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT    # Illumina TruSeq R2
# TGGAATTCTCGGGTGCCAAGG                # Small RNA 3' adapter
```

### Subsampling Large Files

```python
def subsample_analysis(large_fastq, sample_sizes=[10000, 50000, 100000]):
    """Compare quality metrics across different sample sizes."""
    
    results = {}
    
    for size in sample_sizes:
        print(f"Analyzing {size:,} sequences...")
        result = analyze_fastq_quality(large_fastq, subsample_size=size)
        results[size] = result
    
    # Compare basic statistics
    print("\nSample Size Comparison:")
    print("Size\t\tMean Quality\tGC Content\tN Content")
    for size, result in results.items():
        stats = result['basic_stats']
        print(f"{size:,}\t\t{stats['mean_quality']:.2f}\t\t{stats['gc_content']:.1f}%\t\t{stats['n_content']:.1f}%")
    
    return results

# Compare different sample sizes
subsample_results = subsample_analysis("very_large_file.fastq.gz")
```

### Quality Summary Report

```python
def generate_quality_report(results, output_file="quality_report.txt"):
    """Generate a comprehensive text report."""
    
    with open(output_file, 'w') as f:
        f.write("FASTQ Quality Analysis Report\n")
        f.write("=" * 40 + "\n\n")
        
        # Basic statistics
        stats = results['basic_stats']
        f.write("Basic Statistics:\n")
        f.write(f"  Filename: {stats['filename']}\n")
        f.write(f"  Total sequences: {stats['total_sequences']:,}\n")
        f.write(f"  Sequence length: {stats['sequence_length']}\n")
        f.write(f"  GC content: {stats['gc_content']:.1f}%\n")
        f.write(f"  Mean quality: {stats['mean_quality']:.2f}\n\n")
        
        # Quality assessment
        f.write("Quality Assessment:\n")
        
        # Adapter content
        adapter = results['adapter_content']
        f.write(f"  Max adapter content: {adapter['max_adapter_content']:.2f}%\n")
        
        # Duplication
        dup = results['duplication_levels']
        f.write(f"  Sequence duplication: {dup['duplication_percentage']:.2f}%\n")
        
        # N content
        f.write(f"  N content: {stats['n_content']:.2f}%\n\n")
        
        # Recommendations
        f.write("Recommendations:\n")
        if stats['mean_quality'] < 28:
            f.write("  - Consider quality filtering (low mean quality)\n")
        if adapter['max_adapter_content'] > 5:
            f.write("  - Adapter trimming recommended\n")
        if dup['duplication_percentage'] > 20:
            f.write("  - High duplication detected - check for PCR artifacts\n")
        if stats['n_content'] > 5:
            f.write("  - High N content - investigate sequencing issues\n")
        
        if (stats['mean_quality'] >= 28 and 
            adapter['max_adapter_content'] <= 5 and
            dup['duplication_percentage'] <= 20 and
            stats['n_content'] <= 5):
            f.write("  - Data quality appears good for downstream analysis\n")

# Generate report
generate_quality_report(results)
```

## Performance Optimization

### Memory-Efficient Processing

```python
def analyze_large_fastq_efficient(fastq_path, chunk_size=100000):
    """Memory-efficient analysis of very large FASTQ files."""
    
    # Process in chunks to manage memory
    total_sequences = 0
    quality_sum = 0
    gc_sum = 0
    
    # Initialize accumulators
    per_base_quality = {}
    sequence_lengths = {}
    
    # Process file in chunks
    with gzip.open(fastq_path, 'rt') if fastq_path.endswith('.gz') else open(fastq_path, 'r') as f:
        sequences = []
        qualities = []
        
        for line_num, line in enumerate(f):
            if line_num % 4 == 1:  # Sequence line
                sequences.append(line.strip())
            elif line_num % 4 == 3:  # Quality line
                qualities.append(line.strip())
                
                # Process chunk when full
                if len(sequences) >= chunk_size:
                    process_chunk(sequences, qualities, per_base_quality, sequence_lengths)
                    total_sequences += len(sequences)
                    sequences = []
                    qualities = []
        
        # Process final chunk
        if sequences:
            process_chunk(sequences, qualities, per_base_quality, sequence_lengths)
            total_sequences += len(sequences)
    
    print(f"Processed {total_sequences:,} sequences efficiently")

def process_chunk(sequences, qualities, per_base_acc, length_acc):
    """Process a chunk of sequences."""
    # Implementation would update accumulators
    pass
```

## Integration Examples

### Pipeline Integration

```python
def qc_guided_preprocessing(fastq_file, output_dir="output/processed"):
    """Use QC results to guide preprocessing parameters."""
    
    # 1. Run quality control
    qc_results = analyze_fastq_quality(fastq_file)
    
    # 2. Set parameters based on QC
    mean_qual = qc_results['basic_stats']['mean_quality']
    adapter_content = qc_results['adapter_content']['max_adapter_content']
    
    # Quality filtering threshold
    if mean_qual >= 35:
        min_quality = 25
    elif mean_qual >= 30:
        min_quality = 20
    else:
        min_quality = 15
    
    # Adapter trimming
    trim_adapters = adapter_content > 1
    
    print(f"QC-guided parameters:")
    print(f"  Minimum quality: Q{min_quality}")
    print(f"  Adapter trimming: {'Yes' if trim_adapters else 'No'}")
    
    # 3. Apply preprocessing with optimized parameters
    from metainformant.dna.fastq import filter_fastq_quality
    
    output_file = f"{output_dir}/filtered_{Path(fastq_file).name}"
    
    filtered = filter_fastq_quality(
        fastq_file,
        output_file,
        min_quality=min_quality,
        trim_adapters=trim_adapters
    )
    
    return filtered

# Example usage
processed_file = qc_guided_preprocessing("raw_sample.fastq.gz")
```

## Testing

Comprehensive tests are available in `tests/test_quality_fastq.py`:

```bash
# Run FASTQ quality tests
uv run pytest tests/test_quality_fastq.py -v

# Test specific functions
uv run pytest tests/test_quality_fastq.py::test_basic_statistics -v
uv run pytest tests/test_quality_fastq.py::test_per_base_quality -v
```

## Related Documentation

- [Quality Control Overview](./index.md): Module overview and architecture
- [DNA FASTQ Processing](../dna/fastq.md): FASTQ sequence manipulation
- [Core I/O Utilities](../core/io.md): File handling and compression
- [Visualization](../visualization/plots.md): Quality plot generation
