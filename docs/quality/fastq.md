# FASTQ Quality Analysis

The FASTQ quality analysis module provides quality assessment for high-throughput sequencing data, offering FastQC-like functionality with performance and integration.

## Main Analysis Function

### analyze_fastq_quality()

Perform comprehensive quality analysis on FASTQ files:

```python-snippet
from metainformant.quality.io.fastq import analyze_fastq_quality

# Basic analysis
results = analyze_fastq_quality("sample.fastq.gz")

# Analyze only the first 100000 reads of a large file
results = analyze_fastq_quality("sample.fastq.gz", n_reads=100000)

# Print summary
print("Analysis complete for sample.fastq.gz")
print(f"Total sequences: {results['basic_statistics']['total_reads']:,}")
print(f"Mean quality: {results['basic_statistics']['mean_quality']:.2f}")
```

**Parameters:**
- `fastq_path`: Path to FASTQ file (supports .gz compression)
- `n_reads`: Number of reads to analyze from the start of the file (None = all)

**Returns:** Dictionary with one entry per analysis module: `basic_statistics`, `per_base_quality`, `per_sequence_quality`, `sequence_length_distribution`, `gc_content_distribution`, `adapter_content`, `overrepresented_sequences`, `duplication_levels`, `n_content_per_position`, `quality_score_distribution`

## Analysis Modules

Each analysis is also available as a standalone function. These take a list of `FastqRecord` objects (from `read_fastq_records`) rather than a file path, and return the same dictionary that `analyze_fastq_quality()` stores under the matching top-level key.

### Basic Statistics

```python-snippet
from metainformant.quality.io.fastq import basic_statistics, read_fastq_records

# Get basic statistics
records = list(read_fastq_records("sample.fastq.gz"))
stats = basic_statistics(records)
print(stats)
```

**Output:**
```python
{
    'total_reads': 1500000,
    'total_bases': 226500000,
    'min_length': 151,
    'max_length': 151,
    'mean_length': 151,
    'median_length': 151,
    'length_std': 0.0,
    'min_quality': 14.2,
    'max_quality': 38.9,
    'mean_quality': 35.2,
    'median_quality': 35.9,
    'quality_std': 2.8,
    'min_gc': 12.6,
    'max_gc': 78.1,
    'mean_gc': 42.5,
    'median_gc': 42.3,
    'gc_std': 3.9
}
```

### Per-Base Quality Scores

```python-snippet
from metainformant.quality.io.fastq import per_base_quality, read_fastq_records

# Analyze quality scores by position
records = list(read_fastq_records("sample.fastq.gz"))
positions = per_base_quality(records)["positions"]
print(f"Positions analyzed: {len(positions)}")  # one entry per read position

# Each entry: position, mean, median, min, max, count, lower_quartile, upper_quartile
```

### Per-Sequence Quality Scores

```python-snippet
from metainformant.quality.io.fastq import per_sequence_quality, read_fastq_records

# Quality score distribution across sequences
records = list(read_fastq_records("sample.fastq.gz"))
bins = per_sequence_quality(records)["bins"]
print(f"Quality bins: {len(bins)}")

# Bins of mean read quality; each has bin_start, bin_end, count, percentage
```

### Sequence Length Distribution

```python-snippet
from metainformant.quality.io.fastq import read_fastq_records, sequence_length_distribution

# Distribution of sequence lengths
records = list(read_fastq_records("sample.fastq.gz"))
length_dist = sequence_length_distribution(records)["distribution"]
lengths = [entry["length"] for entry in length_dist]
print(f"Length range: {min(lengths)} - {max(lengths)}")

# Each entry: length, count, percentage
```

### GC Content Distribution

```python-snippet
from metainformant.quality.io.fastq import gc_content_distribution, read_fastq_records

# GC content distribution across sequences
records = list(read_fastq_records("sample.fastq.gz"))
gc_dist = gc_content_distribution(records)
print(f"GC content range: 0-100%")

# Returns {"bins": [...]}: 5% GC bins, each with bin_start, bin_end, count, percentage
```

### Adapter Content Analysis

```python-snippet
from metainformant.quality.io.fastq import adapter_content, read_fastq_records

# Detect adapter contamination
records = list(read_fastq_records("sample.fastq.gz"))
adapter_data = adapter_content(
    records,
    adapters=['AGATCGGAAGAG']  # Illumina adapter
)

# Highest adapter percentage at any sampled position
max_adapter = max(
    (p["percentage"] for a in adapter_data["adapters"].values() for p in a["positions"]),
    default=0.0,
)
print(f"Max adapter content: {max_adapter:.2f}%")
print(f"Adapters checked: {len(adapter_data['adapters'])}")

# Position-wise adapter content ("position" is a percentage of read length)
for name, adapter in adapter_data['adapters'].items():
    print(f"{name} ({adapter['adapter_sequence']}):")
    for p in adapter['positions']:
        print(f"  {p['position']}% of read length: {p['percentage']:.2f}%")
```

**Built-in Adapters** (`quality.io.fastq.ILLUMINA_ADAPTERS`, used when `adapters` is omitted):
- TruSeq Universal Adapter: `AGATCGGAAGAG`
- TruSeq Adapter, Read 1: `GATCGGAAGAG`
- TruSeq Adapter, Read 2: `AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG`
- Nextera Transposase Sequence: `CTGTCTCTTAT`

Other adapters, such as the small RNA 3' adapter `TGGAATTCTCGG`, must be passed explicitly.

### Overrepresented Sequences

```python-snippet
from metainformant.quality.io.fastq import overrepresented_sequences, read_fastq_records

# Find highly abundant sequences
records = list(read_fastq_records("sample.fastq.gz"))
overrep = overrepresented_sequences(
    records,
    min_length=20       # Only consider reads at least 20 bp long
)

# Reports up to 10 sequences that each make up more than 1% of reads
for seq_info in overrep["overrepresented"]:
    print(f"Sequence: {seq_info['sequence'][:30]}...")
    print(f"Count: {seq_info['count']}")
    print(f"Percentage: {seq_info['percentage']:.3f}%")
    print()
```

### Duplication Levels

```python-snippet
from metainformant.quality.io.fastq import duplication_levels, read_fastq_records

# Assess sequence duplication
records = list(read_fastq_records("sample.fastq.gz"))
dup_data = duplication_levels(records)

duplicate_reads = dup_data['total_sequences'] - dup_data['unique_sequences']
print(f"Unique sequences: {dup_data['unique_sequences']}")
print(f"Duplicate sequences: {duplicate_reads}")
print(f"Duplication rate: {dup_data['duplication_rate']:.2f}%")

# Duplication level distribution: how many distinct sequences occur once, twice, ...
for level, n_sequences in dup_data['duplication_bins'].items():
    print(f"Occurs {level} time(s): {n_sequences:,} distinct sequences")
```

### N Content per Position

```python-snippet
from metainformant.quality.io.fastq import n_content_per_position, read_fastq_records

# Analyze N content across sequence positions
records = list(read_fastq_records("sample.fastq.gz"))
n_content = n_content_per_position(records)
n_positions = n_content["positions"]

max_n = max((p["n_percentage"] for p in n_positions), default=0.0)
print(f"Positions analyzed: {len(n_positions)}")
print(f"Max N content: {max_n:.2f}%")

# Position-wise N content
for p in n_positions:
    if p["n_percentage"] > 5:  # Report positions with >5% N
        print(f"Position {p['position']}: {p['n_percentage']:.2f}% N")
```

### Quality Score Distribution

```python-snippet
from metainformant.quality.io.fastq import quality_score_distribution, read_fastq_records

# Overall quality score distribution
records = list(read_fastq_records("sample.fastq.gz"))
qual_dist = quality_score_distribution(records)["distribution"]

print("Quality Score Distribution:")
for entry in qual_dist:
    print(f"Q{entry['quality_score']}: {entry['count']:,} ({entry['percentage']:.1f}%)")
```

## Batch Processing

### Multiple File Analysis

```python
from pathlib import Path

def analyze_multiple_fastq(input_dir, pattern="*.fastq.gz"):
    """Analyze multiple FASTQ files."""
    
    fastq_files = list(Path(input_dir).glob(pattern))
    results = {}
    
    for fastq_file in fastq_files:
        print(f"Analyzing {fastq_file.name}...")
        
        # Analyze each file
        result = analyze_fastq_quality(str(fastq_file))
        
        results[fastq_file.name] = result
    
    return results

# Batch analyze all FASTQ files
batch_results = analyze_multiple_fastq("data/raw_reads/")
```

### Paired-End Analysis

```python
def analyze_paired_end(r1_file, r2_file):
    """Analyze paired-end FASTQ files."""
    
    # Analyze R1 and R2 separately
    r1_results = analyze_fastq_quality(r1_file)
    r2_results = analyze_fastq_quality(r2_file)
    
    # Compare basic statistics
    print("Paired-End Comparison:")
    print(f"R1 sequences: {r1_results['basic_statistics']['total_reads']:,}")
    print(f"R2 sequences: {r2_results['basic_statistics']['total_reads']:,}")
    print(f"R1 mean quality: {r1_results['basic_statistics']['mean_quality']:.2f}")
    print(f"R2 mean quality: {r2_results['basic_statistics']['mean_quality']:.2f}")
    
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
    
    mean_quality = results['basic_statistics']['mean_quality']
    
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
    
    max_adapter = max(
        (p['percentage'] for a in results['adapter_content']['adapters'].values() for p in a['positions']),
        default=0.0,
    )
    if max_adapter > 20:
        issues.append("High adapter content detected")
    
    n_positions = results['n_content_per_position']['positions']
    max_n = max((p['n_percentage'] for p in n_positions), default=0.0)
    if max_n > 5:
        issues.append("High N content")
    
    if results['duplication_levels']['duplication_rate'] > 50:
        issues.append("High sequence duplication")
    
    # Check per-base quality
    per_base = results.get('per_base_quality', {}).get('positions', [])
    if per_base and any(pos['mean'] < 25 for pos in per_base):
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

```python-snippet
def assess_gc_content(results, expected_gc=None):
    """Assess GC content against expected values."""
    
    observed_gc = results['basic_statistics']['mean_gc']
    
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
            status = "match" if low <= observed_gc <= high else "     "
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
    
    per_base = results.get('per_base_quality', {}).get('positions', [])
    if not per_base:
        print("No per-base quality data available")
        return
    
    positions = []
    means = []
    q1s = []
    q3s = []
    
    for data in per_base:
        positions.append(data['position'])
        means.append(data['mean'])
        q1s.append(data['lower_quartile'])
        q3s.append(data['upper_quartile'])
    
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
    
    gc_bins = results.get('gc_content_distribution', {}).get('bins', [])
    if not gc_bins:
        print("No GC content distribution data available")
        return
    
    gc_values = [(b['bin_start'] + b['bin_end']) / 2 for b in gc_bins]  # bin centers
    counts = [b['count'] for b in gc_bins]
    
    plt.figure(figsize=(10, 6))
    plt.hist(gc_values, bins=20, range=(0, 100), weights=counts, alpha=0.7, color='skyblue', edgecolor='black')
    
    # Mark mean GC content
    mean_gc = results['basic_statistics']['mean_gc']
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
    records = list(read_fastq_records(fastq_path))
    results = adapter_content(records, adapters=custom_adapters)
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
        result = analyze_fastq_quality(large_fastq, n_reads=size)
        results[size] = result
    
    # Compare basic statistics
    print("\nSample Size Comparison:")
    print("Size\t\tMean Quality\tGC Content\tMax N Content")
    for size, result in results.items():
        stats = result['basic_statistics']
        n_positions = result['n_content_per_position']['positions']
        max_n = max((p['n_percentage'] for p in n_positions), default=0.0)
        print(f"{size:,}\t\t{stats['mean_quality']:.2f}\t\t{stats['mean_gc']:.1f}%\t\t{max_n:.1f}%")
    
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
        stats = results['basic_statistics']
        f.write("Basic Statistics:\n")
        f.write(f"  Total sequences: {stats['total_reads']:,}\n")
        f.write(f"  Mean sequence length: {stats['mean_length']:.1f}\n")
        f.write(f"  GC content: {stats['mean_gc']:.1f}%\n")
        f.write(f"  Mean quality: {stats['mean_quality']:.2f}\n\n")
        
        # Quality assessment
        f.write("Quality Assessment:\n")
        
        # Adapter content
        max_adapter = max(
            (p['percentage'] for a in results['adapter_content']['adapters'].values() for p in a['positions']),
            default=0.0,
        )
        f.write(f"  Max adapter content: {max_adapter:.2f}%\n")
        
        # Duplication
        dup = results['duplication_levels']
        f.write(f"  Sequence duplication: {dup['duplication_rate']:.2f}%\n")
        
        # N content
        n_positions = results['n_content_per_position']['positions']
        max_n = max((p['n_percentage'] for p in n_positions), default=0.0)
        f.write(f"  Max N content: {max_n:.2f}%\n\n")
        
        # Recommendations
        f.write("Recommendations:\n")
        if stats['mean_quality'] < 28:
            f.write("  - Consider quality filtering (low mean quality)\n")
        if max_adapter > 5:
            f.write("  - Adapter trimming recommended\n")
        if dup['duplication_rate'] > 20:
            f.write("  - High duplication detected - check for PCR artifacts\n")
        if max_n > 5:
            f.write("  - High N content - investigate sequencing issues\n")
        
        if (stats['mean_quality'] >= 28 and 
            max_adapter <= 5 and
            dup['duplication_rate'] <= 20 and
            max_n <= 5):
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

```python-snippet
def qc_guided_preprocessing(fastq_file, output_dir="output/processed"):
    """Use QC results to guide preprocessing parameters."""
    
    # 1. Run quality control
    qc_results = analyze_fastq_quality(fastq_file)
    
    # 2. Set parameters based on QC
    mean_qual = qc_results['basic_statistics']['mean_quality']
    
    # Quality filtering threshold
    if mean_qual >= 35:
        min_quality = 25
    elif mean_qual >= 30:
        min_quality = 20
    else:
        min_quality = 15
    
    print(f"QC-guided parameters:")
    print(f"  Minimum quality: Q{min_quality}")
    
    # 3. Apply preprocessing with optimized parameters
    from metainformant.quality.io.fastq import filter_reads
    
    output_file = f"{output_dir}/filtered_{Path(fastq_file).name}"
    
    filter_reads(fastq_file, output_file, min_quality=min_quality)
    
    return output_file

# Example usage
processed_file = qc_guided_preprocessing("raw_sample.fastq.gz")
```

## Testing

Comprehensive tests are available in `tests/quality/test_quality_fastq.py`:

```bash
# Run FASTQ quality tests
uv run pytest tests/quality/test_quality_fastq.py -v

# Test specific functions
uv run pytest tests/quality/test_quality_fastq.py::TestBasicStatistics::test_basic_statistics_calculation -v
uv run pytest tests/quality/test_quality_fastq.py::TestPerBaseQuality::test_per_base_quality_calculation -v
```

## Related Documentation

- [Quality Control Overview](./index.md): Module overview and architecture
- [DNA FASTQ Processing](../dna/fastq.md): FASTQ sequence manipulation
- [Core I/O Utilities](../core/io.md): File handling and compression
- [Visualization](../visualization/plots.md): Quality plot generation
