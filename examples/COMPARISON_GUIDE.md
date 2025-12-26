# Examples Comparison Guide

This guide helps you choose the right METAINFORMANT example for your needs and understand the differences between similar examples.

## Examples vs Scripts

METAINFORMANT provides both examples (for learning) and scripts (for production). Choose based on your use case:

| Aspect | Examples | Scripts |
|--------|----------|---------|
| **Purpose** | Learning concepts | Production workflows |
| **Location** | `examples/` | `scripts/` |
| **Complexity** | Simple, focused | Full-featured, robust |
| **Data** | Sample/synthetic | Real data support |
| **Error Handling** | Basic | Comprehensive |
| **Documentation** | Educational | Production-ready |
| **Testing** | Automated validation | Integration testing |
| **Maintenance** | Educational accuracy | Production reliability |

### When to Use Examples
- Learning new METAINFORMANT concepts
- Understanding API usage patterns
- Testing ideas before production
- Educational purposes and tutorials
- Demonstrating specific algorithms

### When to Use Scripts
- Processing real biological data
- Batch analysis workflows
- Automated pipelines
- Production data analysis
- End-to-end research workflows

## Domain Examples Overview

### Core Examples (`examples/core/`)

Fundamental METAINFORMANT concepts and utilities.

| Example | Purpose | Complexity | Output |
|---------|---------|------------|---------|
| `example_config.py` | Configuration loading with environment overrides | Basic | JSON config |
| `example_io.py` | File I/O patterns (JSON, CSV, JSONL) | Basic | Multiple formats |
| `example_logging.py` | Structured logging setup and usage | Basic | Log files |
| `example_paths.py` | Path validation and containment checking | Basic | Path info |
| `example_workflow.py` | Basic workflow orchestration | Intermediate | Workflow config |

**Choose core examples when**:
- Learning METAINFORMANT fundamentals
- Setting up new projects
- Understanding configuration patterns

### DNA Analysis Examples (`examples/dna/`)

DNA sequence analysis and genomics workflows.

| Example | Analysis Type | Data Size | Key Features |
|---------|---------------|-----------|--------------|
| `example_sequences.py` | Basic sequence operations | Small (sample seqs) | FASTA reading, GC content, reverse complement |
| `example_alignment.py` | Sequence alignment | Small-medium | Global/local alignment, scoring |
| `example_phylogeny.py` | Phylogenetic analysis | Medium | Tree construction, Neighbor-Joining |
| `example_population.py` | Population genetics | Medium | Diversity stats, Fst, Tajima's D |

**Choose DNA examples when**:
- Learning sequence analysis concepts
- Comparing alignment algorithms
- Understanding phylogenetic methods
- Exploring population genetics

**Comparison of DNA examples**:
- **Sequences**: Start here for basic concepts
- **Alignment**: Most computationally intensive
- **Phylogeny**: Good for evolutionary concepts
- **Population**: Best for genetic diversity studies

### RNA Analysis Examples (`examples/rna/`)

Transcriptomic analysis and RNA-seq workflows.

| Example | Analysis Type | Dependencies | Key Features |
|---------|---------------|--------------|--------------|
| `example_amalgkit.py` | Full RNA-seq pipeline | amalgkit CLI | Species detection, quantification, merging |
| `example_quantification.py` | Expression quantification | Basic | Count matrices, normalization |

**Choose RNA examples when**:
- Learning RNA-seq workflows
- Understanding expression analysis
- Working with amalgkit pipelines

### GWAS Examples (`examples/gwas/`)

Genome-wide association studies.

| Example | Analysis Type | Data Type | Key Features |
|---------|---------------|-----------|--------------|
| `example_association.py` | Association testing | Genotype + phenotype | Linear/logistic regression, multiple testing |
| `example_visualization.py` | GWAS visualization | Association results | Manhattan plots, QQ plots |

**Choose GWAS examples when**:
- Learning association testing methods
- Understanding multiple testing correction
- Creating publication-ready plots

**Comparison**:
- **Association**: Focus on statistical methods
- **Visualization**: Focus on result interpretation

### ML Examples (`examples/ml/`)

Machine learning for biological data.

| Example | ML Task | Algorithm | Key Features |
|---------|---------|-----------|--------------|
| `example_pipeline.py` | Classification pipeline | Random Forest | Feature selection, cross-validation, evaluation |

**Choose ML examples when**:
- Learning ML workflows for biology
- Understanding feature selection
- Exploring model evaluation

## Integration Examples (`examples/integration/`)

Cross-domain and multi-omic analysis.

| Example | Integration Type | Domains | Complexity |
|---------|------------------|---------|------------|
| `example_multiomics.py` | Multi-omic correlation | DNA + RNA + Protein | Advanced |
| `example_dna_rna.py` | DNA-RNA integration | DNA + RNA | Intermediate |
| `example_complete_workflow.py` | End-to-end pipeline | Multiple | Advanced |

**Choose integration examples when**:
- Understanding multi-omic analysis
- Learning cross-domain workflows
- Building complex pipelines

## Specialized Domain Examples

### Protein Analysis (`examples/protein/`)
- `example_sequences.py`: Protein sequence analysis and properties

### Epigenome Analysis (`examples/epigenome/`)
- `example_methylation.py`: DNA methylation analysis

### Ontology Analysis (`examples/ontology/`)
- `example_go.py`: Gene Ontology enrichment and analysis

### Phenotype Analysis (`examples/phenotype/`)
- `example_traits.py`: Phenotypic trait analysis

### Ecology Analysis (`examples/ecology/`)
- `example_community.py`: Community diversity and composition

### Information Theory (`examples/information/`)
- `example_entropy.py`: Syntactic and semantic information measures

### Life Events (`examples/life_events/`)
- `example_events.py`: Life course event sequence analysis

### Math Biology (`examples/math/`)
- `example_dynamics.py`: Population dynamics and mathematical models

### Networks (`examples/networks/`)
- `example_networks.py`: Biological network analysis and visualization

### Multi-Omics (`examples/multiomics/`)
- `example_integration.py`: Cross-omics data integration

### Quality Control (`examples/quality/`)
- `example_qc.py`: Data quality assessment and validation

### Simulation (`examples/simulation/`)
- `example_simulation.py`: Synthetic data generation

### Single-Cell (`examples/singlecell/`)
- `example_scrna.py`: Single-cell RNA-seq analysis

### Visualization (`examples/visualization/`)
- `example_plots.py`: Biological data visualization

## Choosing the Right Example

### By Experience Level

#### Beginner (New to METAINFORMANT)
1. Start with `examples/core/example_config.py`
2. Try `examples/core/example_io.py`
3. Explore `examples/dna/example_sequences.py`
4. Learn visualization with `examples/visualization/example_plots.py`

#### Intermediate (Familiar with Basics)
1. Try domain-specific examples in your field
2. Explore integration examples
3. Learn ML with `examples/ml/example_pipeline.py`
4. Understand quality control

#### Advanced (Building Complex Workflows)
1. Study integration examples
2. Learn from complete workflow examples
3. Adapt examples for production use
4. Contribute new examples

### By Analysis Type

#### Sequence Analysis
- **DNA sequences**: `examples/dna/example_sequences.py`
- **Protein sequences**: `examples/protein/example_sequences.py`
- **Alignment**: `examples/dna/example_alignment.py`

#### Population Genetics
- **Diversity statistics**: `examples/dna/example_population.py`
- **GWAS**: `examples/gwas/example_association.py`
- **Phylogenetics**: `examples/dna/example_phylogeny.py`

#### Expression Analysis
- **RNA-seq quantification**: `examples/rna/example_quantification.py`
- **Full RNA-seq pipeline**: `examples/rna/example_amalgkit.py`
- **Single-cell**: `examples/singlecell/example_scrna.py`

#### Systems Biology
- **Network analysis**: `examples/networks/example_networks.py`
- **Multi-omics**: `examples/integration/example_multiomics.py`
- **Pathway analysis**: `examples/ontology/example_go.py`

#### Data Quality
- **General QC**: `examples/quality/example_qc.py`
- **Sequence validation**: `examples/dna/example_sequences.py`
- **Statistical checks**: `examples/math/example_dynamics.py`

### By Computational Requirements

#### Fast Examples (< 5 seconds)
- `examples/core/example_config.py`
- `examples/core/example_paths.py`
- `examples/core/example_logging.py`
- `examples/dna/example_sequences.py`

#### Medium Examples (5-30 seconds)
- `examples/dna/example_alignment.py`
- `examples/gwas/example_association.py`
- `examples/ml/example_pipeline.py`
- `examples/visualization/example_plots.py`

#### Slow Examples (> 30 seconds)
- `examples/dna/example_phylogeny.py`
- `examples/rna/example_amalgkit.py`
- `examples/integration/example_multiomics.py`

## Performance Comparison

### Execution Time by Domain

```
Core Examples     | ████░░░░░░░░░░  (15-30 seconds total)
DNA Examples      | ████████░░░░░░  (30-60 seconds total)
RNA Examples      | ███████████░░░  (60-120 seconds total)
GWAS Examples     | ████████░░░░░░  (30-60 seconds total)
ML Examples       | ███████░░░░░░░  (25-45 seconds total)
Integration       | ██████████████  (90-180 seconds total)
```

### Memory Usage Patterns

| Example Type | Memory Usage | Notes |
|-------------|--------------|--------|
| Core | Low (< 100MB) | Basic data structures |
| DNA | Medium (100-500MB) | Sequence data, alignments |
| RNA | High (500MB+) | Expression matrices |
| GWAS | Medium-High | Genotype matrices |
| ML | Medium | Model training data |
| Networks | Medium | Graph structures |
| Integration | High | Multi-omic datasets |

## Output Format Comparison

All examples follow consistent output patterns:

### Standard JSON Structure
```json
{
  "example": "example_name",
  "domain": "domain_name",
  "description": "Brief description",
  "timestamp": "2024-01-01T12:00:00Z",
  "results": {
    "analysis_type": "specific_analysis",
    "data": {...},
    "statistics": {...}
  }
}
```

### Domain-Specific Variations

#### DNA Examples
- Include sequence metadata
- Report sequence statistics
- May include alignment scores

#### GWAS Examples
- Statistical test results
- P-values and effect sizes
- Multiple testing corrections

#### ML Examples
- Model performance metrics
- Cross-validation results
- Feature importance scores

## Dependencies Comparison

### Minimal Dependencies
- Core examples: Only METAINFORMANT core
- Quality examples: Basic validation libraries
- Ontology examples: NetworkX for graph operations

### Heavy Dependencies
- RNA examples: amalgkit CLI tool
- ML examples: scikit-learn, numpy
- GWAS examples: pandas, matplotlib
- Single-cell examples: scanpy, anndata

### Optional Dependencies
- Visualization examples: matplotlib, seaborn, plotly
- Network examples: networkx, igraph
- Epigenome examples: pyBigWig, pysam

## Testing and Validation

### Automated Testing Coverage

| Domain | Test Status | Notes |
|--------|-------------|--------|
| Core | ✅ Full coverage | All examples tested |
| DNA | ✅ Full coverage | Sequence validation |
| RNA | ⚠️ Partial | Depends on amalgkit |
| GWAS | ✅ Full coverage | Statistical validation |
| ML | ✅ Full coverage | Performance metrics |
| Integration | ✅ Full coverage | Cross-validation |

### Quality Metrics

All examples are validated for:
- ✅ **Syntax correctness**
- ✅ **Import availability**
- ✅ **Output structure**
- ✅ **Data validity**
- ✅ **Performance bounds**
- ✅ **Documentation completeness**

## Migration Guide

### From Old API Versions

If examples fail due to API changes:

1. **Check error messages** for deprecated function names
2. **Update imports** to new module locations
3. **Modify function calls** to match new signatures
4. **Test updated example** with validation script
5. **Update documentation** if parameter names changed

### Adapting Examples for Production

To convert examples to production scripts:

1. **Add configuration files** (YAML/JSON)
2. **Implement proper error handling**
3. **Add logging and progress tracking**
4. **Support real data inputs**
5. **Add validation and quality checks**
6. **Create command-line interfaces**
7. **Add comprehensive documentation**

## Contributing New Examples

### When to Add Examples
- **New analysis methods** not covered
- **Improved implementations** of existing methods
- **Domain-specific workflows** missing
- **Educational content** for complex topics
- **Integration patterns** not demonstrated

### Example Contribution Process
1. **Check existing examples** for similar functionality
2. **Use example generator** for consistent structure
3. **Test thoroughly** with validation scripts
4. **Update documentation** (README, dependencies)
5. **Add to comparison guide** if significantly different
6. **Submit for review** with complete documentation

## Future Development

### Planned Improvements
- **Jupyter notebook versions** for interactive learning
- **Web-based example runner** for online testing
- **Performance comparison dashboard**
- **Automated example updates** for API changes
- **Example difficulty ratings** and learning paths

### Community Requests
- **Domain-specific examples** for specialized fields
- **Real dataset examples** (with appropriate licensing)
- **Performance optimization examples**
- **Cloud deployment examples**
- **Reproducible research examples**

## Quick Reference

### Most Popular Examples
1. `examples/core/example_config.py` - Configuration basics
2. `examples/dna/example_sequences.py` - Sequence analysis
3. `examples/gwas/example_association.py` - GWAS statistics
4. `examples/ml/example_pipeline.py` - ML workflows
5. `examples/visualization/example_plots.py` - Data visualization

### Most Educational Examples
1. `examples/core/example_workflow.py` - Workflow concepts
2. `examples/integration/example_multiomics.py` - Integration patterns
3. `examples/information/example_entropy.py` - Information theory
4. `examples/math/example_dynamics.py` - Mathematical modeling
5. `examples/networks/example_networks.py` - Network analysis

### Most Comprehensive Examples
1. `examples/rna/example_amalgkit.py` - Full RNA-seq pipeline
2. `examples/integration/example_complete_workflow.py` - End-to-end analysis
3. `examples/life_events/example_events.py` - Complex event modeling
4. `examples/singlecell/example_scrna.py` - Single-cell workflows
