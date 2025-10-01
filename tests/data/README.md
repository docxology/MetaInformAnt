# Test Data Directory

This directory contains test fixtures, sample datasets, and reference data used by METAINFORMANT's comprehensive test suite.

## Organization

Test data is organized by biological domain to mirror the source code structure:

```
tests/data/
├── dna/                    # DNA sequence and genomic test data
├── rna/                    # RNA expression and transcriptomic data
├── protein/                # Protein sequence and structure data
├── epigenome/              # Epigenetic modification test data
├── ontology/               # Ontology and annotation test files
├── phenotype/              # Phenotypic trait test data
└── rna/curate/             # Curated RNA datasets for integration testing
    └── Apis_mellifera/     # Species-specific test datasets
        └── tables/         # Tabular data for workflow testing
```

## Data Types

### DNA Test Data (`dna/`)
- **Sequence Files**: Various DNA sequence formats (FASTA, FASTQ)
- **Genomic Data**: Reference genomes and annotations
- **Alignment Files**: Multiple sequence alignments
- **Population Data**: Genetic variation and diversity data

### RNA Test Data (`rna/`)
- **Expression Matrices**: Gene expression count data
- **Transcript Annotations**: GTF/GFF annotation files
- **Quality Metrics**: RNA-seq quality assessment data
- **Workflow Outputs**: Expected results from analysis pipelines

### Protein Test Data (`protein/`)
- **Sequence Databases**: Protein sequence collections
- **Structure Files**: PDB format structure data
- **Annotation Data**: Functional annotation datasets
- **Domain Data**: Protein domain and family classifications

### Epigenome Test Data (`epigenome/`)
- **Methylation Data**: DNA methylation call files
- **Chromatin Data**: ATAC-seq and ChIP-seq peak data
- **Accessibility Data**: Chromatin accessibility measurements

### Ontology Test Data (`ontology/`)
- **GO Annotations**: Gene Ontology term associations
- **OBO Files**: Ontology definition files
- **Mapping Data**: Cross-ontology relationship data

### Phenotype Test Data (`phenotype/`)
- **Trait Measurements**: Morphological and behavioral trait data
- **Species Databases**: Taxonomic and phenotypic databases
- **Measurement Standards**: Trait quantification protocols

## Usage Guidelines

### For Test Development
- **Representative Data**: Use realistic but minimal datasets
- **Edge Cases**: Include boundary conditions and error cases
- **Format Coverage**: Test multiple file formats and versions
- **Size Optimization**: Keep files small but representative

### For Data Addition
1. **Organize by Domain**: Place data in appropriate domain subdirectory
2. **Descriptive Naming**: Use clear, descriptive filenames
3. **Documentation**: Include data source and format information
4. **Attribution**: Credit original data sources where applicable

### For Maintenance
- **Regular Updates**: Refresh test data as formats evolve
- **Version Control**: Track changes to test datasets
- **Compression**: Compress large files to save space
- **Validation**: Ensure data integrity and format correctness

## Integration with Testing

### Automated Testing
```python
# Example: Loading test data in pytest fixtures
@pytest.fixture
def sample_dna_sequences():
    return load_fasta("tests/data/dna/sample_sequences.fasta")

def test_sequence_analysis(sample_dna_sequences):
    # Test implementation using sample data
    pass
```

### Data Validation
- All test data should be validated for format correctness
- Include checksums for large files where appropriate
- Document any preprocessing or normalization applied
- Test data should be stable and reproducible

## Contributing Test Data

### Submission Process
1. Ensure data is appropriately licensed for use
2. Include comprehensive documentation
3. Test with existing test suite
4. Provide rationale for data inclusion

### Quality Standards
- **Relevance**: Data should test specific functionality
- **Correctness**: Data must be accurate and validated
- **Minimalism**: Use smallest representative datasets
- **Documentation**: Clear description of data purpose and format

This test data collection supports comprehensive testing of all METAINFORMANT functionality while maintaining manageable repository size.
