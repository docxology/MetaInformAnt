# Agent Directives: tests/data

## Role
Test fixture data directory containing REAL sample data for test validation.

## Directory Structure
- `dna/` - DNA sequence test fixtures (FASTA, FASTQ samples)
- `epigenome/` - Epigenomic test data (methylation, ChIP-seq)
- `ontology/` - Gene Ontology test files (OBO format)
- `phenotype/` - Phenotype test data
- `protein/` - Protein structure/sequence test files
- `rna/` - RNA-seq and amalgkit test data

## Rules and Constraints

### Real Data Only
- All fixtures contain REAL biological data, not synthetic placeholders
- Sample data is minimal but representative of actual file formats
- Data should be small enough for fast tests but valid for parsing

### File Formats
Test data must match production formats exactly:
- FASTA/FASTQ for sequences
- OBO for ontology
- VCF for variants
- BED/bedGraph for genomic intervals
- JSON/YAML for configurations

### Adding Test Data
When adding new test fixtures:
1. Use real data samples (anonymized if necessary)
2. Keep files minimal (<100KB when possible)
3. Document data source in file comments or README
4. Ensure fixtures exercise edge cases

## Usage in Tests
```python
from pathlib import Path

@pytest.fixture
def test_data_dir():
    return Path(__file__).parent / "data"

def test_with_fixture(test_data_dir):
    fasta_path = test_data_dir / "dna" / "sample.fasta"
    sequences = read_fasta(fasta_path)
```
