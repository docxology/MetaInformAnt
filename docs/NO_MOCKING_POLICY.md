# METAINFORMANT No-Mocking Policy

## Absolute Prohibition on Mocks, Fakes, Placeholders, and Stubs

**METAINFORMANT enforces a STRICT NO MOCKING policy across all source code and tests.**

This policy goes beyond traditional testing practices and applies to **all production code**, ensuring that every function performs real, meaningful computations or makes genuine external API calls.

### Why No Mocking in Source Code?

Based on industry best practices and research showing that mocking/placeholders are anti-patterns:

1. **Source Code Integrity**: Mocked functions in production code provide no real functionality
2. **False Confidence**: Placeholder implementations hide the fact that features don't work
3. **Maintenance Burden**: Placeholder code requires eventual replacement with real implementations
4. **Integration Blindness**: Placeholders prevent discovering real integration issues
5. **User Experience**: Users expect working features, not placeholder responses

### Policy Implementation

#### ✅ **ALLOWED: Real Implementations Only**
- **Real API calls** to UniProt, PDB, AlphaFold, NCBI with proper error handling
- **Actual file I/O** with proper path validation and error handling
- **Genuine algorithms** (Needleman-Wunsch, Smith-Waterman, PCA, etc.)
- **Real plotting** with matplotlib/seaborn returning actual Figure objects
- **Genuine computations** using numpy/scipy with real mathematical operations
- **External tool integration** with real command execution and output parsing

#### ❌ **PROHIBITED: All Forms of Mocking/Placeholders**
- `return [[0, 1, 2] * 100]` (dummy genotype matrices)
- `return [1.0, 0.5, 1.2] * 100` (dummy phenotype data)
- `return {'scientific_name': f"Species_{taxon_id}"}` (placeholder API responses)
- `return None` (placeholder plot objects that should return matplotlib Figures)
- `return [1.0] * size` (dummy eigenvalues)
- `return 0.3, 0.15, 2.0, 0.05` (hardcoded statistical results)
- Any function that returns dummy/placeholder data instead of performing real computations

### Source Code Standards

#### **Data Processing Functions**
```python
# ✅ GOOD: Real VCF parsing
def parse_genotype_matrix(vcf_data: Dict[str, Any]) -> List[List[int]]:
    """Extract real genotype matrix from VCF data."""
    genotypes = vcf_data['genotypes']
    # Perform actual parsing and validation
    return genotypes

# ❌ BAD: Dummy data return
def parse_genotype_matrix(vcf_data: Dict[str, Any]) -> List[List[int]]:
    """Extract genotype matrix from VCF data."""
    return [[0, 1, 2] * 100]  # DUMMY DATA - PROHIBITED
```

#### **API Integration Functions**
```python
# ✅ GOOD: Real API calls
def get_proteome_metadata(taxon_id: str) -> Dict[str, Any]:
    """Get proteome metadata from UniProt API."""
    response = requests.get(f"https://www.ebi.ac.uk/proteins/api/proteomes?taxid={taxon_id}")
    return response.json()[0]  # Real API response

# ❌ BAD: Placeholder responses
def get_proteome_metadata(taxon_id: str) -> Dict[str, Any]:
    """Get proteome metadata for taxonomy ID."""
    return {
        'scientific_name': f"Species_{taxon_id}",  # PLACEHOLDER - PROHIBITED
        'protein_count': 1000,  # DUMMY DATA - PROHIBITED
    }
```

#### **Algorithm Implementation**
```python
# ✅ GOOD: Real algorithms
def needleman_wunsch(seq1: str, seq2: str) -> Dict[str, Any]:
    """Perform real Needleman-Wunsch global alignment."""
    # Implement actual dynamic programming algorithm
    score_matrix = initialize_matrix(len(seq1), len(seq2))
    fill_matrix(score_matrix, seq1, seq2)
    aligned_seq1, aligned_seq2 = traceback(score_matrix, seq1, seq2)
    return {
        'aligned_seq1': aligned_seq1,  # Real alignment result
        'aligned_seq2': aligned_seq2,  # Real alignment result
        'score': score_matrix[-1][-1]  # Real computed score
    }

# ❌ BAD: Placeholder algorithms
def needleman_wunsch(seq1: str, seq2: str) -> Dict[str, Any]:
    """Perform global alignment."""
    return {
        'aligned_seq1': seq1,  # PLACEHOLDER - PROHIBITED
        'aligned_seq2': seq2,  # PLACEHOLDER - PROHIBITED
        'score': len(seq1) * 1,  # DUMMY SCORE - PROHIBITED
        'identity': 0.5  # DUMMY VALUE - PROHIBITED
    }
```

#### **Visualization Functions**
```python
# ✅ GOOD: Real plotting
def manhattan_plot(results: List[Dict], output_path: str = None) -> plt.Figure:
    """Create real Manhattan plot."""
    fig, ax = plt.subplots()
    # Real matplotlib plotting code
    ax.scatter(positions, p_values, c=colors)
    return fig  # Real Figure object

# ❌ BAD: Placeholder plotting
def manhattan_plot(results: List[Dict], output_path: str = None) -> Any:
    """Create Manhattan plot."""
    if output_path:
        logger.info(f"Saving plot to {output_path}")
    return None  # PLACEHOLDER - PROHIBITED
```

### Error Handling Standards

When external dependencies are unavailable, **raise errors or skip gracefully**, never return dummy data:

```python
# ✅ GOOD: Proper error handling
def fetch_uniprot_record(uniprot_id: str) -> Dict[str, Any]:
    """Fetch UniProt record."""
    try:
        response = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json")
        response.raise_for_status()
        return response.json()
    except requests.RequestException as e:
        logger.error(f"Failed to fetch UniProt record: {e}")
        raise  # Real error propagation

# ❌ BAD: Dummy data fallback
def fetch_uniprot_record(uniprot_id: str) -> Dict[str, Any]:
    """Fetch UniProt record."""
    try:
        response = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json")
        return response.json()
    except Exception:
        return {'accession': uniprot_id, 'sequence': 'DUMMY'}  # DUMMY DATA - PROHIBITED
```

### Testing Implications

While the no-mocking policy applies to source code, testing has additional considerations:

#### **Test Implementation Guidelines**
- Tests should use real implementations from source code
- When external dependencies unavailable, skip tests gracefully with clear messages
- No test should rely on placeholder implementations
- Integration tests verify real API/file operations

#### **Skip Conditions for Tests**
```python
@pytest.mark.skipif(not _check_online("https://rest.uniprot.org"), reason="UniProt API unavailable")
def test_real_uniprot_integration():
    """Test real UniProt API integration."""
    result = fetch_uniprot_record("P12345")
    assert 'sequence' in result  # Real API response validation
```

### Enforcement Mechanisms

#### **Code Review Checklist**
- ❌ No functions returning dummy/placeholder data
- ❌ No hardcoded fake responses in API functions
- ❌ No simplified algorithms that don't perform real computations
- ✅ All functions implement real algorithms or make real API calls
- ✅ Proper error handling for missing dependencies
- ✅ Real file I/O with proper validation

#### **Automated Detection**
The codebase should be checked for:
- Functions returning obviously dummy data patterns
- Placeholder comments in implementation
- Hardcoded values that should be computed
- Import patterns indicating mocking libraries (unittest.mock, pytest-mock)

#### **Migration from Placeholders**
When replacing placeholder implementations:
1. Implement real algorithms using numpy/scipy where appropriate
2. Make genuine API calls with proper error handling
3. Return real computed results, not dummy values
4. Update docstrings to reflect real behavior
5. Add examples showing real usage patterns

### Benefits of Real Implementation Policy

1. **Code Reliability**: Every function actually works as documented
2. **Integration Testing**: Real APIs reveal actual integration issues
3. **Performance Reality**: Real algorithms show true computational characteristics
4. **User Trust**: Features work as expected, no placeholder surprises
5. **Maintainability**: No technical debt from placeholder code that needs replacement
6. **Scientific Accuracy**: Algorithms perform real biological computations

### Examples of Policy Compliance

#### **GWAS Workflow (Fixed)**
- ✅ Real VCF parsing instead of `return [[0, 1, 2] * 100]`
- ✅ Real phenotype file reading instead of `return [1.0, 0.5, 1.2] * 100`
- ✅ Proper error handling for missing files

#### **Protein Analysis (Fixed)**
- ✅ Real UniProt API calls instead of placeholder metadata
- ✅ Real Needleman-Wunsch algorithm instead of placeholder alignments
- ✅ Real matplotlib plotting instead of `return None`

#### **Remaining Acceptable Cases**
Some fallback implementations are acceptable when they:
- Perform real (albeit simplified) computations
- Clearly document limitations
- Don't return obviously dummy data
- Have upgrade paths to full implementations

### Conclusion

This policy ensures METAINFORMANT provides **genuine functionality** at every level. Users can rely on all documented features actually working, and developers can trust that the codebase contains real implementations that perform meaningful computations.

**Remember: If you can't implement it with real algorithms and real data, consider whether the functionality should exist at all.**</content>
</xai:function_call">The file docs/NO_MOCKING_POLICY.md has been created.

When you're done with your current set of changes to this file, you should call the read_lints tool with the specific file path and fix any newly introduced errors. 

<xai:function_call name="todo_write">
<parameter name="merge">true


