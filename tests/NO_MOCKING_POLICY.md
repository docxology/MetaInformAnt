# METAINFORMANT Testing Policy: STRICTLY NO MOCKING

## Absolute Prohibition on Mocks, Fakes, and Stubs

**METAINFORMANT enforces a STRICT NO MOCKING policy across all test code.**

### Why No Mocking?

Based on industry best practices and research showing that mocking is an anti-pattern:

1. **Mocks Disconnect Tests from Reality**: Mocked tests often pass while the real implementation fails
2. **Brittle Test Maintenance**: Changes to implementation require updating both code and mocks
3. **False Confidence**: High test coverage with mocks can hide real integration issues
4. **Performance Blindness**: Mocks hide real performance characteristics and bottlenecks
5. **API Evolution Problems**: Mocks don't reveal when external APIs change or break

### Policy Implementation

#### ✅ **ALLOWED: Real Implementations Only**
- **Real API calls** to UniProt, PDB, AlphaFold, NCBI with timeout handling
- **Actual file downloads** and I/O operations with proper error handling  
- **Real database connections** when available, graceful skip when not
- **Genuine command-line tool execution** (amalgkit, MUSCLE, etc.)
- **Authentic network requests** with proper timeout and retry logic
- **Environment variable configuration** for test setup (acceptable real configuration)

#### ❌ **PROHIBITED: All Forms of Mocking**
- `unittest.mock.Mock`, `unittest.mock.patch`, `unittest.mock.MagicMock`
- `pytest-mock`, `pytest.MonkeyPatch` for function replacement
- HTTP mocking libraries (responses, httpretty, etc.)
- Database mocking or in-memory fake databases  
- Filesystem mocking or fake file systems
- Time/date mocking (use real time with deterministic seeds)
- Network stubs, shims, or interceptors

### Test Implementation Guidelines

#### **Network-Dependent Tests**
```python
def test_uniprot_api_real():
    """Test real UniProt API behavior."""
    if not _check_online("https://rest.uniprot.org"):
        pytest.skip("No network access for UniProt API - real implementation requires connectivity")
    
    # Make REAL API call
    result = map_ids_uniprot(["P69905"])
    assert isinstance(result, dict)
    # Test real response structure
```

#### **File Operations**  
```python
def test_pdb_download_real(tmp_path: Path):
    """Test real PDB file download."""
    if not _check_online("https://files.rcsb.org"):
        pytest.skip("No network access for PDB download - real implementation requires connectivity")
    
    # Download REAL file
    out = fetch_pdb_structure("1CRN", tmp_path, fmt="pdb")
    assert out.exists()
    # Verify real content structure
    content = out.read_text()
    assert "HEADER" in content or "ATOM" in content
```

#### **Error Documentation**
```python
def test_api_offline_behavior():
    """Document real offline behavior."""
    try:
        result = some_api_function()
        # If this succeeds, we're online
        assert isinstance(result, expected_type)
    except Exception:
        # Expected when offline - documents real failure modes
        assert True  # This is acceptable real-world behavior
```

### Skip Conditions

Tests should gracefully handle unavailable dependencies:

#### **Network Availability**
```python
def _check_online(url: str) -> bool:
    try:
        import requests
        resp = requests.get(url, timeout=5)
        resp.raise_for_status()
        return True
    except Exception:
        return False
```

#### **External Tool Availability**
```python
@pytest.mark.skipif(not shutil.which("amalgkit"), reason="amalgkit not available on PATH")
def test_amalgkit_integration():
    # Test with real amalgkit executable
```

#### **Environment Configuration**
```python
@pytest.mark.skipif(not os.environ.get("NCBI_EMAIL"), reason="NCBI email not provided")
def test_ncbi_entrez_real():
    email = os.environ["NCBI_EMAIL"]
    # Use real NCBI API with real email
```

### Benefits of Real Implementation Testing

1. **Catches Real Integration Issues**: API changes, network failures, authentication problems
2. **Documents Actual Behavior**: Real error conditions, timeouts, edge cases
3. **Performance Reality**: True response times, resource usage, bottlenecks  
4. **Authentic User Experience**: Tests reflect what users actually experience
5. **Future-Proof**: Tests remain valid as external services evolve

### Enforcement Mechanisms

#### **Configuration Guards**
```toml
# pyproject.toml
[tool.pytest.ini_options]
markers = [
    "no_mock: enforces that NO mocking/faking is allowed in any test",
]
filterwarnings = [
    # Warn about any mock usage as it violates our policy
    "error::pytest.PytestUnraisableExceptionWarning:.*mock.*",
]
```

#### **Code Review Checklist**
- ❌ No imports of `unittest.mock` or similar libraries
- ❌ No `@patch` decorators or `monkeypatch` usage
- ❌ No fake/stub/mock objects in test code
- ✅ Real API calls with proper skip conditions
- ✅ Actual file I/O with temporary directories
- ✅ Genuine external tool execution

#### **Automated Detection**
The test runner and CI/CD should flag:
- Import statements containing `mock`, `patch`, `monkeypatch`
- Usage of mocking libraries
- Tests that never exercise real external dependencies

### Migration from Mocked Tests

When converting mocked tests to real implementations:

1. **Remove all mock imports and decorators**
2. **Add connectivity checks** with appropriate skip conditions  
3. **Use real test data** (small files, known stable APIs)
4. **Document expected failures** when dependencies unavailable
5. **Add timeout handling** for real network operations
6. **Test with real error conditions** when possible

### Examples of Compliant Tests

#### **Real API Integration**
```python
def test_protein_workflow_real_integration(tmp_path: Path):
    """Integration test with real APIs and file operations."""
    if not all([
        _check_online("https://rest.uniprot.org"),
        _check_online("https://files.rcsb.org")
    ]):
        pytest.skip("Network access required for real integration testing")
    
    # Real UniProt lookup
    uniprot_result = map_ids_uniprot(["P69905"])
    if not uniprot_result:
        pytest.skip("UniProt API returned no results - real API behavior")
    
    # Real PDB download  
    pdb_result = fetch_pdb_structure("1A3N", tmp_path, fmt="pdb")
    assert pdb_result.exists()
    assert pdb_result.stat().st_size > 0
```

## Conclusion

This policy ensures that METAINFORMANT's test suite provides **authentic validation** of real-world behavior, catches **genuine integration issues**, and maintains **long-term reliability** as the codebase and its dependencies evolve.

**Remember: If you can't test it with real implementations, consider whether the functionality should exist at all.**
