"""Test configuration and fixtures for the METAINFORMANT test suite.

This module provides shared fixtures, test configuration, and utilities
for the entire test suite following STRICT NO-MOCKING policy.

IMPORTANT: This test suite uses REAL implementations only. No mocks, fakes,
or stubs are allowed. All external APIs use real network calls with graceful
offline handling via skip conditions.
"""

from __future__ import annotations

import os
import random
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Generator, Iterator

import pytest

# Ensure metainformant package is importable by adding src to Python path
# This allows tests to run whether the package is installed or not
_REPO_ROOT = Path(__file__).resolve().parent.parent
_SRC_DIR = _REPO_ROOT / "src"
if _SRC_DIR.exists() and str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

# Test data directory
TEST_DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def test_data_dir() -> Path:
    """Provide path to test data directory."""
    return TEST_DATA_DIR


@pytest.fixture(scope="function")
def isolated_tmp_dir() -> Iterator[Path]:
    """Provide a clean temporary directory for each test function.

    This fixture ensures complete isolation between tests by providing
    a fresh temporary directory that is automatically cleaned up.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield Path(tmp_dir)


@pytest.fixture(scope="function")
def seeded_random() -> Iterator[random.Random]:
    """Provide a seeded random number generator for reproducible tests."""
    rng = random.Random(42)
    yield rng


@pytest.fixture(scope="function")
def sample_dna_sequences() -> dict[str, str]:
    """Provide sample DNA sequences for testing."""
    return {
        "short": "ATCGATCGATCG",
        "medium": "ATCGATCGATCGAATTCCGGAATTCCGG" * 3,
        "long": "ATCGATCGATCGAATTCCGGAATTCCGG" * 10,
        "gc_rich": "GCGCGCGCGCGC",
        "at_rich": "ATATATATATA",
        "with_ambiguous": "ATCGATCGNNCGATCG",
        "empty": "",
    }


@pytest.fixture(scope="function")
def sample_protein_sequences() -> dict[str, str]:
    """Provide sample protein sequences for testing."""
    return {
        "short": "MTEYKLVV",
        "medium": "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTI",
        "long": "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEY" * 3,
        "hydrophobic": "AILVFWYYY",
        "hydrophilic": "KRDENQST",
        "with_stops": "MTEYKL*VVVGA",
        "empty": "",
    }


@pytest.fixture(scope="function")
def clean_environment(monkeypatch) -> Iterator[None]:
    """Provide clean environment for tests by clearing relevant env vars."""
    env_vars_to_clear = [
        "NCBI_EMAIL",
        "PG_HOST",
        "PG_USER",
        "PG_PASSWORD",
        "PG_DATABASE",
        "PG_PORT",
    ]

    # Store original values
    original_values = {}
    for var in env_vars_to_clear:
        original_values[var] = os.environ.get(var)
        monkeypatch.delenv(var, raising=False)

    yield

    # Restore original values
    for var, value in original_values.items():
        if value is not None:
            monkeypatch.setenv(var, value)


@pytest.fixture(scope="session", autouse=True)
def setup_test_environment():
    """Set up test environment including PYTHONPATH and dependency checks.
    
    This fixture runs automatically for all tests and ensures:
    - PYTHONPATH includes src/ so metainformant can be imported
    - User local bin is in PATH for CLI tools
    - Environment variables are set appropriately
    - UV cache and venv directories are configured for filesystem compatibility
    """
    # Ensure src/ is in PYTHONPATH for module imports
    repo_root = Path(__file__).parent.parent
    src_dir = repo_root / "src"
    if src_dir.exists() and str(src_dir) not in sys.path:
        sys.path.insert(0, str(src_dir))
    
    # Ensure user local bin is in PATH (common for --user installs)
    user_bin = Path.home() / ".local" / "bin"
    if user_bin.exists():
        current_path = os.environ.get("PATH", "")
        if str(user_bin) not in current_path:
            os.environ["PATH"] = f"{user_bin}:{current_path}"
    
    # Set PYTHONPATH if not already set
    pythonpath = os.environ.get("PYTHONPATH", "")
    if str(src_dir) not in pythonpath:
        if pythonpath:
            os.environ["PYTHONPATH"] = f"{src_dir}:{pythonpath}"
        else:
            os.environ["PYTHONPATH"] = str(src_dir)
    
    # Auto-configure UV environment variables for filesystem compatibility
    # This handles FAT filesystems that don't support symlinks
    try:
        from metainformant.core.filesystem import get_uv_cache_dir, get_venv_location
        
        # Get UV cache directory (respects existing UV_CACHE_DIR if set)
        uv_cache_dir = get_uv_cache_dir(repo_root)
        if "UV_CACHE_DIR" not in os.environ:
            os.environ["UV_CACHE_DIR"] = str(uv_cache_dir)
        
        # Get venv location (respects existing UV_PROJECT_ENVIRONMENT or METAINFORMANT_VENV if set)
        venv_location = get_venv_location(repo_root)
        if "UV_PROJECT_ENVIRONMENT" not in os.environ:
            os.environ["UV_PROJECT_ENVIRONMENT"] = str(venv_location)
        # Also set METAINFORMANT_VENV for compatibility with existing scripts
        if "METAINFORMANT_VENV" not in os.environ:
            os.environ["METAINFORMANT_VENV"] = str(venv_location)
    except ImportError:
        # If metainformant.core.filesystem is not available (e.g., package not installed),
        # fall back to manual detection using the same logic
        try:
            import subprocess
            import platform
            
            # Detect filesystem type
            fs_type = "unknown"
            system = platform.system().lower()
            if system == "linux" or system == "darwin":
                try:
                    result = subprocess.run(
                        ["df", "-T", str(repo_root)],
                        capture_output=True,
                        text=True,
                        timeout=5,
                    )
                    if result.returncode == 0:
                        lines = result.stdout.strip().split("\n")
                        if len(lines) > 1:
                            parts = lines[1].split()
                            if len(parts) >= 2:
                                fs_type = parts[1].lower()
                except Exception:
                    pass
            
            # Check if FAT filesystem (no symlink support)
            no_symlink_fs = {"exfat", "fat32", "fat", "vfat", "msdos"}
            if fs_type in no_symlink_fs:
                # FAT filesystem - use /tmp locations
                if "UV_CACHE_DIR" not in os.environ:
                    os.environ["UV_CACHE_DIR"] = "/tmp/uv-cache"
                    Path("/tmp/uv-cache").mkdir(parents=True, exist_ok=True)
                if "UV_PROJECT_ENVIRONMENT" not in os.environ:
                    os.environ["UV_PROJECT_ENVIRONMENT"] = "/tmp/metainformant_venv"
                if "METAINFORMANT_VENV" not in os.environ:
                    os.environ["METAINFORMANT_VENV"] = "/tmp/metainformant_venv"
            else:
                # Standard filesystem - use repo locations
                if "UV_CACHE_DIR" not in os.environ:
                    cache_dir = repo_root / ".uv-cache"
                    cache_dir.mkdir(parents=True, exist_ok=True)
                    os.environ["UV_CACHE_DIR"] = str(cache_dir)
                if "UV_PROJECT_ENVIRONMENT" not in os.environ:
                    os.environ["UV_PROJECT_ENVIRONMENT"] = str(repo_root / ".venv")
                if "METAINFORMANT_VENV" not in os.environ:
                    os.environ["METAINFORMANT_VENV"] = str(repo_root / ".venv")
        except Exception:
            # If all detection fails, continue without setting UV vars
            # User can set them manually if needed
            pass


@pytest.fixture(scope="session", autouse=True)
def load_ncbi_config():
    """Load NCBI configuration from config file if NCBI_EMAIL not in environment.
    
    This fixture runs automatically for all tests and ensures NCBI_EMAIL
    is set from config/ncbi.yaml if not already set in environment.
    """
    if "NCBI_EMAIL" not in os.environ:
        try:
            from pathlib import Path
            from metainformant.core.config import load_mapping_from_file
            
            config_path = Path(__file__).parent.parent / "config" / "ncbi.yaml"
            if config_path.exists():
                config = load_mapping_from_file(config_path)
                email = config.get("email", "").strip()
                if email:
                    os.environ["NCBI_EMAIL"] = email
        except Exception:
            # If config loading fails, continue without setting NCBI_EMAIL
            # Tests will skip if they require it
            pass


@pytest.fixture(scope="session", autouse=True)
def ensure_amalgkit_available():
    """Ensure amalgkit CLI is available before running RNA tests.
    
    This fixture runs automatically for all tests and ensures amalgkit
    is available. If not available, it attempts auto-install.
    If installation fails, the test session will fail.
    """
    import os
    from pathlib import Path
    from metainformant.rna.amalgkit import ensure_cli_available
    
    # Ensure user local bin is in PATH (common for --user installs)
    user_bin = Path.home() / ".local" / "bin"
    if user_bin.exists():
        current_path = os.environ.get("PATH", "")
        if str(user_bin) not in current_path:
            os.environ["PATH"] = f"{user_bin}:{current_path}"
    
    ok, msg, install_record = ensure_cli_available(auto_install=True)
    
    if not ok:
        error_msg = f"amalgkit CLI not available: {msg}"
        if install_record and install_record.get("attempted"):
            error_msg += f"\nInstall attempt failed with return code: {install_record.get('return_code')}"
            error_msg += f"\nInstall stderr: {install_record.get('stderr', '')[:500]}"
        pytest.fail(error_msg)
    
    return ok, msg


class TestFileSystem:
    """Helper filesystem for testing file operations using real I/O operations.
    
    This class uses real file operations via tmp_path, following the NO_MOCKING_POLICY.
    It provides convenience methods for creating test files in isolated temporary directories.
    """

    def __init__(self, tmp_path: Path):
        self.tmp_path = tmp_path
        self.files = {}

    def create_file(self, path: str, content: str) -> Path:
        """Create a file with specified content."""
        file_path = self.tmp_path / path
        file_path.parent.mkdir(parents=True, exist_ok=True)
        file_path.write_text(content)
        self.files[path] = content
        return file_path

    def create_binary_file(self, path: str, content: bytes) -> Path:
        """Create a binary file with specified content."""
        file_path = self.tmp_path / path
        file_path.parent.mkdir(parents=True, exist_ok=True)
        file_path.write_bytes(content)
        return file_path


@pytest.fixture(scope="function")
def test_filesystem(isolated_tmp_dir) -> TestFileSystem:
    """Provide a test filesystem helper for testing.
    
    Uses real file operations via tmp_path, following NO_MOCKING_POLICY.
    """
    return TestFileSystem(isolated_tmp_dir)


@pytest.fixture(scope="function", autouse=True)
def cleanup_matplotlib_figures():
    """Automatically close all matplotlib figures after each test.
    
    This prevents RuntimeWarnings about too many open figures and ensures
    proper cleanup of visualization resources.
    """
    yield
    try:
        import matplotlib.pyplot as plt
        plt.close("all")
    except ImportError:
        pass  # matplotlib not available, nothing to clean up


# Pytest hooks for enhanced test reporting
def pytest_configure(config):
    """Configure pytest with custom settings."""
    config.addinivalue_line("markers", "slow: mark test as slow running")
    config.addinivalue_line("markers", "network: mark test as requiring network access")
    config.addinivalue_line("markers", "external_tool: mark test as requiring external tools")


def pytest_collection_modifyitems(config, items):
    """Modify test collection to add automatic markers."""
    # Add 'network' marker to tests that contain network-related keywords
    network_keywords = ["http", "api", "fetch", "download", "request"]

    # Add 'external_tool' marker to tests that require external tools
    external_tool_keywords = ["muscle", "amalgkit", "blast", "ncbi"]

    for item in items:
        # Check if test name contains network-related keywords
        if any(keyword in item.name.lower() for keyword in network_keywords):
            item.add_marker(pytest.mark.network)

        # Check if test name contains external tool keywords
        if any(keyword in item.name.lower() for keyword in external_tool_keywords):
            item.add_marker(pytest.mark.external_tool)

        # Mark long-running tests
        if "integration" in item.name.lower() or "slow" in item.name.lower():
            item.add_marker(pytest.mark.slow)


# Skip conditions for various test scenarios
def pytest_runtest_setup(item):
    """Setup function to handle test skipping based on environment."""
    # Skip network tests if explicitly requested
    if item.get_closest_marker("network") and item.config.getoption("--no-network", False):
        pytest.skip("Network tests skipped (--no-network flag)")

    # Skip external tool tests if tools not available
    if item.get_closest_marker("external_tool"):
        # Check if required tools are available
        # This could be expanded to check specific tools
        pass


def pytest_addoption(parser):
    """Add custom command-line options for pytest."""
    parser.addoption("--no-network", action="store_true", default=False, help="Skip tests that require network access")
    parser.addoption(
        "--no-external-tools", action="store_true", default=False, help="Skip tests that require external tools"
    )
    parser.addoption(
        "--run-slow", action="store_true", default=False, help="Run slow tests (by default they are skipped)"
    )


# Fixtures for common test data files
@pytest.fixture(scope="session")
def sample_fasta_content() -> str:
    """Sample FASTA content for testing."""
    return """>seq1
ATCGATCGATCG
>seq2
GGCCAAGGCCAA
>seq3
TTAACCGGTTAA
"""


@pytest.fixture(scope="session")
def sample_fastq_content() -> str:
    """Sample FASTQ content for testing."""
    return """@read1
ATCGATCGATCG
+
IIIIIIIIIIII
@read2
GGCCAAGGCCAA
+
JJJJJJJJJJJJ
@read3
TTAACCGGTTAA
+
KKKKKKKKKKKK
"""


@pytest.fixture(scope="session")
def sample_vcf_content() -> str:
    """Sample VCF content for testing."""
    return """##fileformat=VCFv4.3
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
chr1	100	.	A	G	60	PASS	DP=30	GT	0/1
chr1	200	.	C	T	80	PASS	DP=40	GT	1/1
"""


@pytest.fixture(scope="session")
def sample_pdb_content() -> str:
    """Sample PDB content for testing."""
    return """HEADER    TEST PROTEIN                            01-JAN-24   TEST
ATOM      1  N   ALA A   1      20.154  16.967  18.250  1.00 10.00           N
ATOM      2  CA  ALA A   1      19.030  17.850  18.740  1.00 10.00           C
ATOM      3  C   ALA A   1      18.540  18.749  17.650  1.00 10.00           C
ATOM      4  O   ALA A   1      18.980  18.889  16.520  1.00 10.00           O
END
"""


# Performance testing utilities
class PerformanceTracker:
    """Utility for tracking test performance metrics."""

    def __init__(self):
        self.metrics = {}

    def track_time(self, test_name: str, duration: float):
        """Track execution time for a test."""
        if test_name not in self.metrics:
            self.metrics[test_name] = []
        self.metrics[test_name].append(duration)

    def get_average_time(self, test_name: str) -> float:
        """Get average execution time for a test."""
        if test_name not in self.metrics:
            return 0.0
        return sum(self.metrics[test_name]) / len(self.metrics[test_name])


@pytest.fixture(scope="session")
def performance_tracker() -> PerformanceTracker:
    """Provide performance tracking for tests."""
    return PerformanceTracker()


# Test result reporting
def pytest_terminal_summary(terminalreporter, exitstatus, config):
    """Provide custom terminal summary after test run."""
    if hasattr(terminalreporter, "stats"):
        # Report coverage information if available
        if hasattr(config.option, "cov") and config.option.cov:
            terminalreporter.write_line("Coverage report generated in output/coverage_html/", yellow=True)

        # Report slow tests
        slow_tests = []
        for item in terminalreporter.stats.get("passed", []):
            if hasattr(item, "duration") and item.duration > 1.0:  # > 1 second
                slow_tests.append((item.nodeid, item.duration))

        if slow_tests:
            terminalreporter.write_line(f"\nSlow tests (>{1.0}s):", yellow=True)
            for test_name, duration in sorted(slow_tests, key=lambda x: x[1], reverse=True)[:5]:
                terminalreporter.write_line(f"  {duration:.2f}s: {test_name}")
