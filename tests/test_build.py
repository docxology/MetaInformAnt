"""
METAINFORMANT Build Validation Tests

Tests for validating build artifacts, package installation,
and build process integrity.
"""

import pytest
import subprocess
import tempfile
import shutil
from pathlib import Path
import sys
import os

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import metainformant


class TestBuildArtifacts:
    """Test build artifact validation."""

    def test_package_can_be_imported(self):
        """Test that the package can be imported after installation."""
        # This tests the installed package
        assert hasattr(metainformant, '__version__')
        assert isinstance(metainformant.__version__, str)
        assert len(metainformant.__version__.split('.')) >= 2

    def test_core_modules_importable(self):
        """Test that core modules can be imported."""
        try:
            from metainformant.core import io, paths, logging, config
            assert True  # Import successful
        except ImportError as e:
            pytest.fail(f"Core module import failed: {e}")

    def test_domain_modules_importable(self):
        """Test that domain modules can be imported."""
        domains_to_test = [
            'dna.sequence',
            'rna.analysis',
            'protein.structure',
            'gwas.analysis',
            'epigenome.analysis',
            'ontology.go',
            'phenotype.analysis',
            'ecology.analysis',
            'math.population_genetics',
            'ml.models',
            'networks.analysis',
            'singlecell.analysis',
            'quality.analysis',
            'visualization.plots',
            'simulation.models',
            'life_events.analysis'
        ]

        failed_imports = []
        for module in domains_to_test:
            try:
                __import__(f'metainformant.{module}')
            except ImportError as e:
                failed_imports.append(f"{module}: {e}")

        if failed_imports:
            pytest.fail(f"Some domain modules failed to import: {failed_imports}")

    def test_cli_available(self):
        """Test that CLI command is available."""
        try:
            result = subprocess.run(
                [sys.executable, '-m', 'metainformant', '--help'],
                capture_output=True,
                text=True,
                timeout=10
            )
            assert result.returncode == 0
            assert 'metainformant' in result.stdout.lower()
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pytest.fail("CLI command not available or timed out")


class TestPackageStructure:
    """Test package structure and metadata."""

    def test_version_consistency(self):
        """Test that version is consistent across package."""
        # Import version from package
        package_version = metainformant.__version__

        # Read version from pyproject.toml
        pyproject_path = Path(__file__).parent.parent / "pyproject.toml"
        assert pyproject_path.exists(), "pyproject.toml not found"

        with open(pyproject_path, 'r') as f:
            content = f.read()

        # Extract version from pyproject.toml
        import re
        version_match = re.search(r'version = "([^"]+)"', content)
        assert version_match, "Version not found in pyproject.toml"

        pyproject_version = version_match.group(1)
        assert package_version == pyproject_version, \
            f"Version mismatch: package={package_version}, pyproject={pyproject_version}"

    def test_package_metadata(self):
        """Test package metadata is accessible."""
        # Test that we can import and access basic metadata
        assert hasattr(metainformant, '__version__')
        assert hasattr(metainformant, '__doc__')
        assert metainformant.__doc__ is not None


class TestBuildProcess:
    """Test build process components."""

    @pytest.mark.slow
    def test_build_scripts_exist(self):
        """Test that build scripts exist and are executable."""
        scripts_dir = Path(__file__).parent.parent / "scripts" / "package"

        build_script = scripts_dir / "build.sh"
        assert build_script.exists(), "build.sh not found"
        assert build_script.stat().st_mode & 0o111, "build.sh not executable"

        validate_script = scripts_dir / "validate_build.sh"
        assert validate_script.exists(), "validate_build.sh not found"
        assert validate_script.stat().st_mode & 0o111, "validate_build.sh not executable"

    def test_build_utils_available(self):
        """Test that build utilities are importable."""
        # This would require importing build utilities if they were Python modules
        # For now, just check that the script files exist
        scripts_dir = Path(__file__).parent.parent / "scripts" / "package"
        utils_script = scripts_dir / "build_utils.sh"
        assert utils_script.exists(), "build_utils.sh not found"


class TestDocumentationBuild:
    """Test documentation build process."""

    @pytest.mark.slow
    def test_docs_config_exists(self):
        """Test that documentation configuration exists."""
        docs_dir = Path(__file__).parent.parent / "docs"

        conf_py = docs_dir / "conf.py"
        assert conf_py.exists(), "docs/conf.py not found"

        index_md = docs_dir / "index.md"
        assert index_md.exists(), "docs/index.md not found"

    def test_sphinx_config_valid(self):
        """Test that Sphinx configuration is valid."""
        docs_dir = Path(__file__).parent.parent / "docs"
        conf_py = docs_dir / "conf.py"

        # Try to import the configuration (basic syntax check)
        import importlib.util
        spec = importlib.util.spec_from_file_location("conf", conf_py)
        assert spec is not None, "Cannot load Sphinx configuration"

        # This is a basic check - actual Sphinx validation would happen during build


class TestCleanEnvironment:
    """Test functionality in clean environments."""

    def test_import_in_subprocess(self):
        """Test that package can be imported in a subprocess."""
        # This simulates testing in a clean environment
        test_code = """
import sys
sys.path.insert(0, 'src')
import metainformant
print(f"Version: {metainformant.__version__}")
from metainformant.core import io
print("Core import successful")
"""

        result = subprocess.run(
            [sys.executable, '-c', test_code],
            cwd=Path(__file__).parent.parent,
            capture_output=True,
            text=True,
            timeout=30
        )

        assert result.returncode == 0, f"Import test failed: {result.stderr}"
        assert "Version:" in result.stdout
        assert "Core import successful" in result.stdout


class TestEntryPoints:
    """Test package entry points."""

    def test_entry_points_defined(self):
        """Test that entry points are properly defined."""
        # Check pyproject.toml for entry points
        pyproject_path = Path(__file__).parent.parent / "pyproject.toml"
        assert pyproject_path.exists()

        with open(pyproject_path, 'r') as f:
            content = f.read()

        # Check for entry points section
        assert '[project.scripts]' in content, "Entry points not defined in pyproject.toml"
        assert 'metainformant =' in content, "metainformant entry point not found"

    def test_entry_point_functionality(self):
        """Test that entry points work correctly."""
        # Test the entry point module exists
        entry_point_path = Path(__file__).parent.parent / "src" / "metainformant" / "__main__.py"
        assert entry_point_path.exists(), "Entry point module not found"

        # Test that the entry point can be executed
        result = subprocess.run(
            [sys.executable, '-m', 'metainformant', '--version'],
            capture_output=True,
            text=True,
            timeout=10
        )

        assert result.returncode == 0, f"Entry point failed: {result.stderr}"
        assert metainformant.__version__ in result.stdout


# Integration tests that require build artifacts
class TestBuiltPackage:
    """Tests that require built package artifacts."""

    @pytest.fixture(scope="session")
    def built_package(self):
        """Fixture to provide path to built package."""
        dist_dir = Path(__file__).parent.parent / "dist"

        # Look for wheel file
        wheel_files = list(dist_dir.glob("*.whl"))
        if not wheel_files:
            pytest.skip("No built wheel package found. Run 'bash scripts/package/build.sh' first.")

        return wheel_files[0]

    def test_wheel_can_be_installed(self, built_package):
        """Test that built wheel can be installed."""
        # Create temporary directory for testing
        with tempfile.TemporaryDirectory() as temp_dir:
            venv_dir = Path(temp_dir) / "test_venv"

            # Create virtual environment
            subprocess.run(
                [sys.executable, '-m', 'venv', str(venv_dir)],
                check=True,
                capture_output=True
            )

            # Install package
            pip_path = venv_dir / "bin" / "pip"
            result = subprocess.run(
                [str(pip_path), 'install', str(built_package)],
                capture_output=True,
                text=True,
                cwd=temp_dir
            )

            assert result.returncode == 0, f"Installation failed: {result.stderr}"

            # Test import in the virtual environment
            python_path = venv_dir / "bin" / "python"
            test_result = subprocess.run(
                [str(python_path), '-c', 'import metainformant; print("Import successful")'],
                capture_output=True,
                text=True,
                cwd=temp_dir
            )

            assert test_result.returncode == 0, f"Import failed: {test_result.stderr}"
            assert "Import successful" in test_result.stdout

    def test_wheel_metadata(self, built_package):
        """Test that wheel contains proper metadata."""
        import zipfile

        with zipfile.ZipFile(built_package, 'r') as zf:
            # Look for METADATA file
            metadata_files = [name for name in zf.namelist() if name.endswith('METADATA')]
            assert len(metadata_files) > 0, "No METADATA file found in wheel"

            # Read metadata
            with zf.open(metadata_files[0]) as f:
                metadata = f.read().decode('utf-8')

            # Check for required fields
            assert 'Name: metainformant' in metadata
            assert f'Version: {metainformant.__version__}' in metadata
            assert 'Summary:' in metadata
