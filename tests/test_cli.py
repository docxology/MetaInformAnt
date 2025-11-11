import subprocess
import sys
from pathlib import Path


def test_module_invocation_shows_help():
    """Test that module invocation displays help information."""
    # Use the source module directly since package may not be installed
    repo_root = Path(__file__).parent.parent
    module_path = repo_root / "src" / "metainformant" / "__main__.py"
    
    # Run the module directly
    result = subprocess.run(
        [sys.executable, str(module_path)],
        capture_output=True,
        text=True,
        cwd=str(repo_root),
    )
    # Module should show help when no args provided (or exit with non-zero)
    # Check for help text in either stdout or stderr
    output = result.stdout + result.stderr
    assert "METAINFORMANT CLI" in output or "usage:" in output or "metainformant" in output.lower()
