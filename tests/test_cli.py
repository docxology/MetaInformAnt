import subprocess
import sys


def test_module_invocation_shows_help():
    """Test that module invocation displays help information."""
    result = subprocess.run([sys.executable, "-m", "metainformant"], capture_output=True, text=True)
    assert result.returncode == 0
    assert "METAINFORMANT CLI" in result.stdout or "usage:" in result.stdout
