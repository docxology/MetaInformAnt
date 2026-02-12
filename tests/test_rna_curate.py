import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path
from metainformant.rna.amalgkit import amalgkit, AmalgkitParams

class TestAmalgkitCurate:
    """Tests for amalgkit curate command wrapper."""

    @pytest.fixture
    def mock_subprocess(self):
        with patch("subprocess.run") as mock:
            yield mock

    @pytest.fixture
    def temp_dir(self, tmp_path):
        return tmp_path

    def test_curate_basic(self, mock_subprocess, temp_dir):
        """Test basic curate command generation."""
        mock_subprocess.return_value = MagicMock(returncode=0, stdout="", stderr="")
        
        params = AmalgkitParams(
            work_dir=temp_dir,
            dist="correlation",
            eval="r2"
        )
        
        amalgkit.curate(params)
        
        # Verify command
        args = mock_subprocess.call_args[0][0]
        assert args[0] == "amalgkit"
        assert args[1] == "curate"
        assert "--out_dir" in args
        assert str(temp_dir) in args
        assert "--dist" in args
        assert "correlation" in args
        assert "--eval" in args
        assert "r2" in args

    def test_curate_with_dict_params(self, mock_subprocess, temp_dir):
        """Test curate with dictionary parameters."""
        mock_subprocess.return_value = MagicMock(returncode=0, stdout="", stderr="")
        
        params = {
            "work_dir": temp_dir,
            "threshold": 0.5,
            "skip_curation": True
        }
        
        amalgkit.curate(params)
        
        args = mock_subprocess.call_args[0][0]
        assert "--threshold" in args
        assert "0.5" in args
        assert "--skip_curation" in args
        assert "yes" in args

    def test_curate_failure(self, mock_subprocess, temp_dir):
        """Test handling of curate failure."""
        mock_subprocess.return_value = MagicMock(returncode=1, stdout="", stderr="Error in R script")
        
        params = AmalgkitParams(work_dir=temp_dir)
        
        # Should not raise by default, but log error
        result = amalgkit.curate(params)
        assert result.returncode == 1
        assert "Error in R script" in result.stderr

    def test_curate_input_validation(self, mock_subprocess, temp_dir):
        """Test that curate accepts specific arguments."""
        mock_subprocess.return_value = MagicMock(returncode=0)
        
        # 'threads' is typically not used in curate (R script usually single-threaded or internal parallel)
        # But our wrapper might pass it if not filtered.
        # amalgkit.py filters 'threads' for 'curate' in build_cli_args.
        
        params = AmalgkitParams(work_dir=temp_dir, threads=4)
        amalgkit.curate(params)
        
        args = mock_subprocess.call_args[0][0]
        # Should NOT contain --threads
        assert "--threads" not in args
