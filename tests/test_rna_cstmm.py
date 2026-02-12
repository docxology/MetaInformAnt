
import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path
from metainformant.rna.amalgkit import amalgkit
from metainformant.rna.amalgkit.amalgkit import AmalgkitParams

class TestAmalgkitCstmm:
    @pytest.fixture
    def mock_subprocess(self):
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="success", stderr="")
            yield mock_run

    @pytest.fixture
    def temp_dir(self, tmp_path):
        return tmp_path

    def test_cstmm_basic(self, mock_subprocess, temp_dir):
        params = AmalgkitParams(
            work_dir=temp_dir, 
            threads=8, 
            species_list=["Species_A"],
            orthogroup_table="groups.tsv"
        )
        # Should NOT pass threads
        amalgkit.cstmm(params)
        
        args = mock_subprocess.call_args[0][0]
        assert args[0] == "amalgkit"
        assert args[1] == "cstmm"
        assert "--out_dir" in args
        assert str(temp_dir) in args
        assert "--threads" not in args  # cstmm does not support threads
        assert "--orthogroup_table" in args
        assert "groups.tsv" in args
        assert "--species" in args
        assert "Species_A" in args

    def test_cstmm_dict_params(self, mock_subprocess, temp_dir):
        params = {
            "work_dir": temp_dir,
            "threads": 4, # Should be ignored
            "gff": "yes",
            "dir_busco": "/path/to/busco"
        }
        amalgkit.cstmm(params)
        args = mock_subprocess.call_args[0][0]
        assert "--threads" not in args
        assert "--gff" in args
        assert "yes" in args
        assert "--dir_busco" in args
        assert "/path/to/busco" in args

    def test_cstmm_failure(self, mock_subprocess, temp_dir):
        mock_subprocess.return_value.returncode = 1
        mock_subprocess.return_value.stderr = "Error: cstmm failed"
        
        params = {"work_dir": temp_dir}
        result = amalgkit.cstmm(params)
        
        assert result.returncode == 1
        assert "Error" in result.stderr
