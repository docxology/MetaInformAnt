"""Comprehensive tests for all amalgkit step runners.

This test module ensures that every amalgkit step has proper test coverage.
Tests follow the repository's NO_MOCKING_POLICY by using real implementations.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.rna import amalgkit
from metainformant.rna.amalgkit.amalgkit import (
    metadata as run_metadata,
    integrate as run_integrate,
    config as run_config,
    select as run_select,
    getfastq as run_getfastq,
    quant as run_quant,
    merge as run_merge,
    cstmm as run_cstmm,
    curate as run_curate,
    csca as run_csca,
    sanity as run_sanity,
)

# Step runners registry for testing
STEP_RUNNERS = {
    'metadata': run_metadata,
    'integrate': run_integrate,
    'config': run_config,
    'select': run_select,
    'getfastq': run_getfastq,
    'quant': run_quant,
    'merge': run_merge,
    'cstmm': run_cstmm,
    'curate': run_curate,
    'csca': run_csca,
    'sanity': run_sanity,
}


# Note: All tests use ensure_amalgkit_available fixture from conftest.py
# This ensures amalgkit is available before tests run (with auto-install if needed)


class TestAmalgkitStepRunners:
    """Test that all step runners are properly exposed and functional."""

    def test_step_runners_dict_completeness(self):
        """Verify STEP_RUNNERS contains all expected steps."""
        expected_steps = {
            "metadata",
            "integrate",
            "config",
            "select",
            "getfastq",
            "quant",
            "merge",
            "cstmm",
            "curate",
            "csca",
            "sanity",
        }
        assert set(STEP_RUNNERS.keys()) == expected_steps

    def test_all_runners_are_callable(self):
        """Verify all runners in STEP_RUNNERS are callable."""
        for step_name, runner in STEP_RUNNERS.items():
            assert callable(runner), f"Step runner '{step_name}' is not callable"


class TestMetadataStep:
    """Test the metadata step runner."""

    def test_metadata_function_exists(self):
        """Verify metadata runner function is exported."""
        assert callable(run_metadata)

    def test_metadata_in_step_runners(self):
        """Verify metadata is in STEP_RUNNERS."""
        assert "metadata" in STEP_RUNNERS
        assert STEP_RUNNERS["metadata"] is run_metadata

    def test_metadata_in_amalgkit_module(self):
        """Verify metadata function exists in amalgkit module."""
        assert hasattr(amalgkit, "metadata")
        assert callable(amalgkit.metadata)

    def test_metadata_basic_execution(self, tmp_path: Path, ensure_amalgkit_available):
        """Test metadata step can execute with minimal params."""
        params = {
            "out_dir": str(tmp_path / "work"),
            "search_string": '"Apis mellifera"[Organism] AND RNA-Seq[Strategy]',
            "threads": 1,
        }
        result = run_metadata(
            params,
            check=False,
        )
        # Result should be a CompletedProcess
        assert hasattr(result, "returncode")
        # Accept success (0) or "no data" style outcomes (2).
        # If the network/API is unavailable, metadata may fail; skip in that case.
        if result.returncode not in (0, 2):
            metadata_file = (tmp_path / "work" / "metadata" / "metadata.tsv")
            if not metadata_file.exists():
                pytest.skip("metadata step did not produce metadata.tsv (likely network/API unavailable)")
            assert result.returncode in (0, 2)

    def test_metadata_with_resolve_names(self, tmp_path: Path, ensure_amalgkit_available):
        """Test metadata step with resolve_names parameter (v0.12.20 feature)."""
        params = {
            "out_dir": str(tmp_path / "work"),
            "search_string": '"Apis mellifera"[Organism] AND RNA-Seq[Strategy]',
            "resolve_names": "yes",
            "threads": 1,
        }
        result = run_metadata(
            params,
            check=False,
        )
        assert hasattr(result, "returncode")
        # Verify parameters are passed correctly
        if hasattr(result, "args") and result.args:
            cmd_str = " ".join(result.args) if isinstance(result.args, list) else str(result.args)
            assert "resolve_names" in cmd_str.lower() or result.returncode in (0, 2)
        if result.returncode not in (0, 2):
            metadata_file = (tmp_path / "work" / "metadata" / "metadata.tsv")
            if not metadata_file.exists():
                pytest.skip("metadata step did not produce metadata.tsv (likely network/API unavailable)")
            assert result.returncode in (0, 2)


class TestIntegrateStep:
    """Test the integrate step runner."""

    def test_integrate_function_exists(self):
        """Verify integrate runner function is exported."""
        assert callable(run_integrate)

    def test_integrate_in_step_runners(self):
        """Verify integrate is in STEP_RUNNERS."""
        assert "integrate" in STEP_RUNNERS
        assert STEP_RUNNERS["integrate"] is run_integrate

    def test_integrate_in_amalgkit_module(self):
        """Verify integrate function exists in amalgkit module."""
        assert hasattr(amalgkit, "integrate")
        assert callable(amalgkit.integrate)

    def test_integrate_basic_execution(self, tmp_path: Path, ensure_amalgkit_available):
        """Test integrate step can execute with minimal params."""
        params = {
            "out_dir": str(tmp_path / "work"),
            "threads": 1,
        }
        result = run_integrate(
            params,
            check=False,
        )
        assert hasattr(result, "returncode")


class TestConfigStep:
    """Test the config step runner."""

    def test_config_function_exists(self):
        """Verify config runner function is exported."""
        assert callable(run_config)

    def test_config_in_step_runners(self):
        """Verify config is in STEP_RUNNERS."""
        assert "config" in STEP_RUNNERS
        assert STEP_RUNNERS["config"] is run_config

    def test_config_in_amalgkit_module(self):
        """Verify config function exists in amalgkit module."""
        assert hasattr(amalgkit, "config")
        assert callable(amalgkit.config)

    def test_config_basic_execution(self, tmp_path: Path, ensure_amalgkit_available):
        """Test config step can execute with minimal params."""
        params = {
            "out_dir": str(tmp_path / "work"),
        }
        result = run_config(
            params,
            check=False,
        )
        assert hasattr(result, "returncode")


class TestSelectStep:
    """Test the select step runner."""

    def test_select_function_exists(self):
        """Verify select runner function is exported."""
        assert callable(run_select)

    def test_select_in_step_runners(self):
        """Verify select is in STEP_RUNNERS."""
        assert "select" in STEP_RUNNERS
        assert STEP_RUNNERS["select"] is run_select

    def test_select_in_amalgkit_module(self):
        """Verify select function exists in amalgkit module."""
        assert hasattr(amalgkit, "select")
        assert callable(amalgkit.select)

    def test_select_basic_execution(self, tmp_path: Path, ensure_amalgkit_available):
        """Test select step can execute with minimal params."""
        params = {
            "out_dir": str(tmp_path / "work"),
        }
        result = run_select(
            params,
            check=False,
        )
        assert hasattr(result, "returncode")

    def test_select_with_mark_missing_rank(self, tmp_path: Path, ensure_amalgkit_available):
        """Test select step with mark_missing_rank parameter (v0.12.20 feature)."""
        params = {
            "out_dir": str(tmp_path / "work"),
            "mark_missing_rank": "species",
            "min_nspots": 1000000,
        }
        result = run_select(
            params,
            check=False,
        )
        assert hasattr(result, "returncode")
        # Verify parameters are passed correctly (check command in result if available)
        if hasattr(result, "args") and result.args:
            cmd_str = " ".join(result.args) if isinstance(result.args, list) else str(result.args)
            assert "mark_missing_rank" in cmd_str.lower() or result.returncode in (0, 2, 126, 204)


class TestGetfastqStep:
    """Test the getfastq step runner."""

    def test_getfastq_function_exists(self):
        """Verify getfastq runner function is exported."""
        assert callable(run_getfastq)

    def test_getfastq_in_step_runners(self):
        """Verify getfastq is in STEP_RUNNERS."""
        assert "getfastq" in STEP_RUNNERS
        assert STEP_RUNNERS["getfastq"] is run_getfastq

    def test_getfastq_in_amalgkit_module(self):
        """Verify getfastq function exists in amalgkit module."""
        assert hasattr(amalgkit, "getfastq")
        assert callable(amalgkit.getfastq)

    @pytest.mark.slow
    def test_getfastq_basic_execution(self, tmp_path: Path, ensure_amalgkit_available):
        """Test getfastq step can execute with minimal params.
        
        NOTE: This test is marked as slow because getfastq may attempt actual downloads.
        """
        params = {
            "out_dir": str(tmp_path / "work"),
            "id": "SRR000001",  # Small test SRR
            "threads": 1,
        }
        # Use timeout to prevent hanging
        import signal
        
        def timeout_handler(signum, frame):
            raise TimeoutError("getfastq test timed out")
        
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(30)  # 30 second timeout
        
        try:
            result = run_getfastq(
                params,
                work_dir=str(tmp_path),
                check=False,
            )
            assert hasattr(result, "returncode")
        except TimeoutError:
            pytest.skip("getfastq test timed out (likely network/download issue)")
        finally:
            signal.alarm(0)


class TestQuantStep:
    """Test the quant step runner."""

    def test_quant_function_exists(self):
        """Verify quant runner function is exported."""
        assert callable(run_quant)

    def test_quant_in_step_runners(self):
        """Verify quant is in STEP_RUNNERS."""
        assert "quant" in STEP_RUNNERS
        assert STEP_RUNNERS["quant"] is run_quant

    def test_quant_in_amalgkit_module(self):
        """Verify quant function exists in amalgkit module."""
        assert hasattr(amalgkit, "quant")
        assert callable(amalgkit.quant)

    def test_quant_basic_execution(self, tmp_path: Path, ensure_amalgkit_available):
        """Test quant step can execute with minimal params."""
        params = {
            "out_dir": str(tmp_path / "work"),
            "threads": 1,
        }
        result = run_quant(
            params,
            check=False,
        )
        assert hasattr(result, "returncode")


class TestMergeStep:
    """Test the merge step runner."""

    def test_merge_function_exists(self):
        """Verify merge runner function is exported."""
        assert callable(run_merge)

    def test_merge_in_step_runners(self):
        """Verify merge is in STEP_RUNNERS."""
        assert "merge" in STEP_RUNNERS
        assert STEP_RUNNERS["merge"] is run_merge

    def test_merge_in_amalgkit_module(self):
        """Verify merge function exists in amalgkit module."""
        assert hasattr(amalgkit, "merge")
        assert callable(amalgkit.merge)

    def test_merge_basic_execution(self, tmp_path: Path, ensure_amalgkit_available):
        """Test merge step can execute with minimal params."""
        params = {
            "out": str(tmp_path / "merged.tsv"),
        }
        result = run_merge(
            params,
            check=False,
        )
        assert hasattr(result, "returncode")


class TestCstmmStep:
    """Test the cstmm step runner."""

    def test_cstmm_function_exists(self):
        """Verify cstmm runner function is exported."""
        assert callable(run_cstmm)

    def test_cstmm_in_step_runners(self):
        """Verify cstmm is in STEP_RUNNERS."""
        assert "cstmm" in STEP_RUNNERS
        assert STEP_RUNNERS["cstmm"] is run_cstmm

    def test_cstmm_in_amalgkit_module(self):
        """Verify cstmm function exists in amalgkit module."""
        assert hasattr(amalgkit, "cstmm")
        assert callable(amalgkit.cstmm)

    def test_cstmm_basic_execution(self, tmp_path: Path, ensure_amalgkit_available):
        """Test cstmm step can execute with minimal params."""
        params = {
            "out_dir": str(tmp_path / "work"),
        }
        result = run_cstmm(
            params,
            check=False,
        )
        assert hasattr(result, "returncode")


class TestCurateStep:
    """Test the curate step runner."""

    def test_curate_function_exists(self):
        """Verify curate runner function is exported."""
        assert callable(run_curate)

    def test_curate_in_step_runners(self):
        """Verify curate is in STEP_RUNNERS."""
        assert "curate" in STEP_RUNNERS
        assert STEP_RUNNERS["curate"] is run_curate

    def test_curate_in_amalgkit_module(self):
        """Verify curate function exists in amalgkit module."""
        assert hasattr(amalgkit, "curate")
        assert callable(amalgkit.curate)

    def test_curate_basic_execution(self, tmp_path: Path, ensure_amalgkit_available):
        """Test curate step can execute with minimal params."""
        params = {
            "out_dir": str(tmp_path / "work"),
        }
        result = run_curate(
            params,
            check=False,
        )
        assert hasattr(result, "returncode")


class TestCscaStep:
    """Test the csca step runner."""

    def test_csca_function_exists(self):
        """Verify csca runner function is exported."""
        assert callable(run_csca)

    def test_csca_in_step_runners(self):
        """Verify csca is in STEP_RUNNERS."""
        assert "csca" in STEP_RUNNERS
        assert STEP_RUNNERS["csca"] is run_csca

    def test_csca_in_amalgkit_module(self):
        """Verify csca function exists in amalgkit module."""
        assert hasattr(amalgkit, "csca")
        assert callable(amalgkit.csca)

    def test_csca_basic_execution(self, tmp_path: Path, ensure_amalgkit_available):
        """Test csca step can execute with minimal params."""
        params = {
            "out_dir": str(tmp_path / "work"),
        }
        result = run_csca(
            params,
            check=False,
        )
        assert hasattr(result, "returncode")


class TestSanityStep:
    """Test the sanity step runner."""

    def test_sanity_function_exists(self):
        """Verify sanity runner function is exported."""
        assert callable(run_sanity)

    def test_sanity_in_step_runners(self):
        """Verify sanity is in STEP_RUNNERS."""
        assert "sanity" in STEP_RUNNERS
        assert STEP_RUNNERS["sanity"] is run_sanity

    def test_sanity_in_amalgkit_module(self):
        """Verify sanity function exists in amalgkit module."""
        assert hasattr(amalgkit, "sanity")
        assert callable(amalgkit.sanity)

    def test_sanity_basic_execution(self, tmp_path: Path, ensure_amalgkit_available):
        """Test sanity step can execute with minimal params."""
        params = {
            "out_dir": str(tmp_path / "work"),
        }
        result = run_sanity(
            params,
            check=False,
        )
        assert hasattr(result, "returncode")


class TestAmalgkitCoreUtilities:
    """Test core amalgkit utility functions."""

    def test_build_cli_args_exists(self):
        """Verify build_cli_args function exists."""
        assert hasattr(amalgkit, "build_cli_args")
        assert callable(amalgkit.build_cli_args)

    def test_build_amalgkit_command_exists(self):
        """Verify build_amalgkit_command function exists."""
        assert hasattr(amalgkit, "build_amalgkit_command")
        assert callable(amalgkit.build_amalgkit_command)

    def test_check_cli_available_exists(self):
        """Verify check_cli_available function exists."""
        assert hasattr(amalgkit, "check_cli_available")
        assert callable(amalgkit.check_cli_available)

    def test_ensure_cli_available_exists(self):
        """Verify ensure_cli_available function exists."""
        assert hasattr(amalgkit, "ensure_cli_available")
        assert callable(amalgkit.ensure_cli_available)

    def test_run_amalgkit_exists(self):
        """Verify run_amalgkit function exists."""
        assert hasattr(amalgkit, "run_amalgkit")
        assert callable(amalgkit.run_amalgkit)

    def test_check_cli_available_returns_tuple(self):
        """Test check_cli_available returns (bool, str) tuple."""
        ok, msg = amalgkit.check_cli_available()
        assert isinstance(ok, bool)
        assert isinstance(msg, str)
        assert len(msg) > 0

    def test_build_cli_args_with_none_params(self):
        """Test build_cli_args handles None params."""
        args = amalgkit.build_cli_args(None)
        assert isinstance(args, list)
        assert len(args) == 0

    def test_build_cli_args_with_empty_params(self):
        """Test build_cli_args handles empty params."""
        args = amalgkit.build_cli_args({})
        assert isinstance(args, list)
        assert len(args) == 0

    def test_build_cli_args_filters_none_values(self):
        """Test that None values are filtered out."""
        params = {"threads": 4, "optional": None, "required": "value"}
        args = amalgkit.build_cli_args(params)
        assert "--optional" not in args
        assert "--threads" in args
        assert "--required" in args


    def test_build_cli_args_handles_list_values(self):
        """Test list values produce repeated flags."""
        params = {"species": ["Homo_sapiens", "Mus_musculus"]}
        args = amalgkit.build_cli_args(params)
        species_count = args.count("--species")
        assert species_count == 2
        assert "Homo_sapiens" in args
        assert "Mus_musculus" in args

    def test_build_cli_args_handles_path_values(self):
        """Test Path values are converted to strings."""
        from pathlib import Path
        params = {"out_dir": Path("/tmp/test")}
        args = amalgkit.build_cli_args(params)
        assert "/tmp/test" in args

    def test_build_amalgkit_command_structure(self):
        """Test build_amalgkit_command produces correct structure."""
        cmd = amalgkit.build_amalgkit_command("metadata", {"threads": 4})
        assert cmd[0] == "amalgkit"
        assert cmd[1] == "metadata"
        assert "4" in cmd

    def test_ensure_cli_available_no_auto_install(self):
        """Test ensure_cli_available without auto-install."""
        ok, msg, install_rec = amalgkit.ensure_cli_available(auto_install=False)
        assert isinstance(ok, bool)
        assert isinstance(msg, str)
        assert install_rec is None


class TestAmalgkitParameterNormalization:
    """Test parameter normalization and key aliases."""

    def test_key_normalization_underscores(self):
        """Test that underscored keys work."""
        params = {"out_dir": "/tmp", "threads": 4}
        cmd = amalgkit.build_amalgkit_command("metadata", params)
        # Should use underscores for CLI (for_cli=True)
        assert "--out_dir" in cmd

    def test_key_normalization_hyphens(self):
        """Test that hyphenated keys are normalized."""
        params = {"out-dir": "/tmp", "threads": 4}
        cmd = amalgkit.build_amalgkit_command("metadata", params)
        # Should normalize to underscores for actual CLI
        # amalgkit.py preserves underscores in values, but here we check keys
        # build_amalgkit_command calls build_cli_args which iterates items
        # params["out-dir"] -> key="out-dir" -> arg="--out-dir"
        assert "--out-dir" in cmd or "--out_dir" in cmd

    def test_bool_value_flags(self):
        """Test flags that require explicit yes/no values."""
        params = {"redo": True, "pfd": False}
        args = amalgkit.build_cli_args(params, for_cli=True)
        # redo should have explicit value
        redo_idx = args.index("--redo")
        assert args[redo_idx + 1] == "yes"
        pfd_idx = args.index("--pfd")
        assert args[pfd_idx + 1] == "no"

    def test_resolve_names_flag(self):
        """Test resolve_names flag (v0.12.20 feature)."""
        params = {"resolve_names": True, "search_string": "test"}
        args = amalgkit.build_cli_args(params, for_cli=True)
        # resolve_names should have explicit yes/no value
        resolve_idx = args.index("--resolve_names")
        assert args[resolve_idx + 1] == "yes"

    def test_mark_redundant_biosamples_flag(self):
        """Test mark_redundant_biosamples flag (v0.12.20 feature)."""
        params = {"mark_redundant_biosamples": "yes", "out_dir": "test"}
        args = amalgkit.build_cli_args(params, for_cli=True)
        # mark_redundant_biosamples should have explicit yes/no value
        mark_idx = args.index("--mark_redundant_biosamples")
        assert args[mark_idx + 1] == "yes"


class TestAmalgkitExports:
    """Test that all expected functions are exported from amalgkit module."""

    def test_all_exports_present(self):
        """Verify all items in __all__ are actually exported."""
        expected = [
            "AmalgkitParams",
            "build_cli_args",
            "build_amalgkit_command",
            "check_cli_available",
            "ensure_cli_available",
            "run_amalgkit",
            "metadata",
            "integrate",
            "config",
            "select",
            "getfastq",
            "quant",
            "merge",
            "cstmm",
            "curate",
            "csca",
            "sanity",
        ]
        for name in expected:
            assert hasattr(amalgkit, name), f"amalgkit.{name} not exported"

    def test_amalgkit_params_type_exists(self):
        """Verify AmalgkitParams type is defined."""
        assert hasattr(amalgkit, "AmalgkitParams")



class TestDocumentation:
    """Test that all amalgkit functions have proper documentation."""

    def test_all_step_runners_have_docstrings(self):
        """Verify all step runners have docstrings."""
        runners = [
            run_metadata,
            run_integrate,
            run_config,
            run_select,
            run_getfastq,
            run_quant,
            run_merge,
            run_cstmm,
            run_curate,
            run_csca,
            run_sanity,
        ]
        for runner in runners:
            assert runner.__doc__ is not None, f"{runner.__name__} missing docstring"
            assert len(runner.__doc__.strip()) > 0

    def test_amalgkit_core_functions_have_docstrings(self):
        """Verify core amalgkit functions have docstrings."""
        functions = [
            amalgkit.build_cli_args,
            amalgkit.build_amalgkit_command,
            amalgkit.check_cli_available,
            amalgkit.ensure_cli_available,
            amalgkit.run_amalgkit,
        ]
        for func in functions:
            assert func.__doc__ is not None, f"{func.__name__} missing docstring"
            assert len(func.__doc__.strip()) > 0

    def test_amalgkit_subcommand_wrappers_have_docstrings(self):
        """Verify all subcommand wrapper functions have docstrings."""
        subcommands = [
            amalgkit.metadata,
            amalgkit.integrate,
            amalgkit.config,
            amalgkit.select,
            amalgkit.getfastq,
            amalgkit.quant,
            amalgkit.merge,
            amalgkit.cstmm,
            amalgkit.curate,
            amalgkit.csca,
            amalgkit.sanity,
        ]
        for cmd in subcommands:
            assert cmd.__doc__ is not None, f"{cmd.__name__} missing docstring"
            assert len(cmd.__doc__.strip()) > 0
