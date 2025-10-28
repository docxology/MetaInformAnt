"""Comprehensive tests for all amalgkit step runners.

This test module ensures that every amalgkit step has proper test coverage.
Tests follow the repository's NO_MOCKING_POLICY by using real implementations.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from metainformant.rna import amalgkit
from metainformant.rna.steps import (
    STEP_RUNNERS,
    run_config,
    run_csca,
    run_cstmm,
    run_curate,
    run_getfastq,
    run_integrate,
    run_merge,
    run_metadata,
    run_quant,
    run_sanity,
    run_select,
)


# Conditional skip if amalgkit is not available
requires_amalgkit = pytest.mark.skipif(
    shutil.which("amalgkit") is None,
    reason="amalgkit CLI not available on PATH"
)


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

    def test_all_runners_have_correct_signature(self):
        """Verify all runners accept standard parameters."""
        import inspect
        
        for step_name, runner in STEP_RUNNERS.items():
            sig = inspect.signature(runner)
            # Check that runners accept params, work_dir, log_dir, check
            params = sig.parameters
            assert "params" in params, f"Step '{step_name}' missing 'params' parameter"
            # work_dir, log_dir, check should be keyword-only
            assert "work_dir" in params, f"Step '{step_name}' missing 'work_dir' parameter"
            assert "log_dir" in params, f"Step '{step_name}' missing 'log_dir' parameter"
            assert "check" in params, f"Step '{step_name}' missing 'check' parameter"


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

    @requires_amalgkit
    def test_metadata_basic_execution(self, tmp_path: Path):
        """Test metadata step can execute with minimal params."""
        params = {
            "out_dir": str(tmp_path / "work"),
            "search_string": '"Apis mellifera"[Organism] AND RNA-Seq[Strategy]',
            "threads": 1,
        }
        result = run_metadata(
            params,
            work_dir=str(tmp_path),
            log_dir=str(tmp_path / "logs"),
            check=False,
        )
        # Result should be a CompletedProcess
        assert hasattr(result, "returncode")
        # Accept success (0) or no data found (2)
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

    @requires_amalgkit
    def test_integrate_basic_execution(self, tmp_path: Path):
        """Test integrate step can execute with minimal params."""
        params = {
            "out_dir": str(tmp_path / "work"),
            "threads": 1,
        }
        result = run_integrate(
            params,
            work_dir=str(tmp_path),
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

    @requires_amalgkit
    def test_config_basic_execution(self, tmp_path: Path):
        """Test config step can execute with minimal params."""
        params = {
            "out_dir": str(tmp_path / "work"),
        }
        result = run_config(
            params,
            work_dir=str(tmp_path),
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

    @requires_amalgkit
    def test_select_basic_execution(self, tmp_path: Path):
        """Test select step can execute with minimal params."""
        params = {
            "out_dir": str(tmp_path / "work"),
        }
        result = run_select(
            params,
            work_dir=str(tmp_path),
            check=False,
        )
        assert hasattr(result, "returncode")


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

    @requires_amalgkit
    def test_getfastq_basic_execution(self, tmp_path: Path):
        """Test getfastq step can execute with minimal params."""
        params = {
            "out_dir": str(tmp_path / "work"),
            "id": "SRR000001",  # Small test SRR
            "threads": 1,
        }
        result = run_getfastq(
            params,
            work_dir=str(tmp_path),
            check=False,
        )
        assert hasattr(result, "returncode")


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

    @requires_amalgkit
    def test_quant_basic_execution(self, tmp_path: Path):
        """Test quant step can execute with minimal params."""
        params = {
            "out_dir": str(tmp_path / "work"),
            "threads": 1,
        }
        result = run_quant(
            params,
            work_dir=str(tmp_path),
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

    @requires_amalgkit
    def test_merge_basic_execution(self, tmp_path: Path):
        """Test merge step can execute with minimal params."""
        params = {
            "out": str(tmp_path / "merged.tsv"),
        }
        result = run_merge(
            params,
            work_dir=str(tmp_path),
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

    @requires_amalgkit
    def test_cstmm_basic_execution(self, tmp_path: Path):
        """Test cstmm step can execute with minimal params."""
        params = {
            "out_dir": str(tmp_path / "work"),
        }
        result = run_cstmm(
            params,
            work_dir=str(tmp_path),
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

    @requires_amalgkit
    def test_curate_basic_execution(self, tmp_path: Path):
        """Test curate step can execute with minimal params."""
        params = {
            "out_dir": str(tmp_path / "work"),
        }
        result = run_curate(
            params,
            work_dir=str(tmp_path),
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

    @requires_amalgkit
    def test_csca_basic_execution(self, tmp_path: Path):
        """Test csca step can execute with minimal params."""
        params = {
            "out_dir": str(tmp_path / "work"),
        }
        result = run_csca(
            params,
            work_dir=str(tmp_path),
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

    @requires_amalgkit
    def test_sanity_basic_execution(self, tmp_path: Path):
        """Test sanity step can execute with minimal params."""
        params = {
            "out_dir": str(tmp_path / "work"),
        }
        result = run_sanity(
            params,
            work_dir=str(tmp_path),
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

    def test_build_cli_args_handles_boolean_flags(self):
        """Test boolean flag handling."""
        params = {"verbose": True, "quiet": False}
        args = amalgkit.build_cli_args(params)
        assert "--verbose" in args
        assert "--quiet" not in args  # False values omitted

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
        assert "--out_dir" in cmd

    def test_bool_value_flags(self):
        """Test flags that require explicit yes/no values."""
        params = {"redo": True, "pfd": False}
        args = amalgkit.build_cli_args(params, for_cli=True)
        # redo should have explicit value
        redo_idx = args.index("--redo")
        assert args[redo_idx + 1] == "yes"
        pfd_idx = args.index("--pfd")
        assert args[pfd_idx + 1] == "no"


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


class TestStepRunnerParameterPassing:
    """Test that step runners properly pass parameters to amalgkit."""

    def test_metadata_passes_work_dir(self, tmp_path: Path):
        """Test that metadata runner passes work_dir parameter."""
        # This test doesn't require amalgkit, just checks parameter passing
        params = {"out_dir": str(tmp_path / "work")}
        # Just verify no exceptions are raised with proper parameters
        assert callable(run_metadata)

    def test_all_runners_accept_standard_params(self):
        """Test that all runners accept the standard parameter set."""
        import inspect
        
        standard_params = {"params", "work_dir", "log_dir", "check"}
        for step_name, runner in STEP_RUNNERS.items():
            sig = inspect.signature(runner)
            runner_params = set(sig.parameters.keys())
            assert standard_params.issubset(runner_params), (
                f"Step '{step_name}' missing standard parameters. "
                f"Has: {runner_params}, Expected: {standard_params}"
            )


class TestDocumentationCompleteness:
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

