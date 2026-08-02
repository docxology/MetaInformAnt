"""Tests for script execution module."""

from __future__ import annotations

import subprocess
import sys
import textwrap
from pathlib import Path

from metainformant.menu.core.executor import (
    execute_bash_script,
    execute_python_script,
    execute_script,
    validate_script_executable,
)


def _run_prompt_child(*, required: list[str], optional: list[str], input_text: str) -> subprocess.CompletedProcess[str]:
    """Exercise ``prompt_for_args`` with real child-process stdin."""

    source = textwrap.dedent(
        f"""
        from pathlib import Path
        from metainformant.menu.core.discovery import ScriptInfo
        from metainformant.menu.core.executor import prompt_for_args

        info = ScriptInfo(
            path=Path("script.py"),
            name="script",
            description="test",
            category="test",
            script_type="python",
            required_args={required!r},
            optional_args={optional!r},
        )
        print(repr(prompt_for_args(info)))
        """
    )
    return subprocess.run(
        [sys.executable, "-c", source],
        input=input_text,
        text=True,
        capture_output=True,
        check=False,
    )


class TestValidateScriptExecutable:
    """Tests for script validation."""

    def test_validate_python_script(self, tmp_path: Path) -> None:
        """Test validating Python script."""
        script_path = tmp_path / "test.py"
        script_path.write_text("print('test')")
        assert validate_script_executable(script_path) is True

    def test_validate_bash_script(self, tmp_path: Path) -> None:
        """Test validating bash script."""
        script_path = tmp_path / "test.sh"
        script_path.write_text("#!/bin/bash\necho 'test'")
        assert validate_script_executable(script_path) is True

    def test_validate_nonexistent_script(self, tmp_path: Path) -> None:
        """Test validating nonexistent script."""
        script_path = tmp_path / "nonexistent.py"
        assert validate_script_executable(script_path) is False


class TestExecutePythonScript:
    """Tests for Python script execution."""

    def test_execute_python_script_success(self, tmp_path: Path) -> None:
        """Test executing successful Python script."""
        script_path = tmp_path / "test.py"
        script_path.write_text("import sys; sys.exit(0)")
        exit_code = execute_python_script(script_path)
        assert exit_code == 0

    def test_execute_python_script_failure(self, tmp_path: Path) -> None:
        """Test executing Python script that fails."""
        script_path = tmp_path / "test.py"
        script_path.write_text("import sys; sys.exit(1)")
        exit_code = execute_python_script(script_path)
        assert exit_code == 1

    def test_execute_python_script_with_args(self, tmp_path: Path) -> None:
        """Test executing Python script with arguments."""
        script_path = tmp_path / "test.py"
        script_path.write_text("import sys; print(sys.argv[1]); sys.exit(0)")
        exit_code = execute_python_script(script_path, ["test_arg"])
        assert exit_code == 0

    def test_execute_python_script_nonexistent(self, tmp_path: Path) -> None:
        """Test executing nonexistent Python script."""
        script_path = tmp_path / "nonexistent.py"
        exit_code = execute_python_script(script_path)
        assert exit_code == 1


class TestExecuteBashScript:
    """Tests for bash script execution."""

    def test_execute_bash_script_success(self, tmp_path: Path) -> None:
        """Test executing successful bash script."""
        script_path = tmp_path / "test.sh"
        script_path.write_text("#!/bin/bash\nexit 0")
        script_path.chmod(0o755)
        exit_code = execute_bash_script(script_path)
        assert exit_code == 0

    def test_execute_bash_script_failure(self, tmp_path: Path) -> None:
        """Test executing bash script that fails."""
        script_path = tmp_path / "test.sh"
        script_path.write_text("#!/bin/bash\nexit 1")
        script_path.chmod(0o755)
        exit_code = execute_bash_script(script_path)
        assert exit_code == 1

    def test_execute_bash_script_with_args(self, tmp_path: Path) -> None:
        """Test executing bash script with arguments."""
        script_path = tmp_path / "test.sh"
        script_path.write_text("#!/bin/bash\necho $1\nexit 0")
        script_path.chmod(0o755)
        exit_code = execute_bash_script(script_path, ["test_arg"])
        assert exit_code == 0

    def test_execute_bash_script_nonexistent(self, tmp_path: Path) -> None:
        """Test executing nonexistent bash script."""
        script_path = tmp_path / "nonexistent.sh"
        exit_code = execute_bash_script(script_path)
        assert exit_code == 1


class TestExecuteScript:
    """Tests for generic script execution."""

    def test_execute_python_via_generic(self, tmp_path: Path) -> None:
        """Test executing Python script via generic function."""
        script_path = tmp_path / "test.py"
        script_path.write_text("import sys; sys.exit(0)")
        exit_code = execute_script(script_path)
        assert exit_code == 0

    def test_execute_bash_via_generic(self, tmp_path: Path) -> None:
        """Test executing bash script via generic function."""
        script_path = tmp_path / "test.sh"
        script_path.write_text("#!/bin/bash\nexit 0")
        script_path.chmod(0o755)
        exit_code = execute_script(script_path)
        assert exit_code == 0

    def test_execute_unsupported_type(self, tmp_path: Path) -> None:
        """Test executing unsupported script type."""
        script_path = tmp_path / "test.txt"
        script_path.write_text("test")
        exit_code = execute_script(script_path)
        assert exit_code == 1


class TestPromptForArgs:
    """Tests for argument prompting through a real child-process stdin."""

    def test_prompt_for_args_no_args(self) -> None:
        """Test prompting when script has no arguments."""
        result = _run_prompt_child(required=[], optional=[], input_text="")
        assert result.returncode == 0
        assert result.stdout.strip() == "[]"

    def test_prompt_for_args_required(self) -> None:
        """Test prompting for required arguments."""
        result = _run_prompt_child(required=["input"], optional=[], input_text="\nvalue\n")
        assert result.returncode == 0
        assert "input is required" in result.stdout
        assert "['--input', 'value']" in result.stdout

    def test_prompt_for_args_optional(self) -> None:
        """Test prompting for optional arguments."""
        result = _run_prompt_child(required=[], optional=["output"], input_text="result.tsv\n")
        assert result.returncode == 0
        assert "['--output', 'result.tsv']" in result.stdout

    def test_prompt_for_args_mixed(self) -> None:
        """Test prompting for both required and optional arguments."""
        result = _run_prompt_child(required=["input"], optional=["output"], input_text="input.tsv\n\n")
        assert result.returncode == 0
        assert "['--input', 'input.tsv']" in result.stdout
