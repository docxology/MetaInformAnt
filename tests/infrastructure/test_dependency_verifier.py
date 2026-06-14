from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load_verifier_module():
    script_path = Path(__file__).resolve().parents[2] / "scripts" / "package" / "verify_dependencies.py"
    spec = importlib.util.spec_from_file_location("verify_dependencies", script_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_run_command_reports_real_python_version() -> None:
    verifier = _load_verifier_module()

    status, detail = verifier.run_command(["python", "--version"])

    assert status == "ok"
    assert "Python" in detail


def test_dependency_report_writes_json_and_markdown(tmp_path: Path) -> None:
    verifier = _load_verifier_module()
    report = verifier.VerificationReport(
        generated_at="2026-05-25T00:00:00+0000",
        python=[verifier.CheckResult("python", "ok", "available")],
        cli=[verifier.CheckResult("python", "ok", "Python available", path="/usr/bin/python")],
    )

    verifier.write_reports(report, tmp_path)

    assert (tmp_path / "verification.json").is_file()
    markdown = (tmp_path / "verification.md").read_text()
    assert "MetaInformAnt Dependency Verification" in markdown
    assert "`python`" in markdown
