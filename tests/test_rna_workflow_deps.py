from __future__ import annotations

from pathlib import Path


def test_workflow_skips_steps_when_missing_deps(monkeypatch, tmp_path: Path):
    from metainformant.rna.workflow import AmalgkitWorkflowConfig, execute_workflow
    from metainformant.rna.workflow import _sanitize_params_for_subcommand
    import metainformant.rna.workflow as wf

    # Make metadata succeed quickly by monkeypatching the runner
    class DummyR:
        def __init__(self, code: int = 0):
            self.returncode = code
            self.stdout = ""
            self.stderr = ""

    def fake_runner(params=None, **kwargs):
        return DummyR(0)

    # Replace runners so we don't call external tools
    monkeypatch.setitem(wf._steps_mod.STEP_RUNNERS, "metadata", fake_runner)
    monkeypatch.setitem(wf._steps_mod.STEP_RUNNERS, "select", fake_runner)
    monkeypatch.setitem(wf._steps_mod.STEP_RUNNERS, "config", fake_runner)
    monkeypatch.setitem(wf._steps_mod.STEP_RUNNERS, "sanity", fake_runner)

    # Force deps check to report OK for metadata/config/select but missing for others
    def fake_check_step_dependencies(step: str):
        class S:
            def __init__(self, step: str):
                self.step = step
                self.required = ["dummy"]
                self.missing = [] if step in {"metadata", "config", "select", "sanity"} else ["dummy"]

            @property
            def ok(self) -> bool:
                return len(self.missing) == 0

        return S(step)

    monkeypatch.setattr(wf, "check_step_dependencies", fake_check_step_dependencies)

    cfg = AmalgkitWorkflowConfig(work_dir=tmp_path / "work", threads=1)
    # Provide a minimal metadata file so sanity runner doesn't choke if called
    (cfg.work_dir / "metadata").mkdir(parents=True, exist_ok=True)
    (cfg.work_dir / "metadata" / "metadata.tsv").write_text("tissue\nbrain\n", encoding="utf-8")
    codes = execute_workflow(cfg, check=False)
    # Expect metadata/config/select/sanity run (0); others skipped (126)
    assert codes.count(0) >= 2
    assert 126 in codes


