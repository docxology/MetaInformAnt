"""The life-events scripts must only offer options that the library implements.

``--outcome-relationship`` offered ``health_focused``, ``education_focused`` and ``complex`` (the
default), none of which the library implements: they silently produced random outcomes. The
orchestrator also offered ``--model-type lstm``, which ``EventSequencePredictor`` cannot fit.
"""

from __future__ import annotations

import contextlib
import importlib.util
import sys
from pathlib import Path
from types import ModuleType
from typing import Iterator, List

import pytest

from metainformant.core import io
from metainformant.life_events.core.utils import generate_synthetic_life_events, load_sequences_from_json

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = REPO_ROOT / "scripts" / "life_events"
IMPLEMENTED_RELATIONSHIPS = ["random", "event_based", "temporal"]
UNIMPLEMENTED_RELATIONSHIPS = ["complex", "health_focused", "education_focused"]
NEGATIVE_EVENTS = {"diagnosis", "divorce", "financial_difficulty"}


def load_script(filename: str) -> ModuleType:
    """Import a script file as a module without leaving its sys.path edit behind."""
    saved_path = list(sys.path)
    try:
        spec = importlib.util.spec_from_file_location(f"life_events_script_{Path(filename).stem}", SCRIPTS / filename)
        assert spec is not None and spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path[:] = saved_path


@contextlib.contextmanager
def command_line(*arguments: str) -> Iterator[None]:
    saved = sys.argv
    sys.argv = ["script", *arguments]
    try:
        yield
    finally:
        sys.argv = saved


@pytest.fixture(scope="module")
def generate_script() -> ModuleType:
    return load_script("generate_synthetic_data.py")


@pytest.fixture(scope="module")
def analysis_script() -> ModuleType:
    return load_script("run_life_events_analysis.py")


@pytest.mark.parametrize("name", UNIMPLEMENTED_RELATIONSHIPS)
def test_unimplemented_outcome_relationships_are_rejected_by_both_scripts(
    generate_script, analysis_script, name, capsys: pytest.CaptureFixture[str]
):
    for script, extra in ((generate_script, []), (analysis_script, ["--synthetic"])):
        with command_line("--outcome-relationship", name, *extra):
            with pytest.raises(SystemExit) as exit_info:
                script.parse_args()
        assert exit_info.value.code == 2
    assert capsys.readouterr().err.count("invalid choice") == 2


@pytest.mark.parametrize("name", IMPLEMENTED_RELATIONSHIPS)
def test_implemented_outcome_relationships_are_accepted_by_both_scripts_and_the_library(
    generate_script, analysis_script, name
):
    with command_line("--outcome-relationship", name):
        assert generate_script.parse_args().outcome_relationship == name
    with command_line("--outcome-relationship", name, "--synthetic"):
        assert analysis_script.parse_args().outcome_relationship == name

    _, outcomes = generate_synthetic_life_events(
        n_sequences=3, generate_outcomes=True, outcome_relationship=name, random_state=0
    )
    assert outcomes is not None and len(outcomes) == 3


def test_default_outcome_relationship_is_implemented_by_the_library(generate_script, analysis_script):
    with command_line():
        generate_default = generate_script.parse_args().outcome_relationship
    with command_line("--synthetic"):
        analysis_default = analysis_script.parse_args().outcome_relationship

    assert generate_default == analysis_default == "random"
    assert generate_default in IMPLEMENTED_RELATIONSHIPS


def test_orchestrator_no_longer_offers_the_unsupported_lstm_model(analysis_script, capsys: pytest.CaptureFixture[str]):
    with command_line("--synthetic", "--model-type", "lstm"):
        with pytest.raises(SystemExit):
            analysis_script.parse_args()
    assert "invalid choice" in capsys.readouterr().err

    for model_type in ("embedding", "simple"):
        with command_line("--synthetic", "--model-type", model_type):
            assert analysis_script.parse_args().model_type == model_type


def assert_outcomes_follow(rule: str, sequences: List, outcomes: List[int]) -> None:
    if rule == "event_based":
        expected = [int(any(event.event_type in NEGATIVE_EVENTS for event in seq.events)) for seq in sequences]
    else:
        expected = [int(len(seq.events) > 12) for seq in sequences]
    assert outcomes == expected
    assert set(outcomes) == {0, 1}


@pytest.mark.parametrize("rule", ["event_based", "temporal"])
def test_generate_script_writes_outcomes_that_follow_the_requested_relationship(generate_script, tmp_path, rule):
    output = tmp_path / "synthetic.json"

    with command_line(
        "--n-sequences",
        "40",
        "--min-events",
        "3",
        "--max-events",
        "20",
        "--generate-outcomes",
        "--outcome-relationship",
        rule,
        "--output",
        str(output),
    ):
        assert generate_script.main() == 0

    sequences = load_sequences_from_json(output)
    outcomes = io.load_json(tmp_path / "synthetic_outcomes.json")["outcomes"]
    assert len(sequences) == len(outcomes) == 40
    assert_outcomes_follow(rule, sequences, outcomes)


@pytest.mark.parametrize("rule", ["event_based", "temporal"])
def test_orchestrator_generates_outcomes_with_the_requested_relationship(analysis_script, rule):
    sequences, outcomes = analysis_script.generate_synthetic_data(
        n_sequences=40,
        min_events=3,
        max_events=20,
        generate_outcomes=True,
        outcome_relationship=rule,
        random_state=7,
    )

    assert_outcomes_follow(rule, sequences, list(outcomes))
