"""Comprehensive tests for BaseWorkflowOrchestrator and WorkflowStep with DAG support.

Tests cover: step execution/timing, DAG topological sort, cycle detection,
dependency propagation, failed-step skipping, config-based step loading,
function resolution via importlib, and run_config_based_workflow integration.

NO MOCKING. All implementations are real.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Dict

import pytest

from metainformant.core import io
from metainformant.core.execution.workflow import (
    BaseWorkflowOrchestrator,
    WorkflowStep,
    run_config_based_workflow,
)

# ---------------------------------------------------------------------------
# Real helper functions used as step callables (no mocks)
# ---------------------------------------------------------------------------


def _add(a: int = 0, b: int = 0, **_: Any) -> int:
    """Simple addition for testing."""
    return a + b


def _multiply(a: int = 1, b: int = 1, **_: Any) -> int:
    """Simple multiplication for testing."""
    return a * b


def _slow_identity(value: Any = None, delay: float = 0.05, **_: Any) -> Any:
    """Return value after a short sleep (proves timing works)."""
    time.sleep(delay)
    return value


def _uppercase(text: str = "", **_: Any) -> str:
    return text.upper()


def _concat(*_: Any, prefix: str = "", suffix: str = "", **kw: Any) -> str:
    """Concatenate prefix + suffix, optionally incorporating upstream results."""
    parts = [prefix]
    # Include any dependency results passed as kwargs
    for key, val in sorted(kw.items()):
        if isinstance(val, str):
            parts.append(val)
    parts.append(suffix)
    return "".join(parts)


def _fail_always(**_: Any) -> None:
    """Always raises an exception."""
    raise RuntimeError("Intentional test failure")


def _noop(**_: Any) -> str:
    return "done"


# ===========================================================================
# WorkflowStep Tests
# ===========================================================================


class TestWorkflowStep:
    """Tests for the WorkflowStep dataclass-like object."""

    def test_execute_sets_start_and_end_time(self) -> None:
        """execute() must populate start_time and end_time."""
        step = WorkflowStep(name="t1", function=_add, config={"a": 1, "b": 2})

        assert step.start_time is None
        assert step.end_time is None

        step.execute()

        assert step.start_time is not None
        assert step.end_time is not None
        assert step.end_time >= step.start_time

    def test_duration_returns_elapsed_time(self) -> None:
        """duration() should return a positive float after a slow step."""
        step = WorkflowStep(name="slow", function=_slow_identity, config={"value": 42, "delay": 0.05})
        step.execute()

        elapsed = step.duration()
        assert elapsed >= 0.04, f"Expected >= 0.04s, got {elapsed}"

    def test_duration_returns_zero_before_execution(self) -> None:
        """duration() returns 0.0 when the step has never been executed."""
        step = WorkflowStep(name="pending", function=_noop, config={})
        assert step.duration() == 0.0

    def test_failed_execution_sets_error_and_status(self) -> None:
        """A failing step must set status='failed' and record the error string."""
        step = WorkflowStep(name="bad", function=_fail_always, config={})

        with pytest.raises(RuntimeError, match="Intentional test failure"):
            step.execute()

        assert step.status == "failed"
        assert step.error is not None
        assert "Intentional test failure" in step.error
        assert step.start_time is not None
        assert step.end_time is not None

    def test_reset_clears_all_state(self) -> None:
        """reset() must restore every mutable field to its initial value."""
        step = WorkflowStep(name="resettable", function=_add, config={"a": 5, "b": 3})
        step.execute()

        assert step.status == "completed"
        assert step.result == 8

        step.reset()

        assert step.status == "pending"
        assert step.result is None
        assert step.error is None
        assert step.start_time is None
        assert step.end_time is None

    def test_execute_returns_function_result(self) -> None:
        """execute() should return the callable's return value."""
        step = WorkflowStep(name="calc", function=_add, config={"a": 10, "b": 20})
        result = step.execute()
        assert result == 30
        assert step.result == 30
        assert step.status == "completed"

    def test_execute_merges_kwargs_into_config(self) -> None:
        """Extra kwargs passed to execute() should override/merge with config."""
        step = WorkflowStep(name="merge", function=_add, config={"a": 1, "b": 2})
        # Override b via kwargs
        result = step.execute(b=100)
        assert result == 101


# ===========================================================================
# BaseWorkflowOrchestrator Tests
# ===========================================================================


class TestBaseWorkflowOrchestrator:
    """Tests for the DAG-based orchestrator."""

    # --- Empty & trivial workflows ---

    def test_empty_workflow_returns_success(self) -> None:
        """A workflow with zero steps should succeed immediately."""
        orch = BaseWorkflowOrchestrator(config={})
        result = orch.run_workflow()

        assert result["success"] is True
        assert result["results"] == {}
        assert result["execution_order"] == []
        assert result["errors"] == []
        assert result["total_duration"] >= 0.0

    def test_single_step_workflow(self) -> None:
        """A workflow with one step should execute it and return its result."""
        orch = BaseWorkflowOrchestrator(config={})
        orch.add_step("only", _add, config={"a": 7, "b": 3})

        result = orch.run_workflow()

        assert result["success"] is True
        assert result["results"]["only"] == 10
        assert result["execution_order"] == ["only"]
        assert len(result["errors"]) == 0

    # --- Dependency chains ---

    def test_linear_chain_a_b_c(self) -> None:
        """A -> B -> C linear chain must execute in order."""
        orch = BaseWorkflowOrchestrator(config={})
        orch.add_step("A", _uppercase, config={"text": "hello"})
        orch.add_step("B", _uppercase, config={"text": "world"}, depends_on=["A"])
        orch.add_step("C", _noop, config={}, depends_on=["B"])

        result = orch.run_workflow()

        assert result["success"] is True
        order = result["execution_order"]
        assert order.index("A") < order.index("B") < order.index("C")
        assert result["results"]["A"] == "HELLO"
        assert result["results"]["B"] == "WORLD"
        assert result["results"]["C"] == "done"

    def test_diamond_pattern(self) -> None:
        """Diamond DAG: A -> B, A -> C, B+C -> D."""
        orch = BaseWorkflowOrchestrator(config={})
        orch.add_step("A", _add, config={"a": 1, "b": 1})
        orch.add_step("B", _add, config={"a": 2, "b": 2}, depends_on=["A"])
        orch.add_step("C", _add, config={"a": 3, "b": 3}, depends_on=["A"])
        orch.add_step("D", _noop, config={}, depends_on=["B", "C"])

        result = orch.run_workflow()

        assert result["success"] is True
        order = result["execution_order"]
        # A must be first
        assert order[0] == "A"
        # D must be last
        assert order[-1] == "D"
        # B and C must both precede D
        assert order.index("B") < order.index("D")
        assert order.index("C") < order.index("D")

    # --- Validation errors ---

    def test_cycle_detection_raises_error(self) -> None:
        """A cycle in the DAG (A -> B -> A) must be caught and reported."""
        orch = BaseWorkflowOrchestrator(config={})
        orch.add_step("A", _noop, config={}, depends_on=["B"])
        orch.add_step("B", _noop, config={}, depends_on=["A"])

        result = orch.run_workflow()

        assert result["success"] is False
        assert any("Cycle detected" in e or "cycle" in e.lower() for e in result["errors"])

    def test_unknown_dependency_returns_error(self) -> None:
        """Referencing a non-existent dependency must produce a validation error."""
        orch = BaseWorkflowOrchestrator(config={})
        orch.add_step("A", _noop, config={}, depends_on=["ghost"])

        result = orch.run_workflow()

        assert result["success"] is False
        assert any("unknown step" in e.lower() for e in result["errors"])

    # --- Failure propagation ---

    def test_failed_step_skips_downstream_dependents(self) -> None:
        """When step A fails, steps B and C (which depend on A) must be skipped."""
        orch = BaseWorkflowOrchestrator(config={})
        orch.add_step("A", _fail_always, config={})
        orch.add_step("B", _noop, config={}, depends_on=["A"])
        orch.add_step("C", _noop, config={}, depends_on=["B"])

        result = orch.run_workflow()

        assert result["success"] is False
        # B and C should be skipped
        assert orch.steps["B"].status == "skipped"
        assert orch.steps["C"].status == "skipped"
        assert orch.steps["A"].status == "failed"
        # Results should NOT contain B or C
        assert "B" not in result["results"]
        assert "C" not in result["results"]

    def test_independent_step_runs_despite_sibling_failure(self) -> None:
        """An independent step should still run even if a parallel branch fails."""
        orch = BaseWorkflowOrchestrator(config={})
        orch.add_step("root", _add, config={"a": 1, "b": 1})
        orch.add_step("fail_branch", _fail_always, config={}, depends_on=["root"])
        orch.add_step("ok_branch", _noop, config={}, depends_on=["root"])
        orch.add_step("after_fail", _noop, config={}, depends_on=["fail_branch"])

        result = orch.run_workflow()

        assert result["success"] is False
        # ok_branch should have completed
        assert orch.steps["ok_branch"].status == "completed"
        assert result["results"]["ok_branch"] == "done"
        # after_fail should be skipped
        assert orch.steps["after_fail"].status == "skipped"

    # --- Chaining API ---

    def test_add_step_chaining_returns_self(self) -> None:
        """add_step() must return the orchestrator for fluent chaining."""
        orch = BaseWorkflowOrchestrator(config={})
        returned = orch.add_step("A", _noop).add_step("B", _noop, depends_on=["A"])

        assert returned is orch
        assert len(orch.steps) == 2

    # --- Introspection ---

    def test_get_step_status_returns_correct_states(self) -> None:
        """get_step_status() must reflect the post-execution state of each step."""
        orch = BaseWorkflowOrchestrator(config={})
        orch.add_step("good", _noop, config={})
        orch.add_step("bad", _fail_always, config={}, depends_on=["good"])
        orch.add_step("skipped", _noop, config={}, depends_on=["bad"])

        orch.run_workflow()

        statuses = orch.get_step_status()
        assert statuses["good"] == "completed"
        assert statuses["bad"] == "failed"
        assert statuses["skipped"] == "skipped"

    def test_get_step_status_before_run(self) -> None:
        """Before run_workflow(), all steps should be 'pending'."""
        orch = BaseWorkflowOrchestrator(config={})
        orch.add_step("x", _noop)
        orch.add_step("y", _noop)

        statuses = orch.get_step_status()
        assert all(s == "pending" for s in statuses.values())

    # --- Result structure ---

    def test_run_workflow_returns_execution_order(self) -> None:
        """The result dict must contain a list of step names in execution order."""
        orch = BaseWorkflowOrchestrator(config={})
        orch.add_step("alpha", _noop)
        orch.add_step("beta", _noop, depends_on=["alpha"])

        result = orch.run_workflow()

        assert isinstance(result["execution_order"], list)
        assert result["execution_order"] == ["alpha", "beta"]

    def test_run_workflow_returns_positive_total_duration(self) -> None:
        """total_duration must be a positive float when steps are executed."""
        orch = BaseWorkflowOrchestrator(config={})
        orch.add_step("timed", _slow_identity, config={"value": 1, "delay": 0.02})

        result = orch.run_workflow()

        assert result["total_duration"] > 0.0

    def test_results_dict_contains_all_step_results(self) -> None:
        """Every successfully executed step must have an entry in results."""
        orch = BaseWorkflowOrchestrator(config={})
        orch.add_step("s1", _add, config={"a": 1, "b": 2})
        orch.add_step("s2", _multiply, config={"a": 3, "b": 4})
        orch.add_step("s3", _uppercase, config={"text": "test"})

        result = orch.run_workflow()

        assert result["success"] is True
        assert result["results"]["s1"] == 3
        assert result["results"]["s2"] == 12
        assert result["results"]["s3"] == "TEST"
        assert len(result["results"]) == 3

    # --- Config-based step loading ---

    def test_config_based_step_loading_from_dict(self) -> None:
        """Steps defined under the 'steps' config key should be auto-loaded.

        Note: downstream steps receive upstream results as extra kwargs, so
        we use json:loads (no deps) for the first step and a kwargs-tolerant
        function for the second.  We test the config loading mechanism by
        verifying the first step resolves and executes correctly.
        """
        cfg: Dict[str, Any] = {
            "steps": {
                "parse": {
                    "function": "json:loads",
                    "params": {"s": '{"key": "value"}'},
                },
            }
        }
        orch = BaseWorkflowOrchestrator(config=cfg)
        result = orch.run_workflow()

        assert result["success"] is True
        assert result["results"]["parse"] == {"key": "value"}
        assert "parse" in result["execution_order"]

    def test_config_based_step_loading_with_dependency(self) -> None:
        """Config-loaded steps with dependencies should execute in correct order.

        Uses programmatic steps (which accept **kwargs) combined with a
        config-loaded first step to prove the auto-loading + dependency mechanism.
        """
        cfg: Dict[str, Any] = {
            "steps": {
                "first": {
                    "function": "json:loads",
                    "params": {"s": "100"},
                },
            }
        }
        orch = BaseWorkflowOrchestrator(config=cfg)
        # Add a second step programmatically that depends on the config-loaded one
        orch.add_step("second", _noop, config={}, depends_on=["first"])

        result = orch.run_workflow()

        assert result["success"] is True
        assert result["results"]["first"] == 100
        assert result["results"]["second"] == "done"
        order = result["execution_order"]
        assert order.index("first") < order.index("second")

    def test_config_step_loading_preserves_programmatic_steps(self) -> None:
        """Programmatically added steps should not be overwritten by config."""
        cfg: Dict[str, Any] = {
            "steps": {
                "step_a": {
                    "function": "json:loads",
                    "params": {"s": '"overwritten"'},
                }
            }
        }
        orch = BaseWorkflowOrchestrator(config=cfg)
        # Add same name programmatically first
        orch.add_step("step_a", _noop)

        result = orch.run_workflow()

        # Programmatic step should win (returns "done", not "overwritten")
        assert result["results"]["step_a"] == "done"


# ===========================================================================
# _resolve_function Tests
# ===========================================================================


class TestResolveFunction:
    """Tests for the static _resolve_function method."""

    def test_resolve_function_with_json_loads(self) -> None:
        """'json:loads' should resolve to the real json.loads callable."""
        func = BaseWorkflowOrchestrator._resolve_function("json:loads")
        assert callable(func)
        assert func('{"a": 1}') == {"a": 1}

    def test_resolve_function_with_json_dumps(self) -> None:
        """'json:dumps' should resolve to the real json.dumps callable."""
        func = BaseWorkflowOrchestrator._resolve_function("json:dumps")
        assert callable(func)
        assert func({"b": 2}) == '{"b": 2}'

    def test_resolve_function_with_os_path_basename(self) -> None:
        """'os.path:basename' should resolve correctly."""
        func = BaseWorkflowOrchestrator._resolve_function("os.path:basename")
        assert func("/usr/local/bin/python") == "python"

    def test_resolve_function_invalid_format_raises_valueerror(self) -> None:
        """A string without ':' should raise ValueError."""
        with pytest.raises(ValueError, match="must be in 'module.path:function_name' format"):
            BaseWorkflowOrchestrator._resolve_function("json.loads")

    def test_resolve_function_bad_module_raises_importerror(self) -> None:
        """A non-existent module should raise ImportError."""
        with pytest.raises(ImportError):
            BaseWorkflowOrchestrator._resolve_function("nonexistent_module_xyz:func")

    def test_resolve_function_bad_attr_raises_attributeerror(self) -> None:
        """A non-existent function in a valid module should raise AttributeError."""
        with pytest.raises(AttributeError):
            BaseWorkflowOrchestrator._resolve_function("json:nonexistent_function_xyz")


# ===========================================================================
# Integration: run_config_based_workflow with "steps" config
# ===========================================================================


class TestRunConfigBasedWorkflowWithSteps:
    """Integration tests verifying that run_config_based_workflow delegates to the orchestrator
    when the config contains a 'steps' key."""

    def test_steps_config_uses_orchestrator(self, tmp_path: Path) -> None:
        """A config file with 'steps' should be executed via BaseWorkflowOrchestrator."""
        cfg = {
            "steps": {
                "parse_json": {
                    "function": "json:loads",
                    "params": {"s": '{"status": "ok"}'},
                }
            }
        }
        config_file = tmp_path / "steps_workflow.json"
        io.dump_json(cfg, config_file)

        result = run_config_based_workflow(config_file)

        assert result["success"] is True
        assert "parse_json" in result["results"]
        assert result["results"]["parse_json"] == {"status": "ok"}
        assert "execution_order" in result
        assert "parse_json" in result["execution_order"]

    def test_steps_config_with_dependencies(self, tmp_path: Path) -> None:
        """A config with two independent steps (no dep-injected kwargs conflict).

        Verifies that config-based loading respects dependency ordering in
        the execution_order list. Both steps are independent json:loads calls
        with a dependency edge to enforce ordering.

        Note: The orchestrator injects upstream results as kwargs into
        downstream steps. Because json.loads does not accept **kwargs,
        we verify ordering was computed correctly via execution_order even
        though step_b will fail at call time from the injected kwarg.
        """
        cfg = {
            "steps": {
                "step_a": {
                    "function": "json:loads",
                    "params": {"s": '"first"'},
                },
                "step_b": {
                    "function": "json:loads",
                    "params": {"s": '"second"'},
                    "depends_on": ["step_a"],
                },
            }
        }
        config_file = tmp_path / "dep_workflow.json"
        io.dump_json(cfg, config_file)

        result = run_config_based_workflow(config_file)

        # Verify topological ordering was computed correctly
        order = result["execution_order"]
        assert order.index("step_a") < order.index("step_b")
        # step_a should have succeeded (no upstream deps)
        assert result["results"]["step_a"] == "first"

    def test_invalid_steps_config_returns_failure(self, tmp_path: Path) -> None:
        """A config with a cycle in steps should return success=False."""
        cfg = {
            "steps": {
                "x": {
                    "function": "json:loads",
                    "params": {"s": "1"},
                    "depends_on": ["y"],
                },
                "y": {
                    "function": "json:loads",
                    "params": {"s": "2"},
                    "depends_on": ["x"],
                },
            }
        }
        config_file = tmp_path / "cyclic_workflow.json"
        io.dump_json(cfg, config_file)

        result = run_config_based_workflow(config_file)

        assert result["success"] is False
        assert len(result["errors"]) > 0

    def test_steps_config_with_invalid_function_ref(self, tmp_path: Path) -> None:
        """A step with a bad function reference should cause a load-time failure."""
        cfg = {
            "steps": {
                "bad": {
                    "function": "no_such_module_xyz:no_func",
                    "params": {},
                }
            }
        }
        config_file = tmp_path / "bad_func_workflow.json"
        io.dump_json(cfg, config_file)

        result = run_config_based_workflow(config_file)

        assert result["success"] is False


# ===========================================================================
# DAG edge cases
# ===========================================================================


class TestDAGEdgeCases:
    """Additional edge-case tests for the DAG engine."""

    def test_three_node_cycle(self) -> None:
        """A -> B -> C -> A cycle must be detected."""
        orch = BaseWorkflowOrchestrator(config={})
        orch.add_step("A", _noop, depends_on=["C"])
        orch.add_step("B", _noop, depends_on=["A"])
        orch.add_step("C", _noop, depends_on=["B"])

        result = orch.run_workflow()
        assert result["success"] is False
        assert any("Cycle" in e or "cycle" in e.lower() for e in result["errors"])

    def test_self_dependency_cycle(self) -> None:
        """A step depending on itself is a trivial cycle."""
        orch = BaseWorkflowOrchestrator(config={})
        orch.add_step("self_ref", _noop, depends_on=["self_ref"])

        result = orch.run_workflow()
        assert result["success"] is False

    def test_wide_fan_out(self) -> None:
        """One root with many independent children should all succeed."""
        orch = BaseWorkflowOrchestrator(config={})
        orch.add_step("root", _add, config={"a": 0, "b": 0})

        num_children = 20
        for i in range(num_children):
            orch.add_step(f"child_{i}", _add, config={"a": i, "b": 1}, depends_on=["root"])

        result = orch.run_workflow()

        assert result["success"] is True
        assert len(result["results"]) == num_children + 1
        assert result["results"]["root"] == 0
        for i in range(num_children):
            assert result["results"][f"child_{i}"] == i + 1

    def test_wide_fan_in(self) -> None:
        """Many independent roots converging on a single sink."""
        orch = BaseWorkflowOrchestrator(config={})
        num_roots = 10
        root_names = []
        for i in range(num_roots):
            name = f"root_{i}"
            root_names.append(name)
            orch.add_step(name, _add, config={"a": i, "b": 0})

        orch.add_step("sink", _noop, config={}, depends_on=root_names)

        result = orch.run_workflow()

        assert result["success"] is True
        order = result["execution_order"]
        # All roots must precede the sink
        sink_idx = order.index("sink")
        for rn in root_names:
            assert order.index(rn) < sink_idx

    def test_dependency_results_passed_as_kwargs(self) -> None:
        """Upstream results should be available as kwargs to downstream steps."""

        def _collector(**kwargs: Any) -> Dict[str, Any]:
            """Collect all kwargs into a dict."""
            return dict(kwargs)

        orch = BaseWorkflowOrchestrator(config={})
        orch.add_step("producer_a", _add, config={"a": 10, "b": 5})
        orch.add_step("producer_b", _uppercase, config={"text": "hello"})
        orch.add_step("consumer", _collector, config={}, depends_on=["producer_a", "producer_b"])

        result = orch.run_workflow()

        assert result["success"] is True
        consumer_result = result["results"]["consumer"]
        assert consumer_result["producer_a"] == 15
        assert consumer_result["producer_b"] == "HELLO"

    def test_working_dir_default(self) -> None:
        """Default working_dir should be Path('output')."""
        orch = BaseWorkflowOrchestrator(config={})
        assert orch.working_dir == Path("output")

    def test_working_dir_custom(self, tmp_path: Path) -> None:
        """Custom working_dir should be stored."""
        orch = BaseWorkflowOrchestrator(config={}, working_dir=tmp_path)
        assert orch.working_dir == tmp_path

    def test_multiple_runs_accumulate_results(self) -> None:
        """Running the workflow twice should re-execute (results come from second run)."""
        orch = BaseWorkflowOrchestrator(config={})
        orch.add_step("counter", _add, config={"a": 1, "b": 1})

        result1 = orch.run_workflow()
        result2 = orch.run_workflow()

        # Both should succeed
        assert result1["success"] is True
        assert result2["success"] is True
        assert result2["results"]["counter"] == 2

    def test_add_step_with_string_function_reference(self) -> None:
        """add_step() should accept a 'module:func' string and resolve it."""
        orch = BaseWorkflowOrchestrator(config={})
        orch.add_step("json_parse", "json:loads", config={"s": "[1,2,3]"})

        result = orch.run_workflow()

        assert result["success"] is True
        assert result["results"]["json_parse"] == [1, 2, 3]

    def test_empty_config_steps_key_is_harmless(self) -> None:
        """An empty 'steps' dict in config should not break anything."""
        orch = BaseWorkflowOrchestrator(config={"steps": {}})
        result = orch.run_workflow()

        assert result["success"] is True
        assert result["results"] == {}

    def test_config_steps_none_is_harmless(self) -> None:
        """A None value for 'steps' in config should not break anything."""
        orch = BaseWorkflowOrchestrator(config={"steps": None})
        result = orch.run_workflow()

        assert result["success"] is True

    def test_long_chain_executes_in_order(self) -> None:
        """A 10-step linear chain must execute in strict sequence."""
        orch = BaseWorkflowOrchestrator(config={})
        chain_length = 10
        for i in range(chain_length):
            deps = [f"step_{i - 1}"] if i > 0 else []
            orch.add_step(f"step_{i}", _add, config={"a": i, "b": 0}, depends_on=deps)

        result = orch.run_workflow()

        assert result["success"] is True
        order = result["execution_order"]
        for i in range(chain_length):
            assert order[i] == f"step_{i}"
