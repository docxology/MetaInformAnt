"""Comprehensive tests for core.engine.workflow_manager module.

Tests cover:
- Stage enum lifecycle values
- PipelineItem creation and defaults
- PipelinePhase creation and filter functions
- BasePipelineManager: add_item, mark_running, mark_done, mark_failed
- BasePipelineManager: run() with custom phases, failing phases, filters
- BasePipelineManager: edge cases (empty pipeline, no phases)
- WorkflowManager backward compatibility: SampleStage, SampleState, add_sample
"""

from __future__ import annotations

from dataclasses import fields
from pathlib import Path
from typing import Any, Dict, List

import pytest

from metainformant.core.engine.workflow_manager import (
    BasePipelineManager,
    PipelineItem,
    PipelinePhase,
    SampleStage,
    SampleState,
    Stage,
    WorkflowManager,
)
from metainformant.core.ui.tui import CYAN

# ---------------------------------------------------------------------------
# Stage enum
# ---------------------------------------------------------------------------


class TestStageEnum:
    """Tests for the Stage lifecycle enum."""

    def test_pending_value(self) -> None:
        """Stage.PENDING has value 'Pending'."""
        assert Stage.PENDING.value == "Pending"

    def test_running_value(self) -> None:
        """Stage.RUNNING has value 'Running'."""
        assert Stage.RUNNING.value == "Running"

    def test_done_value(self) -> None:
        """Stage.DONE has value 'Done'."""
        assert Stage.DONE.value == "Done"

    def test_failed_value(self) -> None:
        """Stage.FAILED has value 'Failed'."""
        assert Stage.FAILED.value == "Failed"

    def test_all_stages_present(self) -> None:
        """Stage enum contains exactly four members."""
        members = list(Stage)
        assert len(members) == 4
        assert set(members) == {Stage.PENDING, Stage.RUNNING, Stage.DONE, Stage.FAILED}


# ---------------------------------------------------------------------------
# PipelineItem dataclass
# ---------------------------------------------------------------------------


class TestPipelineItem:
    """Tests for PipelineItem creation and defaults."""

    def test_creation_defaults(self) -> None:
        """PipelineItem has correct defaults for stage, metadata, error."""
        item = PipelineItem(item_id="sample_001")
        assert item.item_id == "sample_001"
        assert item.metadata == {}
        assert item.stage == Stage.PENDING
        assert item.error == ""

    def test_creation_with_metadata(self) -> None:
        """PipelineItem accepts arbitrary metadata dict."""
        meta = {"species": "Apis mellifera", "tissue": "brain"}
        item = PipelineItem(item_id="SRR123", metadata=meta)
        assert item.metadata["species"] == "Apis mellifera"
        assert item.metadata["tissue"] == "brain"

    def test_creation_with_explicit_stage(self) -> None:
        """PipelineItem can be created with a non-default stage."""
        item = PipelineItem(item_id="x", stage=Stage.RUNNING)
        assert item.stage == Stage.RUNNING

    def test_creation_with_error(self) -> None:
        """PipelineItem can be created with a pre-set error string."""
        item = PipelineItem(item_id="x", stage=Stage.FAILED, error="timeout")
        assert item.error == "timeout"

    def test_metadata_independence(self) -> None:
        """Each PipelineItem gets its own metadata dict (no shared default)."""
        a = PipelineItem(item_id="a")
        b = PipelineItem(item_id="b")
        a.metadata["key"] = "val"
        assert "key" not in b.metadata


# ---------------------------------------------------------------------------
# PipelinePhase dataclass
# ---------------------------------------------------------------------------


class TestPipelinePhase:
    """Tests for PipelinePhase creation."""

    def test_creation_with_handler(self) -> None:
        """PipelinePhase stores name, handler, and color."""

        def dummy_handler(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
            pass

        phase = PipelinePhase(name="Download", handler=dummy_handler, color=CYAN)
        assert phase.name == "Download"
        assert phase.handler is dummy_handler
        assert phase.color == CYAN

    def test_default_filter_selects_pending(self) -> None:
        """Default filter_fn selects items with Stage.PENDING."""

        def noop(m: BasePipelineManager, items: List[PipelineItem]) -> None:
            pass

        phase = PipelinePhase(name="P", handler=noop)
        pending_item = PipelineItem(item_id="a", stage=Stage.PENDING)
        running_item = PipelineItem(item_id="b", stage=Stage.RUNNING)
        done_item = PipelineItem(item_id="c", stage=Stage.DONE)
        failed_item = PipelineItem(item_id="d", stage=Stage.FAILED)

        assert phase.filter_fn(pending_item) is True
        assert phase.filter_fn(running_item) is False
        assert phase.filter_fn(done_item) is False
        assert phase.filter_fn(failed_item) is False

    def test_custom_filter_fn(self) -> None:
        """PipelinePhase accepts a custom filter function."""

        def noop(m: BasePipelineManager, items: List[PipelineItem]) -> None:
            pass

        # Select only items whose metadata has "priority" == "high"
        phase = PipelinePhase(
            name="HighPri",
            handler=noop,
            filter_fn=lambda item: item.metadata.get("priority") == "high",
        )

        high = PipelineItem(item_id="h", metadata={"priority": "high"})
        low = PipelineItem(item_id="l", metadata={"priority": "low"})
        none_meta = PipelineItem(item_id="n")

        assert phase.filter_fn(high) is True
        assert phase.filter_fn(low) is False
        assert phase.filter_fn(none_meta) is False


# ---------------------------------------------------------------------------
# BasePipelineManager
# ---------------------------------------------------------------------------


class TestBasePipelineManagerAddItem:
    """Tests for BasePipelineManager.add_item()."""

    def test_add_item_creates_pipeline_item(self) -> None:
        """add_item registers a PipelineItem in manager.items."""
        manager = BasePipelineManager(phases=[])
        item = manager.add_item("sample_1")
        assert "sample_1" in manager.items
        assert manager.items["sample_1"] is item
        assert item.item_id == "sample_1"
        assert item.stage == Stage.PENDING

    def test_add_item_with_metadata(self) -> None:
        """add_item passes metadata through to the PipelineItem."""
        manager = BasePipelineManager(phases=[])
        meta = {"organism": "honeybee", "reads": 1_000_000}
        item = manager.add_item("SRR999", metadata=meta)
        assert item.metadata["organism"] == "honeybee"
        assert item.metadata["reads"] == 1_000_000

    def test_add_multiple_items(self) -> None:
        """Multiple items can be added and tracked independently."""
        manager = BasePipelineManager(phases=[])
        manager.add_item("a")
        manager.add_item("b")
        manager.add_item("c")
        assert len(manager.items) == 3
        assert set(manager.items.keys()) == {"a", "b", "c"}


class TestBasePipelineManagerStateTransitions:
    """Tests for mark_running, mark_done, mark_failed helpers."""

    def test_mark_running(self) -> None:
        """mark_running transitions item to Stage.RUNNING."""
        manager = BasePipelineManager(phases=[])
        item = manager.add_item("x")
        assert item.stage == Stage.PENDING

        manager.mark_running(item)
        assert item.stage == Stage.RUNNING

    def test_mark_done(self) -> None:
        """mark_done transitions item to Stage.DONE."""
        manager = BasePipelineManager(phases=[])
        item = manager.add_item("x")
        manager.mark_running(item)
        manager.mark_done(item)
        assert item.stage == Stage.DONE

    def test_mark_failed_with_error(self) -> None:
        """mark_failed transitions item to Stage.FAILED and records error."""
        manager = BasePipelineManager(phases=[])
        item = manager.add_item("x")
        manager.mark_running(item)
        manager.mark_failed(item, error="Connection reset")
        assert item.stage == Stage.FAILED
        assert item.error == "Connection reset"

    def test_mark_failed_preserves_metadata(self) -> None:
        """mark_failed does not destroy existing metadata."""
        manager = BasePipelineManager(phases=[])
        item = manager.add_item("x", metadata={"tissue": "wing"})
        manager.mark_failed(item, error="disk full")
        assert item.metadata["tissue"] == "wing"
        assert item.error == "disk full"

    def test_full_lifecycle_pending_to_done(self) -> None:
        """Item can go through PENDING -> RUNNING -> DONE."""
        manager = BasePipelineManager(phases=[])
        item = manager.add_item("lifecycle")
        assert item.stage == Stage.PENDING
        manager.mark_running(item)
        assert item.stage == Stage.RUNNING
        manager.mark_done(item)
        assert item.stage == Stage.DONE

    def test_full_lifecycle_pending_to_failed(self) -> None:
        """Item can go through PENDING -> RUNNING -> FAILED."""
        manager = BasePipelineManager(phases=[])
        item = manager.add_item("lifecycle")
        manager.mark_running(item)
        manager.mark_failed(item, error="bad data")
        assert item.stage == Stage.FAILED


class TestBasePipelineManagerRun:
    """Tests for BasePipelineManager.run() with various phase configurations."""

    def test_run_all_items_reach_done(self) -> None:
        """Pipeline with a single phase that marks all items DONE."""

        def mark_all_done(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
            for item in items:
                manager.mark_running(item)
                manager.mark_done(item)

        phase = PipelinePhase(name="Process", handler=mark_all_done)
        manager = BasePipelineManager(phases=[phase])
        manager.add_item("s1")
        manager.add_item("s2")
        manager.add_item("s3")

        results = manager.run()

        assert results == {"s1": True, "s2": True, "s3": True}
        for item in manager.items.values():
            assert item.stage == Stage.DONE

    def test_run_with_failing_phase_handler(self) -> None:
        """Pipeline with a phase that marks items as FAILED."""

        def fail_all(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
            for item in items:
                manager.mark_running(item)
                manager.mark_failed(item, error="intentional failure")

        phase = PipelinePhase(name="FailPhase", handler=fail_all)
        manager = BasePipelineManager(phases=[phase])
        manager.add_item("f1")
        manager.add_item("f2")

        results = manager.run()

        assert results == {"f1": False, "f2": False}
        for item in manager.items.values():
            assert item.stage == Stage.FAILED
            assert item.error == "intentional failure"

    def test_run_multiple_phases_in_sequence(self) -> None:
        """Items flow through multiple phases sequentially."""
        execution_order: List[str] = []

        def phase_a_handler(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
            execution_order.append("A")
            for item in items:
                item.stage = Stage.RUNNING
                item.metadata["phase_a"] = True
                # Move to a "ready for B" state -- we use DONE for simplicity.
                # Phase B uses a custom filter to pick up items with phase_a=True.
                item.stage = Stage.DONE

        def phase_b_handler(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
            execution_order.append("B")
            for item in items:
                item.metadata["phase_b"] = True

        phase_a = PipelinePhase(
            name="A",
            handler=phase_a_handler,
            filter_fn=lambda item: item.stage == Stage.PENDING,
        )
        phase_b = PipelinePhase(
            name="B",
            handler=phase_b_handler,
            filter_fn=lambda item: item.stage == Stage.DONE and item.metadata.get("phase_a"),
        )

        manager = BasePipelineManager(phases=[phase_a, phase_b])
        manager.add_item("seq1")

        manager.run()

        assert execution_order == ["A", "B"]
        item = manager.items["seq1"]
        assert item.metadata.get("phase_a") is True
        assert item.metadata.get("phase_b") is True

    def test_run_with_phase_filter_partial(self) -> None:
        """Phase filter selects only some items; unselected items stay unchanged."""

        def process_high(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
            for item in items:
                manager.mark_running(item)
                manager.mark_done(item)

        phase = PipelinePhase(
            name="HighPri",
            handler=process_high,
            filter_fn=lambda item: item.metadata.get("priority") == "high",
        )

        manager = BasePipelineManager(phases=[phase])
        manager.add_item("important", metadata={"priority": "high"})
        manager.add_item("regular", metadata={"priority": "low"})
        manager.add_item("unknown")

        results = manager.run()

        assert results["important"] is True
        assert results["regular"] is False  # Never processed, still PENDING
        assert results["unknown"] is False  # Never processed, still PENDING
        assert manager.items["important"].stage == Stage.DONE
        assert manager.items["regular"].stage == Stage.PENDING
        assert manager.items["unknown"].stage == Stage.PENDING

    def test_run_empty_pipeline_no_items(self) -> None:
        """Pipeline with phases but no items runs without error."""

        def should_not_run(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
            raise AssertionError("Handler should not be called with no eligible items")

        phase = PipelinePhase(name="Ghost", handler=should_not_run)
        manager = BasePipelineManager(phases=[phase])

        results = manager.run()
        assert results == {}

    def test_run_no_phases(self) -> None:
        """Pipeline with items but no phases runs without error."""
        manager = BasePipelineManager(phases=[])
        manager.add_item("orphan")

        results = manager.run()

        # Item was never processed, so result is False (not DONE).
        assert results == {"orphan": False}
        assert manager.items["orphan"].stage == Stage.PENDING

    def test_run_empty_pipeline_no_items_no_phases(self) -> None:
        """Completely empty pipeline runs without error."""
        manager = BasePipelineManager(phases=[])
        results = manager.run()
        assert results == {}

    def test_run_mixed_success_and_failure(self) -> None:
        """Phase can mark some items done and others failed."""

        def mixed_handler(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
            for item in items:
                manager.mark_running(item)
                if item.metadata.get("will_fail"):
                    manager.mark_failed(item, error="bad quality")
                else:
                    manager.mark_done(item)

        phase = PipelinePhase(name="QC", handler=mixed_handler)
        manager = BasePipelineManager(phases=[phase])
        manager.add_item("good_1")
        manager.add_item("bad_1", metadata={"will_fail": True})
        manager.add_item("good_2")
        manager.add_item("bad_2", metadata={"will_fail": True})

        results = manager.run()

        assert results["good_1"] is True
        assert results["good_2"] is True
        assert results["bad_1"] is False
        assert results["bad_2"] is False

    def test_run_phase_handler_exception_propagates(self) -> None:
        """If a phase handler raises, the exception propagates but TUI still stops."""

        def exploding_handler(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
            raise RuntimeError("Kaboom")

        phase = PipelinePhase(name="Boom", handler=exploding_handler)
        manager = BasePipelineManager(phases=[phase])
        manager.add_item("victim")

        with pytest.raises(RuntimeError, match="Kaboom"):
            manager.run()

        # TUI should have been stopped despite the exception (finally block).
        assert manager._running is False

    def test_config_passed_to_manager(self) -> None:
        """Config dict is accessible from the manager inside phase handlers."""
        captured_config: Dict[str, Any] = {}

        def config_reader(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
            captured_config.update(manager.config)
            for item in items:
                manager.mark_done(item)

        phase = PipelinePhase(name="ReadCfg", handler=config_reader)
        cfg = {"threads": 16, "species": "Apis mellifera"}
        manager = BasePipelineManager(phases=[phase], config=cfg)
        manager.add_item("x")

        manager.run()

        assert captured_config["threads"] == 16
        assert captured_config["species"] == "Apis mellifera"


class TestBasePipelineManagerTUIOptional:
    """Test that TUI integration is optional -- pipeline works without starting TUI explicitly."""

    def test_pipeline_works_without_explicit_tui_start(self) -> None:
        """Pipeline run() handles the TUI lifecycle internally; callers never call ui.start()."""

        def simple(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
            for item in items:
                manager.mark_done(item)

        phase = PipelinePhase(name="Auto", handler=simple)
        manager = BasePipelineManager(phases=[phase])
        manager.add_item("auto1")

        # run() internally calls ui.start() and ui.stop() -- we do not.
        results = manager.run()
        assert results["auto1"] is True


# ---------------------------------------------------------------------------
# WorkflowManager backward compatibility
# ---------------------------------------------------------------------------


class TestSampleStageEnum:
    """Tests for SampleStage enum backward compatibility."""

    def test_pending_value(self) -> None:
        assert SampleStage.PENDING.value == "Pending"

    def test_downloading_value(self) -> None:
        assert SampleStage.DOWNLOADING.value == "Downloading"

    def test_downloaded_value(self) -> None:
        assert SampleStage.DOWNLOADED.value == "Downloaded"

    def test_extracting_value(self) -> None:
        assert SampleStage.EXTRACTING.value == "Extracting"

    def test_extracted_value(self) -> None:
        assert SampleStage.EXTRACTED.value == "Extracted"

    def test_quantifying_value(self) -> None:
        assert SampleStage.QUANTIFYING.value == "Quantifying"

    def test_done_value(self) -> None:
        assert SampleStage.DONE.value == "Done"

    def test_failed_value(self) -> None:
        assert SampleStage.FAILED.value == "Failed"

    def test_truncated_value(self) -> None:
        assert SampleStage.TRUNCATED.value == "Truncated"

    def test_all_members_count(self) -> None:
        """SampleStage has exactly nine members."""
        assert len(list(SampleStage)) == 9


class TestSampleStateDataclass:
    """Tests for SampleState dataclass backward compatibility."""

    def test_required_fields(self) -> None:
        """SampleState requires sample_id, sra_url, dest_path."""
        state = SampleState(
            sample_id="SRR123",
            sra_url="https://example.com/SRR123",
            dest_path=Path("/tmp/SRR123.sra"),
        )
        assert state.sample_id == "SRR123"
        assert state.sra_url == "https://example.com/SRR123"
        assert state.dest_path == Path("/tmp/SRR123.sra")

    def test_default_values(self) -> None:
        """SampleState defaults: PENDING stage, 0 bytes, empty error."""
        state = SampleState(
            sample_id="SRR999",
            sra_url="https://example.com/SRR999",
            dest_path=Path("/tmp/SRR999.sra"),
        )
        assert state.stage == SampleStage.PENDING
        assert state.current_bytes == 0
        assert state.total_bytes == 0
        assert state.error == ""

    def test_field_names_unchanged(self) -> None:
        """SampleState has the expected set of field names."""
        field_names = {f.name for f in fields(SampleState)}
        expected = {"sample_id", "sra_url", "dest_path", "stage", "current_bytes", "total_bytes", "error"}
        assert field_names == expected


class TestWorkflowManagerInstantiation:
    """Tests for WorkflowManager creation with a real config file."""

    def test_instantiation_with_minimal_config(self, tmp_path: Path) -> None:
        """WorkflowManager loads a YAML config and sets work_dir/species."""
        config_file = tmp_path / "workflow.yaml"
        config_file.write_text("work_dir: output/test_wf\nspecies: apis_mellifera\nthreads: 4\n")

        wm = WorkflowManager(config_path=config_file, max_threads=2)

        assert wm.config["species"] == "apis_mellifera"
        assert wm.config["threads"] == 4
        assert wm.species == "apis_mellifera"
        assert wm.work_dir == Path("output/test_wf")
        assert wm.max_threads == 2

    def test_instantiation_default_work_dir(self, tmp_path: Path) -> None:
        """WorkflowManager uses default work_dir when config omits it."""
        config_file = tmp_path / "minimal.yaml"
        config_file.write_text("species: drosophila\n")

        wm = WorkflowManager(config_path=config_file)

        assert wm.work_dir == Path("output/amalgkit")
        assert wm.species == "drosophila"

    def test_instantiation_default_species(self, tmp_path: Path) -> None:
        """WorkflowManager uses 'unknown' species when config omits it."""
        config_file = tmp_path / "bare.yaml"
        config_file.write_text("threads: 1\n")

        wm = WorkflowManager(config_path=config_file)

        assert wm.species == "unknown"

    def test_has_three_phases(self, tmp_path: Path) -> None:
        """WorkflowManager creates exactly three RNA-seq phases."""
        config_file = tmp_path / "wf.yaml"
        config_file.write_text("species: test\n")

        wm = WorkflowManager(config_path=config_file)

        assert len(wm.phases) == 3
        phase_names = [p.name for p in wm.phases]
        assert phase_names == ["Download", "getfastq", "quant"]


class TestWorkflowManagerAddSample:
    """Tests for WorkflowManager.add_sample() backward-compatible method."""

    def test_add_sample_creates_sample_state(self, tmp_path: Path) -> None:
        """add_sample registers a SampleState in the samples dict."""
        config_file = tmp_path / "wf.yaml"
        config_file.write_text("species: honeybee\n")

        wm = WorkflowManager(config_path=config_file)
        dest = tmp_path / "SRR001.sra"
        wm.add_sample("SRR001", "https://sra.example.com/SRR001", dest)

        assert "SRR001" in wm.samples
        state = wm.samples["SRR001"]
        assert state.sample_id == "SRR001"
        assert state.sra_url == "https://sra.example.com/SRR001"
        assert state.dest_path == dest
        assert state.stage == SampleStage.PENDING

    def test_add_sample_also_creates_pipeline_item(self, tmp_path: Path) -> None:
        """add_sample registers both a SampleState AND a PipelineItem."""
        config_file = tmp_path / "wf.yaml"
        config_file.write_text("species: test\n")

        wm = WorkflowManager(config_path=config_file)
        dest = tmp_path / "SRR002.sra"
        wm.add_sample("SRR002", "https://sra.example.com/SRR002", dest)

        # Pipeline item should exist in the generic items dict.
        assert "SRR002" in wm.items
        pi = wm.items["SRR002"]
        assert pi.item_id == "SRR002"
        assert pi.metadata["sra_url"] == "https://sra.example.com/SRR002"
        assert pi.metadata["dest_path"] == str(dest)

    def test_add_multiple_samples(self, tmp_path: Path) -> None:
        """Multiple samples can be added and tracked."""
        config_file = tmp_path / "wf.yaml"
        config_file.write_text("species: test\n")

        wm = WorkflowManager(config_path=config_file)
        for i in range(5):
            sid = f"SRR{i:03d}"
            wm.add_sample(sid, f"https://sra.example.com/{sid}", tmp_path / f"{sid}.sra")

        assert len(wm.samples) == 5
        assert len(wm.items) == 5


class TestWorkflowManagerStageColors:
    """Tests for WorkflowManager.STAGE_COLORS backward compatibility."""

    def test_stage_colors_mapping_exists(self, tmp_path: Path) -> None:
        """STAGE_COLORS maps all SampleStage values."""
        config_file = tmp_path / "wf.yaml"
        config_file.write_text("species: test\n")

        wm = WorkflowManager(config_path=config_file)

        for stage in SampleStage:
            assert stage in wm.STAGE_COLORS, f"Missing color for {stage}"

    def test_stage_colors_are_strings(self, tmp_path: Path) -> None:
        """Every value in STAGE_COLORS is a string (ANSI code)."""
        config_file = tmp_path / "wf.yaml"
        config_file.write_text("species: test\n")

        wm = WorkflowManager(config_path=config_file)

        for stage, color in wm.STAGE_COLORS.items():
            assert isinstance(color, str), f"Color for {stage} is not a string: {type(color)}"


# ---------------------------------------------------------------------------
# Integration-style: multi-phase pipeline with realistic workflow
# ---------------------------------------------------------------------------


class TestMultiPhasePipeline:
    """Integration-style tests simulating a realistic multi-phase pipeline."""

    def test_three_phase_pipeline_all_succeed(self) -> None:
        """Simulate download -> extract -> quantify with all items succeeding."""
        log: List[str] = []

        def download(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
            for item in items:
                manager.mark_running(item, status="Downloading")
                item.metadata["downloaded"] = True
                # Transition to a custom intermediate stage via metadata.
                manager.mark_done(item, status="Downloaded")
            log.append("download_done")

        def extract(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
            for item in items:
                item.stage = Stage.RUNNING
                item.metadata["extracted"] = True
                item.stage = Stage.DONE
            log.append("extract_done")

        def quantify(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
            for item in items:
                item.stage = Stage.RUNNING
                item.metadata["quantified"] = True
                item.stage = Stage.DONE
            log.append("quantify_done")

        phases = [
            PipelinePhase(
                name="Download",
                handler=download,
                filter_fn=lambda i: i.stage == Stage.PENDING,
            ),
            PipelinePhase(
                name="Extract",
                handler=extract,
                filter_fn=lambda i: i.stage == Stage.DONE
                and i.metadata.get("downloaded")
                and not i.metadata.get("extracted"),
            ),
            PipelinePhase(
                name="Quantify",
                handler=quantify,
                filter_fn=lambda i: i.stage == Stage.DONE
                and i.metadata.get("extracted")
                and not i.metadata.get("quantified"),
            ),
        ]

        manager = BasePipelineManager(phases=phases, config={"threads": 2})
        manager.add_item("S1")
        manager.add_item("S2")

        results = manager.run()

        assert results == {"S1": True, "S2": True}
        assert log == ["download_done", "extract_done", "quantify_done"]
        for item in manager.items.values():
            assert item.metadata.get("downloaded") is True
            assert item.metadata.get("extracted") is True
            assert item.metadata.get("quantified") is True

    def test_failure_in_middle_phase_stops_downstream(self) -> None:
        """If items fail in phase 2, phase 3 filter excludes them."""

        def phase_1(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
            for item in items:
                manager.mark_done(item)
                item.metadata["p1"] = True

        def phase_2(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
            for item in items:
                # Fail every other item.
                if int(item.item_id[-1]) % 2 == 0:
                    manager.mark_failed(item, error="corrupted")
                else:
                    manager.mark_done(item)
                    item.metadata["p2"] = True

        def phase_3(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
            for item in items:
                manager.mark_done(item)
                item.metadata["p3"] = True

        phases = [
            PipelinePhase(name="P1", handler=phase_1, filter_fn=lambda i: i.stage == Stage.PENDING),
            PipelinePhase(
                name="P2", handler=phase_2, filter_fn=lambda i: i.stage == Stage.DONE and i.metadata.get("p1")
            ),
            PipelinePhase(
                name="P3",
                handler=phase_3,
                filter_fn=lambda i: i.stage == Stage.DONE and i.metadata.get("p2"),
            ),
        ]

        manager = BasePipelineManager(phases=phases)
        manager.add_item("item0")  # Will fail in P2 (even)
        manager.add_item("item1")  # Will succeed
        manager.add_item("item2")  # Will fail in P2 (even)
        manager.add_item("item3")  # Will succeed

        results = manager.run()

        assert results["item0"] is False
        assert results["item1"] is True
        assert results["item2"] is False
        assert results["item3"] is True

        # Failed items should NOT have p3 metadata.
        assert "p3" not in manager.items["item0"].metadata
        assert "p3" not in manager.items["item2"].metadata
        # Successful items should have p3.
        assert manager.items["item1"].metadata.get("p3") is True
        assert manager.items["item3"].metadata.get("p3") is True

    def test_phase_sees_only_filtered_items(self) -> None:
        """Phase handler receives only items that pass the filter."""
        received_ids: List[List[str]] = []

        def spy_handler(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
            received_ids.append([i.item_id for i in items])
            for item in items:
                manager.mark_done(item)

        phase = PipelinePhase(
            name="Spy",
            handler=spy_handler,
            filter_fn=lambda i: i.metadata.get("eligible", False),
        )

        manager = BasePipelineManager(phases=[phase])
        manager.add_item("yes1", metadata={"eligible": True})
        manager.add_item("no1", metadata={"eligible": False})
        manager.add_item("yes2", metadata={"eligible": True})

        manager.run()

        assert len(received_ids) == 1
        assert sorted(received_ids[0]) == ["yes1", "yes2"]
