"""Pipeline Manager for orchestrating multi-stage bioinformatics workflows with TUI.

This module provides:
- ``BasePipelineManager``: A generic, domain-agnostic pipeline manager that
  coordinates sequential phases of work items with real-time TUI visualization.
- ``WorkflowManager``: An RNA-seq-specific subclass (backward compatible) that
  wires up download, getfastq, and quant phases using amalgkit.

The generic pipeline is built around three primitives:
- ``Stage``: Minimal lifecycle enum (PENDING, RUNNING, DONE, FAILED).
- ``PipelineItem``: A tracked work item flowing through the pipeline.
- ``PipelinePhase``: A named phase with a handler callable and item filter.

Backward compatibility
----------------------
``SampleStage``, ``SampleState``, and ``WorkflowManager`` retain their original
public API.  Existing callers (e.g. ``run_workflow_tui.py``) require zero changes.
"""

from __future__ import annotations

import logging
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence

from metainformant.core.ui.tui import (
    BLUE,
    CYAN,
    GREEN,
    MAGENTA,
    RED,
    YELLOW,
    TerminalInterface,
)
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Generic pipeline primitives
# ---------------------------------------------------------------------------


class Stage(Enum):
    """Minimal lifecycle stages for any pipeline item."""

    PENDING = "Pending"
    RUNNING = "Running"
    DONE = "Done"
    FAILED = "Failed"


@dataclass
class PipelineItem:
    """A single work item flowing through the pipeline.

    Parameters
    ----------
    item_id:
        Unique identifier (e.g. SRA accession, sample barcode).
    metadata:
        Arbitrary key/value bag carried alongside the item.
    stage:
        Current lifecycle stage.
    error:
        Human-readable error message when ``stage`` is ``FAILED``.
    """

    item_id: str
    metadata: Dict[str, Any] = field(default_factory=dict)
    stage: Stage = Stage.PENDING
    error: str = ""


@dataclass
class PipelinePhase:
    """Describes a single phase in a multi-phase pipeline.

    Parameters
    ----------
    name:
        Short display name (shown in TUI footer, e.g. ``"Download"``).
    handler:
        Callable that executes the phase.  Signature::

            handler(manager: BasePipelineManager, items: List[PipelineItem]) -> None

        The handler is responsible for mutating each item's ``stage`` and
        calling ``manager.ui`` for progress updates.
    filter_fn:
        Predicate selecting which items enter this phase.  Receives a
        ``PipelineItem`` and returns ``True`` if the item should be processed.
        Defaults to selecting all ``Stage.PENDING`` items.
    color:
        ANSI color code used in the TUI for this phase.
    """

    name: str
    handler: Callable[[BasePipelineManager, List[PipelineItem]], None]
    filter_fn: Callable[[PipelineItem], bool] = lambda item: item.stage == Stage.PENDING
    color: str = CYAN


# ---------------------------------------------------------------------------
# BasePipelineManager -- generic, domain-agnostic
# ---------------------------------------------------------------------------


class BasePipelineManager:
    """Generic multi-phase pipeline manager with TUI visualization.

    Subclass this and supply ``phases`` to build domain-specific pipelines
    (RNA-seq, GWAS, proteomics, etc.).

    Parameters
    ----------
    phases:
        Ordered sequence of ``PipelinePhase`` objects executed left-to-right.
    max_threads:
        Thread pool size for phases that use ``self.executor``.
    config:
        Optional configuration dict available to phase handlers.
    """

    def __init__(
        self,
        phases: Sequence[PipelinePhase],
        max_threads: int = 5,
        config: Optional[Dict[str, Any]] = None,
    ) -> None:
        self.phases: List[PipelinePhase] = list(phases)
        self.max_threads = max_threads
        self.config: Dict[str, Any] = config or {}
        self.ui = TerminalInterface()
        self.items: Dict[str, PipelineItem] = {}
        self.executor = ThreadPoolExecutor(max_workers=max_threads)
        self._running = False

    # -- Public API ---------------------------------------------------------

    def add_item(self, item_id: str, metadata: Optional[Dict[str, Any]] = None) -> PipelineItem:
        """Register a new work item in the pipeline.

        Parameters
        ----------
        item_id:
            Unique identifier for the item.
        metadata:
            Optional key/value bag associated with this item.

        Returns
        -------
        The newly created ``PipelineItem``.
        """
        item = PipelineItem(item_id=item_id, metadata=metadata or {})
        self.items[item_id] = item
        self.ui.add_bar(item_id, item_id[:20], total=0.0, unit="MB")
        self.ui.update(item_id, stage="Pending", status="Queued", color=BLUE)
        return item

    def run(self) -> Dict[str, bool]:
        """Execute all phases sequentially.

        Returns
        -------
        Dict mapping ``item_id`` to ``True`` (done) / ``False`` (failed).
        """
        self._running = True
        self.ui.start()

        try:
            for phase in self.phases:
                eligible = [item for item in self.items.values() if phase.filter_fn(item)]
                if not eligible:
                    continue
                self.ui.set_footer(f"Phase: {phase.name} | Processing {len(eligible)} items...")
                phase.handler(self, eligible)
        finally:
            self.ui.stop()
            self.executor.shutdown(wait=True)
            self._running = False

        return {item_id: item.stage == Stage.DONE for item_id, item in self.items.items()}

    # -- Convenience helpers for phase handlers -----------------------------

    def mark_running(self, item: PipelineItem, *, status: str = "Running...", color: str = CYAN) -> None:
        """Transition an item to RUNNING and update the TUI."""
        item.stage = Stage.RUNNING
        self.ui.update(item.item_id, status=status, color=color)

    def mark_done(self, item: PipelineItem, *, status: str = "Done", color: str = GREEN) -> None:
        """Transition an item to DONE and update the TUI."""
        item.stage = Stage.DONE
        self.ui.update(item.item_id, status=status, color=color)

    def mark_failed(self, item: PipelineItem, error: str, *, status: str = "Failed", color: str = RED) -> None:
        """Transition an item to FAILED, record the error, and update the TUI."""
        item.stage = Stage.FAILED
        item.error = error
        self.ui.update(item.item_id, status=status, color=color)


# ---------------------------------------------------------------------------
# RNA-seq backward-compatible types
# ---------------------------------------------------------------------------


class SampleStage(Enum):
    """Stages a sample goes through in the RNA-seq workflow.

    Kept for backward compatibility -- existing callers import this directly.
    """

    PENDING = "Pending"
    DOWNLOADING = "Downloading"
    DOWNLOADED = "Downloaded"
    EXTRACTING = "Extracting"
    EXTRACTED = "Extracted"
    QUANTIFYING = "Quantifying"
    DONE = "Done"
    FAILED = "Failed"
    TRUNCATED = "Truncated"


@dataclass
class SampleState:
    """Tracks a sample's state through the RNA-seq workflow.

    Kept for backward compatibility -- existing callers import this directly.
    """

    sample_id: str
    sra_url: str
    dest_path: Path
    stage: SampleStage = SampleStage.PENDING
    current_bytes: int = 0
    total_bytes: int = 0
    error: str = ""


# ---------------------------------------------------------------------------
# RNA-seq phase handlers (free functions used as PipelinePhase.handler)
# ---------------------------------------------------------------------------

# Color map shared by the RNA-seq WorkflowManager and its phase handlers.
_RNA_STAGE_COLORS: Dict[SampleStage, str] = {
    SampleStage.PENDING: BLUE,
    SampleStage.DOWNLOADING: CYAN,
    SampleStage.DOWNLOADED: YELLOW,
    SampleStage.EXTRACTING: MAGENTA,
    SampleStage.EXTRACTED: YELLOW,
    SampleStage.QUANTIFYING: CYAN,
    SampleStage.DONE: GREEN,
    SampleStage.FAILED: RED,
    SampleStage.TRUNCATED: RED,
}


def _download_single_sample(
    manager: WorkflowManager,
    sample_id: str,
    url: str,
    dest: Path,
    state: SampleState,
) -> bool:
    """Download one SRA file, updating the TUI as we go."""
    from metainformant.core.io.download_robust import get_remote_file_size, robust_download_url

    tid = threading.get_native_id()
    manager.ui.update(sample_id, status=f"Init (TID:{tid})", color=CYAN)

    total_bytes = get_remote_file_size(url)
    total_mb = total_bytes / 1024 / 1024 if total_bytes else 0
    state.total_bytes = total_bytes

    # Cache hit -- file already present with correct size.
    if dest.exists():
        existing_size = dest.stat().st_size
        existing_mb = existing_size / 1024 / 1024
        if total_bytes > 0 and existing_size >= (total_bytes * 0.99):
            manager.ui.update(sample_id, current=existing_mb, total=existing_mb, status="Cached", color=GREEN)
            return True
        elif total_bytes == 0:
            manager.ui.update(sample_id, current=existing_mb, total=existing_mb, status="Cached", color=GREEN)
            return True

    if total_mb > 0:
        manager.ui.update(sample_id, total=total_mb, status=f"Downloading (TID:{tid})")
    else:
        manager.ui.update(sample_id, status=f"Downloading (TID:{tid})")

    success = robust_download_url(url, dest)

    if success and dest.exists():
        final_size = dest.stat().st_size
        final_mb = final_size / 1024 / 1024

        if total_bytes > 0 and final_size < (total_bytes * 0.99):
            manager.ui.update(sample_id, status="Truncated", color=RED)
            state.stage = SampleStage.TRUNCATED
            return False

        manager.ui.update(sample_id, current=final_mb, total=max(total_mb, final_mb))
        return True

    manager.ui.update(sample_id, status="Failed", color=RED)
    return False


def _rna_download_phase(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
    """Phase handler: download all SRA files in parallel."""
    # ``manager`` is actually a WorkflowManager here (the only user of this handler).
    assert isinstance(manager, WorkflowManager)

    futures: Dict[Any, str] = {}

    for item in items:
        sid = item.item_id
        state = manager.samples[sid]
        state.stage = SampleStage.DOWNLOADING
        manager.ui.update(sid, stage="Download", status="Starting...", color=CYAN)

        future = manager.executor.submit(
            _download_single_sample,
            manager,
            sid,
            state.sra_url,
            state.dest_path,
            state,
        )
        futures[future] = sid

    # Poll file sizes while downloads are in flight.
    while any(not f.done() for f in futures):
        _update_download_progress(manager)
        time.sleep(0.5)

    # Collect results.
    for f in as_completed(futures):
        sid = futures[f]
        try:
            success = f.result()
            if success:
                manager.samples[sid].stage = SampleStage.DOWNLOADED
                manager.ui.update(sid, stage="Download", status="Downloaded", color=YELLOW)
            else:
                manager.samples[sid].stage = SampleStage.FAILED
        except Exception as exc:
            manager.samples[sid].stage = SampleStage.FAILED
            manager.samples[sid].error = str(exc)
            logger.error("Download failed for %s: %s", sid, exc)


def _update_download_progress(manager: WorkflowManager) -> None:
    """Poll file sizes during download phase and push updates to TUI."""
    active = 0
    done = 0
    for sid, state in manager.samples.items():
        if state.stage == SampleStage.DOWNLOADING:
            dest = state.dest_path
            temp = dest.with_suffix(dest.suffix + ".part")
            current_size = 0
            if dest.exists():
                current_size = dest.stat().st_size
                done += 1
            elif temp.exists():
                current_size = temp.stat().st_size
                active += 1
            manager.ui.update(sid, current=current_size / 1024 / 1024, speed=f"{current_size / 1024 / 1024:.1f} MB")
        elif state.stage == SampleStage.DOWNLOADED:
            done += 1
    manager.ui.set_footer(f"Phase: Download | Active: {active} | Done: {done} | Total: {len(manager.samples)}")


def _rna_getfastq_phase(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
    """Phase handler: run amalgkit getfastq on downloaded samples."""
    assert isinstance(manager, WorkflowManager)
    from metainformant.rna.amalgkit.amalgkit import run_amalgkit

    downloaded = [manager.samples[item.item_id] for item in items]

    for state in downloaded:
        state.stage = SampleStage.EXTRACTING
        manager.ui.update(
            state.sample_id, stage="Extract", status="Extracting FASTQ...", color=MAGENTA, current=0.0, total=100.0
        )

    manager.ui.set_footer(f"Phase: getfastq | Processing {len(downloaded)} samples...")

    try:
        metadata_path = manager.work_dir / "metadata" / "metadata_selected.tsv"
        params: Dict[str, Any] = {
            "out_dir": str(manager.work_dir),
            "metadata": str(metadata_path),
            "threads": manager.config.get("threads", 8),
        }
        result = run_amalgkit("getfastq", params)

        if result.returncode == 0:
            for state in downloaded:
                state.stage = SampleStage.EXTRACTED
                manager.ui.update(
                    state.sample_id,
                    current=100.0,
                    total=100.0,
                    stage="Extract",
                    status="Extracted",
                    color=YELLOW,
                )
        else:
            for state in downloaded:
                state.stage = SampleStage.FAILED
                state.error = "getfastq failed"
                manager.ui.update(state.sample_id, status="Extract Failed", color=RED)
    except Exception as exc:
        logger.error("getfastq phase failed: %s", exc)
        for state in downloaded:
            state.stage = SampleStage.FAILED
            state.error = str(exc)
            manager.ui.update(state.sample_id, status="Error", color=RED)


def _rna_quant_phase(manager: BasePipelineManager, items: List[PipelineItem]) -> None:
    """Phase handler: run amalgkit quant on extracted samples."""
    assert isinstance(manager, WorkflowManager)
    from metainformant.rna.amalgkit.amalgkit import run_amalgkit

    extracted = [manager.samples[item.item_id] for item in items]

    for state in extracted:
        state.stage = SampleStage.QUANTIFYING
        manager.ui.update(state.sample_id, stage="Quant", status="Quantifying...", color=CYAN, current=0.0, total=100.0)

    manager.ui.set_footer(f"Phase: quant | Processing {len(extracted)} samples...")

    try:
        metadata_path = manager.work_dir / "metadata" / "metadata_selected.tsv"
        params: Dict[str, Any] = {
            "out_dir": str(manager.work_dir),
            "metadata": str(metadata_path),
            "threads": manager.config.get("threads", 8),
        }
        result = run_amalgkit("quant", params)

        if result.returncode == 0:
            for state in extracted:
                state.stage = SampleStage.DONE
                manager.ui.update(
                    state.sample_id,
                    current=100.0,
                    total=100.0,
                    stage="Done",
                    status="Quantified",
                    color=GREEN,
                )
        else:
            for state in extracted:
                state.stage = SampleStage.FAILED
                state.error = "quant failed"
                manager.ui.update(state.sample_id, status="Quant Failed", color=RED)
    except Exception as exc:
        logger.error("quant phase failed: %s", exc)
        for state in extracted:
            state.stage = SampleStage.FAILED
            state.error = str(exc)
            manager.ui.update(state.sample_id, status="Error", color=RED)

    # Final summary.
    done = sum(1 for s in manager.samples.values() if s.stage == SampleStage.DONE)
    failed = sum(1 for s in manager.samples.values() if s.stage in (SampleStage.FAILED, SampleStage.TRUNCATED))
    manager.ui.set_footer(f"Complete | Done: {done} | Failed: {failed} | Total: {len(manager.samples)}")


# ---------------------------------------------------------------------------
# WorkflowManager -- RNA-seq specific, backward compatible
# ---------------------------------------------------------------------------


class WorkflowManager(BasePipelineManager):
    """RNA-seq workflow manager (download -> getfastq -> quant) with TUI.

    This is the original public API preserved for backward compatibility.
    Internally it delegates to :class:`BasePipelineManager` with RNA-specific
    :class:`PipelinePhase` definitions.

    Parameters
    ----------
    config_path:
        Path to a YAML config containing at least ``work_dir`` and ``species``.
    max_threads:
        Thread pool size for the download phase.
    """

    STAGE_COLORS = _RNA_STAGE_COLORS

    def __init__(self, config_path: Path, max_threads: int = 5) -> None:
        self.config_path = Path(config_path)
        self.samples: Dict[str, SampleState] = {}

        # Load YAML config first so we can feed it to the base class.
        self._load_config()

        # Define the three RNA-seq phases with their filter predicates.
        phases = [
            PipelinePhase(
                name="Download",
                handler=_rna_download_phase,
                filter_fn=lambda item: item.item_id in self.samples
                and self.samples[item.item_id].stage == SampleStage.PENDING,
                color=CYAN,
            ),
            PipelinePhase(
                name="getfastq",
                handler=_rna_getfastq_phase,
                filter_fn=lambda item: item.item_id in self.samples
                and self.samples[item.item_id].stage == SampleStage.DOWNLOADED,
                color=MAGENTA,
            ),
            PipelinePhase(
                name="quant",
                handler=_rna_quant_phase,
                filter_fn=lambda item: item.item_id in self.samples
                and self.samples[item.item_id].stage == SampleStage.EXTRACTED,
                color=CYAN,
            ),
        ]

        super().__init__(phases=phases, max_threads=max_threads, config=self.config)

    # -- Config loading -----------------------------------------------------

    def _load_config(self) -> None:
        """Load workflow configuration from YAML."""
        import yaml

        with open(self.config_path) as fh:
            self.config: Dict[str, Any] = yaml.safe_load(fh)  # type: ignore[no-redef]

        self.work_dir = Path(self.config.get("work_dir", "output/amalgkit"))
        self.species: str = self.config.get("species", "unknown")

    # -- Backward-compatible public API -------------------------------------

    def add_sample(self, sample_id: str, sra_url: str, dest_path: Path) -> None:
        """Add a sample to track through the workflow.

        This is the original public method.  It creates both a
        :class:`SampleState` (RNA-specific) and a :class:`PipelineItem`
        (generic) so the base-class ``run()`` loop sees the item.
        """
        self.samples[sample_id] = SampleState(
            sample_id=sample_id,
            sra_url=sra_url,
            dest_path=dest_path,
        )
        # Register in the generic items dict as well so BasePipelineManager.run()
        # iterates over it and phase filter_fn sees it.
        self.add_item(sample_id, metadata={"sra_url": sra_url, "dest_path": str(dest_path)})

    def run(self) -> Dict[str, bool]:
        """Execute the full RNA-seq workflow.

        Returns
        -------
        Dict mapping ``sample_id`` to ``True`` (done) / ``False`` (failed).
        """
        # Silence the download logger during TUI rendering.
        robust_logger = logging.getLogger("metainformant.core.io.download_robust")
        original_level = robust_logger.getEffectiveLevel()
        robust_logger.setLevel(logging.ERROR)

        try:
            return super().run()
        finally:
            robust_logger.setLevel(original_level)
