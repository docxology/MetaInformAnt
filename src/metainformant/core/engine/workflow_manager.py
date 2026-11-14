"""Pipeline manager for orchestrating multi-phase workflows with a TUI.

``BasePipelineManager`` is a generic, domain-agnostic pipeline manager that
coordinates sequential phases of work items with real-time TUI visualization.
It is built around three primitives:

- ``Stage``: Minimal lifecycle enum (PENDING, RUNNING, DONE, FAILED).
- ``PipelineItem``: A tracked work item flowing through the pipeline.
- ``PipelinePhase``: A named phase with a handler callable and item filter.

Domain packages subclass ``BasePipelineManager`` and supply their own phases;
this module imports no domain code.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable, Dict, List, Optional, Sequence

from metainformant.core.ui.tui import BLUE, CYAN, GREEN, RED, TerminalInterface

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
