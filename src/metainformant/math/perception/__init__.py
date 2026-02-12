"""Perception and Psychophysics submodule.

This module provides tools for modeling sensory perception and decision making,
including Psychophysics (Weber-Fechner, Stevens) and Signal Detection Theory."""
from __future__ import annotations

from . import psychophysics, signal_detection

__all__ = ['psychophysics', 'signal_detection']
