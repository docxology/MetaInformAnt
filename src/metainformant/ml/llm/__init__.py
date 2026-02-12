"""LLM module for METAINFORMANT.

This module provides local Large Language Model capabilities through
various backends, with Ollama as the primary supported provider.

Subpackages:
    ollama: Ollama backend for local LLM inference"""
from __future__ import annotations

from . import ollama

__all__ = ['ollama']
