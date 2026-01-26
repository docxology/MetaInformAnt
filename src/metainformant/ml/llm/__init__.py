"""LLM module for METAINFORMANT.

This module provides local Large Language Model capabilities through
various backends, with Ollama as the primary supported provider.

Subpackages:
    ollama: Ollama backend for local LLM inference
"""

from __future__ import annotations

# Import subpackages
from . import ollama

# Convenience exports from ollama
from .ollama import (
    OllamaClient,
    OllamaConfig,
    PromptTemplate,
    ChatMessage,
    Chain,
    SequentialChain,
)

__all__ = [
    # Subpackages
    "ollama",
    # Core classes
    "OllamaClient",
    "OllamaConfig",
    # Prompts
    "PromptTemplate",
    "ChatMessage",
    # Chains
    "Chain",
    "SequentialChain",
]
