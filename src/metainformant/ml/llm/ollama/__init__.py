"""Ollama backend for local LLM inference.

This subpackage provides a composable interface to the Ollama server
for local large language model inference. It supports text generation,
chat conversations, streaming, and composable chain patterns.

Example:
    >>> from metainformant.ml.llm.ollama import OllamaClient, OllamaConfig
    >>> config = OllamaConfig(model="smollm2:135m-instruct-q4_K_S")
    >>> client = OllamaClient(config)
    >>> response = client.generate("Explain DNA in one sentence.")
    >>> print(response.text)"""
from __future__ import annotations

from . import chains, client, config, prompts

__all__ = ['chains', 'client', 'config', 'prompts']
