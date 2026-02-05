"""Ollama backend for local LLM inference.

This subpackage provides a composable interface to the Ollama server
for local large language model inference. It supports text generation,
chat conversations, streaming, and composable chain patterns.

Example:
    >>> from metainformant.ml.llm.ollama import OllamaClient, OllamaConfig
    >>> config = OllamaConfig(model="smollm2:135m-instruct-q4_K_S")
    >>> client = OllamaClient(config)
    >>> response = client.generate("Explain DNA in one sentence.")
    >>> print(response.text)
"""

from __future__ import annotations

from .chains import (
    Chain,
    ChainResult,
    ConversationChain,
    MapReduceChain,
    PromptChain,
    RouterChain,
    SequentialChain,
    TransformChain,
)
from .client import OllamaClient
from .config import OllamaConfig
from .prompts import (
    ChatMessage,
    PromptTemplate,
    SystemPrompt,
    build_bioinformatics_prompt,
)

__all__ = [
    # Core
    "OllamaClient",
    "OllamaConfig",
    # Prompts
    "PromptTemplate",
    "SystemPrompt",
    "ChatMessage",
    "build_bioinformatics_prompt",
    # Chains
    "Chain",
    "ChainResult",
    "PromptChain",
    "SequentialChain",
    "MapReduceChain",
    "RouterChain",
    "TransformChain",
    "ConversationChain",
]
