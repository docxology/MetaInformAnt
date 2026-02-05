"""Configuration for Ollama client.

This module provides the OllamaConfig dataclass for configuring
connections and inference parameters for the Ollama server.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class OllamaConfig:
    """Configuration for Ollama client.

    Attributes:
        host: Ollama server URL. Defaults to localhost.
        model: Default model name for inference.
        temperature: Sampling temperature (0.0-2.0). Higher = more creative.
        top_p: Nucleus sampling threshold (0.0-1.0).
        top_k: Top-k sampling limit.
        timeout: Request timeout in seconds.
        max_retries: Number of retry attempts on failure.
        stream: Whether to stream responses by default.
        context_length: Maximum context window size.

    Example:
        >>> config = OllamaConfig(
        ...     model="llama3:latest",
        ...     temperature=0.7,
        ...     timeout=60.0
        ... )
    """

    host: str = "http://localhost:11434"
    model: str = "smollm2:135m-instruct-q4_K_S"
    temperature: float = 0.7
    top_p: float = 0.9
    top_k: int = 40
    timeout: float = 120.0
    max_retries: int = 3
    stream: bool = False
    context_length: int = 4096

    # Optional system prompt for all requests
    system_prompt: Optional[str] = None

    # Generation stop sequences
    stop_sequences: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        """Validate configuration values."""
        if not self.host.startswith(("http://", "https://")):
            self.host = f"http://{self.host}"

        if not 0.0 <= self.temperature <= 2.0:
            raise ValueError(f"temperature must be 0.0-2.0, got {self.temperature}")

        if not 0.0 <= self.top_p <= 1.0:
            raise ValueError(f"top_p must be 0.0-1.0, got {self.top_p}")

        if self.top_k < 1:
            raise ValueError(f"top_k must be >= 1, got {self.top_k}")

        if self.timeout <= 0:
            raise ValueError(f"timeout must be > 0, got {self.timeout}")

    @property
    def api_url(self) -> str:
        """Get the base API URL."""
        return f"{self.host}/api"

    @property
    def generate_url(self) -> str:
        """Get the generate endpoint URL."""
        return f"{self.api_url}/generate"

    @property
    def chat_url(self) -> str:
        """Get the chat endpoint URL."""
        return f"{self.api_url}/chat"

    @property
    def tags_url(self) -> str:
        """Get the tags (list models) endpoint URL."""
        return f"{self.api_url}/tags"

    def to_options_dict(self) -> dict:
        """Convert sampling parameters to Ollama options format."""
        return {
            "temperature": self.temperature,
            "top_p": self.top_p,
            "top_k": self.top_k,
            "num_ctx": self.context_length,
        }

    def with_model(self, model: str) -> "OllamaConfig":
        """Return a new config with a different model."""
        return OllamaConfig(
            host=self.host,
            model=model,
            temperature=self.temperature,
            top_p=self.top_p,
            top_k=self.top_k,
            timeout=self.timeout,
            max_retries=self.max_retries,
            stream=self.stream,
            context_length=self.context_length,
            system_prompt=self.system_prompt,
            stop_sequences=self.stop_sequences.copy(),
        )
