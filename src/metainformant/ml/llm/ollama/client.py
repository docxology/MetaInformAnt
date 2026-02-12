"""Ollama client for local LLM inference.

This module provides the OllamaClient class for interacting with
a local Ollama server for text generation and chat completions.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from typing import Any, Callable, Iterator, Optional

from metainformant.core.utils import logging

from .config import OllamaConfig

logger = logging.get_logger(__name__)


@dataclass
class GenerateResponse:
    """Response from a generate request.

    Attributes:
        text: The generated text content.
        model: Model used for generation.
        done: Whether generation is complete.
        total_duration: Total time in nanoseconds.
        prompt_eval_count: Number of tokens in prompt.
        eval_count: Number of tokens generated.
        context: Context tokens for continuation.
    """

    text: str
    model: str
    done: bool = True
    total_duration: int = 0
    prompt_eval_count: int = 0
    eval_count: int = 0
    context: Optional[list[int]] = None

    @property
    def duration_seconds(self) -> float:
        """Get duration in seconds."""
        return self.total_duration / 1e9

    @property
    def tokens_per_second(self) -> float:
        """Calculate tokens per second."""
        if self.duration_seconds > 0 and self.eval_count > 0:
            return self.eval_count / self.duration_seconds
        return 0.0

    @classmethod
    def from_dict(cls, data: dict) -> "GenerateResponse":
        """Create response from API dict."""
        return cls(
            text=data.get("response", ""),
            model=data.get("model", ""),
            done=data.get("done", True),
            total_duration=data.get("total_duration", 0),
            prompt_eval_count=data.get("prompt_eval_count", 0),
            eval_count=data.get("eval_count", 0),
            context=data.get("context"),
        )


@dataclass
class ChatMessage:
    """A chat message with role and content.

    Attributes:
        role: Message role - 'system', 'user', or 'assistant'.
        content: Message text content.
    """

    role: str
    content: str

    def to_dict(self) -> dict:
        """Convert to API format."""
        return {"role": self.role, "content": self.content}


@dataclass
class ChatResponse:
    """Response from a chat request.

    Attributes:
        message: The assistant's response message.
        model: Model used for generation.
        done: Whether generation is complete.
        total_duration: Total time in nanoseconds.
        prompt_eval_count: Number of tokens in prompt.
        eval_count: Number of tokens generated.
    """

    message: ChatMessage
    model: str
    done: bool = True
    total_duration: int = 0
    prompt_eval_count: int = 0
    eval_count: int = 0

    @property
    def text(self) -> str:
        """Get the response text content."""
        return self.message.content

    @property
    def duration_seconds(self) -> float:
        """Get duration in seconds."""
        return self.total_duration / 1e9

    @classmethod
    def from_dict(cls, data: dict) -> "ChatResponse":
        """Create response from API dict."""
        msg_data = data.get("message", {})
        return cls(
            message=ChatMessage(
                role=msg_data.get("role", "assistant"),
                content=msg_data.get("content", ""),
            ),
            model=data.get("model", ""),
            done=data.get("done", True),
            total_duration=data.get("total_duration", 0),
            prompt_eval_count=data.get("prompt_eval_count", 0),
            eval_count=data.get("eval_count", 0),
        )


@dataclass
class ModelInfo:
    """Information about an available model.

    Attributes:
        name: Model name/tag.
        size: Model size in bytes.
        modified_at: Last modification timestamp.
        digest: Model digest hash.
    """

    name: str
    size: int = 0
    modified_at: str = ""
    digest: str = ""

    @property
    def size_gb(self) -> float:
        """Get size in gigabytes."""
        return self.size / (1024**3)

    @classmethod
    def from_dict(cls, data: dict) -> "ModelInfo":
        """Create from API dict."""
        return cls(
            name=data.get("name", ""),
            size=data.get("size", 0),
            modified_at=data.get("modified_at", ""),
            digest=data.get("digest", ""),
        )


class OllamaClient:
    """Client for interacting with Ollama server.

    Provides methods for text generation, chat, and model management.
    Uses the Ollama REST API with pure Python urllib (no external deps).

    Example:
        >>> client = OllamaClient()
        >>> response = client.generate("Hello, world!")
        >>> print(response.text)

        >>> # Chat conversation
        >>> messages = [ChatMessage("user", "What is DNA?")]
        >>> response = client.chat(messages)
        >>> print(response.text)
    """

    def __init__(self, config: Optional[OllamaConfig] = None) -> None:
        """Initialize the Ollama client.

        Args:
            config: Optional configuration. Uses defaults if not provided.
        """
        self.config = config or OllamaConfig()
        self._conversation_history: list[ChatMessage] = []

    def _request(
        self,
        url: str,
        data: Optional[dict] = None,
        method: str = "POST",
        stream: bool = False,
    ) -> Any:
        """Make an HTTP request to the Ollama API.

        Args:
            url: API endpoint URL.
            data: Request body data.
            method: HTTP method.
            stream: Whether to stream the response.

        Returns:
            Parsed JSON response or iterator for streaming.

        Raises:
            ConnectionError: If server is unreachable.
            RuntimeError: If request fails.
        """
        headers = {"Content-Type": "application/json"}

        body = json.dumps(data).encode("utf-8") if data else None
        request = urllib.request.Request(url, data=body, headers=headers, method=method)

        for attempt in range(self.config.max_retries):
            try:
                response = urllib.request.urlopen(request, timeout=self.config.timeout)

                if stream:
                    return self._stream_response(response)

                return json.loads(response.read().decode("utf-8"))

            except urllib.error.URLError as e:
                if attempt < self.config.max_retries - 1:
                    wait_time = 2**attempt
                    logger.warning(f"Request failed, retrying in {wait_time}s: {e}")
                    time.sleep(wait_time)
                else:
                    raise ConnectionError(f"Failed to connect to Ollama at {self.config.host}: {e}") from e
            except urllib.error.HTTPError as e:
                raise RuntimeError(f"Ollama API error {e.code}: {e.read().decode('utf-8')}") from e

    def _stream_response(self, response) -> Iterator[dict]:
        """Stream NDJSON response from Ollama.

        Args:
            response: HTTP response object.

        Yields:
            Parsed JSON objects from each line.
        """
        for line in response:
            if line.strip():
                yield json.loads(line.decode("utf-8"))

    def is_available(self) -> bool:
        """Check if Ollama server is available.

        Returns:
            True if server responds, False otherwise.
        """
        try:
            self.list_models()
            return True
        except Exception:
            return False

    def list_models(self) -> list[ModelInfo]:
        """List available models on the server.

        Returns:
            List of ModelInfo objects.

        Example:
            >>> client = OllamaClient()
            >>> models = client.list_models()
            >>> for m in models:
            ...     print(f"{m.name}: {m.size_gb:.1f}GB")
        """
        request = urllib.request.Request(self.config.tags_url, method="GET")

        try:
            response = urllib.request.urlopen(request, timeout=self.config.timeout)
            data = json.loads(response.read().decode("utf-8"))
            return [ModelInfo.from_dict(m) for m in data.get("models", [])]
        except Exception as e:
            raise ConnectionError(f"Failed to list models: {e}") from e

    def generate(
        self,
        prompt: str,
        model: Optional[str] = None,
        system: Optional[str] = None,
        context: Optional[list[int]] = None,
        stream: bool = False,
        stream_callback: Optional[Callable[[str], None]] = None,
        **kwargs,
    ) -> GenerateResponse:
        """Generate text completion.

        Args:
            prompt: The prompt text.
            model: Model to use (overrides config).
            system: System prompt (overrides config).
            context: Previous context for continuation.
            stream: Whether to stream the response.
            stream_callback: Function called with each text chunk when streaming.
            **kwargs: Additional Ollama options.

        Returns:
            GenerateResponse with generated text.

        Example:
            >>> response = client.generate(
            ...     "Explain photosynthesis in one sentence.",
            ...     temperature=0.3
            ... )
            >>> print(response.text)
        """
        data = {
            "model": model or self.config.model,
            "prompt": prompt,
            "stream": stream or (stream_callback is not None),
            "options": {**self.config.to_options_dict(), **kwargs},
        }

        if system or self.config.system_prompt:
            data["system"] = system or self.config.system_prompt

        if context:
            data["context"] = context

        if self.config.stop_sequences:
            data["stop"] = self.config.stop_sequences

        if data["stream"]:
            return self._generate_stream(data, stream_callback)

        response_data = self._request(self.config.generate_url, data)
        return GenerateResponse.from_dict(response_data)

    def _generate_stream(
        self,
        data: dict,
        callback: Optional[Callable[[str], None]] = None,
    ) -> GenerateResponse:
        """Handle streaming generation.

        Args:
            data: Request data.
            callback: Function to call with each chunk.

        Returns:
            Complete GenerateResponse.
        """
        chunks: list[str] = []
        final_response: dict = {}

        for chunk in self._request(self.config.generate_url, data, stream=True):
            text = chunk.get("response", "")
            chunks.append(text)

            if callback:
                callback(text)

            if chunk.get("done"):
                final_response = chunk
                break

        return GenerateResponse(
            text="".join(chunks),
            model=data["model"],
            done=True,
            total_duration=final_response.get("total_duration", 0),
            prompt_eval_count=final_response.get("prompt_eval_count", 0),
            eval_count=final_response.get("eval_count", 0),
            context=final_response.get("context"),
        )

    def chat(
        self,
        messages: list[ChatMessage],
        model: Optional[str] = None,
        stream: bool = False,
        stream_callback: Optional[Callable[[str], None]] = None,
        **kwargs,
    ) -> ChatResponse:
        """Conduct a chat conversation.

        Args:
            messages: List of ChatMessage objects.
            model: Model to use (overrides config).
            stream: Whether to stream the response.
            stream_callback: Function called with each text chunk.
            **kwargs: Additional Ollama options.

        Returns:
            ChatResponse with assistant message.

        Example:
            >>> messages = [
            ...     ChatMessage("system", "You are a biology expert."),
            ...     ChatMessage("user", "What is mitosis?"),
            ... ]
            >>> response = client.chat(messages)
            >>> print(response.text)
        """
        data = {
            "model": model or self.config.model,
            "messages": [m.to_dict() for m in messages],
            "stream": stream or (stream_callback is not None),
            "options": {**self.config.to_options_dict(), **kwargs},
        }

        if data["stream"]:
            return self._chat_stream(data, stream_callback)

        response_data = self._request(self.config.chat_url, data)
        return ChatResponse.from_dict(response_data)

    def _chat_stream(
        self,
        data: dict,
        callback: Optional[Callable[[str], None]] = None,
    ) -> ChatResponse:
        """Handle streaming chat.

        Args:
            data: Request data.
            callback: Function to call with each chunk.

        Returns:
            Complete ChatResponse.
        """
        chunks: list[str] = []
        final_response: dict = {}

        for chunk in self._request(self.config.chat_url, data, stream=True):
            msg = chunk.get("message", {})
            text = msg.get("content", "")
            chunks.append(text)

            if callback:
                callback(text)

            if chunk.get("done"):
                final_response = chunk
                break

        return ChatResponse(
            message=ChatMessage("assistant", "".join(chunks)),
            model=data["model"],
            done=True,
            total_duration=final_response.get("total_duration", 0),
            prompt_eval_count=final_response.get("prompt_eval_count", 0),
            eval_count=final_response.get("eval_count", 0),
        )

    def conversation(
        self,
        message: str,
        system: Optional[str] = None,
        reset: bool = False,
    ) -> ChatResponse:
        """Continue or start a conversation with history.

        Maintains conversation history internally for multi-turn chat.

        Args:
            message: User message.
            system: System prompt (only used at start).
            reset: Whether to reset conversation history.

        Returns:
            ChatResponse with assistant reply.

        Example:
            >>> client.conversation("What is RNA?")
            >>> client.conversation("How is it different from DNA?")
        """
        if reset:
            self._conversation_history.clear()

        # Add system message at start if provided
        if system and not self._conversation_history:
            self._conversation_history.append(ChatMessage("system", system))

        # Add user message
        self._conversation_history.append(ChatMessage("user", message))

        # Get response
        response = self.chat(self._conversation_history)

        # Add assistant response to history
        self._conversation_history.append(response.message)

        return response

    def reset_conversation(self) -> None:
        """Clear conversation history."""
        self._conversation_history.clear()

    @property
    def conversation_history(self) -> list[ChatMessage]:
        """Get current conversation history."""
        return self._conversation_history.copy()
