"""Composable chain patterns for LLM workflows.

This module provides chain abstractions for building complex
LLM workflows by composing simple operations.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Callable, Optional, TypeVar, Generic

from .client import OllamaClient, ChatMessage
from .prompts import PromptTemplate

logger = logging.getLogger(__name__)

T = TypeVar("T")
R = TypeVar("R")


@dataclass
class ChainResult(Generic[T]):
    """Result from a chain execution.
    
    Attributes:
        value: The output value.
        success: Whether execution succeeded.
        error: Error message if failed.
        metadata: Additional execution metadata.
    """
    
    value: Optional[T]
    success: bool = True
    error: Optional[str] = None
    metadata: dict[str, Any] = field(default_factory=dict)
    
    @classmethod
    def ok(cls, value: T, **metadata) -> "ChainResult[T]":
        """Create a successful result."""
        return cls(value=value, success=True, metadata=metadata)
    
    @classmethod
    def fail(cls, error: str, **metadata) -> "ChainResult[T]":
        """Create a failed result."""
        return cls(value=None, success=False, error=error, metadata=metadata)


class Chain(ABC, Generic[T, R]):
    """Base class for composable chains.
    
    A chain takes an input and produces an output, potentially
    using an LLM for processing.
    """
    
    def __init__(self, client: Optional[OllamaClient] = None) -> None:
        """Initialize the chain.
        
        Args:
            client: Ollama client to use. Creates default if not provided.
        """
        self.client = client or OllamaClient()
    
    @abstractmethod
    def run(self, input_data: T) -> ChainResult[R]:
        """Execute the chain.
        
        Args:
            input_data: Input to process.
            
        Returns:
            ChainResult with output or error.
        """
        pass
    
    def __call__(self, input_data: T) -> ChainResult[R]:
        """Execute the chain (callable interface)."""
        return self.run(input_data)
    
    def then(self, next_chain: "Chain[R, Any]") -> "SequentialChain":
        """Chain this with another chain.
        
        Args:
            next_chain: Chain to run after this one.
            
        Returns:
            SequentialChain combining both.
        """
        return SequentialChain([self, next_chain], client=self.client)


class PromptChain(Chain[dict, str]):
    """Chain that applies a prompt template and generates a response.
    
    Example:
        >>> template = PromptTemplate("Summarize: {text}")
        >>> chain = PromptChain(template)
        >>> result = chain.run({"text": "Long document..."})
        >>> print(result.value)
    """
    
    def __init__(
        self,
        template: PromptTemplate,
        client: Optional[OllamaClient] = None,
        **generate_kwargs,
    ) -> None:
        """Initialize the prompt chain.
        
        Args:
            template: Prompt template to use.
            client: Ollama client.
            **generate_kwargs: Additional args for generate().
        """
        super().__init__(client)
        self.template = template
        self.generate_kwargs = generate_kwargs
    
    def run(self, input_data: dict) -> ChainResult[str]:
        """Execute the prompt chain.
        
        Args:
            input_data: Variables for the template.
            
        Returns:
            ChainResult with generated text.
        """
        try:
            prompt = self.template.format(**input_data)
            response = self.client.generate(prompt, **self.generate_kwargs)
            return ChainResult.ok(
                response.text,
                tokens=response.eval_count,
                duration=response.duration_seconds,
            )
        except Exception as e:
            logger.error(f"PromptChain failed: {e}")
            return ChainResult.fail(str(e))


class SequentialChain(Chain[Any, Any]):
    """Chain that runs multiple chains in sequence.
    
    Each chain's output becomes the next chain's input.
    
    Example:
        >>> chain = SequentialChain([
        ...     ExtractChain(),
        ...     SummarizeChain(),
        ...     FormatChain(),
        ... ])
        >>> result = chain.run(document)
    """
    
    def __init__(
        self,
        chains: list[Chain],
        client: Optional[OllamaClient] = None,
    ) -> None:
        """Initialize the sequential chain.
        
        Args:
            chains: List of chains to run in order.
            client: Shared Ollama client.
        """
        super().__init__(client)
        self.chains = chains
        
        # Share client with all chains
        for chain in self.chains:
            chain.client = self.client
    
    def run(self, input_data: Any) -> ChainResult[Any]:
        """Execute all chains in sequence.
        
        Args:
            input_data: Initial input.
            
        Returns:
            ChainResult from the final chain.
        """
        current = input_data
        total_tokens = 0
        total_duration = 0.0
        
        for i, chain in enumerate(self.chains):
            result = chain.run(current)
            
            if not result.success:
                return ChainResult.fail(
                    f"Chain {i} failed: {result.error}",
                    failed_at=i,
                )
            
            current = result.value
            total_tokens += result.metadata.get("tokens", 0)
            total_duration += result.metadata.get("duration", 0.0)
        
        return ChainResult.ok(
            current,
            tokens=total_tokens,
            duration=total_duration,
            chain_count=len(self.chains),
        )
    
    def add(self, chain: Chain) -> "SequentialChain":
        """Add a chain to the sequence.
        
        Args:
            chain: Chain to add.
            
        Returns:
            Self for chaining.
        """
        chain.client = self.client
        self.chains.append(chain)
        return self


class MapReduceChain(Chain[list, str]):
    """Chain that maps over items and reduces results.
    
    Useful for processing large documents by chunks or
    analyzing multiple items and aggregating.
    
    Example:
        >>> chain = MapReduceChain(
        ...     map_template=PromptTemplate("Summarize: {item}"),
        ...     reduce_template=PromptTemplate("Combine: {summaries}"),
        ... )
        >>> result = chain.run(["doc1", "doc2", "doc3"])
    """
    
    def __init__(
        self,
        map_template: PromptTemplate,
        reduce_template: PromptTemplate,
        client: Optional[OllamaClient] = None,
        separator: str = "\n\n",
    ) -> None:
        """Initialize map-reduce chain.
        
        Args:
            map_template: Template for mapping each item.
            reduce_template: Template for reducing mapped results.
            client: Ollama client.
            separator: Separator for joining mapped results.
        """
        super().__init__(client)
        self.map_template = map_template
        self.reduce_template = reduce_template
        self.separator = separator
    
    def run(self, input_data: list) -> ChainResult[str]:
        """Execute map-reduce.
        
        Args:
            input_data: List of items to process.
            
        Returns:
            ChainResult with reduced output.
        """
        if not input_data:
            return ChainResult.fail("Empty input list")
        
        # Map phase
        mapped_results: list[str] = []
        total_tokens = 0
        
        for i, item in enumerate(input_data):
            try:
                prompt = self.map_template.format(item=item)
                response = self.client.generate(prompt)
                mapped_results.append(response.text)
                total_tokens += response.eval_count
            except Exception as e:
                logger.warning(f"Map failed for item {i}: {e}")
                continue
        
        if not mapped_results:
            return ChainResult.fail("All map operations failed")
        
        # Reduce phase
        try:
            combined = self.separator.join(mapped_results)
            reduce_prompt = self.reduce_template.format(summaries=combined)
            response = self.client.generate(reduce_prompt)
            total_tokens += response.eval_count
            
            return ChainResult.ok(
                response.text,
                tokens=total_tokens,
                mapped_count=len(mapped_results),
            )
        except Exception as e:
            return ChainResult.fail(f"Reduce failed: {e}")


class RouterChain(Chain[str, str]):
    """Chain that routes input to different chains based on classification.
    
    Example:
        >>> router = RouterChain({
        ...     "question": QuestionChain(),
        ...     "task": TaskChain(),
        ...     "chat": ChatChain(),
        ... })
        >>> result = router.run("What is the weather?")
    """
    
    def __init__(
        self,
        routes: dict[str, Chain],
        client: Optional[OllamaClient] = None,
        default_route: Optional[str] = None,
    ) -> None:
        """Initialize router chain.
        
        Args:
            routes: Mapping of route names to chains.
            client: Ollama client.
            default_route: Route to use if classification fails.
        """
        super().__init__(client)
        self.routes = routes
        self.default_route = default_route or list(routes.keys())[0]
        
        # Share client
        for chain in routes.values():
            chain.client = self.client
    
    def classify(self, input_text: str) -> str:
        """Classify input to determine route.
        
        Args:
            input_text: Input to classify.
            
        Returns:
            Route name.
        """
        route_names = ", ".join(self.routes.keys())
        prompt = f"""Classify the following input into exactly one category.
Categories: {route_names}

Input: {input_text}

Respond with only the category name, nothing else."""
        
        try:
            response = self.client.generate(prompt, temperature=0.1)
            route = response.text.strip().lower()
            
            # Find matching route
            for name in self.routes.keys():
                if name.lower() in route:
                    return name
            
            return self.default_route
        except Exception:
            return self.default_route
    
    def run(self, input_data: str) -> ChainResult[str]:
        """Route and execute appropriate chain.
        
        Args:
            input_data: Input text.
            
        Returns:
            ChainResult from the selected chain.
        """
        route = self.classify(input_data)
        chain = self.routes.get(route)
        
        if chain is None:
            return ChainResult.fail(f"Unknown route: {route}")
        
        result = chain.run(input_data)
        result.metadata["route"] = route
        return result


class TransformChain(Chain[T, R]):
    """Chain that applies a Python function transformation.
    
    Useful for pre/post-processing without LLM calls.
    
    Example:
        >>> chain = TransformChain(lambda x: x.upper())
        >>> result = chain.run("hello")
        >>> print(result.value)  # "HELLO"
    """
    
    def __init__(
        self,
        transform_fn: Callable[[T], R],
        client: Optional[OllamaClient] = None,
    ) -> None:
        """Initialize transform chain.
        
        Args:
            transform_fn: Function to apply.
            client: Ollama client (unused but kept for interface).
        """
        super().__init__(client)
        self.transform_fn = transform_fn
    
    def run(self, input_data: T) -> ChainResult[R]:
        """Apply the transformation.
        
        Args:
            input_data: Input to transform.
            
        Returns:
            ChainResult with transformed output.
        """
        try:
            result = self.transform_fn(input_data)
            return ChainResult.ok(result)
        except Exception as e:
            return ChainResult.fail(str(e))


class ConversationChain(Chain[str, str]):
    """Chain for multi-turn conversations with context.
    
    Maintains conversation history and supports system prompts.
    
    Example:
        >>> chain = ConversationChain(
        ...     system="You are a helpful biology tutor."
        ... )
        >>> chain.run("What is DNA?")
        >>> chain.run("How is it replicated?")
    """
    
    def __init__(
        self,
        system: Optional[str] = None,
        client: Optional[OllamaClient] = None,
    ) -> None:
        """Initialize conversation chain.
        
        Args:
            system: System prompt.
            client: Ollama client.
        """
        super().__init__(client)
        self.system = system
        self.messages: list[ChatMessage] = []
        
        if system:
            self.messages.append(ChatMessage("system", system))
    
    def run(self, input_data: str) -> ChainResult[str]:
        """Continue the conversation.
        
        Args:
            input_data: User message.
            
        Returns:
            ChainResult with assistant response.
        """
        self.messages.append(ChatMessage("user", input_data))
        
        try:
            response = self.client.chat(self.messages)
            self.messages.append(response.message)
            
            return ChainResult.ok(
                response.text,
                tokens=response.eval_count,
                duration=response.duration_seconds,
            )
        except Exception as e:
            # Remove failed user message
            self.messages.pop()
            return ChainResult.fail(str(e))
    
    def reset(self) -> None:
        """Reset conversation history."""
        self.messages.clear()
        if self.system:
            self.messages.append(ChatMessage("system", self.system))
