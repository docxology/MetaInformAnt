"""Tests for Ollama chains functionality.

These tests use real Ollama API calls - no mocks.
"""

import pytest

from metainformant.ml.llm.ollama import (
    OllamaClient,
    OllamaConfig,
    PromptTemplate,
    Chain,
    SequentialChain,
    MapReduceChain,
    TransformChain,
    ConversationChain,
)
from metainformant.ml.llm.ollama.chains import PromptChain, ChainResult


@pytest.fixture
def client() -> OllamaClient:
    """Create an Ollama client with fast model."""
    config = OllamaConfig(
        model="smollm2:135m-instruct-q4_K_S",
        temperature=0.1,
        timeout=60.0,
    )
    return OllamaClient(config)


@pytest.fixture
def check_ollama(client: OllamaClient):
    """Skip tests if Ollama is not available or model is not installed."""
    if not client.is_available():
        pytest.skip("Ollama server not available")
    # Also check that the specific model is installed
    try:
        models = client.list_models()
        model_names = [m.name for m in models]
        if not any(client.config.model in name for name in model_names):
            pytest.skip(f"Model {client.config.model} not installed on Ollama server")
    except Exception:
        pytest.skip("Could not verify model availability")


class TestChainResult:
    """Tests for ChainResult class."""

    def test_ok_result(self):
        """Test successful result creation."""
        result = ChainResult.ok("test value", tokens=10)

        assert result.success is True
        assert result.value == "test value"
        assert result.error is None
        assert result.metadata["tokens"] == 10

    def test_fail_result(self):
        """Test failed result creation."""
        result = ChainResult.fail("something went wrong")

        assert result.success is False
        assert result.value is None
        assert result.error == "something went wrong"


class TestPromptTemplate:
    """Tests for PromptTemplate class."""

    def test_simple_format(self):
        """Test simple variable substitution."""
        template = PromptTemplate("Hello, {name}!")
        result = template.format(name="Alice")

        assert result == "Hello, Alice!"

    def test_with_defaults(self):
        """Test template with default values."""
        template = PromptTemplate(
            "Hello, {name}! You are {age} years old.",
            {"name": "Unknown", "age": "unknown"},
        )

        # Use defaults
        result1 = template.format()
        assert "Unknown" in result1

        # Override one default
        result2 = template.format(name="Bob")
        assert "Bob" in result2
        assert "unknown" in result2

    def test_partial(self):
        """Test partial template filling."""
        template = PromptTemplate("Hello, {name}! Age: {age}")
        partial = template.partial(name="Alice")

        result = partial.format(age=30)
        assert result == "Hello, Alice! Age: 30"

    def test_variable_names(self):
        """Test extracting variable names."""
        template = PromptTemplate("{greeting}, {name}! {message}")

        names = template.variable_names
        assert set(names) == {"greeting", "name", "message"}

    def test_validate(self):
        """Test validation of variables."""
        template = PromptTemplate("{a}, {b}, {c}")

        assert template.validate(a=1, b=2, c=3) is True
        assert template.validate(a=1, b=2) is False

        # With defaults
        template2 = PromptTemplate("{a}, {b}", {"a": "default"})
        assert template2.validate(b=2) is True


class TestTransformChain:
    """Tests for TransformChain class."""

    def test_simple_transform(self):
        """Test simple transformation."""
        chain = TransformChain(lambda x: x.upper())
        result = chain.run("hello")

        assert result.success is True
        assert result.value == "HELLO"

    def test_transform_failure(self):
        """Test transformation that raises exception."""
        chain = TransformChain(lambda x: x / 0)
        result = chain.run(10)

        assert result.success is False
        assert "division" in result.error.lower()

    def test_chain_composition(self):
        """Test chaining transforms."""
        upper = TransformChain(lambda x: x.upper())
        add_exclaim = TransformChain(lambda x: x + "!")

        pipeline = SequentialChain([upper, add_exclaim])
        result = pipeline.run("hello")

        assert result.success is True
        assert result.value == "HELLO!"


class TestPromptChain:
    """Tests for PromptChain with real Ollama."""

    def test_simple_prompt_chain(self, client: OllamaClient, check_ollama):
        """Test simple prompt chain execution."""
        template = PromptTemplate("Say the word: {word}")
        chain = PromptChain(template, client=client)

        result = chain.run({"word": "hello"})

        assert result.success is True
        assert result.value  # Non-empty response
        assert "tokens" in result.metadata


class TestSequentialChain:
    """Tests for SequentialChain."""

    def test_transform_sequence(self):
        """Test sequence of transforms."""
        step1 = TransformChain(lambda x: x + 1)
        step2 = TransformChain(lambda x: x * 2)
        step3 = TransformChain(lambda x: x - 3)

        chain = SequentialChain([step1, step2, step3])
        result = chain.run(5)

        # (5 + 1) * 2 - 3 = 9
        assert result.success is True
        assert result.value == 9
        assert result.metadata["chain_count"] == 3

    def test_sequence_with_failure(self):
        """Test sequence that fails midway."""
        step1 = TransformChain(lambda x: x + 1)
        step2 = TransformChain(lambda x: 1 / 0)  # Will fail
        step3 = TransformChain(lambda x: x - 3)

        chain = SequentialChain([step1, step2, step3])
        result = chain.run(5)

        assert result.success is False
        assert result.metadata["failed_at"] == 1

    def test_add_chain(self):
        """Test adding chain to sequence."""
        chain = SequentialChain([TransformChain(lambda x: x + 1)])
        chain.add(TransformChain(lambda x: x * 2))

        result = chain.run(5)
        assert result.value == 12


class TestConversationChain:
    """Tests for ConversationChain with real Ollama."""

    def test_conversation_without_system(self, client: OllamaClient, check_ollama):
        """Test conversation chain without system prompt."""
        chain = ConversationChain(client=client)

        result = chain.run("Say hello briefly.")

        assert result.success is True
        assert result.value
        assert len(chain.messages) == 2  # user + assistant

    def test_conversation_with_system(self, client: OllamaClient, check_ollama):
        """Test conversation chain with system prompt."""
        chain = ConversationChain(
            system="You are a helpful assistant. Be very brief.",
            client=client,
        )

        result = chain.run("What are you?")

        assert result.success is True
        assert result.value
        assert len(chain.messages) == 3  # system + user + assistant

    def test_conversation_multi_turn(self, client: OllamaClient, check_ollama):
        """Test multi-turn conversation."""
        chain = ConversationChain(
            system="Remember what I tell you. Be brief.",
            client=client,
        )

        chain.run("My favorite color is blue.")
        chain.run("What is my favorite color?")

        # Should have: system, user, assistant, user, assistant
        assert len(chain.messages) == 5

    def test_conversation_reset(self, client: OllamaClient, check_ollama):
        """Test conversation reset."""
        chain = ConversationChain(
            system="Test system prompt.",
            client=client,
        )

        chain.run("Hello")
        assert len(chain.messages) == 3

        chain.reset()
        # After reset, should only have system message
        assert len(chain.messages) == 1
        assert chain.messages[0].role == "system"


class TestMapReduceChain:
    """Tests for MapReduceChain."""

    def test_map_reduce_with_transforms(self):
        """Test map-reduce using transform chains."""
        # Simple numeric map-reduce
        items = [1, 2, 3, 4, 5]

        # Map: square each
        # Reduce: sum (manual since we're using transforms)
        map_chain = TransformChain(lambda x: [i**2 for i in x])
        reduce_chain = TransformChain(lambda x: sum(x))

        pipeline = SequentialChain([map_chain, reduce_chain])
        result = pipeline.run(items)

        assert result.success is True
        assert result.value == 55  # 1 + 4 + 9 + 16 + 25
