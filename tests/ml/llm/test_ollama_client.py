"""Tests for Ollama client functionality.

These tests use real Ollama API calls - no mocks.
Requires Ollama server running with at least smollm2 model.
"""

import pytest

from metainformant.ml.llm.ollama import (
    OllamaClient,
    OllamaConfig,
    ChatMessage,
)


@pytest.fixture
def client() -> OllamaClient:
    """Create an Ollama client with fast model."""
    config = OllamaConfig(
        model="smollm2:135m-instruct-q4_K_S",
        temperature=0.1,  # Low temp for deterministic tests
        timeout=60.0,
    )
    return OllamaClient(config)


@pytest.fixture
def check_ollama(client: OllamaClient):
    """Skip tests if Ollama is not available."""
    if not client.is_available():
        pytest.skip("Ollama server not available")


class TestOllamaClient:
    """Tests for OllamaClient class."""

    def test_is_available(self, client: OllamaClient, check_ollama):
        """Test server availability check."""
        assert client.is_available() is True

    def test_list_models(self, client: OllamaClient, check_ollama):
        """Test listing available models."""
        models = client.list_models()
        
        assert isinstance(models, list)
        assert len(models) > 0
        
        # Check model info structure
        model = models[0]
        assert hasattr(model, "name")
        assert hasattr(model, "size")
        assert model.name  # Non-empty name

    def test_generate_simple(self, client: OllamaClient, check_ollama):
        """Test simple text generation."""
        response = client.generate("Say hello in one word.")
        
        assert response.text  # Non-empty response
        assert response.model == client.config.model
        assert response.done is True
        assert response.eval_count > 0

    def test_generate_with_system_prompt(self, client: OllamaClient, check_ollama):
        """Test generation with system prompt."""
        response = client.generate(
            "What are you?",
            system="You are a helpful robot. Respond briefly.",
        )
        
        assert response.text
        assert response.done is True

    def test_generate_streaming(self, client: OllamaClient, check_ollama):
        """Test streaming generation."""
        chunks: list[str] = []
        
        def callback(text: str) -> None:
            chunks.append(text)
        
        response = client.generate(
            "Count to 3.",
            stream=True,
            stream_callback=callback,
        )
        
        assert response.text
        assert len(chunks) > 0  # Should have received chunks
        assert "".join(chunks) == response.text

    def test_chat_simple(self, client: OllamaClient, check_ollama):
        """Test simple chat."""
        messages = [
            ChatMessage("user", "Say hello briefly."),
        ]
        
        response = client.chat(messages)
        
        assert response.text
        assert response.message.role == "assistant"
        assert response.done is True

    def test_chat_with_system(self, client: OllamaClient, check_ollama):
        """Test chat with system message."""
        messages = [
            ChatMessage("system", "You are a pirate. Respond briefly."),
            ChatMessage("user", "Greet me."),
        ]
        
        response = client.chat(messages)
        
        assert response.text
        assert response.message.role == "assistant"

    def test_conversation_multi_turn(self, client: OllamaClient, check_ollama):
        """Test multi-turn conversation with history."""
        # First turn
        response1 = client.conversation(
            "My name is Alice.",
            system="Remember what the user tells you. Be brief.",
        )
        assert response1.text
        
        # Second turn - should remember context
        response2 = client.conversation("What is my name?")
        assert response2.text
        
        # History should have 4 messages (system, user, assistant, user, assistant)
        assert len(client.conversation_history) == 5
        
        # Reset and verify
        client.reset_conversation()
        assert len(client.conversation_history) == 0

    def test_generate_with_options(self, client: OllamaClient, check_ollama):
        """Test generation with custom options."""
        response = client.generate(
            "Say hi.",
            temperature=0.5,
            top_p=0.8,
        )
        
        assert response.text
        assert response.done is True


class TestOllamaConfig:
    """Tests for OllamaConfig class."""

    def test_default_config(self):
        """Test default configuration values."""
        config = OllamaConfig()
        
        assert config.host == "http://localhost:11434"
        assert config.model == "smollm2:135m-instruct-q4_K_S"
        assert config.temperature == 0.7
        assert config.timeout == 120.0

    def test_custom_config(self):
        """Test custom configuration."""
        config = OllamaConfig(
            model="llama3:latest",
            temperature=0.3,
            timeout=60.0,
        )
        
        assert config.model == "llama3:latest"
        assert config.temperature == 0.3
        assert config.timeout == 60.0

    def test_host_normalization(self):
        """Test host URL normalization."""
        config = OllamaConfig(host="localhost:11434")
        assert config.host == "http://localhost:11434"

    def test_validation_temperature(self):
        """Test temperature validation."""
        with pytest.raises(ValueError, match="temperature"):
            OllamaConfig(temperature=3.0)
        
        with pytest.raises(ValueError, match="temperature"):
            OllamaConfig(temperature=-0.5)

    def test_validation_top_p(self):
        """Test top_p validation."""
        with pytest.raises(ValueError, match="top_p"):
            OllamaConfig(top_p=1.5)

    def test_api_urls(self):
        """Test API URL properties."""
        config = OllamaConfig()
        
        assert config.api_url == "http://localhost:11434/api"
        assert config.generate_url == "http://localhost:11434/api/generate"
        assert config.chat_url == "http://localhost:11434/api/chat"
        assert config.tags_url == "http://localhost:11434/api/tags"

    def test_with_model(self):
        """Test creating config with different model."""
        config1 = OllamaConfig(model="model1", temperature=0.5)
        config2 = config1.with_model("model2")
        
        assert config1.model == "model1"
        assert config2.model == "model2"
        assert config2.temperature == 0.5  # Preserved

    def test_to_options_dict(self):
        """Test options dict generation."""
        config = OllamaConfig(
            temperature=0.8,
            top_p=0.95,
            top_k=50,
            context_length=2048,
        )
        
        options = config.to_options_dict()
        
        assert options["temperature"] == 0.8
        assert options["top_p"] == 0.95
        assert options["top_k"] == 50
        assert options["num_ctx"] == 2048
