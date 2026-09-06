"""Tests for Ollama prompt builders and response dataclasses.

Pure offline tests - real implementations, no network required.
"""

from __future__ import annotations

import pytest

from metainformant.ml.llm.ollama.client import (
    ChatResponse,
    GenerateResponse,
    ModelInfo,
)
from metainformant.ml.llm.ollama.prompts import (
    ChatMessage,
    SystemPrompt,
    build_bioinformatics_prompt,
    build_conversation_messages,
)


class TestSystemPrompt:
    def test_empty_build(self):
        assert SystemPrompt().build() == ""

    def test_fluent_builder(self):
        prompt = (
            SystemPrompt()
            .role("bioinformatics expert")
            .context("analyzing genomic data")
            .instruction("be concise")
            .constraint("cite sources")
            .build()
        )

        assert "You are a bioinformatics expert." in prompt
        assert "Context: analyzing genomic data" in prompt
        assert "Instructions:\n- be concise" in prompt
        assert "Constraints:\n- cite sources" in prompt

    def test_to_message_is_system_role(self):
        message = SystemPrompt().role("tutor").to_message()

        assert isinstance(message, ChatMessage)
        assert message.role == "system"
        assert "tutor" in message.content


class TestPromptBuilders:
    def test_build_bioinformatics_prompt(self):
        prompt = build_bioinformatics_prompt(
            task="identify functional domains",
            data="MTEYKLVVVG",
            organism="Mus musculus",
            output_format="JSON",
        )

        assert "Mus musculus" in prompt
        assert "identify functional domains" in prompt
        assert "MTEYKLVVVG" in prompt
        assert "JSON" in prompt

    def test_build_conversation_messages_alternates(self):
        messages = build_conversation_messages(
            "You are a helpful assistant.",
            ["Hello!", "What is DNA?"],
            ["Hi there!"],
        )

        roles = [m.role for m in messages]
        assert roles == ["system", "user", "assistant", "user"]
        assert messages[0].content == "You are a helpful assistant."
        assert messages[-1].content == "What is DNA?"

    def test_build_conversation_messages_no_assistant(self):
        messages = build_conversation_messages("system", ["q1", "q2"])

        assert [m.role for m in messages] == ["system", "user", "user"]


class TestResponseDataclasses:
    def test_generate_response_from_dict(self):
        response = GenerateResponse.from_dict(
            {
                "response": "text",
                "model": "stub:latest",
                "done": True,
                "total_duration": 2_000_000_000,
                "eval_count": 10,
                "prompt_eval_count": 5,
                "context": [7, 8],
            }
        )

        assert response.text == "text"
        assert response.duration_seconds == pytest.approx(2.0)
        assert response.tokens_per_second == pytest.approx(5.0)
        assert response.context == [7, 8]

    def test_generate_response_defaults(self):
        response = GenerateResponse.from_dict({})

        assert response.text == ""
        assert response.done is True
        assert response.context is None
        assert response.tokens_per_second == 0.0

    def test_chat_response_text_property(self):
        response = ChatResponse.from_dict(
            {"message": {"role": "assistant", "content": "answer"}, "model": "stub"}
        )

        assert response.text == "answer"
        assert response.message.role == "assistant"

    def test_model_info_size_gb(self):
        info = ModelInfo.from_dict({"name": "m", "size": 3 * 1024**3, "digest": "abc"})

        assert info.size_gb == pytest.approx(3.0)
        assert info.digest == "abc"
        assert info.modified_at == ""

    def test_chat_message_to_dict(self):
        assert ChatMessage("user", "hi").to_dict() == {"role": "user", "content": "hi"}


def test_generate_response_zero_tokens_per_second():
    """No eval data means tokens_per_second must be exactly 0."""
    response = GenerateResponse(text="x", model="m", eval_count=0)

    assert response.tokens_per_second == 0.0
