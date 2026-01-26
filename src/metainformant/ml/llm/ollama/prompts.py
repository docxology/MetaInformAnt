"""Prompt templates and utilities for Ollama.

This module provides composable prompt templates and message builders
for constructing prompts for LLM inference.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Optional


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
    
    @classmethod
    def system(cls, content: str) -> "ChatMessage":
        """Create a system message."""
        return cls(role="system", content=content)
    
    @classmethod
    def user(cls, content: str) -> "ChatMessage":
        """Create a user message."""
        return cls(role="user", content=content)
    
    @classmethod
    def assistant(cls, content: str) -> "ChatMessage":
        """Create an assistant message."""
        return cls(role="assistant", content=content)


@dataclass
class PromptTemplate:
    """A template for generating prompts with variable substitution.
    
    Uses {variable_name} syntax for placeholders.
    
    Attributes:
        template: The template string with placeholders.
        variables: Default values for variables.
        
    Example:
        >>> template = PromptTemplate(
        ...     "Analyze the {organism} gene {gene_name}.",
        ...     {"organism": "human"}
        ... )
        >>> prompt = template.format(gene_name="BRCA1")
        >>> print(prompt)
        Analyze the human gene BRCA1.
    """
    
    template: str
    variables: dict[str, str] = field(default_factory=dict)
    
    def format(self, **kwargs) -> str:
        """Format the template with provided variables.
        
        Args:
            **kwargs: Variable values to substitute.
            
        Returns:
            Formatted prompt string.
            
        Raises:
            KeyError: If a required variable is missing.
        """
        all_vars = {**self.variables, **kwargs}
        return self.template.format(**all_vars)
    
    def partial(self, **kwargs) -> "PromptTemplate":
        """Create a new template with some variables filled in.
        
        Args:
            **kwargs: Variable values to pre-fill.
            
        Returns:
            New PromptTemplate with updated defaults.
        """
        new_vars = {**self.variables, **kwargs}
        return PromptTemplate(self.template, new_vars)
    
    @property
    def variable_names(self) -> list[str]:
        """Get list of variable names in the template."""
        return re.findall(r"\{(\w+)\}", self.template)
    
    def validate(self, **kwargs) -> bool:
        """Check if all required variables are provided.
        
        Args:
            **kwargs: Variable values to check.
            
        Returns:
            True if all variables are available.
        """
        all_vars = {**self.variables, **kwargs}
        return all(name in all_vars for name in self.variable_names)


@dataclass
class SystemPrompt:
    """Builder for system prompts with common patterns.
    
    Provides a fluent interface for constructing system prompts
    with role, context, and instruction components.
    
    Example:
        >>> prompt = (SystemPrompt()
        ...     .role("bioinformatics expert")
        ...     .context("analyzing genomic data")
        ...     .instruction("be concise and cite sources")
        ...     .build())
    """
    
    _role: Optional[str] = None
    _context: Optional[str] = None
    _instructions: list[str] = field(default_factory=list)
    _constraints: list[str] = field(default_factory=list)
    
    def role(self, description: str) -> "SystemPrompt":
        """Set the assistant's role.
        
        Args:
            description: Role description.
            
        Returns:
            Self for chaining.
        """
        self._role = description
        return self
    
    def context(self, description: str) -> "SystemPrompt":
        """Set the context for the conversation.
        
        Args:
            description: Context description.
            
        Returns:
            Self for chaining.
        """
        self._context = description
        return self
    
    def instruction(self, instruction: str) -> "SystemPrompt":
        """Add an instruction.
        
        Args:
            instruction: Instruction text.
            
        Returns:
            Self for chaining.
        """
        self._instructions.append(instruction)
        return self
    
    def constraint(self, constraint: str) -> "SystemPrompt":
        """Add a constraint.
        
        Args:
            constraint: Constraint text.
            
        Returns:
            Self for chaining.
        """
        self._constraints.append(constraint)
        return self
    
    def build(self) -> str:
        """Build the complete system prompt.
        
        Returns:
            Formatted system prompt string.
        """
        parts: list[str] = []
        
        if self._role:
            parts.append(f"You are a {self._role}.")
        
        if self._context:
            parts.append(f"Context: {self._context}")
        
        if self._instructions:
            parts.append("Instructions:")
            for inst in self._instructions:
                parts.append(f"- {inst}")
        
        if self._constraints:
            parts.append("Constraints:")
            for const in self._constraints:
                parts.append(f"- {const}")
        
        return "\n".join(parts)
    
    def to_message(self) -> ChatMessage:
        """Convert to a ChatMessage.
        
        Returns:
            ChatMessage with system role.
        """
        return ChatMessage.system(self.build())


# Pre-built templates for bioinformatics tasks

GENE_ANALYSIS_TEMPLATE = PromptTemplate(
    """Analyze the following gene information:

Gene: {gene_name}
Organism: {organism}

Provide a summary including:
1. Gene function
2. Associated pathways
3. Known disease associations
4. Relevant literature references
""",
    {"organism": "Homo sapiens"},
)


SEQUENCE_ANALYSIS_TEMPLATE = PromptTemplate(
    """Analyze the following {sequence_type} sequence:

```
{sequence}
```

Identify:
1. Sequence characteristics
2. Potential functional elements
3. Notable patterns or motifs
""",
    {"sequence_type": "DNA"},
)


LITERATURE_SUMMARY_TEMPLATE = PromptTemplate(
    """Summarize the following scientific abstract in {style} style:

Title: {title}

Abstract:
{abstract}

Provide:
1. Key findings
2. Methods used
3. Implications
""",
    {"style": "concise"},
)


def build_bioinformatics_prompt(
    task: str,
    data: str,
    organism: str = "Homo sapiens",
    output_format: str = "structured",
) -> str:
    """Build a bioinformatics-focused prompt.
    
    Args:
        task: The analysis task to perform.
        data: The data to analyze.
        organism: Target organism.
        output_format: Desired output format.
        
    Returns:
        Formatted prompt string.
        
    Example:
        >>> prompt = build_bioinformatics_prompt(
        ...     task="identify functional domains",
        ...     data="MTEYKLVVVG...",
        ...     organism="Homo sapiens"
        ... )
    """
    return f"""You are a bioinformatics expert specializing in {organism} genomics.

Task: {task}

Data:
{data}

Please provide a {output_format} response with:
1. Analysis results
2. Confidence assessment
3. Recommendations for follow-up analysis
"""


def build_conversation_messages(
    system_prompt: str,
    user_messages: list[str],
    assistant_messages: Optional[list[str]] = None,
) -> list[ChatMessage]:
    """Build a conversation message list.
    
    Args:
        system_prompt: System message content.
        user_messages: List of user messages.
        assistant_messages: Optional list of assistant responses.
        
    Returns:
        List of ChatMessage objects.
        
    Example:
        >>> messages = build_conversation_messages(
        ...     "You are a helpful assistant.",
        ...     ["Hello!", "What is DNA?"],
        ...     ["Hi there!", None]  # None means no prior response
        ... )
    """
    messages = [ChatMessage.system(system_prompt)]
    
    assistant_messages = assistant_messages or []
    
    for i, user_msg in enumerate(user_messages):
        messages.append(ChatMessage.user(user_msg))
        
        if i < len(assistant_messages) and assistant_messages[i]:
            messages.append(ChatMessage.assistant(assistant_messages[i]))
    
    return messages
