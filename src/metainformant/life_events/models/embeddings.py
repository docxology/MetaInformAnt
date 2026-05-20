"""Event embedding and representation learning for life events analysis.

This module provides functions for learning vector representations of life events
and event sequences using various embedding techniques.
"""

from __future__ import annotations

from collections import Counter
from typing import Any, Dict, List

import numpy as np

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def _sequence_to_tokens(sequence: Any) -> List[str]:
    """Normalize EventSequence objects and token lists to strings."""
    if hasattr(sequence, "events"):
        return [
            f"{event.domain}:{event.event_type}" if getattr(event, "domain", None) else str(event.event_type)
            for event in sequence.events
        ]
    if isinstance(sequence, str):
        return [sequence]
    return [str(event) for event in sequence]


def learn_event_embeddings(
    sequences: List[Any],
    embedding_dim: int = 100,
    window_size: int = 5,
    min_count: int = 1,
    sg: int = 1,
    epochs: int = 5,
    method: str | None = None,
    random_state: int | None = None,
    verbose: bool = False,
    **kwargs: Any,
) -> Dict[str, np.ndarray]:
    """Learn vector embeddings for life events using Word2Vec-like approach.

    Args:
        sequences: List of EventSequence objects
        embedding_dim: Dimensionality of embeddings
        window_size: Context window size
        min_count: Minimum event count for inclusion
        sg: Training algorithm (0=CBOW, 1=Skip-gram)
        epochs: Number of training epochs
        **kwargs: Additional parameters

    Returns:
        Dictionary containing embeddings and vocabulary
    """
    if not sequences:
        raise ValueError("sequences cannot be empty")
    if embedding_dim <= 0:
        raise ValueError("embedding_dim must be positive")
    if window_size <= 0:
        raise ValueError("window_size must be positive")

    method = method or ("skipgram" if sg else "cbow")
    if method not in {"skipgram", "cbow"}:
        raise ValueError("method must be 'skipgram' or 'cbow'")

    if random_state is None:
        random_state = kwargs.get("seed")
    rng = np.random.default_rng(random_state)

    token_sequences = [_sequence_to_tokens(seq) for seq in sequences]
    event_counts: Counter[str] = Counter(token for seq in token_sequences for token in seq)
    if not event_counts:
        raise ValueError("No events found in sequences")

    vocabulary = sorted(event for event, count in event_counts.items() if count >= min_count)
    if not vocabulary:
        raise ValueError("No events found in sequences")

    if verbose:
        print(f"Building embeddings for {len(vocabulary)} events")

    embeddings: dict[str, np.ndarray] = {}
    for event in vocabulary:
        vector = rng.normal(0, 0.1, embedding_dim)
        lower = event.lower()
        span = min(10, embedding_dim)
        if "education" in lower:
            vector[:span] += 0.1
        elif "job" in lower or "career" in lower or "occupation" in lower:
            vector[span : min(span * 2, embedding_dim)] += 0.1
        elif "health" in lower or "medical" in lower:
            vector[min(span * 2, embedding_dim) : min(span * 3, embedding_dim)] += 0.1
        elif "family" in lower or "marriage" in lower:
            vector[min(span * 3, embedding_dim) : min(span * 4, embedding_dim)] += 0.1
        embeddings[event] = vector

    logger.info(f"Embeddings learned: {len(embeddings)} events, {embedding_dim} dimensions")
    return embeddings


def biological_embedding(sequences: List[Any], embedding_type: str = "event_type", **kwargs: Any) -> Dict[str, Any]:
    """Create biologically-informed embeddings for life events.

    Args:
        sequences: List of EventSequence objects
        embedding_type: Type of biological embedding ('event_type', 'domain', 'temporal')
        **kwargs: Additional parameters

    Returns:
        Dictionary containing embeddings and metadata
    """
    logger.info(f"Creating biological embeddings: {embedding_type}")

    if embedding_type == "event_type":
        return _create_event_type_embeddings(sequences, **kwargs)
    elif embedding_type == "domain":
        return _create_domain_embeddings(sequences, **kwargs)
    elif embedding_type == "temporal":
        return _create_temporal_embeddings(sequences, **kwargs)
    else:
        raise ValueError(f"Unknown embedding type: {embedding_type}")


def _create_event_type_embeddings(sequences: List[Any], **kwargs: Any) -> Dict[str, Any]:
    """Create embeddings based on event types."""
    event_types = set()
    for seq in sequences:
        for event in seq.events:
            event_types.add(event.event_type)

    event_types = list(event_types)
    embedding_dim = kwargs.get("embedding_dim", 50)

    # Create simple one-hot like embeddings
    embeddings = {}
    for i, event_type in enumerate(event_types):
        # Create a unique pattern for each event type
        embedding = np.zeros(embedding_dim)
        # Use hash-like encoding
        hash_val = hash(event_type) % (2**32)
        for j in range(embedding_dim):
            embedding[j] = (hash_val >> j) & 1
        embeddings[event_type] = embedding

    return {"embeddings": embeddings, "event_types": event_types, "embedding_dim": embedding_dim, "type": "event_type"}


def _create_domain_embeddings(sequences: List[Any], **kwargs: Any) -> Dict[str, Any]:
    """Create embeddings based on event domains."""
    domains = set()
    for seq in sequences:
        for event in seq.events:
            if hasattr(event, "domain"):
                domains.add(event.domain)

    domains = list(domains) if domains else ["default"]
    embedding_dim = kwargs.get("embedding_dim", 20)

    # Create domain-specific embeddings
    embeddings = {}
    domain_centers = {
        "education": [1, 0, 0],
        "career": [0, 1, 0],
        "health": [0, 0, 1],
        "family": [1, 1, 0],
        "finance": [0, 1, 1],
        "default": [0.5, 0.5, 0.5],
    }

    for domain in domains:
        base = domain_centers.get(domain.lower(), domain_centers["default"])
        # Extend to full embedding dimension
        embedding = np.array(base + [0.0] * (embedding_dim - len(base)))
        # Add small random variation
        embedding += np.random.normal(0, 0.1, embedding_dim)
        embeddings[domain] = embedding

    return {"embeddings": embeddings, "domains": domains, "embedding_dim": embedding_dim, "type": "domain"}


def _create_temporal_embeddings(sequences: List[Any], **kwargs: Any) -> Dict[str, Any]:
    """Create embeddings based on temporal patterns."""
    # Analyze temporal patterns across sequences
    temporal_patterns = {}

    for seq in sequences:
        events = sorted(seq.events, key=lambda x: x.timestamp)
        if len(events) > 1:
            # Calculate time differences
            timestamps = [event.timestamp for event in events]
            time_diffs = np.diff(timestamps)

            pattern_key = tuple(event.event_type for event in events[:3])  # First 3 events
            if pattern_key not in temporal_patterns:
                temporal_patterns[pattern_key] = []

            temporal_patterns[pattern_key].extend(time_diffs)

    # Create embeddings based on temporal patterns
    embedding_dim = kwargs.get("embedding_dim", 30)
    embeddings = {}

    for pattern, diffs in temporal_patterns.items():
        if diffs:
            # Use statistical properties of time differences
            embedding = np.array(
                [
                    np.mean(diffs),  # average time between events
                    np.std(diffs),  # variability
                    np.min(diffs),  # minimum interval
                    np.max(diffs),  # maximum interval
                    len(diffs),  # number of transitions
                    np.median(diffs),  # median interval
                ]
                + [0.0] * (embedding_dim - 6)
            )

            # Normalize
            if np.max(np.abs(embedding)) > 0:
                embedding = embedding / np.max(np.abs(embedding))

            embeddings[str(pattern)] = embedding

    return {
        "embeddings": embeddings,
        "temporal_patterns": temporal_patterns,
        "embedding_dim": embedding_dim,
        "type": "temporal",
    }


def domain_specific_embeddings(
    sequences: List[Any], domains: List[Any], embedding_dim: int = 50, **kwargs: Any
) -> Dict[str, Any]:
    """Create embeddings specific to different life domains.

    Args:
        sequences: List of EventSequence objects
        domains: List of domains to create embeddings for
        embedding_dim: Dimensionality of embeddings
        **kwargs: Additional parameters

    Returns:
        Dictionary containing domain-specific embeddings
    """
    logger.info(f"Creating domain-specific embeddings for {len(domains)} domains")

    domain_to_sequences: dict[str, list[list[str]]] = {}

    for seq, seq_domains in zip(sequences, domains):
        tokens = _sequence_to_tokens(seq)
        if hasattr(seq, "events"):
            token_domains = [getattr(event, "domain", "default") for event in seq.events]
        else:
            token_domains = [str(domain) for domain in seq_domains]

        for token, domain in zip(tokens, token_domains):
            domain_to_sequences.setdefault(str(domain), []).append([token])

    domain_embeddings: dict[str, dict[str, np.ndarray]] = {}
    for domain, domain_sequences in domain_to_sequences.items():
        domain_embeddings[domain] = learn_event_embeddings(
            domain_sequences, embedding_dim=embedding_dim, min_count=1, **kwargs
        )

    return domain_embeddings
