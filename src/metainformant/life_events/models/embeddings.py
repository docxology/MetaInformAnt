"""Event embedding and representation learning for life events analysis.

This module provides functions for learning vector representations of life events
and event sequences using various embedding techniques.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)

# Optional imports for ML functionality
try:
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE

    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    logger.warning("scikit-learn not available, some embedding methods disabled")


def learn_event_embeddings(
    sequences: List[Any],
    embedding_dim: int = 100,
    window_size: int = 5,
    min_count: int = 1,
    sg: int = 1,
    epochs: int = 5,
    **kwargs: Any,
) -> Dict[str, Any]:
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
    logger.info(f"Learning event embeddings with dim={embedding_dim}, window={window_size}")

    # Build vocabulary
    event_counts = {}
    for seq in sequences:
        for event in seq.events:
            event_type = event.event_type
            event_counts[event_type] = event_counts.get(event_type, 0) + 1

    # Filter by min_count
    vocabulary = {event: count for event, count in event_counts.items() if count >= min_count}

    if not vocabulary:
        raise ValueError("No events meet minimum count threshold")

    logger.info(f"Vocabulary size: {len(vocabulary)} events")

    # Initialize embeddings randomly
    vocab_size = len(vocabulary)
    embeddings = np.random.normal(0, 0.1, (vocab_size, embedding_dim))

    # Create event-to-index mapping
    event_to_idx = {event: i for i, event in enumerate(vocabulary.keys())}

    # Simple training (placeholder - would implement full Word2Vec in production)
    # For now, just return random embeddings
    logger.warning("Using random embeddings - full Word2Vec implementation needed for production")

    # Add some structure based on event types
    for event, idx in event_to_idx.items():
        # Add small bias based on event type
        if "education" in event.lower():
            embeddings[idx, :10] += 0.1
        elif "job" in event.lower() or "career" in event.lower():
            embeddings[idx, 10:20] += 0.1
        elif "health" in event.lower() or "medical" in event.lower():
            embeddings[idx, 20:30] += 0.1
        elif "family" in event.lower() or "marriage" in event.lower():
            embeddings[idx, 30:40] += 0.1

    result = {
        "embeddings": embeddings,
        "vocabulary": vocabulary,
        "event_to_idx": event_to_idx,
        "embedding_dim": embedding_dim,
        "vocab_size": vocab_size,
        "training_params": {"window_size": window_size, "min_count": min_count, "sg": sg, "epochs": epochs},
    }

    logger.info(f"Embeddings learned: {vocab_size} events, {embedding_dim} dimensions")
    return result


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
    sequences: List[Any], domains: List[str], embedding_dim: int = 50, **kwargs: Any
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

    domain_embeddings = {}

    for domain in domains:
        # Filter sequences to this domain
        domain_sequences = []
        for seq in sequences:
            domain_events = [event for event in seq.events if hasattr(event, "domain") and event.domain == domain]
            if domain_events:
                # Create a sequence with only this domain's events
                domain_seq = type(seq)(person_id=seq.person_id, events=domain_events)
                domain_sequences.append(domain_seq)

        if domain_sequences:
            # Learn embeddings for this domain
            domain_result = learn_event_embeddings(domain_sequences, embedding_dim=embedding_dim, **kwargs)
            domain_embeddings[domain] = domain_result

    return {"domain_embeddings": domain_embeddings, "domains": domains, "embedding_dim": embedding_dim}
