"""Event embedding methods for learning dense vector representations.

This module provides methods for learning embeddings of events and event sequences,
inspired by NLP word embedding techniques.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray


def _build_event_contexts(
    sequences: List[List[str]],
    window_size: int = 5
) -> List[Tuple[str, str]]:
    """Build context pairs from event sequences using sliding window.
    
    Args:
        sequences: List of event sequences (each sequence is list of event tokens)
        window_size: Size of context window
        
    Returns:
        List of (target, context) tuples
    """
    contexts = []
    for seq in sequences:
        for i, target in enumerate(seq):
            # Get context events within window
            start = max(0, i - window_size)
            end = min(len(seq), i + window_size + 1)
            for j in range(start, end):
                if j != i:
                    contexts.append((target, seq[j]))
    return contexts


def learn_event_embeddings(
    sequences: List[List[str]],
    embedding_dim: int = 100,
    window_size: int = 5,
    method: str = "skipgram",
    epochs: int = 10,
    learning_rate: float = 0.01,
    random_state: Optional[int] = None,
    verbose: bool = False
) -> Dict[str, NDArray[np.float64]]:
    """Learn dense vector representations of events using Word2Vec-style methods.
    
    Implements Skip-gram or CBOW to learn embeddings where similar events
    (those that appear in similar contexts) have similar vectors.
    
    Args:
        sequences: List of event sequences, where each sequence is a list of
            event tokens (e.g., "health:diagnosis" or "occupation:job_change")
        embedding_dim: Dimension of embedding vectors
        window_size: Size of context window for training
        method: Embedding method ("skipgram" or "cbow")
        epochs: Number of training epochs
        learning_rate: Learning rate for gradient descent
        random_state: Random seed for reproducibility
        verbose: If True, print progress information for large datasets
        
    Returns:
        Dictionary mapping event tokens to embedding vectors
        
    Examples:
        >>> sequences = [
        ...     ["health:diagnosis", "occupation:job_change", "income:raise"],
        ...     ["education:degree", "occupation:job_change", "address:move"],
        ... ]
        >>> embeddings = learn_event_embeddings(sequences, embedding_dim=50)
        >>> "health:diagnosis" in embeddings
        True
        >>> embeddings["health:diagnosis"].shape
        (50,)
    """
    if not sequences:
        raise ValueError("sequences list cannot be empty")
    
    if embedding_dim <= 0:
        raise ValueError(f"embedding_dim must be positive, got {embedding_dim}")
    
    if window_size <= 0:
        raise ValueError(f"window_size must be positive, got {window_size}")
    
    if method not in ["skipgram", "cbow"]:
        raise ValueError(f"method must be 'skipgram' or 'cbow', got {method}")
    
    if random_state is not None:
        np.random.seed(random_state)
    
    # Build vocabulary
    vocab = set()
    for seq in sequences:
        if seq:  # Skip empty sequences
            vocab.update(seq)
    
    if not vocab:
        raise ValueError("No events found in sequences (all sequences may be empty)")
    
    vocab = sorted(list(vocab))
    vocab_size = len(vocab)
    
    if verbose:
        print(f"Building embeddings for {vocab_size} unique events from {len(sequences)} sequences")
    
    # Create token to index mapping
    token_to_idx = {token: i for i, token in enumerate(vocab)}
    idx_to_token = {i: token for token, i in token_to_idx.items()}
    
    # Initialize embeddings (input and output)
    W_in = np.random.normal(0, 0.01, (vocab_size, embedding_dim))
    W_out = np.random.normal(0, 0.01, (vocab_size, embedding_dim))
    
    # Build context pairs
    contexts = _build_event_contexts(sequences, window_size=window_size)
    
    if not contexts:
        # Return random embeddings if no contexts
        return {token: np.random.normal(0, 0.01, embedding_dim) for token in vocab}
    
    # Training loop
    for epoch in range(epochs):
        total_loss = 0.0
        np.random.shuffle(contexts)
        
        if verbose and epochs > 1:
            print(f"Epoch {epoch + 1}/{epochs} ({len(contexts)} context pairs)")
        
        for target, context in contexts:
            target_idx = token_to_idx[target]
            context_idx = token_to_idx[context]
            
            if method == "skipgram":
                # Skip-gram: predict context from target
                # Forward pass
                h = W_in[target_idx]
                u = W_out[context_idx]
                score = np.dot(h, u)
                
                # Simple gradient update (SGD)
                # For simplicity, using basic update rule
                error = learning_rate * (1.0 - score)
                W_in[target_idx] += error * u
                W_out[context_idx] += error * h
                
                total_loss += (1.0 - score) ** 2
            else:  # CBOW
                # CBOW: predict target from context
                # For simplicity, using single context
                h = W_out[context_idx]
                u = W_in[target_idx]
                score = np.dot(h, u)
                
                error = learning_rate * (1.0 - score)
                W_out[context_idx] += error * u
                W_in[target_idx] += error * h
                
                total_loss += (1.0 - score) ** 2
    
    # Return input embeddings as final representations
    # (Can also use average of input and output, or just output)
    embeddings = {}
    for token in vocab:
        idx = token_to_idx[token]
        # Use input embedding as final representation
        embeddings[token] = W_in[idx].copy()
    
    return embeddings


def sequence_embeddings(
    sequences: List[List[str]],
    event_embeddings: Dict[str, NDArray[np.float64]],
    method: str = "mean",
    temporal_weighting: bool = False
) -> NDArray[np.float64]:
    """Aggregate event sequences into fixed-size vector representations.
    
    Args:
        sequences: List of event sequences (each sequence is list of event tokens)
        event_embeddings: Dictionary mapping event tokens to embedding vectors
        method: Aggregation method ("mean", "sum", "max", "attention")
        temporal_weighting: If True, apply temporal weighting (more recent events weighted higher)
        
    Returns:
        Array of shape (n_sequences, embedding_dim) with sequence embeddings
        
    Examples:
        >>> event_embeddings = {
        ...     "health:diagnosis": np.array([0.1, 0.2]),
        ...     "occupation:job_change": np.array([0.3, 0.4]),
        ... }
        >>> sequences = [
        ...     ["health:diagnosis", "occupation:job_change"],
        ...     ["health:diagnosis"],
        ... ]
        >>> seq_embeddings = sequence_embeddings(sequences, event_embeddings, method="mean")
        >>> seq_embeddings.shape
        (2, 2)
    """
    if not sequences:
        return np.array([])
    
    # Get embedding dimension
    sample_emb = next(iter(event_embeddings.values()))
    embedding_dim = sample_emb.shape[0]
    
    result = []
    
    for seq in sequences:
        if not seq:
            # Empty sequence -> zero vector
            result.append(np.zeros(embedding_dim))
            continue
        
        # Get embeddings for events in sequence
        seq_embeddings = []
        for event_token in seq:
            if event_token in event_embeddings:
                seq_embeddings.append(event_embeddings[event_token])
            else:
                # Unknown event -> zero vector
                seq_embeddings.append(np.zeros(embedding_dim))
        
        seq_embeddings_array = np.array(seq_embeddings)
        
        if method == "mean":
            if temporal_weighting:
                # Weight by position (more recent = higher weight)
                weights = np.linspace(0.5, 1.0, len(seq_embeddings))
                weights = weights / weights.sum()
                embedding = np.average(seq_embeddings_array, axis=0, weights=weights)
            else:
                embedding = np.mean(seq_embeddings_array, axis=0)
        
        elif method == "sum":
            embedding = np.sum(seq_embeddings_array, axis=0)
        
        elif method == "max":
            embedding = np.max(seq_embeddings_array, axis=0)
        
        elif method == "attention":
            # Simple attention: uniform weights for now
            # (Can be enhanced with learned attention weights)
            embedding = np.mean(seq_embeddings_array, axis=0)
        
        else:
            raise ValueError(f"Unknown aggregation method: {method}")
        
        result.append(embedding)
    
    return np.array(result)


def domain_specific_embeddings(
    sequences: List[List[str]],
    domain_labels: List[List[str]],
    embedding_dim: int = 100,
    window_size: int = 5,
    **kwargs: Any
) -> Dict[str, Dict[str, NDArray[np.float64]]]:
    """Learn separate embeddings per domain.
    
    Creates separate embedding spaces for each domain (health, education, etc.),
    allowing domain-specific semantic relationships to be captured.
    
    Args:
        sequences: List of event sequences
        domain_labels: List of domain labels for each sequence (same structure as sequences)
        embedding_dim: Dimension of embedding vectors
        window_size: Size of context window
        **kwargs: Additional arguments passed to learn_event_embeddings
        
    Returns:
        Dictionary mapping domain names to event embedding dictionaries
        
    Examples:
        >>> sequences = [
        ...     ["diagnosis", "job_change", "raise"],
        ...     ["degree", "job_change", "move"],
        ... ]
        >>> domains = [
        ...     ["health", "occupation", "income"],
        ...     ["education", "occupation", "address"],
        ... ]
        >>> embeddings = domain_specific_embeddings(sequences, domains, embedding_dim=50)
        >>> "health" in embeddings
        True
        >>> "diagnosis" in embeddings["health"]
        True
    """
    # Group events by domain
    domain_sequences: Dict[str, List[List[str]]] = defaultdict(list)
    
    for seq, domain_seq in zip(sequences, domain_labels):
        domain_to_events: Dict[str, List[str]] = defaultdict(list)
        for event, domain in zip(seq, domain_seq):
            domain_to_events[domain].append(event)
        
        for domain, events in domain_to_events.items():
            domain_sequences[domain].append(events)
    
    # Learn embeddings per domain
    domain_embeddings = {}
    for domain, domain_seqs in domain_sequences.items():
        if domain_seqs:
            embeddings = learn_event_embeddings(
                domain_seqs,
                embedding_dim=embedding_dim,
                window_size=window_size,
                **kwargs
            )
            domain_embeddings[domain] = embeddings
    
    return domain_embeddings

