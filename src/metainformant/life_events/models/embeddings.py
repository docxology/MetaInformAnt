"""Event embedding and representation learning for life events analysis.

This module provides functions for learning vector representations of life events
and event sequences using various embedding techniques.
"""

from __future__ import annotations

import zlib
from collections import Counter
from typing import Any, Dict, List, Tuple

import numpy as np

from metainformant.core.utils import logging
from metainformant.life_events.core.utils import _sequence_to_tokens

logger = logging.get_logger(__name__)


# Logits are clipped before the sigmoid, as in the reference word2vec implementation; exp(35) is far
# below the float64 overflow threshold, so no overflow warnings are raised.
_MAX_EXP = 35.0

_DEFAULT_NEGATIVE_SAMPLES = 5
_DEFAULT_LEARNING_RATE = 0.025
_DEFAULT_MIN_LEARNING_RATE = 1e-4
_DEFAULT_NS_EXPONENT = 0.75
_EMBEDDING_COMBINATIONS = ("sum", "input")

# Domain embeddings are a fixed centre per domain plus a small variation that keeps domains with the same
# centre (every domain without a dedicated centre uses the generic one) apart. The variation is seeded.
_DEFAULT_DOMAIN_SEED = 0
_DOMAIN_VARIATION_SCALE = 0.1


def _sigmoid(values: np.ndarray) -> np.ndarray:
    """Numerically stable logistic function."""
    return np.asarray(1.0 / (1.0 + np.exp(-np.clip(values, -_MAX_EXP, _MAX_EXP))))


def _take_option(options: Dict[str, Any], names: Tuple[str, ...], default: Any) -> Any:
    """Remove ``names`` from ``options`` and return the first non-None value, or ``default``.

    Earlier names take precedence over later aliases.
    """
    found = default
    for name in reversed(names):
        value = options.pop(name, None)
        if value is not None:
            found = value
    return found


def _train_word2vec(
    sentences: List[np.ndarray],
    counts: np.ndarray,
    *,
    embedding_dim: int,
    window_size: int,
    skipgram: bool,
    epochs: int,
    negative: int,
    learning_rate: float,
    min_learning_rate: float,
    ns_exponent: float,
    cbow_mean: bool,
    rng: np.random.Generator,
    verbose: bool = False,
) -> Tuple[np.ndarray, np.ndarray, List[float]]:
    """Train word2vec (skip-gram or CBOW) with negative sampling on integer-coded sequences.

    This follows Mikolov et al. (2013): vectors are updated by stochastic gradient ascent on the
    negative-sampling objective, the effective window of every position is drawn uniformly from
    ``1..window_size`` (so nearer events are used more often), negatives are drawn from the unigram
    distribution raised to ``ns_exponent``, and the learning rate decays linearly to
    ``min_learning_rate`` over the whole run.

    Args:
        sentences: One integer array of vocabulary indices per sequence.
        counts: Occurrence count of every vocabulary entry (defines the negative-sampling distribution).
        embedding_dim: Vector dimensionality.
        window_size: Maximum distance between an event and its context events.
        skipgram: Skip-gram when True (predict each context event from the centre event); otherwise
            CBOW (predict the centre event from the mean of its context events).
        epochs: Passes over all sequences.
        negative: Negative samples per positive example.
        learning_rate: Initial learning rate.
        min_learning_rate: Final learning rate.
        ns_exponent: Exponent applied to the counts for the noise distribution.
        cbow_mean: CBOW only; average (True) or sum (False) the context vectors.
        rng: Random generator driving initialisation, sequence order, window sizes and negatives.
        verbose: Print the mean loss after every epoch.

    Returns:
        ``(input_vectors, output_vectors, epoch_losses)`` where ``epoch_losses`` holds the mean
        negative-sampling loss per positive example for every epoch.
    """
    vocab_size = int(counts.shape[0])
    input_vectors = (rng.random((vocab_size, embedding_dim)) - 0.5) / embedding_dim
    output_vectors = np.zeros((vocab_size, embedding_dim))

    noise = np.asarray(counts, dtype=float) ** ns_exponent
    noise_cdf = np.cumsum(noise / noise.sum())

    labels = np.zeros(negative + 1)
    labels[0] = 1.0

    total_steps = max(1, epochs * sum(int(sentence.shape[0]) for sentence in sentences))
    step = 0
    epoch_losses: List[float] = []

    for epoch in range(epochs):
        loss_sum = 0.0
        loss_terms = 0
        for sentence_index in rng.permutation(len(sentences)):
            sentence = sentences[int(sentence_index)]
            length = int(sentence.shape[0])
            if length < 2:
                step += length
                continue

            reach = rng.integers(1, window_size + 1, size=length)
            max_context = min(2 * window_size, length - 1)
            noise_ids = np.minimum(
                np.searchsorted(noise_cdf, rng.random((length, max_context, negative)), side="right"),
                vocab_size - 1,
            )

            for position in range(length):
                alpha = max(min_learning_rate, learning_rate - (learning_rate - min_learning_rate) * step / total_steps)
                step += 1
                half = int(reach[position])
                context = np.concatenate(
                    (sentence[max(0, position - half) : position], sentence[position + 1 : position + 1 + half])
                )
                n_context = int(context.shape[0])
                if n_context == 0:
                    continue
                centre = int(sentence[position])

                if skipgram:
                    # Input: the centre vector.  Targets: every context event (label 1) and
                    # ``negative`` noise events per context event (label 0).
                    targets = np.concatenate((context[:, None], noise_ids[position, :n_context]), axis=1)
                    active = targets != context[:, None]
                    active[:, 0] = True
                    hidden = input_vectors[centre]
                    target_vectors = output_vectors[targets]
                    scores = _sigmoid(target_vectors @ hidden)
                    gradient = (labels - scores) * active * alpha
                    input_update = np.einsum("ck,ckd->d", gradient, target_vectors)
                    np.add.at(output_vectors, targets, gradient[:, :, None] * hidden)
                    input_vectors[centre] += input_update
                    loss_sum -= float(
                        np.log(np.maximum(scores[:, 0], 1e-12)).sum()
                        + (active[:, 1:] * np.log(np.maximum(1.0 - scores[:, 1:], 1e-12))).sum()
                    )
                    loss_terms += n_context
                else:
                    # Input: the mean (or sum) of the context vectors.  Targets: the centre event
                    # (label 1) and ``negative`` noise events (label 0).
                    hidden = input_vectors[context].mean(axis=0) if cbow_mean else input_vectors[context].sum(axis=0)
                    targets = np.concatenate(([centre], noise_ids[position, 0]))
                    active = targets != centre
                    active[0] = True
                    target_vectors = output_vectors[targets]
                    scores = _sigmoid(target_vectors @ hidden)
                    gradient = (labels - scores) * active * alpha
                    input_update = gradient @ target_vectors
                    if not cbow_mean:
                        input_update = input_update / n_context
                    np.add.at(output_vectors, targets, gradient[:, None] * hidden)
                    np.add.at(input_vectors, context, input_update)
                    loss_sum -= float(
                        np.log(max(float(scores[0]), 1e-12))
                        + (active[1:] * np.log(np.maximum(1.0 - scores[1:], 1e-12))).sum()
                    )
                    loss_terms += 1

        epoch_loss = loss_sum / loss_terms if loss_terms else 0.0
        epoch_losses.append(epoch_loss)
        logger.debug(f"Embedding epoch {epoch + 1}/{epochs}: loss {epoch_loss:.4f}")
        if verbose:
            print(f"Epoch {epoch + 1}/{epochs}: loss {epoch_loss:.4f}")

    return input_vectors, output_vectors, epoch_losses


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
    """Learn vector embeddings for life events with word2vec (negative sampling).

    Events that occur close together in the same sequence receive similar vectors. Training is real
    word2vec in NumPy (Mikolov et al., 2013): skip-gram or CBOW, trained by stochastic gradient
    ascent on the negative-sampling objective, with a learning rate decaying linearly over the run.

    The vector returned for an event is the sum of its input and output (context) vectors, the
    ``w + c`` combination used for GloVe vectors and studied by Levy et al. (2015). Cosine similarity
    then reflects both direct co-occurrence inside the window and shared contexts, whereas the input
    vectors alone only capture shared contexts and are close to arbitrary for events that merely
    occur next to each other. Pass ``combine="input"`` for the plain word2vec input vectors.

    Args:
        sequences: List of EventSequence objects or token lists. Co-occurrence is counted inside each
            sequence only.
        embedding_dim: Dimensionality of embeddings.
        window_size: Maximum distance (in positions) between an event and its context events.
        min_count: Events with fewer occurrences are dropped before training and are not returned.
        sg: Training algorithm when ``method`` is not given (0=CBOW, 1=Skip-gram).
        epochs: Number of passes over the sequences.
        method: ``"skipgram"`` or ``"cbow"``; takes precedence over ``sg``.
        random_state: Seed; identical seeds give identical embeddings.
        verbose: Print progress to stdout.
        **kwargs: ``seed`` (alias for ``random_state``), ``negative`` (negative samples per positive
            example, default 5), ``learning_rate`` / ``alpha`` (initial rate, default 0.025),
            ``min_learning_rate`` / ``min_alpha`` (final rate, default 1e-4), ``ns_exponent``
            (noise distribution exponent, default 0.75), ``cbow_mean`` (average rather than sum the
            context vectors in CBOW, default True) and ``combine`` (``"sum"`` or ``"input"``,
            default ``"sum"``). Other keywords are ignored with a warning.

    Returns:
        Dictionary mapping each retained event to its ``(embedding_dim,)`` float vector.

    Raises:
        ValueError: If the inputs are empty or any numeric option is out of range, or no event
            reaches ``min_count``.
    """
    if not sequences:
        raise ValueError("sequences cannot be empty")
    if embedding_dim <= 0:
        raise ValueError("embedding_dim must be positive")
    if window_size <= 0:
        raise ValueError("window_size must be positive")
    if epochs <= 0:
        raise ValueError("epochs must be positive")

    method = method or ("skipgram" if sg else "cbow")
    if method not in {"skipgram", "cbow"}:
        raise ValueError("method must be 'skipgram' or 'cbow'")

    options = dict(kwargs)
    seed = _take_option(options, ("seed",), None)
    if random_state is None:
        random_state = seed
    negative = int(_take_option(options, ("negative",), _DEFAULT_NEGATIVE_SAMPLES))
    learning_rate = float(_take_option(options, ("learning_rate", "alpha"), _DEFAULT_LEARNING_RATE))
    min_learning_rate = float(_take_option(options, ("min_learning_rate", "min_alpha"), _DEFAULT_MIN_LEARNING_RATE))
    ns_exponent = float(_take_option(options, ("ns_exponent",), _DEFAULT_NS_EXPONENT))
    cbow_mean = bool(_take_option(options, ("cbow_mean",), True))
    combine = str(_take_option(options, ("combine",), "sum"))
    if options:
        logger.warning(f"learn_event_embeddings ignores unsupported options: {sorted(options)}")
    if negative < 1:
        raise ValueError("negative must be at least 1")
    if learning_rate <= 0:
        raise ValueError("learning_rate must be positive")
    if min_learning_rate < 0:
        raise ValueError("min_learning_rate must not be negative")
    if combine not in _EMBEDDING_COMBINATIONS:
        raise ValueError(f"combine must be one of {_EMBEDDING_COMBINATIONS}")
    min_learning_rate = min(min_learning_rate, learning_rate)

    token_sequences = [_sequence_to_tokens(seq) for seq in sequences]
    event_counts: Counter[str] = Counter(token for seq in token_sequences for token in seq)
    if not event_counts:
        raise ValueError("No events found in sequences")

    min_count = int(min_count)
    vocabulary = sorted(event for event, count in event_counts.items() if count >= min_count)
    if not vocabulary:
        raise ValueError(f"No events found in sequences with at least min_count={min_count} occurrences")

    if verbose:
        print(f"Building embeddings for {len(vocabulary)} events")

    index = {event: i for i, event in enumerate(vocabulary)}
    encoded = [
        np.fromiter((index[token] for token in seq if token in index), dtype=np.int64) for seq in token_sequences
    ]
    sentences = [sentence for sentence in encoded if sentence.size > 0]
    counts = np.array([event_counts[event] for event in vocabulary], dtype=float)

    input_vectors, output_vectors, epoch_losses = _train_word2vec(
        sentences,
        counts,
        embedding_dim=embedding_dim,
        window_size=window_size,
        skipgram=method == "skipgram",
        epochs=epochs,
        negative=negative,
        learning_rate=learning_rate,
        min_learning_rate=min_learning_rate,
        ns_exponent=ns_exponent,
        cbow_mean=cbow_mean,
        rng=np.random.default_rng(random_state),
        verbose=verbose,
    )

    vectors = input_vectors + output_vectors if combine == "sum" else input_vectors
    embeddings = {event: vectors[i].copy() for i, event in enumerate(vocabulary)}

    logger.info(
        f"Embeddings learned: {len(embeddings)} events, {embedding_dim} dimensions, "
        f"{method}, final loss {epoch_losses[-1]:.4f}"
    )
    return embeddings


def biological_embedding(
    sequences: List[Any], embedding_type: str = "event_type", random_state: int | None = None, **kwargs: Any
) -> Dict[str, Any]:
    """Create biologically-informed embeddings for life events.

    The result is reproducible: the same sequences and options always give the same embeddings. Only the
    ``'domain'`` embeddings contain a random component (a small variation around each domain's centre);
    it is drawn from ``random_state`` and, when none is given, from a fixed default seed, so repeated calls
    and separate processes agree. The ``'event_type'`` and ``'temporal'`` embeddings involve no randomness.

    Args:
        sequences: List of EventSequence objects
        embedding_type: Type of biological embedding ('event_type', 'domain', 'temporal')
        random_state: Seed for the random variation of the ``'domain'`` embeddings; different seeds give
            different variations. ``seed`` is accepted as an alias. Defaults to a fixed seed.
        **kwargs: Additional parameters (``embedding_dim`` for every type)

    Returns:
        Dictionary containing embeddings and metadata. The ``'domain'`` result also records the
        ``random_state`` that was used.
    """
    logger.info(f"Creating biological embeddings: {embedding_type}")

    seed = kwargs.pop("seed", None)
    if random_state is None:
        random_state = seed

    if embedding_type == "event_type":
        return _create_event_type_embeddings(sequences, **kwargs)
    elif embedding_type == "domain":
        return _create_domain_embeddings(sequences, random_state=random_state, **kwargs)
    elif embedding_type == "temporal":
        return _create_temporal_embeddings(sequences, **kwargs)
    else:
        raise ValueError(f"Unknown embedding type: {embedding_type}")


def _create_event_type_embeddings(sequences: List[Any], **kwargs: Any) -> Dict[str, Any]:
    """Create embeddings based on event types."""
    event_types_set: set[str] = set()
    for seq in sequences:
        for event in seq.events:
            event_types_set.add(event.event_type)

    event_types = list(event_types_set)
    embedding_dim = kwargs.get("embedding_dim", 50)

    # Create simple one-hot like embeddings
    embeddings = {}
    for i, event_type in enumerate(event_types):
        # Create a unique pattern for each event type
        embedding = np.zeros(embedding_dim)
        # Use hash-like encoding. zlib.crc32 is stable across interpreter runs, unlike the
        # per-process salted built-in hash(), so the same event type always gets the same bits.
        hash_val = zlib.crc32(event_type.encode("utf-8"))
        for j in range(embedding_dim):
            embedding[j] = (hash_val >> j) & 1
        embeddings[event_type] = embedding

    return {"embeddings": embeddings, "event_types": event_types, "embedding_dim": embedding_dim, "type": "event_type"}


def _create_domain_embeddings(sequences: List[Any], random_state: int | None = None, **kwargs: Any) -> Dict[str, Any]:
    """Create embeddings based on event domains.

    Every domain draws its variation from its own generator, seeded from ``random_state`` and a stable
    hash of the domain name. A domain's embedding therefore does not depend on the iteration order of the
    domain set (which changes between processes with the string hash seed), on which other domains are
    present, or on global NumPy state.
    """
    domains_set: set[str] = set()
    for seq in sequences:
        for event in seq.events:
            if hasattr(event, "domain"):
                domains_set.add(event.domain)

    domains = sorted(domains_set) if domains_set else ["default"]
    embedding_dim = kwargs.get("embedding_dim", 20)
    seed = _DEFAULT_DOMAIN_SEED if random_state is None else int(random_state)

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
        base: Any = domain_centers.get(domain.lower(), domain_centers["default"])
        # Extend to full embedding dimension
        embedding = np.array(list(base) + [0.0] * (embedding_dim - len(base)))
        # Add small seeded variation (reproducible: see the docstring)
        variation = np.random.default_rng([seed, zlib.crc32(str(domain).encode("utf-8"))])
        embedding += variation.normal(0.0, _DOMAIN_VARIATION_SCALE, embedding_dim)
        embeddings[domain] = embedding

    return {
        "embeddings": embeddings,
        "domains": domains,
        "embedding_dim": embedding_dim,
        "type": "domain",
        "random_state": seed,
    }


def _create_temporal_embeddings(sequences: List[Any], **kwargs: Any) -> Dict[str, Any]:
    """Create embeddings based on temporal patterns."""
    # Analyze temporal patterns across sequences
    temporal_patterns: Dict[tuple[str, ...], list[Any]] = {}

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

        # Keep the events of one domain together, in order: events of the same domain that follow
        # each other inside a sequence are the training contexts for that domain's embeddings.
        per_domain: dict[str, list[str]] = {}
        for token, domain in zip(tokens, token_domains):
            per_domain.setdefault(str(domain), []).append(token)
        for domain, domain_tokens in per_domain.items():
            domain_to_sequences.setdefault(domain, []).append(domain_tokens)

    domain_embeddings: dict[str, dict[str, np.ndarray]] = {}
    for domain, domain_sequences in domain_to_sequences.items():
        domain_embeddings[domain] = learn_event_embeddings(
            domain_sequences, embedding_dim=embedding_dim, min_count=1, **kwargs
        )

    return domain_embeddings
