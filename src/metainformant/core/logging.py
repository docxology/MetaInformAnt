import logging


def get_logger(name: str) -> logging.Logger:
    """Get or create a logger with default console handler.
    
    Args:
        name: Logger name (typically __name__)
        
    Returns:
        Configured logger with console handler if none exists
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            fmt="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    return logger


def setup_logger(name: str, log_file: str | None = None, level: str = "INFO") -> logging.Logger:
    """Set up a logger with file and/or console output.

    Args:
        name: Logger name
        log_file: Optional log file path
        level: Logging level

    Returns:
        Configured logger
    """
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, level.upper()))

    # Clear existing handlers
    logger.handlers.clear()

    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # File handler if specified
    if log_file:
        from pathlib import Path

        Path(log_file).parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def log_with_metadata(logger: logging.Logger, message: str, metadata: dict) -> None:
    """Log message with structured metadata.

    Args:
        logger: Logger instance
        message: Log message
        metadata: Dictionary of metadata to include
    """
    import json

    # Format metadata as JSON for structured logging
    metadata_str = json.dumps(metadata, default=str)
    full_message = f"{message} | {metadata_str}"

    logger.info(full_message)
