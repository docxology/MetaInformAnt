from __future__ import annotations

import logging
from typing import Tuple

import psycopg2

from .config import load_postgres_config_from_env


def get_db_client() -> Tuple["psycopg2.extensions.connection", "psycopg2.extensions.cursor"]:
    config = load_postgres_config_from_env()
    if config is None:
        raise RuntimeError("Postgres configuration not found in environment variables")

    logging.debug("Initializing connection to PostgreSQL database")
    conn = psycopg2.connect(
        host=config.host,
        port=config.port,
        database=config.database,
        user=config.user,
        password=config.password,
    )
    cur = conn.cursor()
    cur.execute("SELECT 1")
    return conn, cur


