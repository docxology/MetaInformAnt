import psycopg2
import os
import logging


def get_db_client():
    try:
        logging.debug("Fetching env var")
        host = os.environ.get("PG_HOST", "localhost")
        db_name =os.environ.get("DB_NAME", "postgres")
        user = os.environ.get("DB_USER", "genome_manager")
        password = os.environ.get("DB_PASSWORD", "")

        logging.debug("Init connection to pg db")
        conn = psycopg2.connect(
        host=host,
        database= db_name,
        user=user,
        password=password)

        cur = conn.cursor()
        cur.execute('SELECT version()')
    except Exception as e:
        logging.error("Problem initializing db")
        raise e
    return cur