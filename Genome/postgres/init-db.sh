#!/bin/bash
set -e

psql -v ON_ERROR_STOP=1 --username "$POSTGRES_USER" --dbname "$POSTGRES_DB" <<-EOSQL
    CREATE USER genome_manager;
    CREATE DATABASE metainformant;
    GRANT ALL PRIVILEGES ON DATABASE metainformant TO genome-manager;
    ALTER DEFAULT PRIVILEGES IN SCHEMA metainformant GRANT ALL ON TABLES TO genome_manager;
    CREATE TABLE metainformant.genomes (
        genome_id serial PRIMARY KEY,
        ncbi_id VARCHAR(255) NOT NULL,
        file_path TEXT NOT NULL,
        created_on TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
    );
EOSQL