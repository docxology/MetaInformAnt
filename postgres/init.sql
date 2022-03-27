CREATE USER genome_manager;
CREATE schema metainformant;
grant usage on schema metainformant to genome_manager;
GRANT ALL PRIVILEGES ON DATABASE postgres TO genome_manager;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA metainformant TO genome_manager;
ALTER DEFAULT PRIVILEGES FOR USER genome_manager IN SCHEMA metainformant GRANT SELECT, INSERT, UPDATE, DELETE ON TABLES to genome_manager;
CREATE TABLE metainformant.genomes (
    genome_id serial PRIMARY KEY,
    ncbi_id VARCHAR(255) NOT NULL,
    file_path TEXT NOT NULL,
    created_on TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);
