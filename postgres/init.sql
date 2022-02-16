CREATE USER genome_manager;
CREATE schema metainformant;
GRANT ALL PRIVILEGES ON DATABASE postgres TO genome_manager;
CREATE TABLE metainformant.genomes (
    genome_id serial PRIMARY KEY,
    ncbi_id VARCHAR(255) NOT NULL,
    file_path TEXT NOT NULL,
    created_on TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);
