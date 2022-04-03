CREATE USER genome_manager;
CREATE schema metainformant;
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";
CREATE TABLE metainformant.taxonomic_names (
    id UUID NOT NULL DEFAULT uuid_generate_v4() , 
    genus VARCHAR(255) NOT NULL,
    epithet VARCHAR(255) NOT NULL,
    created_on TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    UNIQUE(genus, epithet),
    CONSTRAINT taxonomic_names_pkey_ PRIMARY KEY (id)
);

CREATE TABLE metainformant.genome (
    assembly_accession_id VARCHAR(255) PRIMARY KEY,
    taxon_name UUID NOT NULL REFERENCES metainformant.taxonomic_names,
    infraspecies_connecting_term TEXT NULL,
    infraspecies_name TEXT NULL,
    file_path TEXT NOT NULL,
    created_on TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);


grant usage on schema metainformant to genome_manager;
GRANT ALL PRIVILEGES ON DATABASE postgres TO genome_manager;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA metainformant TO genome_manager;
ALTER DEFAULT PRIVILEGES FOR USER genome_manager IN SCHEMA metainformant GRANT SELECT, INSERT, UPDATE, DELETE ON TABLES to genome_manager;