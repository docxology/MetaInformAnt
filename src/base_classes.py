from dataclasses import dataclass
from ncbi_datasets_helper import (
    download_genome_data_package,
    get_metadata_by_single_accesion,
)


@dataclass
class OrganismGenome:
    """Class for tracking the genome object"""

    genome_id: str
    file_path: str
    file_contents: str = None


class NCBIGenome:

    search_term: str
    search_term_type: str

    def __init__(self, search_term, search_term_type) -> None:
        self.search_term = search_term
        self.search_type = search_term_type

    def search_router(self):
        if self.search_type == "assembly_accession_id":
            self.accesion_id = self.search_term
            self.assembly_metadata = get_metadata_by_single_accesion([self.accesion_id])
            self._set_taxon_name()
            self._set_filename()
            data = download_genome_data_package([self.accesion_id], self.filename)
        # implement other search types here

    def _set_taxon_name(self):
        self.taxon_name = self.assembly_metadata["assembly"]["org"]["sci_name"]

    def _set_intraspecies_details():
        pass

    def _set_filename(self):
        self.filename = f"{self.taxon_name, self.accesion_id}"
