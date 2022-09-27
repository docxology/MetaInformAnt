from dataclasses import dataclass
from ncbi_datasets_helper import (
    download_genome_data_package,
    get_metadata_by_single_accesion,
)
from utils import save_file_to_local_dir
import logging

import os
from pathlib import Path


@dataclass
class OrganismGenome:
    """Class for tracking the genome object"""

    genome_id: str
    file_path: str
    file_contents: str = None


class NCBIGenome:

    search_term: str
    search_term_type: str
    mode: str

    @staticmethod
    def _normalize_scientitfic_name(text: str):
        return text.strip().lower()

    def __init__(self, search_term, search_term_type) -> None:
        self.search_term = search_term
        self.search_type = search_term_type
        self.mode = "local"

    def search_router(self):
        if self.search_type == "assembly_accession_id":
            self.accesion_id = self.search_term
            self.assembly_metadata = get_metadata_by_single_accesion([self.accesion_id])
            logging.debug(self.assembly_metadata)
            self._set_taxon_name()
            self._set_epithet()
            self._set_tx_id()
            self._set_filename()
        else:
            raise Exception("Only search by assembly_accession_id implemented.")
            # download_genome_data_package([self.accesion_id], filename=self.filename)
        # implement other search types here

    def _set_taxon_name(self):
        self.taxon_name = self._normalize_scientitfic_name(
            self.assembly_metadata["assembly"]["org"]["sci_name"].split()[0]
        )

    def _set_epithet(self):
        self.epithet = self._normalize_scientitfic_name(
            self.assembly_metadata["assembly"]["org"]["sci_name"].split()[-1]
        )

    def _set_intraspecies_details():
        pass

    def _set_tx_id(self):
        self.tx_id = self.assembly_metadata["assembly"]["org"]["tax_id"]

    def _set_filename(self):
        taxon_underscore = f"{self.taxon_name}_{self.epithet}"
        if self.mode == "local":
            dir_path = os.path.dirname(os.path.realpath(__file__))
            dir_path = dir_path.replace("src", "test_data")
            Path(f"{dir_path}/{taxon_underscore}").mkdir(parents=True, exist_ok=True)

        self.filename = f"{dir_path}/{taxon_underscore}/{self.accesion_id}.zip"
        logging.debug(self.filename)
