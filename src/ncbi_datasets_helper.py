import sys
from typing import List

from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets import GenomeApi as DatasetsGenomeApi

from ncbi.datasets.metadata.genome import print_assembly_metadata_by_fields
from ncbi.datasets.metadata.genome import get_assembly_metadata_by_asm_accessions
from ncbi.datasets.package import dataset
import logging
from utils import save_file_to_local_dir

# accessions: List[str] = ["GCF_000001405.39"]
# zipfile_name = "human_reference.zip"


def download_genome_data_package(
    accessions: List[str],
    filename: str,
    annotation_types_include=["RNA_FASTA", "PROT_FASTA"],
):
    # download an NCBI Datasets Genome Data Package given a list of NCBI Assembly accessions
    with DatasetsApiClient() as api_client:
        genome_api = DatasetsGenomeApi(api_client)
        try:
            logging.info("Begin download of genome data package ...")
            genome_ds_download = genome_api.download_assembly_package(
                accessions,
                # include_annotation_type=annotation_types_include,
                # _preload_content=False,
                filename=filename,
            )
            logging.info(genome_ds_download)
            logging.info(f"Download completed")
        except DatasetsApiException as e:
            raise (f"Exception when calling download_assembly_package: {e}\n")


def get_metadata_by_single_accesion(genome_assembly_accessions: List[str]):
    with DatasetsApiClient() as api_client:
        genome_api = DatasetsGenomeApi(api_client)
        genome_metadata = genome_api.assembly_descriptors_by_accessions(
            genome_assembly_accessions
        )
        # raise error if empty
        if not genome_metadata:
            logging.info("Genome not found for: %s", genome_assembly_accessions)
            return None
        return genome_metadata["assemblies"][0]


# # open the package zip archive so we can retrieve files from it
# package = dataset.AssemblyDataset(zipfile_name)
# # print the names and types of all files in the downloaded zip file
# print(package.get_catalog())

# # search by file type to get the names of all the genomic fasta files in the package
# for file_name in package.get_file_names_by_type("GENOMIC_NUCLEOTIDE_FASTA"):
#     print(file_name)

# # get the data report and print the organism name and assembly level for each genome
# for report in package.get_data_reports():
#     print(f"{report.organism_name}\t{report.assembly_info.assembly_level}")
