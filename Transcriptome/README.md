# Tissue-specific gene expression meta-analysis pipeline with amalgkit

This bioinformatics pipeline consists of six scripts designed to facilitate the process of downloading, processing, and analyzing transcriptomic data. 

## Pipeline Scripts

The pipeline includes the following scripts:

1. `0_setEnv.sh`: Set up the environment by installing required dependencies and tools.
2. `1_download_genome.py`: Download the reference genome sequence.
3. `2_get_metadata.py`: Generate metadata for the genomic data and organize the output.
4. `3_parallel_download.py`: Download SRA files using a specified number of parallel threads and a `fastp` option.
5. `4_quant.py`: Perform transcript quantification using kallisto and amalgkit.
6. `5_curate.py`: Merge and curate the expression data using amalgkit.

## Requirements

To run the pipeline, you need Python 3.6 or higher installed on your system, as well as the required bioinformatics tools, which will be installed by the `0_setEnv.sh` script.

## Usage

1. **Set up the environment**: Before running the pipeline, execute the `0_setEnv.sh` script to install required dependencies and tools:
    ```bash
    chmod +x 0_setEnv.sh 
    ./0_setEnv.sh
    ```

2. **Download the reference genome**: Run the `1_download_genome.py` script to download the reference genome sequence:
    ```bash
    python 1_download_genome.py
    ```

3. **Generate and organize metadata**: Use the `2_get_metadata.py` script to generate metadata for the genomic data and organize the output:
    ```bash
    python 2_get_metadata.py
    python 2.5_update_metadata.py
    ```
    The `update_metadata.py` script updates some curate groups, and performs other updates that aim to include the maximum number of samples.

4. **Download SRA files in parallel**: To download SRA files using a specified number of parallel threads and a `fastp` option, run the `3_parallel_download.py` script:
    ```bash
    python 3_parallel_download.py <number_threads> --fastp <fastp_option> --metadata
    ```
    Replace `<number_threads>` with the desired number of parallel threads and `<fastp_option>` with either "yes" or "no" (default is yes).
    
    Example with normal metadata location, and 8 threads:
    ```bash
    python 3_parallel_download.py 8 ./metadata/metadata.tsv
    ```

5. **Perform transcript quantification**: Execute the `4_quant.py` script to perform transcript quantification using kallisto and amalgkit:
    ```bash
    python 4_quant.py
    ```

6. **Merge and curate expression data**: Finally, run the `5_curate.py` script to merge and curate the expression data using amalgkit:
    ```bash
    python 5_curate.py
    ```

7. **Analysis**: In the `Analysis` folder there are some scripts for descriptive and visual analysis. 

## Troubleshooting and Contact

If you encounter issues while running this pipeline, please use Github issues so we can resolve it. 

For further help or queries, please contact the repository owner. 

