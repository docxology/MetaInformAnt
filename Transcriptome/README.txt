Bioinformatics Pipeline README

This bioinformatics pipeline consists of six scripts designed to facilitate the process of downloading, processing, and analyzing transcriptomic data. The pipeline includes the following scripts:

0_setEnv.sh: Set up the environment by installing required dependencies and tools.
1_download_genome.py: Download the reference genome sequence.
2_get_metadata.py: Generate metadata for the genomic data and organize the output.
3_parallel_download.py: Download SRA files using a specified number of parallel threads and a fastp option.
4_quant.py: Perform transcript quantification using kallisto and amalgkit.
5_curate.py: Merge and curate the expression data using amalgkit.
Requirements

To run the pipeline, you need Python 3.6 or higher installed on your system, as well as the required bioinformatics tools, which will be installed by the 0_setEnv.sh script.

Getting Started

1. Set up the environment
Before running the pipeline, execute the 0_setEnv.sh script to install required dependencies and tools:

chmod +x 0_setEnv.sh ./0_setEnv.sh

2. Download the reference genome
Run the 1_download_genome.py script to download the reference genome sequence:

python3 1_download_genome.py

3. Generate and organize metadata
Use the 2_get_metadata.py script to generate metadata for the genomic data and organize the output:

python3 2_get_metadata.py

4. Download SRA files in parallel
To download SRA files using a specified number of parallel threads and a fastp option, run the 3_parallel_download.py script:

python3 3_parallel_download.py <number_threads> --fastp <fastp_option>

Replace <number_threads> with the desired number of parallel threads and <fastp_option> with either "yes" or "no" (default is yes).

5. Perform transcript quantification
Execute the 4_quant.py script to perform transcript quantification using kallisto and amalgkit:

python3 4_quant.py

6. Merge and curate expression data
Finally, run the 5_curate.py script to merge and curate the expression data using amalgkit:

python3 5_curate.py