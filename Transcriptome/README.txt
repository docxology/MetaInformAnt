# Amalgkit RNA-Seq Analysis

This repository contains a series of shell scripts for running RNA-seq analysis using the Amalgkit pipeline. For detailed information on function use, config files, tissue/genome selection, etc., please refer to the [Amalgkit repository and wiki](https://github.com/kfuku52/amalgkit/wiki/).

## Getting Started

1. Download the shell script files from [this repository](https://github.com/docxology/apis-seq) into a new folder.

2. Run `0_setEnv.sh` to set up the environment and install dependencies:

```
bash 0_setEnv.sh
```

During the installation process, you may need to press enter or write "y" multiple times. This step may take several minutes. `0_setEnv.sh` will also create config files for later steps.

3. Run `1_metadata.sh` to create a metadata file with selected SRA and reference genome downloaded:

```
bash 1_metadata.sh
```

This exact script is also provided as 1_metadata.py , a python script.

Note there is an optional script 1_1_split_metadata.py , this splits the big metadata.tsv file into sections, which can be useful for speeding up the downloading of files in a federated environment.

4. Inspect the metadata file in `/metadata/metadata/metadata.tsv`. Check all the columns, especially the "exclusion" column. If exclusion=no, the sample is included in the next analysis.

5. Download the raw FASTQ RNA-seq reads for the target included libraries:

```
bash 2_getfastq.sh
```

This step may take a while, especially if it's your first time downloading the SRA. Each SRA can be several GB in size, so downloading dozens or hundreds of SRAs may take a long time depending on your connection. After each SRA is downloaded, `parallel-fastq-dump` performs filtering (for low-quality and duplicate reads, adapter trimming), which may also take some time depending on your processor.

6. Run `3_quant.sh`:

```
bash 3_quant.sh
```

This uses kallisto to quantify gene expression for each SRA.


7. Run `curate.sh`:

```
bash curate.sh
```

This script uses amalgkit curate (default output is "log2p1-fpkm" normalized) to harmonize expression across SRA, preparing the dataset for downstream analysis.

**Note:** If you encounter a "fastq-dump error! exit code: 3", ensure that your disk is in NTFS format to handle files over 4 GB.
