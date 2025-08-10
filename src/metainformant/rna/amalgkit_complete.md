Repository: kfuku52/amalgkit
Files analyzed: 42

Estimated tokens: 108.0k

Directory structure:
└── kfuku52-amalgkit/
    ├── README.md
    ├── LICENSE
    ├── MANIFEST.in
    ├── setup.py
    ├── amalgkit/
    │   ├── __init__.py
    │   ├── amalgkit
    │   ├── config.py
    │   ├── csca.py
    │   ├── csca.r
    │   ├── cstmm.py
    │   ├── cstmm.r
    │   ├── curate.py
    │   ├── curate.r
    │   ├── getfastq.py
    │   ├── integrate.py
    │   ├── merge.py
    │   ├── merge.r
    │   ├── metadata.py
    │   ├── quant.py
    │   ├── sanity.py
    │   ├── select.py
    │   ├── util.py
    │   ├── util.r
    │   └── config_dir/
    │       ├── __init__.py
    │       ├── base/
    │       │   ├── __init__.py
    │       │   ├── control_term.config
    │       │   ├── exclude_keyword.config
    │       │   └── group_attribute.config
    │       ├── plantae/
    │       │   ├── __init__.py
    │       │   ├── control_term.config
    │       │   ├── exclude_keyword.config
    │       │   └── group_attribute.config
    │       ├── test/
    │       │   ├── __init__.py
    │       │   ├── control_term.config
    │       │   ├── exclude_keyword.config
    │       │   └── group_attribute.config
    │       └── vertebrate/
    │           ├── __init__.py
    │           ├── control_term.config
    │           ├── exclude_keyword.config
    │           └── group_attribute.config
    ├── logo/
    ├── util/
    │   └── batch_curate.sh
    └── .github/
        └── workflows/
            └── tag_from_version.yml


================================================
FILE: README.md
================================================
![](logo/logo_amalgkit_large.png)

## Overview
**AMALGKIT** ([/əm`ælgkit/](http://ipa-reader.xyz/?text=%C9%99m%60%C3%A6lgkit&voice=Joanna)) is a toolkit to integrate RNA-seq data from [the NCBI SRA database](https://www.ncbi.nlm.nih.gov/sra) and from private fastq files to generate unbiased cross-species transcript abundance dataset for a large-scale evolutionary gene expression analysis.

![](logo/flowchart_00.png)

## Installation
```
# Installation with pip
pip install git+https://github.com/kfuku52/amalgkit

# This should show complete options
amalgkit -h
```

## Functions
See [Wiki](https://github.com/kfuku52/amalgkit/wiki) for details.

- [`amalgkit metadata`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-metadata): NCBI SRA metadata retrieval

- [`amalgkit integrate`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-integrate): Appending local fastq info to a metadata table

- [`amalgkit config`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-config): Creating a series of config files for the metadata selection

- [`amalgkit select`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-select): Selecting SRA entries for analysis

- [`amalgkit getfastq`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-getfastq): Generating fastq files

- [`amalgkit quant`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-quant): Transcript abundance estimation

- [`amalgkit merge`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-merge): Generating transcript abundance tables

- [`amalgkit cstmm`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-cstmm): Cross-species TMM normalization using single-copy genes

- [`amalgkit curate`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-curate): Automatic removal of outlier samples and unwanted biases

- [`amalgkit csca`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-csca): Generating plots with cross-species correlation analysis

- [`amalgkit sanity`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-sanity): Checking the integrity of AMALGKIT input and output files

## Citation
Although **AMALGKIT** supports novel unpublished functions, some functionalities including metadata curation, expression level quantification, and further curation steps have been described in this paper, in which we reported the transcriptome amalgamation of 21 vertebrate species.

Fukushima K*, Pollock DD*. 2020. Amalgamated cross-species transcriptomes reveal organ-specific propensity in gene expression evolution. Nature Communications 11: 4459 (DOI: 10.1038/s41467-020-18090-8) [open access](https://www.nature.com/articles/s41467-020-18090-8)

## Licensing
**amalgkit** is BSD-licensed (3 clause). See [LICENSE](LICENSE) for details.



================================================
FILE: LICENSE
================================================
BSD 3-Clause License

Copyright (c) 2019, Kenji Fukushima
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



================================================
FILE: MANIFEST.in
================================================
recursive-include amalgkit *
#include README.md
#recursive-include /amalgkit/transcriptome_curation.R
include *.md



================================================
FILE: setup.py
================================================
import os
import re
import ast

from setuptools import setup, find_packages

with open(os.path.join('amalgkit', '__init__.py')) as f:
    match = re.search(r'__version__\s+=\s+(.*)', f.read())
version = str(ast.literal_eval(match.group(1)))

setup(
    name                    = 'amalgkit',
    version                 = version,
    description             = 'Tools for transcriptome amalgamation',
    license                 = "BSD 3-clause License",
    author                  = "Kenji Fukushima, Matthias Freund",
    author_email            = 'kfuku52@gmail.com, matthias_freund@outlook.com',
    url                     = 'https://github.com/kfuku52/amalgkit.git',
    keywords                = 'transcriptome amalgamation',
    packages                = find_packages(),
    install_requires        = ['numpy','pandas','biopython','lxml'],
    scripts                 = ['amalgkit/amalgkit',],
    include_package_data    = True
)


================================================
FILE: amalgkit/__init__.py
================================================
__version__ = '0.12.19'


================================================
FILE: amalgkit/amalgkit
================================================
#! /usr/bin/env python

import argparse
import sys
import time

from amalgkit.__init__ import __version__

print('AMALGKIT version: {}'.format(__version__))
print('AMALGKIT command: {}'.format(' '.join(sys.argv)))
print('AMALGKIT bug report: https://github.com/kfuku52/amalgkit/issues')

def strtobool(val):
    val = val.lower()
    if val in ("y", "yes", "t", "true", "on", "1"):
        return True
    elif val in ("n", "no", "f", "false", "off", "0"):
        return False
    else:
        raise ValueError(f"invalid truth value {val!r}")

def command_metadata(args):
    sys.stdout.write('amalgkit metadata: start\n')
    start = time.time()
    from amalgkit.metadata import metadata_main
    metadata_main(args)
    print('Time elapsed: {:,} sec'.format(int(time.time() - start)))
    sys.stdout.write('amalgkit metadata: end\n')

def command_select(args):
    sys.stdout.write('amalgkit select: start\n')
    start = time.time()
    from amalgkit.select import select_main
    select_main(args)
    print('Time elapsed: {:,} sec'.format(int(time.time() - start)))
    sys.stdout.write('amalgkit select: end\n')

def command_getfastq(args):
    sys.stdout.write('amalgkit getfastq: start\n')
    start = time.time()
    from amalgkit.getfastq import getfastq_main
    getfastq_main(args)
    print('Time elapsed: {:,} sec'.format(int(time.time() - start)))
    sys.stdout.write('amalgkit getfastq: end\n')

def command_quant(args):
    sys.stdout.write('amalgkit quant: start\n')
    start = time.time()
    from amalgkit.quant import quant_main
    try:
        quant_main(args)
    except ValueError as err:
        print("ERROR: ", err)
        sys.exit(1)
    print('Time elapsed: {:,} sec'.format(int(time.time() - start)))
    sys.stdout.write('amalgkit quant: end\n')

def command_cstmm(args):
    sys.stdout.write('amalgkit cstmm: start\n')
    start = time.time()
    from amalgkit.cstmm import cstmm_main
    cstmm_main(args)
    print('Time elapsed: {:,} sec'.format(int(time.time() - start)))
    sys.stdout.write('amalgkit cstmm: end\n')

def command_curate(args):
    sys.stdout.write('amalgkit curate: start\n')
    start = time.time()
    from amalgkit.curate import curate_main
    curate_main(args)
    print('Time elapsed: {:,} sec'.format(int(time.time() - start)))
    sys.stdout.write('amalgkit curate: end\n')

def command_merge(args):
    sys.stdout.write('amalgkit merge: start\n')
    start = time.time()
    from amalgkit.merge import merge_main
    merge_main(args)
    print('Time elapsed: {:,} sec'.format(int(time.time() - start)))
    sys.stdout.write('amalgkit merge: end\n')

def command_sanity(args):
    sys.stdout.write('amalgkit sanity: start\n')
    start = time.time()
    from amalgkit.sanity import sanity_main
    sanity_main(args)
    print('Time elapsed: {:,} sec'.format(int(time.time() - start)))
    sys.stdout.write('amalgkit sanity: end\n')

def command_integrate(args):
    sys.stdout.write('amalgkit integrate: start\n')
    start = time.time()
    from amalgkit.integrate import integrate_main
    integrate_main(args)
    print('Time elapsed: {:,} sec'.format(int(time.time() - start)))
    sys.stdout.write('amalgkit integrate: end\n')

def command_csca(args):
    sys.stdout.write('amalgkit csca: start\n')
    start = time.time()
    from amalgkit.csca import csca_main
    csca_main(args)
    print('Time elapsed: {:,} sec'.format(int(time.time() - start)))
    sys.stdout.write('amalgkit csca: end\n')

def command_config(args):
    sys.stdout.write('amalgkit config: start\n')
    start = time.time()
    from amalgkit.config import config_main
    config_main(args)
    print('Time elapsed: {:,} sec'.format(int(time.time() - start)))
    sys.stdout.write('amalgkit config: end\n')

def command_help(args):
    print(parser.parse_args([args.command, '--help']))


# Main parser
parser = argparse.ArgumentParser(description='A toolkit for cross-species transcriptome amalgamation')
parser.add_argument('--version', action='version', version='amalgkit version ' + __version__)
subparsers = parser.add_subparsers()

# Parent parsers
pp_meta = argparse.ArgumentParser(add_help=False)
pp_meta.add_argument('--metadata', metavar='PATH', default='inferred', type=str, required=False, action='store',
                 help='default=%(default)s: "inferred" = out_dir/metadata/metadata.tsv. '
                      'PATH to metadata table, the output file of `amalgkit metadata`.')
pp_out = argparse.ArgumentParser(add_help=False)
pp_out.add_argument('--out_dir', metavar='PATH', default='./', type=str, required=False, action='store',
                 help='default=%(default)s: PATH to the directory where intermediate and output files are generated.')
pp_batch = argparse.ArgumentParser(add_help=False)
pp_batch.add_argument('--batch', metavar='INT', default=None, type=int, required=False, action='store',
                 help='default=%(default)s: One-based index of metadata table (--metadata). '
                      'If set, process only one SRA record. This function is intended for array job processing.')
pp_redo = argparse.ArgumentParser(add_help=False)
pp_redo.add_argument('--redo', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                 help='default=%(default)s: Redo the analysis even if previous output files are detected.')
pp_threads = argparse.ArgumentParser(add_help=False)
pp_threads.add_argument('--threads', metavar='INT', default=1, type=int, required=False, action='store',
                 help='default=%(default)s: Number of threads.')
pp_sg = argparse.ArgumentParser(add_help=False)
pp_sg.add_argument('--sample_group', metavar='tissueA,tissueB,tissueC,...', default=None, type=str, required=False, action='store',
                 help='default=%(default)s: "Comma separated list of sample groups. '
                      'By default, all sample_group values in metadata.tsv are passed.')
pp_sgc = argparse.ArgumentParser(add_help=False)
pp_sgc.add_argument('--sample_group_color', metavar='#d95f02ff,#1b9e77ff,#7570b3ff,...', default='DEFAULT', type=str, required=False, action='store',
                 help='default=%(default)s: "Comma separated list of sample groups colors. '
                      'The order should be the same as --sample_group. '
                      'The number of colors should be the same as the number of selected sample groups. '
                      'By default, all colors are automatically assigned.')

# Sub parser: metadata
pme_help = 'NCBI SRA metadata retrieval and curation. See `amalgkit metadata -h`'
pme = subparsers.add_parser('metadata', help=pme_help, parents=[pp_out, pp_redo])
pme.add_argument('--search_string', metavar='PATH', default=None, type=str, required=True, action='store',
                 help='default=%(default)s: Entrez search string. See https://www.ncbi.nlm.nih.gov/books/NBK25499/ for details. '
                      'The search string is used to identify SRA entries that can be found at https://www.ncbi.nlm.nih.gov/sra/ using the same string. '
                      'Example: "Cephalotus follicularis"[Organism] AND "Illumina"[Platform] AND "RNA-seq"[Strategy]')
pme.add_argument('--entrez_email', metavar='aaa@bbb.com', default='', type=str, required=False, action='store',
                 help='default=%(default)s: Your email address. See https://www.ncbi.nlm.nih.gov/books/NBK25497/')
pme.set_defaults(handler=command_metadata)

# Sub parser: select
pse_help = 'Selecting SRA entries for analysis. See `amalgkit select -h`'
pse = subparsers.add_parser('select', help=pse_help, parents=[pp_out, pp_meta, pp_sg])
pse.add_argument('--config_dir', metavar='PATH', default='inferred', type=str, required=False, action='store',
                 help='default=%(default)s: PATH to the config directory. "inferred" = out_dir/config')
pse.add_argument('--min_nspots', metavar='INT', default=5000000, type=int, required=False, action='store',
                 help='default=%(default)s: Minimum number of RNA-seq reads per sample.')
pse.add_argument('--max_sample', metavar='INT', default=99999, type=int, required=False, action='store',
                 help='default=%(default)s: Maximum number of RNA-seq data to retain for one sample group in a species.')
pse.add_argument('--mark_redundant_biosamples', metavar='no|yes', default='no', type=strtobool,
                 required=False, action='store',
                 help='default=%(default)s: Whether to label SRAs with the same BioSample ID as unqualified.')
pse.set_defaults(handler=command_select)

# Sub parser: getfastq
pge_help = 'Retrieving fastq files. See `amalgkit getfastq -h`'
pge = subparsers.add_parser('getfastq', help=pge_help, parents=[pp_out, pp_meta, pp_threads, pp_redo, pp_batch])
pge.add_argument('--entrez_email', metavar='aaa@bbb.com', default='', type=str, required=False, action='store',
                 help='default=%(default)s: Your email address. See https://www.ncbi.nlm.nih.gov/books/NBK25497/')
pge.add_argument('--id', metavar='XXXXX0000', default=None, type=str, required=False, action='store',
                 help='default=%(default)s: BioProject/BioSample/SRR ID. This option can be used to directly specify '
                      'an ID to start FASTQ generation without running `amalgkit metadata` beforehand.')
pge.add_argument('--id_list', metavar='PATH', default=None, type=str, required=False, action='store',
                 help='default=%(default)s: Location of file containing a list of SRA IDs. Otherwise works like --id')
pge.add_argument('--layout', metavar='single|paired|auto', default='auto', type=str, required=False, action='store',
                 choices=['single', 'paired', 'auto'],
                 help='default=%(default)s: Library layout of RNA-seq data to be dumped. '
                      '"auto" prioritizes paird-end libraries if both types are available.')
pge.add_argument('--max_bp', metavar='INT', default='999,999,999,999,999', type=str, required=False, action='store',
                 help='default=%(default)s: Target sequence size (bp) to be dumped.')
pge.add_argument('--min_read_length', metavar='INT', default=25, type=int, required=False, action='store',
                 help='default=%(default)s: Minimum read length.')
pge.add_argument('--pfd', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                 help='default=%(default)s: Run parallel-fastq-dump.')
pge.add_argument('--pfd_exe', metavar='PATH', default='parallel-fastq-dump', type=str, required=False, action='store',
                 help='default=%(default)s: PATH to parallel-fastq-dump executable.')
pge.add_argument('--prefetch_exe', metavar='PATH', default='prefetch', type=str, required=False, action='store',
                 help='default=%(default)s: PATH to prefetch executable.')
pge.add_argument('--fastp', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                 help='default=%(default)s: Run fastp.')
pge.add_argument('--fastp_exe', metavar='PATH', default='fastp', type=str, required=False, action='store',
                 help='default=%(default)s: PATH to fastp executable.')
pge.add_argument('--fastp_option', metavar='STR', default='-j /dev/null -h /dev/null', type=str, required=False,
                 action='store',
                 help='default=%(default)s: Options to be passed to fastp. Do not include --length_required option here. '
                      'It can be specified throught --min_read_length in amalgkit. ')
pge.add_argument('--remove_sra', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                 help='default=%(default)s: Remove downloaded SRA files after fastq extraction.')
pge.add_argument('--remove_tmp', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                 help='default=%(default)s: Remove temporary files.')
pge.add_argument('--pfd_print', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                 help='default=%(default)s: Show parallel-fastq-dump stdout and stderr.')
pge.add_argument('--fastp_print', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                 help='default=%(default)s: Show fastp stdout and stderr.')
pge.add_argument('--sci_name', metavar='STR', default=None, type=str, required=False, action='store',
                 help='default=%(default)s: Species name in case the BioProject covers multiple species. Example: "Homo sapiens"')
pge.add_argument('--ncbi', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                 help='default=%(default)s: Download SRA files using wget from NCBI cloud, if available.')
pge.add_argument('--aws', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                 help='default=%(default)s: Download SRA files from Amazon Cloud (AWS), if available.')
pge.add_argument('--gcp', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                 help='default=%(default)s: Download SRA files from Google Cloud (GCP), if available.')
pge.add_argument('--read_name', metavar='default|trinity', default='default', type=str, required=False, action='store',
                 choices=['default', 'trinity'],
                 help='default=%(default)s: read name formatting for downstream analysis.')
pge.add_argument('--entrez_additional_search_term', metavar='STR',
                 default=None,
                # default='"platform illumina"[Properties] AND "type rnaseq"[Filter] AND "sra biosample"[Filter]',
                 type=str, required=False, action='store',
                 help='default=%(default)s: Entrez search terms in addition to --id option to further restrict the SRA entry.')
pge.add_argument('--tol', metavar='FLOAT', default=1, type=float, required=False, action='store',
                 help='default=%(default)s: Acceptable percentage loss of reads relative to --max_bp. If the 1st-round sequence '
                      'generation could not produce enough reads, the 2nd-round sequence generation is activated to '
                      'compensate the loss.')
pge.set_defaults(handler=command_getfastq)

# Sub parser: quant
pqu_help = 'Estimating transcript abundance with kallisto. See `amalgkit quant -h`'
pqu = subparsers.add_parser('quant', help=pqu_help, parents=[pp_out, pp_meta, pp_threads, pp_redo, pp_batch])

pqu.add_argument('--index_dir', metavar='PATH', default=None, type=str, required=False, action='store',
                 help='default=%(default)s: PATH to index directory. Only required if index directory is not '
                      'out_dir/index/')
pqu.add_argument('--clean_fastq', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                 help='default=%(default)s: Remove getfastq-processed fastq files when quant is successfully completed.')
pqu.add_argument('--fasta_dir', metavar='PATH', default='inferred', type=str, required=False, action='store',
                 help='default=%(default)s: "inferred" = out_dir/fasta. '
                      'PATH to directory containing reference transcriptome fasta files required for kallisto index building (see --build_index). '
                      'In this directory, file names of fasta files are expected to start with the string '
                      'in the "scientific_name" column of the metadata table, with a space replaced with an underbar. '
                      'Example: Arabidopsis_thaliana_v1.fasta for Arabidopsis thaliana.')
pqu.add_argument('--build_index', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                 help='default=%(default)s: Allows AMALGKIT to build kallisto index from reference fasta files. '
                      'Will only do this for a species if an index file is not already present. One fasta file per species should be put in --fasta_dir. '
                      'AMALGKIT will read the species from the metadata, try to find the fasta file (.fa or .fasta) and build the index for further use.')
pqu.set_defaults(handler=command_quant)

# Sub parser: merge
pmg_help = 'Generating transcript abundance tables. See `amalgkit merge -h`'
pmg = subparsers.add_parser('merge', help=pmg_help, parents=[pp_out, pp_meta])
pmg.set_defaults(handler=command_merge)

# Sub parser: cstmm
pcs_help = 'Applying cross-species TMM normalization using single-copy genes. See `amalgkit cstmm -h`'
pcs = subparsers.add_parser('cstmm', help=pcs_help, parents=[pp_out, pp_meta])
pcs.add_argument('--orthogroup_table', metavar='PATH', default=None, type=str, required=False, action='store',
                 help='default=%(default)s: PATH to orthogroup table, which is, for example, Orthogroups.tsv and N0.tsv in OrthoFinder.'
                      'Specify `--orthogroup_table ""` for single-species TMM normalization.')
pcs.add_argument('--dir_busco', metavar='PATH', default=None, type=str, required=False, action='store',
                 help='default=%(default)s: PATH to the directory where per-species BUSCO full tables are stored. '
                      'File names in this directory are expected to be GENUS_SPECIES_MISC.tsv: e.g., Arabidopsis_thaliana_full_table.tsv')
pcs.add_argument('--dir_count', metavar='PATH', default='inferred', type=str, required=False, action='store',
                 help='default=%(default)s: AMALGKIT subfolder PATH to per-species transcript abundance data as produced by `amalgkit merge`. '
                      '"inferred" = out_dir/merge')
pcs.set_defaults(handler=command_cstmm)

# Sub parser: csca
pca_help = 'Generating plots with cross-species correlation analysis. See `amalgkit csca -h`'
pca = subparsers.add_parser('csca', help=pca_help, parents=[pp_out, pp_meta, pp_sg, pp_sgc])
pca.add_argument('--orthogroup_table', metavar='PATH', default=None, type=str, required=False, action='store',
                 help='default=%(default)s: PATH to orthogroup table, which is, for example, Orthogroups.tsv and N0.tsv in OrthoFinder.'
                      'Specify `--orthogroup_table ""` for single-species TMM normalization.')
pca.add_argument('--dir_busco', metavar='PATH', default=None, type=str, required=False, action='store',
                 help='default=%(default)s: PATH to the directory where per-species BUSCO full tables are stored. '
                      'File names in this directory are expected to be GENUS_SPECIES_MISC.tsv: e.g., Arabidopsis_thaliana_full_table.tsv')
pca.add_argument('--batch_effect_alg', metavar='(no|sva|ruvseq|combatseq)', choices=['no', 'sva', 'ruvseq', 'combatseq'],
                 default='sva', type=str, required=False, action='store',
                 help='default=%(default)s: Batch-effect removal algorithm used in `amalgkit curate`.')
pca.set_defaults(handler=command_csca)

# Sub parser: curate
pcu_help = 'Automatic removal of outlier samples and SVA-based unwanted biases. See `amalgkit curate -h`'
pcu = subparsers.add_parser('curate', help=pcu_help, parents=[pp_out, pp_meta, pp_batch, pp_sg, pp_sgc, pp_redo])
pcu.add_argument('--input_dir', metavar='PATH', default='inferred', type=str, required=False, action='store',
                 help='default=%(default)s: PATH to `amalgkit merge` or `amalgkit cstmm` output folder. '
                      '"inferred" = out_dir/cstmm if exist, else out_dir/merge.')
pcu.add_argument('--dist_method', metavar='STR', default='pearson', type=str, required=False, action='store',
                 help='default=%(default)s: Method for calculating distance.')
pcu.add_argument('--mapping_rate', metavar='FLOAT', default=0.20, type=float, required=False, action='store',
                 help='default=%(default)s: Cutoff for mapping rate.')
pcu.add_argument('--correlation_threshold', metavar='FLOAT', default=0.30, type=float, required=False, action='store',
                 help='default=%(default)s: Lower cutoff for pearson r during outlier removal.')
pcu.add_argument('--plot_intermediate', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                 help='default=%(default)s: If yes, calculates and plots SVA correction after each iteration of outlier removal. Drastically increases computing times!')
pcu.add_argument('--one_outlier_per_iter', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                 help='default=%(default)s: If yes, allows curate to remove only 1 sample per same-sample-group or same-bioproject. Increases computing times!')
pcu.add_argument('--norm', metavar='(logn|log2|lognp1|log2p1|none)-(fpkm|tpm|none)',
                 default='log2p1-fpkm', type=str, required=False, action='store',
                 help='default=%(default)s: Expression level transformation before the batch effect removal. '
                      'SVA is best performed with log-transformed values. '
                      'logn: log_n normalization after FPKM/TPM transformation. '
                      'log2: log_2 normalization after FPKM/TPM transformation. '
                      'lognp1: log_n(x+1) normalization after FPKM/TPM transformation. '
                      'log2p1: log_2(x+1) normalization after FPKM/TPM transformation. '
                      'fpkm/tpm/none: FPKM, TPM, or no transformation. ')
pcu.add_argument('--batch_effect_alg', metavar='(no|sva|ruvseq|combatseq)', choices=['no', 'sva', 'ruvseq', 'combatseq'],
                 default='sva', type=str, required=False, action='store',
                 help='default=%(default)s: Batch-effect removal algorithm. '
                 'no: No batch-effect removal. '
                 'sva: Surrogate variable analysis. Use with log-transformed values. '
                 'ruvseq: Experimental. Batch effect removal based on control genes. Control genes are obtained from the residuals of a GLM. '
                 'combatseq: Experimental. Batch effect removal based on BioProject IDs. '
                 'If log-fpkm/tpm is set in combination with ruvseq or combatseq, transformation will be applied after batch-effect removal. ')
pcu.add_argument('--clip_negative', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                 help='default=%(default)s: Negative values will be clipped to 0 after the batch effect removal, '
                      'if the log*p1-* transformation is applied before the batch effect removal.')
pcu.add_argument('--maintain_zero', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                 help='default=%(default)s: Any instances of zero expression levels in the input will remain as '
                      'zero-values in the output tables, even if the process of batch effect removal causes deviation.')
pcu.add_argument('--skip_curation', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                 help='default=%(default)s: stops curate before mapping-rate based sample removal step. Outputs only '
                      'uncorrected + transformed (see --norm) count-table and the corresponding mean count-table.')
pcu.set_defaults(handler=command_curate)

# Sub parser: sanity
psa_help = 'Checking the integrity of AMALGKIT input and output files. See `amalgkit sanity -h`'
psa = subparsers.add_parser('sanity', help=psa_help, parents=[pp_out, pp_meta])
psa.add_argument('--index', required=False, action='store_true',
                 help='set this option if you want to check for availability of index files '
                      'based on species name in metadata file.')
psa.add_argument('--quant', required=False, action='store_true',
                 help='set this option if you want to check for availability quant output files '
                      'based on SRA IDs in metadata file.')
psa.add_argument('--getfastq', required=False, action='store_true',
                 help='set this option if you want to check for availability of getfastq output files '
                      'based on SRA IDs in metadata file.')
psa.add_argument('--all', required=False, action='store_true',
                 help='setting this option runs amalgkit sanity as if --index, --quant, --getfastq were set')
psa.add_argument('--index_dir', metavar='PATH', default=None, type=str, required=False, action='store',
                 help='default=%(default)s: PATH to index directory. Only required if index directory is not '
                      'out_dir/index/')
psa.add_argument('--quant_dir', metavar='PATH', default=None, type=str, required=False, action='store',
                 help='default=%(default)s: PATH to quant directory. Only required if quant directory is not '
                      'out_dir/quant/')
psa.add_argument('--getfastq_dir', metavar='PATH', default=None, type=str, required=False, action='store',
                 help='default=%(default)s: PATH to index directory. Only required if getfastq directory is not '
                      'out_dir/getfastq/')
psa.set_defaults(handler=command_sanity)

# Sub parser: integrate
pin_help = 'Appending local fastq info to a metadata table. See `amalgkit integrate -h`'
pin = subparsers.add_parser('integrate', help=pin_help, parents=[pp_out, pp_meta, pp_threads])
pin.add_argument('--fastq_dir', metavar='PATH', default=None, type=str, required=False, action='store',
                 help='default=%(default)s: PATH to input directory where fastq files are stored.')
pin.add_argument('--getfastq_dir', metavar='PATH', default=None, type=str, required=False, action='store',
                 help='default=%(default)s: PATH to index directory. Only required if getfastq directory is not '
                      'out_dir/getfastq/')
pin.add_argument('--remove_tmp', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                 help='default=%(default)s: Remove temporary files.')
pin.add_argument('--accurate_size', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                 help='default=%(default)s: ONLY APPLIES TO .gz COMPRESSED FASTQ FILES. If no, runs seqkit only on the first 1000 sequences in the fastq file to get an estimate for information like average read length. '
                      'If yes, runs seqkit on the whole fastq file. More accurate, but comes with much higher runtime.')
pin.set_defaults(handler=command_integrate)

# Sub parser: config
pco_help = 'Creating a series of config files for the metadata search. See `amalgkit config -h`'
pco = subparsers.add_parser('config', help=pco_help, parents=[pp_out])
pco.add_argument('--config', metavar='base|test|plantae|vertebrate', default='base', type=str, required=False, action='store',
                 help='default=%(default)s: Name of config dataset to be exported. Options: '
                      '"base": a minimal set of .config files for the purpose of creating custom config files. '
                      '"base_all": a complete set of near-empty .config files. '
                      '"test": short animal set for testing amalgkit metadtata. '
                      '"vertebrate" preconfigured set of config files for vertebrate animal data.'
                      '"plantae": preconfigured set of config files for plant data.')
pco.add_argument('--overwrite', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                 help='default=%(default)s: allow to overwrite config files in out_dir/config/config_name/ .')
pco.set_defaults(handler=command_config)

# Sub parser: help
parser_help = subparsers.add_parser('help', help='Printing help messages')
parser_help.set_defaults(handler=command_help)

# Handler
args = parser.parse_args()
if hasattr(args, 'handler'):
    args.handler(args)
else:
    parser.print_help()



================================================
FILE: amalgkit/config.py
================================================
import os
import sys
from glob import glob

try:
    import importlib.resources as ir
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as ir

def check_directory(args):
    path_config = os.path.join(args.out_dir, 'config_'+args.config)
    if os.path.exists(path_config):
        print('Output directory already exists: {}'.format(path_config))
        if args.overwrite:
            print('--overwrite is set to "yes". Any config files will be overwritten in: {}'.format(path_config))
        else:
            print('--overwrite is set to "no". Exiting.')
            sys.exit()
    else:
        os.makedirs(path_config)
    return

def create_config_from_package(args):
    path_config = os.path.join(args.out_dir, 'config_'+args.config)
    config_base = 'amalgkit.config_dir' + '.' + args.config
    config_files = ir.files(config_base).rglob('*.config')
    for config_file in config_files:
        file_content = ir.files(config_base).joinpath(os.path.basename(config_file)).read_bytes()
        print('Copying from {} to {}'.format(config_file, path_config))
        with open(os.path.join(path_config, os.path.basename(config_file)), mode='wb') as f:
            f.write(file_content)

def config_main(args):
    check_directory(args)
    create_config_from_package(args)



================================================
FILE: amalgkit/csca.py
================================================
from amalgkit.util import *

import re
import subprocess
import os
import shutil
import sys

def get_sample_group_string(args):
    if args.sample_group is None:
        metadata = load_metadata(args)
        sample_group = metadata.df.loc[:, 'sample_group'].dropna().unique()
    else:
        sample_group = re.findall(r"[\w]+", args.sample_group)
    print('sample_groups to be included: {}'.format(', '.join(sample_group)))
    sample_group_string = '|'.join(sample_group)
    return sample_group_string

def get_spp_from_dir(dir_curate):
    files = os.listdir(dir_curate)
    spp = [ f for f in files if (not f.startswith('.')) & (not f.startswith('tmp.')) ]
    return spp

def generate_csca_input_symlinks(dir_csca_input_table, dir_curate, spp):
    if os.path.exists(dir_csca_input_table):
        shutil.rmtree(dir_csca_input_table)
    os.makedirs(dir_csca_input_table)
    for sp in spp:
        files = os.listdir(os.path.join(dir_curate, sp, 'tables'))
        files = [ f for f in files if f.endswith('.tsv') ]
        for file in files:
            path_src = os.path.join(dir_curate, sp, 'tables', file)
            path_dst = os.path.join(dir_csca_input_table, file)
            if os.path.exists(path_dst):
                os.remove(path_dst)
            os.symlink(path_src, path_dst)
    return None

def csca_main(args):
    check_rscript()
    check_ortholog_parameter_compatibility(args)
    dir_out = os.path.realpath(args.out_dir)
    dir_curate = os.path.join(dir_out, 'curate')
    dir_csca = os.path.join(dir_out, 'csca')
    dir_csca_input_table = os.path.join(dir_csca, 'csca_input_symlinks')
    spp = get_spp_from_dir(dir_curate)
    generate_csca_input_symlinks(dir_csca_input_table, dir_curate, spp)
    sample_group_string = get_sample_group_string(args)
    if not os.path.exists(dir_csca):
        os.makedirs(dir_csca)
    if args.dir_busco is not None:
        file_orthogroup_table = os.path.join(dir_csca, 'multispecies_busco_table.tsv')
        generate_multisp_busco_table(dir_busco=args.dir_busco, outfile=file_orthogroup_table)
    elif args.orthogroup_table is not None:
        file_orthogroup_table = os.path.realpath(args.orthogroup_table)
    dir_amalgkit_script = os.path.dirname(os.path.realpath(__file__))
    csca_r_script_path = os.path.join(dir_amalgkit_script, 'csca.r')
    r_util_path = os.path.join(dir_amalgkit_script, 'util.r')
    file_genecount = os.path.join(dir_csca, 'multispecies_genecount.tsv')
    orthogroup2genecount(file_orthogroup=file_orthogroup_table, file_genecount=file_genecount, spp=spp)
    call_list = ['Rscript',
                 csca_r_script_path,
                 sample_group_string,
                 args.sample_group_color,
                 dir_out,
                 dir_csca_input_table,
                 file_orthogroup_table,
                 file_genecount,
                 r_util_path,
                 dir_csca,
                 args.batch_effect_alg,
                 ]
    print(f"Rscript command: {' '.join(call_list)}")
    subprocess.call(call_list)
    for f in glob.glob("tmp.amalgkit.*"):
        os.remove(f)



================================================
FILE: amalgkit/csca.r
================================================
#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(amap, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(colorspace, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(dendextend, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(NMF, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(MASS, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(pvclust, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(Rtsne, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(pcaMethods, quietly = TRUE)))
options(stringsAsFactors = FALSE)

debug_mode = ifelse(length(commandArgs(trailingOnly = TRUE)) == 1, "debug", "batch")
font_size = 8

if (debug_mode == "debug") {
    developer = 'kf'
    if (developer == 'mf') {
        selected_sample_groups = c('root', 'flower', 'leaf')
        dir_work = '/home/s229181/projects/amalgkit_paper/Plant_set'
        #selected_sample_groups = c('Amicoumacin_low','Amicoumacin_high','anaerobic','Azithromycin_low','Azithromycin_high','control','ciprofloxacin_low','ciprofloxacin_high','colistin_low','colistin_high','H2O2','heat','meropenem_low','meropenem_high','NaCl','Anaerobic','IPTG','biofilm_medium','acid','H202','butyrate','aerobic','Ampicillin','Vancomycin','ciprofloxacin','colistin','glucose','Glucose','rifampicin','probiotic','cleaner','ampicillin','tetracycline','pediocin','glycerol','Pyruvate','Ca2+','Glycerol','H2o2','anhydrotetracycline','TB47','treated','iron','Lactate','rion','phage','Ag','biofilm','MIC','AZM','citrate','NaNO2','Acetate','sucrose','coumermycin','copper','mitomycin','arabinose','Cefotaxime','Cellulose','vancomycin','mupirocin','galactose','macrophages','tobramycin')
        sample_group_colors = 'DEFAULT'
        #dir_work = '/home/s229181/projects/amalgkit_paper/Prokaryote_set'
        dir_csca_input_table = file.path(dir_work, 'csca/csca_input_symlinks')
        file_orthogroup = file.path(dir_work, 'csca/multispecies_busco_table.tsv')
        file_genecount = file.path(dir_work, 'csca/multispecies_genecount.tsv')
        r_util_path = '/home/s229181/projects/amalgkit_paper/amalgkit/amalgkit/util.r'
        dir_csca = file.path(dir_work, 'csca')
        batch_effect_alg = 'sva'
    } else if (developer == 'kf') {
        selected_sample_groups = c('root', 'flower', 'leaf')
        sample_group_colors = c('#d95f02ff', '#1b9e77ff', '#7570b3ff')
        dir_work = '/Users/kf/Library/CloudStorage/GoogleDrive-kenji.fukushima@nig.ac.jp/My Drive/psnl/data/evolutionary_transcriptomics/20230527_amalgkit/amalgkit_out'
        dir_csca_input_table = file.path(dir_work, 'csca/csca_input_symlinks')
        file_orthogroup = file.path(dir_work, 'csca/multispecies_busco_table.tsv')
        file_genecount = file.path(dir_work, 'csca/multispecies_genecount.tsv')
        r_util_path = '/Users/kf/Library/CloudStorage/GoogleDrive-kenji.fukushima@nig.ac.jp/My Drive/psnl/repos/amalgkit/amalgkit/util.r'
        dir_csca = file.path(dir_work, 'csca')
        batch_effect_alg = 'sva'
    }
} else if (debug_mode == "batch") {
    args = commandArgs(trailingOnly = TRUE)
    selected_sample_groups = strsplit(args[1], "\\|")[[1]]
    sample_group_colors = strsplit(args[2], ",")[[1]]
    dir_work = args[3]
    dir_csca_input_table = args[4]
    file_orthogroup = args[5]
    file_genecount = args[6]
    r_util_path = args[7]
    dir_csca = args[8]
    batch_effect_alg = args[9]
}
source(r_util_path)
setwd(dir_csca)
cat('selected_sample_groups:', selected_sample_groups, "\n")
cat('selected_sample_group_colors:', sample_group_colors, "\n")
cat('dir_work:', dir_work, "\n")
cat('dir_csca_input_table:', dir_csca_input_table, "\n")
cat('file_orthogroup:', file_orthogroup, "\n")
cat('file_genecount:', file_genecount, "\n")
cat('r_util_path:', r_util_path, "\n")
cat('dir_csca:', dir_csca, "\n")
cat('batch_effect_alg:', batch_effect_alg, "\n")

sort_labels = function(df_label, label_orders) {
    df_tmp = data.frame()
    for (lo in label_orders) {
        splits = strsplit(lo, '_')[[1]]
        scientific_name = paste(splits[1], splits[2])
        sample_group = paste(splits[3:length(splits)], collapse = '_')
        df_tmp = rbind(df_tmp, df_label[(df_label[['scientific_name']] == scientific_name) & (df_label[['sample_group']] == sample_group),])
    }
    return(df_tmp)
}

sort_averaged_tc = function(tc) {
    split_colnames = strsplit(colnames(tc), "_")
    genus_names = c()
    specific_names = c()
    sample_group_names = c()
    for (i in 1:length(split_colnames)) {
        genus_names = c(genus_names, split_colnames[[i]][1])
        specific_names = c(specific_names, split_colnames[[i]][2])
        sample_group_names = c(sample_group_names, paste0(split_colnames[[i]][3:length(split_colnames[[i]])], collapse = '_'))
    }
    colname_order = order(sample_group_names, genus_names, specific_names)
    tc = tc[, colname_order]
    return(tc)
}

color_children2parent = function(node) {
    if (length(node) != 2) {
        return(node)
    }
    if (!is.null(attributes(node[[1]])) && !is.null(attributes(node[[1]])$edgePar)) {
        child1_color = attributes(node[[1]])$edgePar[['col']]
    } else {
        child1_color = NULL
    }
    if (!is.null(attributes(node[[2]])) && !is.null(attributes(node[[2]])$edgePar)) {
        child2_color = attributes(node[[2]])$edgePar[['col']]
    } else {
        child2_color = NULL
    }
    if (is.null(child1_color) | is.null(child2_color)) {
        return(node)
    }
    if (is.na(child1_color) | is.na(child2_color)) {
        return(node)
    }
    if (child1_color == child2_color) {
        attributes(node)$edgePar[['col']] = child1_color
    }
    return(node)
}

map_color = function(redundant_variables, c) {
    uniq_var = unique(redundant_variables)
    uniq_col = rainbow_hcl(length(uniq_var), c = c)
    df_unique = data.frame(var = uniq_var, col = uniq_col, stringsAsFactors = FALSE)
    df_redundant = data.frame(var = redundant_variables, order = seq(1, length(redundant_variables)), stringsAsFactors = FALSE)
    df_redundant = merge(df_redundant, df_unique, by = "var", all.x = TRUE, stringsAsFactors = FALSE)
    df_redundant = df_redundant[order(df_redundant[['order']]),]
    return(df_redundant[['col']])
}

draw_multisp_heatmap = function(tc, df_label) {
    tc_dist_matrix = cor(tc, method = 'pearson')
    tc_dist_matrix[is.na(tc_dist_matrix)] = 0
    ann_label = df_label[, c('scientific_name', 'sample_group')]
    colnames(ann_label) = c('species', 'sample_group')
    is_sp_first_appearance = (!duplicated(df_label[['sp_color']]))
    sp_color = df_label[is_sp_first_appearance, 'sp_color']
    sp_color = sp_color[order(df_label[is_sp_first_appearance, 'scientific_name'])]
    is_cg_first_appearance = (!duplicated(df_label[['sample_group_color']]))
    sample_group_color = df_label[is_cg_first_appearance, 'sample_group_color']
    sample_group_color = sample_group_color[order(df_label[is_cg_first_appearance, 'sample_group'])]
    ann_color = list(species = sp_color, sample_group = sample_group_color)
    breaks = c(0, seq(0.3, 1, 0.01))
    NMF::aheatmap(tc_dist_matrix, color = "-RdYlBu2:71", Rowv = NA, Colv = NA, revC = TRUE, legend = TRUE, breaks = breaks,
                  annCol = ann_label, annRow = ann_label, annColors = ann_color, annLegend = FALSE, labRow = NA, labCol = NA)
}

draw_multisp_dendrogram = function(tc, df_label, df_metadata, nboot, cex.xlab, cex.yaxis, pvclust_file = 'pvclust.RData') {
    colnames(tc) = sub("_.*", "", sub('_', ' ', colnames(tc)))
    dist_fun = function(x) { Dist(t(x), method = 'pearson') }
    if (file.exists(pvclust_file)) {
        if (file.info(pvclust_file)$size) {
            cat('pvclust file found.\n')
            load(pvclust_file)
        }
    } else {
        cat('No pvclust file found. Start bootstrapping.\n')
        result = pvclust(tc, method.dist = dist_fun, method.hclust = "average", nboot = nboot, parallel = FALSE) # UPGMA
        save(result, file = pvclust_file)
    }
    dend = as.dendrogram(result)
    dend_colors = df_label[order.dendrogram(dend), 'sample_group_color']
    label_colors = df_label[order.dendrogram(dend), 'sp_color']
    labels_colors(dend) = label_colors
    dend_labels <- df_metadata[order.dendrogram(dend), 'run']
    dend <- color_branches(dend, labels = dend_labels, col = dend_colors)
    dend <- set(dend, "branches_lwd", 2)
    for (i in 1:ncol(tc)) {
        dend = dendrapply(dend, color_children2parent)
    }
    par(cex = cex.xlab)
    plot(dend, las = 1, yaxt = 'n')
    text(result, print.num = FALSE, cex = 1, col.pv = 'black')
    par(cex = cex.yaxis)
    axis(2, las = 1)
    mtext(text = 'Distance', side = 2, line = 4, cex = cex.yaxis)

    n = ncol(tc)
    f = 100
    sample_group_unique = unique(df_metadata['sample_group'])
    sp_unique = unique(df_metadata[['scientific_name']])
    bp_unique = unique(df_metadata[['bioproject']])
    sample_group_color_unique = unique(df_metadata[['sample_group_color']])
    sp_color_unique = unique(df_metadata[['sp_color']])
    bp_color_unique = unique(df_metadata[['bp_color']])
    legend_text = c(as.character(sample_group_unique), "", as.character(sp_unique), "", as.character(bp_unique))
    legend_bg = c(sample_group_color_unique, "white", sp_color_unique, "white", bp_color_unique)
    legend_fg = c(rep("black", length(sample_group_color_unique)), "white", rep("black", length(sp_color_unique)), "white", rep("black", length(bp_color_unique)))
    #plot.new() ; par(mar=c(0,0,0,0))
    #legend("center", legend=legend_text, cex=1, pch=22, lty=0, lwd=1, pt.bg=legend_bg, col=legend_fg)
}

draw_multisp_pca = function(tc, df_label) {
    tc_dist_matrix = cor(tc, method = 'pearson')
    tc_dist_matrix[is.na(tc_dist_matrix)] = 0
    set.seed(1)
    pca = prcomp(tc_dist_matrix)
    xlabel = paste0("PC 1 (", round(summary(pca)$importance[2, 1] * 100, digits = 1), "%)")
    ylabel = paste0("PC 2 (", round(summary(pca)$importance[2, 2] * 100, digits = 1), "%)")
    plot(
        pca[['x']][, 1],
        pca[['x']][, 2],
        pch = 21,
        cex = 2,
        lwd = 1,
        bg = df_label[['sample_group_color']],
        col = df_label[['sp_color']],
        xlab = xlabel,
        ylab = ylabel,
        las = 1
    )
}

draw_multisp_mds = function(tc, df_label) {
    tc_dist_dist = Dist(t(tc), method = 'pearson') + 0.000000001
    tc_dist_dist[is.na(tc_dist_dist)] = 1
    set.seed(1)
    try_out = tryCatch(
    { isoMDS(tc_dist_dist, k = 2, maxit = 100) },
    error = function(a) { return("MDS failed.") }
    )
    if (mode(try_out) == "character") {
        cat('MDS failed.\n')
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    } else {
        mds <- try_out
        plot(mds$points[, 1], mds$points[, 2], pch = 21, cex = 2, lwd = 1, bg = df_label$sample_group_color, col = df_label$sp_color, xlab = "MDS dimension 1", ylab = "MDS dimension 2", las = 1)
    }
}

draw_multisp_tsne = function(tc, df_label) {
    perplexity = min(30, floor(ncol(tc) / 4), na.rm = TRUE)
    set.seed(1)
    out_tsne = Rtsne(as.matrix(t(tc)), theta = 0, check_duplicates = FALSE, verbose = FALSE, perplexity = perplexity, dims = 2)
    try_out = tryCatch(
    {
        plot(out_tsne$Y[, 1], out_tsne$Y[, 2], pch = 21, cex = 2, lwd = 1, bg = df_label$sample_group_color, col = df_label$sp_color,
             xlab = "t-SNE dimension 1", ylab = "t-SNE dimension 2", las = 1)
    },
        error = function(a) { return("t-SNE plot failed.") }
    )
    if (mode(try_out) == "character") {
        cat('t-SNE failed.\n')
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    }
}

draw_multisp_legend = function(df_label) {
    cex_axis = 0.7
    sample_group_unique = df_label$sample_group[!duplicated(df_label$sample_group)]
    sp_unique = df_label$scientific_name[!duplicated(df_label$scientific_name)]
    sample_group_color_unique = df_label$sample_group_color[!duplicated(df_label$sample_group_color)]
    sp_color_unique = df_label$sp_color[!duplicated(df_label$sp_color)]
    toumei = rgb(1, 1, 1, 0)
    legend_text = c('Tissue', as.character(sample_group_unique), "", 'Species', as.character(sp_unique))
    legend_bg = c(toumei, sample_group_color_unique, toumei, toumei, rep(toumei, length(sp_color_unique)))
    legend_fg = c(toumei, rep(toumei, length(sample_group_color_unique)), toumei, toumei, sp_color_unique)
    legend_pch = c(1, rep(21, length(sample_group_color_unique)), 1, 1, rep(1, length(sp_color_unique)))
    legend_font = c(2, rep(1, length(sample_group_color_unique)), 1, 2, rep(3, length(sp_color_unique)))
    plot.new()
    legend("right", legend = legend_text, pt.cex = 1, pch = legend_pch, lty = 0, lwd = 2, pt.bg = legend_bg, col = legend_fg, cex = cex_axis, text.font = legend_font)
}

prepare_metadata_table = function(dir_csca_input_table, selected_sample_groups, spp) {
    files = list.files(dir_csca_input_table, pattern = ".*metadata.*")
    df_metadata = data.frame()
    for (file in files) {
        metadata_path = file.path(dir_csca_input_table, file)
        tmp_metadata = read.table(metadata_path, header = TRUE, sep = '\t', quote = '', comment.char = '', check.names = FALSE)
        df_metadata = rbind(df_metadata, tmp_metadata)
    }
    df_metadata = df_metadata[(df_metadata[['sample_group']] %in% selected_sample_groups) & (df_metadata[['scientific_name']] %in% spp),]
    df_metadata = df_metadata[, !startsWith(colnames(df_metadata), 'Unnamed')]
    return(df_metadata)
}

get_label_orders = function(df_metadata) {
    order_cg = order(df_metadata[['sample_group']])
    label_orders = unique(paste(df_metadata[order_cg, 'scientific_name'], df_metadata[order_cg, 'sample_group'], sep = '_'))
    label_orders = sub(' ', '_', label_orders)
    return(label_orders)
}

extract_ortholog_mean_expression_table = function(df_singleog, averaged_tcs, label_orders) {
    averaged_orthologs = list()
    rowname_order = rownames(df_singleog)
    for (d in c('uncorrected', 'corrected')) {
        averaged_orthologs[[d]] = df_singleog
        averaged_orthologs[[d]][, 'busco_id'] = rownames(averaged_orthologs[[d]])
        for (sp_filled in colnames(df_singleog)) {
            tc = averaged_tcs[[d]][[sp_filled]]
            averaged_orthologs[[d]] = merge(averaged_orthologs[[d]], tc, by.x = sp_filled, by.y = "row.names", all.x = TRUE, all.y = FALSE, sort = FALSE)
        }
        num_remove_col = ncol(df_singleog) + 1
        rownames(averaged_orthologs[[d]]) = averaged_orthologs[[d]][['busco_id']]
        averaged_orthologs[[d]] = averaged_orthologs[[d]][rowname_order,]
        averaged_orthologs[[d]] = averaged_orthologs[[d]][, -(1:num_remove_col)]
        averaged_orthologs[[d]] = sort_averaged_tc(averaged_orthologs[[d]])
        available_label_orders = label_orders[label_orders %in% colnames(averaged_orthologs[[d]])]
        averaged_orthologs[[d]] = averaged_orthologs[[d]][, available_label_orders]
    }
    cat(nrow(averaged_orthologs[[d]]), 'orthologs were found before filtering.\n')
    return(averaged_orthologs)
}

load_unaveraged_expression_tables = function(dir_csca_input_table, spp_filled, batch_effect_alg) {
    unaveraged_tcs = list()
    unaveraged_tcs[['uncorrected']] = list()
    unaveraged_tcs[['corrected']] = list()
    all_files = list.files(dir_csca_input_table, pattern = "*.tc.tsv")
    uncorrected_files = all_files[grepl("uncorrected", all_files)]
    corrected_files = all_files[((!grepl("uncorrected", all_files)) & (grepl(batch_effect_alg, all_files)))]
    for (sp in spp_filled) {
        uncorrected_file = uncorrected_files[startsWith(uncorrected_files, sp)]
        uncorrected_path = file.path(dir_csca_input_table, uncorrected_file)
        corrected_file = corrected_files[startsWith(corrected_files, sp)]
        corrected_path = file.path(dir_csca_input_table, corrected_file)
        if ((length(uncorrected_path) == 0) | (length(corrected_path) == 0)) {
            cat(paste("Skipping. `amalgkit curate` output(s) not found:", sp, "\n"), file = stderr())
            next
        }
        unaveraged_tcs[['uncorrected']][[sp]] = read.delim(uncorrected_path, header = TRUE, row.names = 1, sep = '\t', check.names = FALSE)
        unaveraged_tcs[['corrected']][[sp]] = read.delim(corrected_path, header = TRUE, row.names = 1, sep = '\t', check.names = FALSE)
    }
    return(unaveraged_tcs)
}

extract_ortholog_unaveraged_expression_table = function(df_singleog, unaveraged_tcs) {
    unaveraged_orthologs = list()
    rowname_order = rownames(df_singleog)
    for (d in c('uncorrected', 'corrected')) {
        unaveraged_orthologs[[d]] = df_singleog
        unaveraged_orthologs[[d]][, 'busco_id'] = rownames(unaveraged_orthologs[[d]])
        for (sp_filled in colnames(df_singleog)) {
            tc = unaveraged_tcs[[d]][[sp_filled]]
            colnames(tc) = paste(sp_filled, colnames(tc), sep = '_')
            unaveraged_orthologs[[d]] = merge(unaveraged_orthologs[[d]], tc, by.x = sp_filled, by.y = "row.names", all.x = TRUE, all.y = FALSE, sort = FALSE)
        }
        num_remove_col = length(spp) + 1
        rownames(unaveraged_orthologs[[d]]) = unaveraged_orthologs[[d]][['busco_id']]
        unaveraged_orthologs[[d]] = unaveraged_orthologs[[d]][rowname_order,]
        unaveraged_orthologs[[d]] = unaveraged_orthologs[[d]][, -(1:num_remove_col)]
        tc_order = order(sub('.*_', '', colnames(unaveraged_orthologs[[d]])))
        unaveraged_orthologs[[d]] = unaveraged_orthologs[[d]][, tc_order]
    }
    return(unaveraged_orthologs)
}

get_df_labels_averaged = function(df_metadata, label_orders, selected_sample_groups, sample_group_colors) {
    metadata_tmp = df_metadata[(df_metadata[['exclusion']] == 'no'),]
    df_label = unique(metadata_tmp[, c('scientific_name', 'sample_group')])
    categories = list(scientific_name = metadata_tmp[['scientific_name']], sample_group = metadata_tmp[['sample_group']])
    df_bp = aggregate(metadata_tmp[['bioproject']], by = categories, function(x) { length(unique(x)) })
    colnames(df_bp) = c('scientific_name', 'sample_group', 'num_bp')
    df_label = merge(df_label, df_bp, all.x = TRUE, all.y = FALSE)
    df_run = aggregate(metadata_tmp[['run']], by = categories, function(x) { length(unique(x)) })
    colnames(df_run) = c('scientific_name', 'sample_group', 'num_run')
    df_label = merge(df_label, df_run, all.x = TRUE, all.y = FALSE)
    df_label = df_label[order(df_label[['sample_group']], df_label[['scientific_name']]),]
    df_label = sort_labels(df_label, label_orders)
    df_label = add_color_to_metadata(df_label, selected_sample_groups, sample_group_colors)
    df_label = sort_labels(df_label, label_orders)
    rownames(df_label) = NULL
    write.table(df_label, paste0('csca_color_averaged.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
    return(df_label)
}

get_df_labels_unaveraged = function(df_metadata, selected_sample_groups, sample_group_colors) {
    cols = c('run', 'bioproject', 'sample_group', 'scientific_name', 'sp_color', 'sample_group_color', 'bp_color')
    metadata_tmp = df_metadata[(df_metadata[['exclusion']] == 'no'),]
    df_color = add_color_to_metadata(metadata_tmp, selected_sample_groups, sample_group_colors)
    df_color = df_color[, cols]
    label_order = order(df_color[['run']])
    df_color = df_color[label_order,]
    write.table(df_color, paste0('csca_color_unaveraged.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
    return(df_color)
}

add_color_to_metadata = function(df, selected_sample_groups, sample_group_colors) {
    df = df[, (!colnames(df) %in% c('bp_color', 'sp_color', 'sample_group_color'))]
    scientific_name = as.character(df[['scientific_name']])
    sample_group = as.character(df[['sample_group']])
    scientific_name_unique = sort(scientific_name[!duplicated(scientific_name)])
    if (length(sample_group_colors) == 1 && sample_group_colors == 'DEFAULT') {
        if (length(selected_sample_groups) <= 8) {
            sample_group_color = brewer.pal(length(unique(sample_group)), "Dark2")
            sp_color = rainbow_hcl(length(unique(scientific_name)), c = 100)
        } else if (length(selected_sample_groups) <= 12) {
            sample_group_color = brewer.pal(length(unique(sample_group)), "Paired")
            sp_color = rainbow_hcl(length(unique(scientific_name)), c = 100)
        } else {
            sample_group_color = rainbow_hcl(length(selected_sample_groups), c = 100)
            sp_color = rainbow_hcl(length(unique(scientific_name)), c = 150)
        }
    } else {
        if (length(sample_group_colors) != length(selected_sample_groups)) {
            stop("Length of sample_group_colors must match length of selected_sample_groups")
        }
        sample_group_color = sample_group_colors
        sp_color = rainbow_hcl(length(unique(scientific_name)), c = 100)
    }
    sample_group_unique = selected_sample_groups
    df_sample_group = data.frame(sample_group = sample_group_unique, sample_group_color = sample_group_color[1:length(sample_group_unique)], stringsAsFactors = FALSE)
    df_sp = data.frame(scientific_name = scientific_name_unique, sp_color = sp_color[1:length(scientific_name_unique)], stringsAsFactors = FALSE)
    df = merge(df, df_sp, sort = FALSE, all.y = FALSE)
    df = merge(df, df_sample_group, sort = FALSE, all.y = FALSE)
    if ('bioproject' %in% colnames(df)) {
        bioproject = as.character(df$bioproject)
        bp_color = rainbow_hcl(length(unique(bioproject)), c = 50)
        df_bp = data.frame(bioproject = unique(bioproject), bp_color = bp_color[1:length(unique(bioproject))], stringsAsFactors = FALSE)
        df = merge(df, df_bp, sort = FALSE, all.y = FALSE)
    }
    return(df)
}

save_averaged_tsne_plot = function(tc, df_label) {
    cat('Generating averaged t-SNE plot.\n')
    perplexity = min(30, floor(ncol(tc) / 4), na.rm = TRUE)
    set.seed(1)
    out_tsne = try(Rtsne(as.matrix(t(tc)), theta = 0, check_duplicates = FALSE, verbose = FALSE, perplexity = perplexity, dims = 2))
    if ("try-error" %in% class(out_tsne)) {
        flag_tsne_success = FALSE
        print(out_tsne)
    } else {
        flag_tsne_success = TRUE
    }
    if (!flag_tsne_success) {
        return()
    }
    try_out = tryCatch(
    {
        plot(out_tsne$Y[, 1], out_tsne$Y[, 2], pch = 21, cex = 2, lwd = 1, bg = df_label[['sample_group_color']],
             col = df_label[['sp_color']], xlab = "t-SNE dimension 1", ylab = "t-SNE dimension 2", las = 1)
    },
        error = function(a) { return("t-SNE plot failed.") }
    )
    if (mode(try_out) == "character") {
        cat('t-SNE failed.\n')
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    }
}

get_pca_coordinates = function(tc, df_label, by = 'species_sample_group') {
    tc_dist_matrix = cor(tc, method = 'pearson')
    tc_dist_matrix[is.na(tc_dist_matrix)] = 0
    #set.seed(1)
    pca = prcomp(tc_dist_matrix)
    labels = c()
    for (i in 1:5) {
        labels = c(labels, paste0("Principal component ", i, " (", round(summary(pca)$importance[2, i] * 100, digits = 1), "%)"))
    }
    PC1 = pca[['x']][, 'PC1']
    PC2 = pca[['x']][, 'PC2']
    PC3 = pca[['x']][, 'PC3']
    PC4 = pca[['x']][, 'PC4']
    PC5 = pca[['x']][, 'PC5']
    tmp = data.frame(PC1, PC2, PC3, PC4, PC5)
    if (by == 'species_sample_group') {
        df_label[by] = paste0(sub(' ', '_', df_label[['scientific_name']]), '_', df_label[['sample_group']])
        tmp[by] = rownames(tmp)
    } else if (by == 'run') {
        tmp[by] = sub('.*_', '', rownames(tmp))
    } else {
        tmp[by] = rownames(tmp)
    }
    tmp = merge(df_label, tmp, by = by)
    return(list(tmp, labels))
}

save_unaveraged_pca_plot = function(unaveraged_orthologs, df_color_unaveraged, df_metadata) {
    cat('Generating unaveraged PCA plot.\n')
    for (d in c('uncorrected', 'corrected')) {
        out = get_pca_coordinates(tc = unaveraged_orthologs[[d]], df_label = df_color_unaveraged, by = 'run')
        tmp = out[[1]]
        pc_contributions = out[[2]]
        pc_cols = c('PC1', 'PC2', 'PC3', 'PC4', 'PC5')
        pc_cols2 = paste(pc_cols, d, sep = '_')
        sorted_cols = c(colnames(df_metadata), pc_cols2)
        tmp2 = tmp[, c('run', pc_cols)]
        colnames(tmp2) = c('run', pc_cols2)
        df_metadata = merge(df_metadata, tmp2, all.x = TRUE, by = 'run', sort = FALSE)
        df_metadata = df_metadata[, sorted_cols]

        df_color_uniq = unique(df_color_unaveraged[, c('sample_group', 'sample_group_color')])
        sample_group_colors = df_color_uniq[['sample_group_color']]
        names(sample_group_colors) = df_color_uniq[['sample_group']]

        for (pcxy in list(c(1, 2), c(3, 4))) {
            pcx = pcxy[1]
            pcy = pcxy[2]
            colx = paste0('PC', pcx)
            coly = paste0('PC', pcy)

            # Check if data is finite
            xvals = tmp[[colx]]
            yvals = tmp[[coly]]
            finite_x = xvals[is.finite(xvals)]
            finite_y = yvals[is.finite(yvals)]

            if (length(finite_x) == 0 || length(finite_y) == 0) {
                message(sprintf("Skipping PC%d vs PC%d plot for '%s' - no finite data.", pcx, pcy, d))
                next
            }

            xmin = min(finite_x)
            xmax = max(finite_x)
            ymin = min(finite_y)
            ymax = max(finite_y)

            # Check if there's any range to plot
            if (xmin == xmax || ymin == ymax) {
                message(sprintf("Skipping PC%d vs PC%d plot for '%s' - no variation in data.", pcx, pcy, d))
                next
            }

            # Add a small buffer around the min/max
            xunit = (xmax - xmin) * 0.01
            yunit = (ymax - ymin) * 0.01
            xmin = xmin - xunit
            xmax = xmax + xunit
            ymin = ymin - yunit
            ymax = ymax + yunit

            g = ggplot(tmp, aes(x = !!rlang::sym(colx), y = !!rlang::sym(coly), color = sample_group)) +
                theme_bw() +
                geom_point(size = 0.5, alpha = 0.3) +
                # Add density lines only if there is more than one unique point
                # This reduces the chance of zero contour issues.
                geom_density_2d(mapping = aes(color = sample_group), bins = 12, linewidth = 0.25) +
                scale_color_manual(values = sample_group_colors) +
                xlab(pc_contributions[pcx]) +
                ylab(pc_contributions[pcy]) +
                xlim(xmin, xmax) +
                ylim(ymin, ymax) +
                theme(
                    axis.text = element_text(size = font_size),
                    axis.title = element_text(size = font_size),
                    legend.text = element_text(size = font_size),
                    legend.title = element_text(size = font_size)
                )

            filename = paste0('csca_unaveraged_pca_PC', pcx, pcy, '.', d, '.pdf')
            tryCatch(
                {
                    ggsave(file = filename, g, height = 2.15, width = 4.25)
                },
                error = function(cond) {
                    message(paste("PCA could not be computed for file", filename))
                }
            )
        }
    }
    return(df_metadata)
}

get_tsne_coordinates = function(tc, df_label, by = 'run') {
    perplexity = min(30, floor(ncol(tc) / 4), na.rm = TRUE)
    set.seed(1)
    out_tsne = Rtsne(as.matrix(t(tc)), theta = 0, check_duplicates = FALSE, verbose = FALSE, perplexity = perplexity, dims = 2)
    tmp = data.frame(tsne1 = out_tsne$Y[, 1], tsne2 = out_tsne$Y[, 2])
    tmp[[by]] = sub('.*_', '', colnames(tc))
    tmp = merge(df_label, tmp, by = by)
    return(tmp)
}

save_unaveraged_tsne_plot = function(unaveraged_orthologs, df_color_unaveraged) {
    cat('Generating unaveraged t-SNE plot.\n')
    df_color_uniq = unique(df_color_unaveraged[, c('sample_group', 'sample_group_color')])
    sample_group_colors = df_color_uniq[['sample_group_color']]
    names(sample_group_colors) = df_color_uniq[['sample_group']]
    for (d in c('uncorrected', 'corrected')) {
        tmp = get_tsne_coordinates(tc = unaveraged_orthologs[[d]], df_label = df_color_unaveraged)
        pcx = 1
        pcy = 2
        colx = paste0('tsne', pcx)
        coly = paste0('tsne', pcy)
        xmin = min(tmp[[colx]], na.rm = TRUE)
        xmax = max(tmp[[colx]], na.rm = TRUE)
        xunit = (xmax - xmin) * 0.01
        xmin = xmin - xunit
        xmax = xmax + xunit

        ymin = min(tmp[[coly]], na.rm = TRUE)
        ymax = max(tmp[[coly]], na.rm = TRUE)
        yunit = (ymax - ymin) * 0.01
        ymin = ymin - yunit
        ymax = ymax + yunit

        g = ggplot(tmp, aes(x = !!rlang::sym(colx), !!rlang::sym(coly), color = sample_group))
        g = g + theme_bw()
        g = g + geom_point(size = 0.5)
        g = g + geom_density_2d(mapping = aes(color = sample_group), bins = 12, linewidth = 0.25)
        g = g + scale_color_manual(values = sample_group_colors)
        g = g + xlab('t-SNE dimension 1')
        g = g + ylab('t-SNE dimension 2')
        g = g + xlim(xmin, xmax)
        g = g + ylim(ymin, ymax)
        g = g + theme(
            axis.text = element_text(size = font_size),
            axis.title = element_text(size = font_size),
            legend.text = element_text(size = font_size),
            legend.title = element_text(size = font_size)
        )
        filename = paste0('csca_unaveraged_tsne.', d, '.pdf')
        tryCatch(
        {
            ggsave(file = filename, g, height = 2.15, width = 4.25)
        },
            error = function(cond) {
                message(paste("t-SNE could not be computed for file", filename))
                #message("original error message:")
                #message(conditionMessage(cond))
            }
        )
    }
}

save_averaged_heatmap_plot = function(averaged_orthologs, df_color_averaged) {
    cat('Generating averaged heatmap.\n')
    file_name = 'csca_SVA_heatmap.pdf'
    pdf(file_name, height = 3.3, width = 7.2) # full figure size = 9.7 x 7.2
    layout_matrix = matrix(c(
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        1, 3, 3, 3, 3, 3, 3, 3, 3, 3),
                           2, 10, byrow = TRUE)
    layout(t(layout_matrix))
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(0.27, 0.5, 'Uncorrected', srt = 0, cex = 2)
    text(0.80, 0.5, 'Corrected', srt = 0, cex = 2)
    for (d in c('uncorrected', 'corrected')) {
        tc = averaged_orthologs[[d]]
        df_label = df_color_averaged
        par(mar = c(0, 0, 0, 0))
        draw_multisp_heatmap(tc = tc, df_label = df_label)
    }
    graphics.off()
}

save_averaged_dendrogram_plot = function(averaged_orthologs, df_color_averaged) {
    cat('Generating averaged dendrogram.\n')
    file_name = 'csca_SVA_dendrogram.pdf'
    pdf(file_name, height = 2.5, width = 7.2) # full figure size = 9.7 x 7.2
    layout_matrix = matrix(c(1, 2), 2, 1, byrow = TRUE)
    layout(t(layout_matrix))
    for (d in c('uncorrected', 'corrected')) {
        tc = averaged_orthologs[[d]]
        df_label = df_color_averaged
        par(cex = 0.5, mar = c(10, 5.5, 0, 0), mgp = c(4, 0.7, 0))
        pvclust_file = paste0('csca_pvclust_', d, '.RData')
        draw_multisp_dendrogram(tc = tc, df_label = df_label, df_metadata = df_metadata, pvclust_file = pvclust_file,
                                nboot = 1, cex.xlab = 0.3, cex.yaxis = 0.5)
    }
    graphics.off()
}

save_averaged_dimensionality_reduction_summary = function(averaged_orthologs, df_color_averaged) {
    cat('Generating averaged dimensionality reduction summary.\n')
    par(cex = 1)
    file_name = 'csca_averaged_summary.pdf'
    pdf(file_name, height = 7.2, width = 7.2) # full figure size = 9.7 x 7.2
    layout_matrix = matrix(c(1, 1, 1, 4, 4, 4, 7, 7, 2, 2, 2, 5, 5, 5, 7, 7, 3, 3, 3, 6, 6, 6, 7, 7), 3, 8, byrow = TRUE)
    layout(layout_matrix)
    for (d in c('uncorrected', 'corrected')) {
        tc = averaged_orthologs[[d]]
        df_label = df_color_averaged
        par(mar = c(4, 4, 0.1, 1)); draw_multisp_pca(tc = tc, df_label = df_label)
        par(mar = c(4, 4, 0.1, 1)); draw_multisp_tsne(tc = tc, df_label = df_label)
        par(mar = c(4, 4, 0.1, 1)); draw_multisp_mds(tc = tc, df_label = df_label)
    }
    par(mar = c(0, 0, 0, 0)); draw_multisp_legend(df_label)
    graphics.off()
}

draw_multisp_boxplot = function(df_metadata, tc_dist_matrix, fontsize = 8) {
    is_same_sp = outer(df_metadata[['scientific_name']], df_metadata[['scientific_name']], function(x, y) { x == y })
    is_same_sample_group = outer(df_metadata[['sample_group']], df_metadata[['sample_group']], function(x, y) { x == y })
    plot(c(0.5, 4.5), c(0, 1), type = 'n', xlab = '', ylab = "Pearson's correlation\ncoefficient", las = 1, xaxt = 'n')
    boxplot(tc_dist_matrix[(!is_same_sp) & (!is_same_sample_group)], at = 1, add = TRUE, col = 'gray', yaxt = 'n')
    boxplot(tc_dist_matrix[(is_same_sp) & (!is_same_sample_group)], at = 2, add = TRUE, col = 'gray', yaxt = 'n')
    boxplot(tc_dist_matrix[(!is_same_sp) & (is_same_sample_group)], at = 3, add = TRUE, col = 'gray', yaxt = 'n')
    boxplot(tc_dist_matrix[(is_same_sp) & (is_same_sample_group)], at = 4, add = TRUE, col = 'gray', yaxt = 'n')
    labels = c('bw\nbw', 'bw\nwi', 'wi\nbw', 'wi\nwi')
    axis(side = 1, at = c(1, 2, 3, 4), labels = labels, padj = 0.5)
    axis(side = 1, at = 0.35, labels = 'Group\nSpecies', padj = 0.5, hadj = 1, tick = FALSE)
}

save_averaged_box_plot = function(averaged_orthologs, df_color_averaged) {
    cat('Generating averaged boxplot.\n')
    file_name = 'csca_boxplot.pdf'
    pdf(file_name, height = 3.6, width = 7.2) # full figure size = 9.7 x 7.2
    par(mfrow = c(1, 2))
    for (d in c('uncorrected', 'corrected')) {
        tc = averaged_orthologs[[d]]
        tc[tc < 0] = 0
        tc_dist_matrix = cor(tc, method = 'pearson')
        tc_dist_matrix[is.na(tc_dist_matrix)] = 0
        draw_multisp_boxplot(df_color_averaged, tc_dist_matrix, fontsize = 8)
    }
    graphics.off()
}

calculate_correlation_within_group = function(unaveraged_orthologs, averaged_orthologs, df_metadata, selected_sample_groups, dist_method = 'pearson') {
    my_fun1 = function(x) { median(x, na.rm = TRUE) }
    for (d in c('uncorrected', 'corrected')) {
        ortholog_med = data.frame(matrix(NA, nrow(averaged_orthologs[[d]]), length(selected_sample_groups)))
        colnames(ortholog_med) = selected_sample_groups
        rownames(ortholog_med) = rownames(averaged_orthologs[[d]])
        for (sample_group in selected_sample_groups) {
            is_sample_group = endsWith(colnames(averaged_orthologs[[d]]), sample_group)
            ortholog_med[, sample_group] = apply(as.data.frame(averaged_orthologs[[d]][, is_sample_group]), 1, my_fun1)
        }
        stopifnot(all(rownames(unaveraged_orthologs[[d]]) == rownames(ortholog_med)))
        target_col = paste0('within_group_cor_', d)
        nongroup_col = paste0('max_nongroup_cor_', d)
        df_metadata[, target_col] = NA
        df_metadata[, nongroup_col] = NA
        for (sp_and_run in colnames(unaveraged_orthologs[[d]])) {
            split_string = strsplit(sp_and_run, "_")[[1]]
            sp = paste(split_string[1:2], collapse = " ")
            sra_run = paste(split_string[3:length(split_string)], collapse = '_')
            is_sra = (df_metadata[['run']] == sra_run) & (df_metadata[['scientific_name']] == sp)
            if (sum(is_sra) == 0) {
                warning(paste('Sample skipped:', sp_and_run))
                next
            }
            sample_cg = df_metadata[is_sra, 'sample_group']
            sample_values = unaveraged_orthologs[[d]][, sp_and_run]
            for (sample_group in colnames(ortholog_med)) {
                med_values = ortholog_med[, sample_group]
                is_na = (is.na(sample_values) | is.na(med_values))
                sample_values2 = sample_values[!is_na]
                med_values2 = med_values[!is_na]
                cor_coef = cor(sample_values2, med_values2, method = dist_method)
                if (sample_cg == sample_group) {
                    df_metadata[is_sra, target_col] = cor_coef
                } else {
                    max_value = max(cor_coef, df_metadata[is_sra, nongroup_col], na.rm = TRUE)
                    df_metadata[is_sra, nongroup_col] = max_value
                }
            }
        }
    }
    return(df_metadata)
}

save_group_cor_histogram = function(df_metadata, font_size = 8) {
    cat('Generating unaveraged group correlation histogram.\n')

    # Columns we will iterate over
    cor_cols = c('within_group_cor_uncorrected', 'within_group_cor_corrected')
    fill_by_vars = c('sample_group', 'scientific_name')

    # Pre-filter the data to remove rows with NA/Non-finite values in these columns
    df_clean = df_metadata
    for (col in cor_cols) {
        df_clean = df_clean[is.finite(df_clean[[col]]), ]
    }

    # Limit values to avoid out-of-range warnings if needed, e.g. to [0, 1]
    # Adjust or remove this if not required.
    for (col in cor_cols) {
        df_clean = df_clean[!is.na(df_clean[[col]]) & df_clean[[col]] >= 0 & df_clean[[col]] <= 1, ]
    }

    max_count <- 0
    plot_list <- list()

    for (col in cor_cols) {
        for (fill_by in fill_by_vars) {
            tmp = df_clean[!is.na(df_clean[[col]]), ]
            g = ggplot2::ggplot(tmp) +
                geom_histogram(aes(x = !!rlang::sym(col), fill = !!rlang::sym(fill_by)),
                               position = "stack", alpha = 0.7, bins = 40, na.rm = TRUE) +
                theme_bw(base_size = font_size) +
                xlim(c(0, 1)) +
                labs(x = col, y = 'Count') +
                theme(
                    axis.text = element_text(size = font_size, color = 'black'),
                    axis.title = element_text(size = font_size, color = 'black'),
                    panel.grid.major.x = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),
                    rect = element_rect(fill = "transparent"),
                    plot.margin = unit(rep(0.1, 4), "cm")
                )
            plot_list[[paste0(col, "_", fill_by)]] <- g
        }
    }

    final_plot <- (plot_list[['within_group_cor_uncorrected_scientific_name']] +
        plot_list[['within_group_cor_corrected_scientific_name']]) /
        (plot_list[['within_group_cor_uncorrected_sample_group']] +
            plot_list[['within_group_cor_corrected_sample_group']])

    ggsave(filename = "csca_within_group_cor.pdf", plot = final_plot, width = 7.2, height = 6.0)
}

extract_selected_tc_only = function(unaveraged_tcs, df_metadata) {
    selected_runs = df_metadata[(df_metadata[['exclusion']] == 'no'), 'run']
    for (d in c('uncorrected', 'corrected')) {
        scientific_names = names(unaveraged_tcs[[d]])
        for (sci_name in scientific_names) {
            is_selected = colnames(unaveraged_tcs[[d]][[sci_name]]) %in% selected_runs
            if (sum(is_selected) == 0) {
                warning(paste('No', d, 'samples were selected:', sci_name))
                next
            }
            unaveraged_tcs[[d]][[sci_name]] = unaveraged_tcs[[d]][[sci_name]][, is_selected]
        }
    }
    return(unaveraged_tcs)
}

unaveraged2averaged = function(unaveraged_tcs, df_metadata, selected_sample_groups) {
    is_sample_groups = list()
    for (sample_group in selected_sample_groups) {
        is_sample_groups[[sample_group]] = (df_metadata[['sample_group']] == sample_group)
    }
    is_sci_names = list()
    for (sci_name in unique(df_metadata[['scientific_name']])) {
        sci_name_ub = sub(' ', '_', sci_name)
        is_sci_names[[sci_name_ub]] = (df_metadata[['scientific_name']] == sci_name)
    }
    is_not_excluded = (df_metadata[['exclusion']] == 'no')
    averaged_tcs = list()
    for (d in c('uncorrected', 'corrected')) {
        averaged_tcs[[d]] = list()
        scientific_names = names(unaveraged_tcs[[d]])
        for (sci_name in scientific_names) {
            n_rows = nrow(unaveraged_tcs[[d]][[sci_name]])
            averaged_tcs[[d]][[sci_name]] = data.frame(matrix(ncol = 0, nrow = n_rows))
            rownames(averaged_tcs[[d]][[sci_name]]) = rownames(unaveraged_tcs[[d]][[sci_name]])
            for (sample_group in selected_sample_groups) {
                is_target = (is_sample_groups[[sample_group]] &
                    is_sci_names[[sci_name]] &
                    is_not_excluded)
                target_runs = df_metadata[is_target, 'run']
                target_runs = target_runs[target_runs %in% colnames(unaveraged_tcs[[d]][[sci_name]])]
                if (length(target_runs) == 0) {
                    next
                }
                label = paste(sub(' ', '_', sci_name), sample_group, sep = '_')
                if (sum(is_target) == 1) {
                    averaged_tcs[[d]][[sci_name]][, label] = unaveraged_tcs[[d]][[sci_name]][, target_runs]
                } else {
                    averaged_tcs[[d]][[sci_name]][, label] = apply(unaveraged_tcs[[d]][[sci_name]][, target_runs], 1, mean)
                }
            }
            if (ncol(averaged_tcs[[d]][[sci_name]]) == 0) {
                averaged_tcs[[d]][[sci_name]] = NULL
            }
        }
    }
    return(averaged_tcs)
}

save_group_cor_scatter = function(df_metadata, font_size = 8) {
    cat('Generating unaveraged group correlation scatter plot.\n')
    alpha_value = 0.2
    improvement_xymin = 0.5
    improvement_xymax = 2.0

    # Compute new columns
    df_metadata[['corrected_per_uncorrected_group_cor']] = df_metadata[['within_group_cor_corrected']] / df_metadata[['within_group_cor_uncorrected']]
    df_metadata[['corrected_per_uncorrected_max_nongroup_cor']] = df_metadata[['max_nongroup_cor_corrected']] / df_metadata[['max_nongroup_cor_uncorrected']]

    # Limit values to the specified range
    df_metadata[['corrected_per_uncorrected_group_cor']] = pmax(pmin(df_metadata[['corrected_per_uncorrected_group_cor']], improvement_xymax), improvement_xymin)
    df_metadata[['corrected_per_uncorrected_max_nongroup_cor']] = pmax(pmin(df_metadata[['corrected_per_uncorrected_max_nongroup_cor']], improvement_xymax), improvement_xymin)

    # Remove rows with NA, NaN, or Inf
    df_clean <- df_metadata[is.finite(df_metadata$within_group_cor_uncorrected) &
                            is.finite(df_metadata$within_group_cor_corrected) &
                            is.finite(df_metadata$max_nongroup_cor_uncorrected) &
                            is.finite(df_metadata$max_nongroup_cor_corrected) &
                            is.finite(df_metadata$corrected_per_uncorrected_group_cor) &
                            is.finite(df_metadata$corrected_per_uncorrected_max_nongroup_cor), ]

    # Plotting
    ps = list()
    ps[[1]] = ggplot(df_clean, aes(x = max_nongroup_cor_uncorrected, y = within_group_cor_uncorrected)) +
        xlim(c(0, 1)) + ylim(c(0, 1))
    ps[[2]] = ggplot(df_clean, aes(x = max_nongroup_cor_corrected, y = within_group_cor_corrected)) +
        xlim(c(0, 1)) + ylim(c(0, 1))
    ps[[3]] = ggplot(df_clean, aes(x = within_group_cor_uncorrected, y = within_group_cor_corrected)) +
        xlim(c(0, 1)) + ylim(c(0, 1))
    ps[[4]] = ggplot(df_clean, aes(x = max_nongroup_cor_uncorrected, y = max_nongroup_cor_corrected)) +
        xlim(c(0, 1)) + ylim(c(0, 1))
    ps[[5]] = ggplot(df_clean, aes(x = corrected_per_uncorrected_max_nongroup_cor, y = corrected_per_uncorrected_group_cor)) +
        xlim(c(improvement_xymin, improvement_xymax)) + ylim(c(improvement_xymin, improvement_xymax))

    for (i in seq_along(ps)) {
        ps[[i]] = ps[[i]] +
            geom_point(alpha = alpha_value, na.rm = TRUE) +
            geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'blue') +
            theme_bw() +
            stat_density_2d(bins = 12, linewidth = 0.25, color = 'gray', na.rm = TRUE) +
            theme(
                text = element_text(size = font_size),
                axis.text = element_text(size = font_size),
                axis.title = element_text(size = font_size),
                legend.text = element_text(size = font_size),
                legend.title = element_text(size = font_size)
            )
    }

    ps[[6]] = ggplot() + theme_void()
    combined_plot = wrap_plots(ps)
    ggsave(filename = "csca_group_cor_scatter.pdf", plot = combined_plot, width = 7.2, height = 4.8)
}

write_pivot_table = function(df_metadata, unaveraged_tcs, selected_sample_groups) {
    d = 'corrected'
    is_selected_cg = (df_metadata[['sample_group']] %in% selected_sample_groups)
    is_loaded_run = FALSE
    for (sci_name_ub in names(unaveraged_tcs[[d]])) {
        sci_name = sub('_', ' ', sci_name_ub)
        is_loaded_run_sp = ((df_metadata[['run']] %in% colnames(unaveraged_tcs[[d]][[sci_name_ub]])) & (df_metadata[['scientific_name']] == sci_name))
        is_loaded_run = (is_loaded_run | is_loaded_run_sp)
    }
    is_selected = (is_selected_cg & is_loaded_run)
    tmp = df_metadata[is_selected, c('scientific_name', 'sample_group')]
    pivot_table = as.data.frame.matrix(table(tmp[['scientific_name']], tmp[['sample_group']]))
    pivot_table = cbind(data.frame(scientific_name = rownames(pivot_table)), pivot_table)
    write.table(pivot_table, 'csca_pivot_selected_samples.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
}

save_delta_pcc_plot = function(directory, plot_title) {

    first_row_means <- list()
    last_row_means <- list()
    all_data_points <- list()

    file_list <- list.files(directory, pattern = "\\.correlation_statistics\\.tsv$")
    for (filename in file_list) {
        filepath <- file.path(directory, filename)

        df <- read.table(filepath, sep = '\t', header = TRUE)

        first_row_means[[filename]] <- df[1, c('bwbw_mean', 'bwwi_mean', 'wiwi_mean', 'wibw_mean')]
        last_row_means[[filename]] <- df[nrow(df), c('bwbw_mean', 'bwwi_mean', 'wiwi_mean', 'wibw_mean')]

        all_data_points[[filename]] <- as.numeric(unlist(df[c('bwbw_mean', 'bwwi_mean', 'wiwi_mean', 'wibw_mean')]))
    }

    # Concatenate the mean values for the first and last rows into DataFrames
    first_row_means_df <- do.call(rbind, first_row_means)
    colnames(first_row_means_df) <- c('bwbw_uncorrected', 'bwwi_uncorrected', 'wiwi_uncorrected', 'wibw_uncorrected')

    last_row_means_df <- do.call(rbind, last_row_means)
    colnames(last_row_means_df) <- c('bwbw_corrected', 'bwwi_corrected', 'wiwi_corrected', 'wibw_corrected')

    combined_means_df <- cbind(first_row_means_df, last_row_means_df)

    #cat("Combined Means DataFrame:\n")
    #print(combined_means_df)

    # Calculate deltas
    combined_means_df$delta_bwbw_wibw_uncorrected <- combined_means_df$bwbw_uncorrected - combined_means_df$wibw_uncorrected
    combined_means_df$delta_bwbw_wibw_corrected <- combined_means_df$bwbw_corrected - combined_means_df$wibw_corrected
    combined_means_df$delta_wiwi_bwwi_uncorrected <- combined_means_df$bwwi_uncorrected - combined_means_df$wiwi_uncorrected
    combined_means_df$delta_wiwi_bwwi_corrected <- combined_means_df$bwwi_corrected - combined_means_df$wiwi_corrected

    # Select only the delta columns and take the absolute values
    delta_means_df <- combined_means_df[, c('delta_bwbw_wibw_uncorrected', 'delta_bwbw_wibw_corrected',
                                            'delta_wiwi_bwwi_uncorrected', 'delta_wiwi_bwwi_corrected')]
    delta_means_df <- abs(delta_means_df)
    delta_means_df <- na.omit(delta_means_df)

    #cat("\nDelta Means DataFrame:\n")
    #print(delta_means_df)

    if (nrow(delta_means_df) > 2) {
        cat("\nShapiro-Wilk Test Results:\n")
        cat("delta_bwbw_wibw_uncorrected: ", shapiro.test(delta_means_df$delta_bwbw_wibw_uncorrected)$p.value, "\n")
        cat("delta_bwbw_wibw_corrected: ", shapiro.test(delta_means_df$delta_bwbw_wibw_corrected)$p.value, "\n")
        cat("delta_wiwi_bwwi_uncorrected: ", shapiro.test(delta_means_df$delta_wiwi_bwwi_uncorrected)$p.value, "\n")
        cat("delta_wiwi_bwwi_corrected: ", shapiro.test(delta_means_df$delta_wiwi_bwwi_corrected)$p.value, "\n")

        cat("\nDependent T-Test Results:\n")
        t_test_result_bw = t.test(delta_means_df$delta_bwbw_wibw_uncorrected, delta_means_df$delta_bwbw_wibw_corrected, paired = TRUE)$p.value
        t_test_result_wi = t.test(delta_means_df$delta_wiwi_bwwi_uncorrected, delta_means_df$delta_wiwi_bwwi_corrected, paired = TRUE)$p.value
        cat("delta_bwbw_wibw_uncorrected vs delta_bwbw_wibw_corrected: ", t_test_result_bw, "\n")
        cat("delta_wiwi_bwwi_uncorrected vs delta_wiwi_bwwi_corrected: ", t_test_result_wi, "\n")

        p_label1 <- paste("p =", signif(t_test_result_bw, 3))
        p_label2 <- paste("p =", signif(t_test_result_wi, 3))
    } else {
        cat("Not enough data points to perform statistical tests. P values will not be shown.\n")
    }


    pdf(plot_title, height = 4.5, width = 4.5, fonts = "Helvetica", pointsize = 7)

    plot(c(0.5, 4.5), c(0, 0.45), type = 'n', xlab = '', ylab = expression(Delta ~ "mean PCC"), las = 1, xaxt = 'n')
    boxplot(delta_means_df$delta_bwbw_wibw_uncorrected, at = 1, add = TRUE, col = 'gray', yaxt = 'n')
    boxplot(delta_means_df$delta_bwbw_wibw_corrected, at = 2, add = TRUE, col = 'gray', yaxt = 'n')
    boxplot(delta_means_df$delta_wiwi_bwwi_uncorrected, at = 3, add = TRUE, col = 'gray', yaxt = 'n')
    boxplot(delta_means_df$delta_wiwi_bwwi_corrected, at = 4, add = TRUE, col = 'gray', yaxt = 'n')

    if (nrow(delta_means_df) > 2) {
        segments(x0 = 1, y0 = mean(delta_means_df$delta_bwbw_wibw_uncorrected) + 0.2, x1 = 2, y1 = mean(delta_means_df$delta_bwbw_wibw_uncorrected) + 0.2, col = "black", lwd = 0.5)
        segments(x0 = 3, y0 = mean(delta_means_df$delta_wiwi_bwwi_uncorrected) + 0.2, x1 = 4, y1 = mean(delta_means_df$delta_wiwi_bwwi_uncorrected) + 0.2, col = "black", lwd = 0.5)
        text(x = 1.5, y = mean(delta_means_df$delta_bwbw_wibw_uncorrected) + 0.22, labels = p_label1, xpd = TRUE, srt = 0)
        text(x = 3.5, y = mean(delta_means_df$delta_wiwi_bwwi_uncorrected) + 0.22, labels = p_label2, xpd = TRUE, srt = 0)
    }

    axis(side = 1, at = c(1, 2, 3, 4), labels = c("uncorr.\n", "corr.\n", "uncorr.\n", "corr.\n"), padj = 0.5, tick = FALSE)
    axis(side = 1, at = 0.35, labels = 'Correction\nSample Group', padj = 0.5, hadj = 1, tick = FALSE)
    axis(side = 1, at = c(1.5, 3.5), labels = c("\nbetween group", "\nwithin group"), padj = 0.5, tick = FALSE)
    graphics.off()
}

save_sample_number_heatmap <- function(df_metadata, font_size = 8, dpi = 300) {
    sampled_data <- df_metadata[df_metadata$exclusion == 'no', c('scientific_name', 'sample_group')]
    freq_table <- table(sampled_data$scientific_name, sampled_data$sample_group)
    df_sample_count <- as.data.frame(freq_table)
    names(df_sample_count) <- c('scientific_name', 'sample_group', 'Freq')
    all_scientific_names <- unique(df_metadata$scientific_name)
    all_sample_groups <- unique(df_metadata$sample_group)
    complete_grid <- expand.grid(scientific_name = all_scientific_names, sample_group = all_sample_groups, stringsAsFactors = FALSE)
    df_sample_count <- merge(complete_grid, df_sample_count, by = c("scientific_name", "sample_group"), all.x = TRUE)
    df_sample_count$Freq[is.na(df_sample_count$Freq)] <- 0
    df_sample_count$log2_Freq <- ifelse(df_sample_count$Freq > 0, log2(df_sample_count$Freq + 1), 0)
    fill_max <- max(df_sample_count$log2_Freq)
    freq_breaks <- c(0, 1, 2, 4, 8, 16, 32, 64, 128)
    freq_breaks <- freq_breaks[freq_breaks <= max(df_sample_count$Freq)]
    if (!0 %in% freq_breaks) { freq_breaks <- c(0, freq_breaks) }
    if (!1 %in% freq_breaks) { freq_breaks <- c(freq_breaks, 1) }
    freq_breaks <- sort(unique(freq_breaks))
    log2_breaks <- log2(freq_breaks + 1)
    log2_breaks <- log2_breaks[log2_breaks <= fill_max]
    n_sample_groups <- length(all_sample_groups)
    n_scientific_names <- length(all_scientific_names)
    base_width <- 5
    base_height <- 5
    per_sample_group_width <- 0.1
    per_scientific_name_height <- 0.1
    width <- base_width + (per_sample_group_width * n_sample_groups)
    height <- base_height + (per_scientific_name_height * n_scientific_names)
    p <- ggplot(data = df_sample_count, aes(x = sample_group, y = scientific_name, fill = log2_Freq)) +
        geom_tile(color = "grey80") +
        geom_text(aes(label = Freq), size = 3, color = "black") +
        scale_fill_gradientn(
            colors = c("white", "#440154", "#482878", "#3E4989", "#31688E", "#35B779", "#4DCD63", "#76D730", "#B8DE29", "#FDE725"),
            values = sort(unique(c(0, log2_breaks / fill_max, 1))),
            limits = c(0, fill_max),
            name = "# of samples",
            breaks = log2_breaks,
            labels = freq_breaks
        ) +
        xlab('') +
        ylab('') +
        scale_y_discrete(
            limits = rev(unique(df_sample_count$scientific_name))
        ) +
        theme_minimal() +
        theme(
            axis.text.y = element_text(size = font_size),
            axis.text.x = element_text(size = font_size, angle = 90, hjust = 1),
            legend.title = element_text(size = font_size),
            legend.text = element_text(size = font_size)
        )
    ggsave('csca_sample_number_heatmap.pdf', plot = p, width = width, height = height, dpi = dpi)
}

df_og = read.table(file_orthogroup, header = TRUE, sep = '\t', row.names = 1, quote = '', check.names = FALSE)
df_gc = read.table(file_genecount, header = TRUE, sep = '\t', quote = '', check.names = FALSE)
rownames(df_gc) = df_gc[['orthogroup_id']]
df_gc[, 'orthogroup_id'] = NULL
spp_filled = colnames(df_gc)

is_singlecopy = get_singlecopy_bool_index(df_gc, spp_filled)
df_singleog = df_og[is_singlecopy, spp_filled]
spp = sub('_', ' ', spp_filled)
df_metadata = prepare_metadata_table(dir_csca_input_table, selected_sample_groups, spp)
label_orders = get_label_orders(df_metadata)
df_color_averaged = get_df_labels_averaged(df_metadata, label_orders, selected_sample_groups, sample_group_colors)
df_color_unaveraged = get_df_labels_unaveraged(df_metadata, selected_sample_groups, sample_group_colors)
cat('Number of orthologs in input table:', nrow(df_og), '\n')
cat('Number of selected single-copy orthologs:', nrow(df_singleog), '\n')
cat('Number of selected species:', length(spp), '\n')
save_sample_number_heatmap(df_metadata, font_size = font_size, dpi = 300)

unaveraged_tcs = load_unaveraged_expression_tables(dir_csca_input_table, spp_filled, batch_effect_alg)
unaveraged_tcs = extract_selected_tc_only(unaveraged_tcs, df_metadata)

# if a species was skipped during load_unaveraged_expression_tables(), it will cause indexing issues down the line, due to df_singleog having more species entries than the tcs
if (length(unaveraged_tcs[['corrected']]) < length(colnames(df_singleog))) {
    df_singleog = df_singleog[, names(unaveraged_tcs[['corrected']])]
}

unaveraged_orthologs = extract_ortholog_unaveraged_expression_table(df_singleog, unaveraged_tcs)
averaged_tcs = unaveraged2averaged(unaveraged_tcs, df_metadata, selected_sample_groups)
averaged_orthologs = extract_ortholog_mean_expression_table(df_singleog, averaged_tcs, label_orders)
write_pivot_table(df_metadata, unaveraged_tcs, selected_sample_groups)

cat('Applying expression level imputation for missing orthologs.\n')
imputed_averaged_orthologs = list()
imputed_unaveraged_orthologs = list()
for (d in c('uncorrected', 'corrected')) {
    imputed_averaged_orthologs[[d]] = impute_expression(averaged_orthologs[[d]])
    imputed_unaveraged_orthologs[[d]] = impute_expression(unaveraged_orthologs[[d]])
    write_table_with_index_name(
        df = averaged_orthologs[[d]],
        file_path = paste0('csca_ortholog_averaged.', d, '.tsv'),
        index_name = 'target_id',
        sort = FALSE
    )
    write_table_with_index_name(
        df = unaveraged_orthologs[[d]],
        file_path = paste0('csca_ortholog_unaveraged.', d, '.tsv'),
        index_name = 'target_id',
        sort = FALSE
    )
    write_table_with_index_name(
        df = imputed_averaged_orthologs[[d]],
        file_path = paste0('csca_ortholog_averaged.imputed.', d, '.tsv'),
        index_name = 'target_id',
        sort = FALSE
    )
    write_table_with_index_name(
        df = imputed_unaveraged_orthologs[[d]],
        file_path = paste0('csca_ortholog_unaveraged.imputed.', d, '.tsv'),
        index_name = 'target_id',
        sort = FALSE
    )
}
cat(nrow(imputed_unaveraged_orthologs[[d]]), 'orthologs were found after filtering and imputation.\n')
df_metadata = calculate_correlation_within_group(unaveraged_orthologs, averaged_orthologs, df_metadata, selected_sample_groups)
save_group_cor_scatter(df_metadata, font_size = 8)
save_group_cor_histogram(df_metadata, font_size = 8)
save_averaged_tsne_plot(tc = imputed_unaveraged_orthologs[['corrected']], df_label = df_color_unaveraged)
save_averaged_heatmap_plot(imputed_averaged_orthologs, df_color_averaged)
save_averaged_dendrogram_plot(imputed_averaged_orthologs, df_color_averaged)
save_averaged_dimensionality_reduction_summary(imputed_averaged_orthologs, df_color_averaged)
save_averaged_box_plot(imputed_averaged_orthologs, df_color_averaged)
df_metadata = save_unaveraged_pca_plot(imputed_unaveraged_orthologs, df_color_unaveraged, df_metadata)
save_unaveraged_tsne_plot(imputed_unaveraged_orthologs, df_color_unaveraged)
save_delta_pcc_plot(directory = dir_csca_input_table, plot_title = 'csca_delta_pcc_boxplot.pdf')

file_metadata_out = file.path(dir_csca, 'metadata.tsv')
write.table(df_metadata, file_metadata_out, row.names = FALSE, sep = '\t', quote = FALSE)

cat(sprintf('Number of SRA samples for exclusion potting: %s\n', formatC(nrow(df_metadata), format = 'd', big.mark = ',')))
out_path = file.path(dir_csca, 'csca_exclusion.pdf')
save_exclusion_plot(df = df_metadata, out_path = out_path, font_size = 8)

if (file.exists('Rplots.pdf')) {
    file.remove('Rplots.pdf')
}
cat('csca.r completed!\n')



================================================
FILE: amalgkit/cstmm.py
================================================
import pandas

import glob
import subprocess
import os
import re
import sys

from amalgkit.util import *

def get_count_files(dir_count):
    sciname_dirs = os.listdir(dir_count)
    count_files = list()
    for sciname_dir in sciname_dirs:
        sciname_path = os.path.join(dir_count, sciname_dir)
        if not os.path.isdir(sciname_path):
            continue
        files = os.listdir(sciname_path)
        for file in files:
            sciname_file = os.path.join(sciname_path, file)
            if not os.path.isfile(sciname_file):
                continue
            if not file.endswith('est_counts.tsv'):
                continue
            count_files.append(sciname_file)
        sciname_count_files = [ f for f in count_files if '/'+sciname_dir+'/' in f ]
        num_sciname_count_file = len(sciname_count_files)
        if (num_sciname_count_file==0):
            sys.stderr.write('No est_counts.tsv file found in: {}\n'.format(sciname_path))
        elif (num_sciname_count_file==1):
            continue # good to go
        elif (num_sciname_count_file>=2):
            raise Exception('Multiple est_counts.tsv files found in: {}\n'.format(sciname_path))
    if (len(count_files)==0):
        raise Exception('No est_counts.tsv file was detected.')
    return count_files

def filepath2spp(file_paths):
    spp = [ os.path.basename(cf) for cf in file_paths]
    spp = [ re.sub('_', 'PLACEHOLDER', sp, count=1) for sp in spp ]
    spp = [ re.sub('_.*', '', sp) for sp in spp ]
    spp = [ re.sub('PLACEHOLDER', '_', sp) for sp in spp ]
    return spp

def cstmm_main(args):
    check_rscript()
    check_ortholog_parameter_compatibility(args)
    dir_out = os.path.realpath(args.out_dir)
    dir_cstmm = os.path.join(dir_out, 'cstmm')
    if not os.path.exists(dir_cstmm):
        os.makedirs(dir_cstmm)
    if args.dir_count=='inferred':
        dir_count = os.path.join(dir_out, 'merge')
    else:
        dir_count = os.path.realpath(args.dir_count)
    if args.dir_busco is not None:
        file_orthogroup_table = os.path.join(dir_cstmm, 'cstmm_multispecies_busco_table.tsv')
        generate_multisp_busco_table(dir_busco=args.dir_busco, outfile=file_orthogroup_table)
    elif args.orthogroup_table is not None:
        file_orthogroup_table = os.path.realpath(args.orthogroup_table)
    count_files = get_count_files(dir_count)
    if (len(count_files)==1):
        txt = 'Only one species was detected. Standard TMM normalization will be applied.'
        print(txt, flush=True)
        mode_tmm = 'single_species'
    else:
        txt = 'Multiple species were detected. ' \
              'Cross-species TMM normalization will be applied with single-copy orthologs.'
        print(txt, flush=True)
        mode_tmm = 'multi_species'
    file_genecount = os.path.join(dir_cstmm, 'cstmm_orthogroup_genecount.tsv')
    spp = filepath2spp(count_files)
    orthogroup2genecount(file_orthogroup=file_orthogroup_table, file_genecount=file_genecount, spp=spp)
    dir_amalgkit_script = os.path.dirname(os.path.realpath(__file__))
    r_cstmm_path = os.path.join(dir_amalgkit_script, 'cstmm.r')
    r_util_path = os.path.join(dir_amalgkit_script, 'util.r')
    r_command = ['Rscript', r_cstmm_path, dir_count, file_orthogroup_table, file_genecount, dir_cstmm, mode_tmm, r_util_path]
    print('')
    print('Starting R script: {}'.format(' '.join(r_command)), flush=True)
    subprocess.check_call(r_command)
    for f in glob.glob("tmp.amalgkit.*"):
        os.remove(f)



================================================
FILE: amalgkit/cstmm.r
================================================
#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(edgeR, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork, quietly = TRUE)))

mode = ifelse(length(commandArgs(trailingOnly = TRUE)) == 1, 'debug', 'batch')

if (mode == "debug") {
    dir_work = '/Users/kf/Dropbox/data/evolutionary_transcriptomics/20230527_gfe_pipeline/amalgkit_out'
    dir_count = file.path(dir_work, "merge")
    file_orthogroup_table = file.path(dir_work, 'cstmm', 'cstmm_multispecies_busco_table.tsv')
    file_genecount = file.path(dir_work, 'cstmm', 'cstmm_orthogroup_genecount.tsv')
    dir_cstmm = file.path(dir_work, "cstmm")
    mode_tmm = 'multi_species'
    r_util_path = '/Users/kf/Dropbox/repos/amalgkit/amalgkit/util.r'
    setwd(dir_work)
} else if (mode == "batch") {
    args = commandArgs(trailingOnly = TRUE)
    dir_count = args[1]
    file_orthogroup_table = args[2]
    file_genecount = args[3]
    dir_cstmm = args[4]
    mode_tmm = args[5]
    r_util_path = args[6]
}
source(r_util_path)

get_spp_filled = function(dir_count, df_gc = NA) {
    sciname_dirs = list.dirs(dir_count, full.names = FALSE, recursive = FALSE)
    spp_filled = c()
    for (sciname_dir in sciname_dirs) {
        count_files = list.files(path = file.path(dir_count, sciname_dir), pattern = ".*est_counts\\.tsv")
        if (length(count_files) == 1) {
            spp_filled = c(spp_filled, count_files)
        } else {
            warning(paste0('Multiple or no est_counts files were detected for ', sciname_dir, ': ', paste(count_files, collapse = ', ')))
        }
    }
    spp_filled = sub('_', '|', spp_filled)
    spp_filled = sub('_.*', '', spp_filled)
    spp_filled = sub('\\|', '_', spp_filled)
    if ('data.frame' %in% class(df_gc)) {
        is_missing_in_genecount = (!spp_filled %in% colnames(df_gc))
        if (sum(is_missing_in_genecount)) {
            for (sp in spp_filled[is_missing_in_genecount]) {
                warning(paste0('Species excluded. Not found in the orthogroup table: ', sp))
            }
        }
        spp_filled = spp_filled[!is_missing_in_genecount]
    }
    return(spp_filled)
}

read_est_counts = function(dir_count, sp) {
    sciname_path = file.path(dir_count, sp)
    infile = list.files(path = sciname_path, pattern = ".*est_counts\\.tsv")
    if (length(infile) > 1) {
        stop(paste0("Multiple *count.tsv files found: ", sp, "\n"))
    } else if (length(infile) == 0) {
        warning(paste0("Skipping. No *est_counts.tsv files found: ", sp, "\n"))
        return(NULL)
    }
    infile_path = file.path(sciname_path, infile[1])
    cat('Input file found, reading:', infile[1], '\n')
    dat = read.delim(infile_path, header = TRUE, row.names = 1, sep = '\t', check.names = FALSE)
    dat = dat[, (colnames(dat) != 'length'), drop = FALSE]
    colnames(dat) = paste(sp, colnames(dat), sep = '_')
    return(dat)
}

get_uncorrected = function(dir_count, file_genecount = NA) {
    if (is.na(file_genecount)) {
        df_gc = NA
    } else {
        df_gc = read.table(file_genecount, header = TRUE, sep = '\t', check.names = FALSE, quote = '', comment.char = '')
        rownames(df_gc) = df_gc[['orthogroup_id']]
        df_gc[, 'orthogroup_id'] = NULL
    }
    spp_filled = get_spp_filled(dir_count, df_gc)
    uncorrected = list()
    for (sp in spp_filled) {
        dat = read_est_counts(dir_count, sp)
        if (is.null(dat)) {
            next
        }
        uncorrected[[sp]] = dat
    }
    return(uncorrected)
}

get_df_exp_single_copy_ortholog = function(file_genecount, file_orthogroup_table, dir_count, uncorrected) {
    df_gc = read.table(file_genecount, header = TRUE, sep = '\t', check.names = FALSE, quote = '', comment.char = '')
    rownames(df_gc) = df_gc[['orthogroup_id']]
    df_gc[, 'orthogroup_id'] = NULL
    df_og = read.table(file_orthogroup_table, header = TRUE, sep = '\t', row.names = 1, check.names = FALSE, quote = '', comment.char = '')
    spp_filled = get_spp_filled(dir_count, df_gc)
    is_singlecopy = get_singlecopy_bool_index(df_gc, spp_filled)
    df_singleog = df_og[is_singlecopy, spp_filled, drop = FALSE]
    df_sog = df_singleog
    for (sp in spp_filled) {
        if (!sp %in% names(uncorrected)) {
            next
        }
        df_sog = merge(df_sog, uncorrected[[sp]], by.x = sp, by.y = "row.names", all.x = TRUE, all.y = FALSE, sort = FALSE)
    }
    df_sog = df_sog[, -(1:length(spp_filled))]
    rownames(df_sog) = rownames(df_singleog)
    return(df_sog)
}

get_df_nonzero = function(df_sog, imputation = TRUE) {
    is_no_count_col = apply(df_sog, 2, function(x) { sum(x, na.rm = TRUE) == 0 })
    txt = 'Removing %s out of %s samples whose read mapping values are all zero.\n'
    cat(sprintf(txt, formatC(sum(is_no_count_col), big.mark = ','), formatC(ncol(df_sog), big.mark = ',')))
    df_nonzero = df_sog[, !is_no_count_col]
    if (imputation) {
        df_nonzero = impute_expression(df_nonzero)
    } else {
        is_na_containing_row = apply(df_sog, 1, function(x) { any(is.na(x)) })
        txt = 'Removing %s out of %s orthogroups because missing values are observed in at least one species.\n'
        cat(sprintf(txt, formatC(sum(is_na_containing_row), big.mark = ','), formatC(nrow(df_sog), big.mark = ',')))
        df_nonzero = df_sog[!is_na_containing_row,]
    }
    return(df_nonzero)
}

create_eff_length_symlink = function(dir_count, dir_cstmm, sp) {
    path_sp = file.path(dir_count, sp)
    eff_length_files = list.files(path = path_sp, pattern = ".*eff_length\\.tsv")
    if (length(eff_length_files) == 1) {
        path_target = file.path(path_sp, eff_length_files[1])
        path_link = file.path(dir_cstmm, sp, eff_length_files[1])
        cat('Copying file from', path_target, 'to', path_link, '\n')
        file.copy(from = path_target, to = path_link, overwrite = TRUE)
    } else {
        warning(paste0('No eff_length.tsv file found: ', path_sp))
    }
}

append_tmm_stats_to_metadata = function(df_metadata, cnf_out2) {

    my_fun = function(x) {
        split_sample_name = strsplit(x, '_')[[1]]
        run_name = paste(split_sample_name[3:length(split_sample_name)], collapse = '_')
        return(run_name)
    }

    df_nf = cnf_out2[[2]]
    df_nf[['sample']] = rownames(df_nf)
    df_nf[['scientific_name']] = df_nf[['sample']]
    df_nf[['scientific_name']] = sub('_', 'PLACEHOLDER', df_nf[['scientific_name']])
    df_nf[['scientific_name']] = sub('_.*', '', df_nf[['scientific_name']])
    df_nf[['scientific_name']] = sub('PLACEHOLDER', ' ', df_nf[['scientific_name']])
    df_nf[['run']] = sapply(df_nf[['sample']], function(x) { my_fun(x) })
    df_nf = df_nf[, c('scientific_name', 'run', 'lib.size', 'norm.factors')]
    colnames(df_nf) = c('scientific_name', 'run', 'tmm_library_size', 'tmm_normalization_factor')
    out_cols = c(colnames(df_metadata), colnames(df_nf)[3:ncol(df_nf)])
    df_metadata = merge(df_metadata, df_nf, by = c('scientific_name', 'run'), sort = FALSE, all.x = TRUE, all.y = FALSE)
    df_metadata = df_metadata[, out_cols]
    filled_mapping_rate = df_metadata[['mapping_rate']]
    filled_mapping_rate[is.na(filled_mapping_rate)] = -999
    df_metadata[((filled_mapping_rate == 0) & (df_metadata[['exclusion']] == 'no')), 'exclusion'] = 'no_mapping'
    df_metadata[((!df_metadata[['run']] %in% df_nf[['run']]) & (df_metadata[['exclusion']] == 'no')), 'exclusion'] = 'no_cstmm_output'
    df_metadata[(is.na(df_metadata[['tmm_normalization_factor']]) & (df_metadata[['exclusion']] == 'no')), 'exclusion'] = 'cstmm_failed'
    return(df_metadata)
}

plot_norm_factor_histogram = function(df_metadata, font_size = 8) {
    tmp = df_metadata[(!is.na(df_metadata[['tmm_normalization_factor']])),]
    x_limit = max(abs(log2(tmp[['tmm_normalization_factor']])), na.rm = TRUE)
    for (fill_by in c('scientific_name', 'sample_group')) {
        g = ggplot2::ggplot(tmp) +
            geom_histogram(aes(x = log2(tmm_normalization_factor), fill = !!rlang::sym(fill_by)), position = "stack", alpha = 0.7, bins = 40) +
            theme_bw(base_size = font_size) +
            xlim(c(-x_limit, x_limit)) +
            labs(x = 'log2(TMM normalization factor)', y = 'Count') +
            guides(fill = guide_legend(ncol = 1)) +
            theme(
                axis.text = element_text(size = font_size, color = 'black'),
                axis.title = element_text(size = font_size, color = 'black'),
                #panel.grid.major.y=element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.grid.minor.x = element_blank(),
                legend.title = element_blank(),
                legend.text = element_text(size = font_size, color = 'black'),
                rect = element_rect(fill = "transparent"),
                plot.margin = unit(rep(0.1, 4), "cm")
            )
        out_path = file.path(dir_cstmm, paste0('cstmm_normalization_factor_histogram.', fill_by, '.pdf'))
        ggsave(out_path, plot = g, width = 4.8, height = 2.4, units = 'in')
    }
}

plot_norm_factor_scatter = function(df_metadata, font_size = 8) {
    tmp = df_metadata[(!is.na(df_metadata[['tmm_normalization_factor']])),]
    x_limit = max(abs(log2(tmp[['tmm_normalization_factor']])), na.rm = TRUE)
    g = ggplot2::ggplot(tmp, aes(x = log10(tmm_library_size), y = log2(tmm_normalization_factor), fill = scientific_name, color = sample_group)) +
        geom_point(shape = 21, alpha = 0.7) +
        scale_fill_hue(l = 65) +
        scale_color_hue(l = 45) +
        theme_bw(base_size = font_size) +
        ylim(c(-x_limit, x_limit)) +
        labs(x = 'log10(Library size)', y = 'log2(TMM normalization factor)') +
        guides(fill = guide_legend(ncol = 1), color = guide_legend(ncol = 1)) +
        theme(
            axis.text = element_text(size = font_size, color = 'black'),
            axis.title = element_text(size = font_size, color = 'black'),
            #panel.grid.major.y=element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            legend.title = element_text(size = font_size, color = 'black'),
            legend.text = element_text(size = font_size, color = 'black'),
            rect = element_rect(fill = "transparent"),
            plot.margin = unit(rep(0.1, 4), "cm")
        )
    ggsave(file.path(dir_cstmm, 'cstmm_normalization_factor_scatter.pdf'), plot = g, width = 4.8, height = 2.0, units = 'in')
}

get_library_sizes = function(df_nonzero, uncorrected) {
    df_libsize = data.frame(library_size = rep(NA, ncol(df_nonzero)))
    rownames(df_libsize) = colnames(df_nonzero)
    for (sp in names(uncorrected)) {
        for (sample in colnames(uncorrected[[sp]])) {
            if (sample %in% rownames(df_libsize)) {
                df_libsize[sample, 'library_size'] = sum(uncorrected[[sp]][[sample]], na.rm = TRUE)
            }
        }
    }
    library_sizes = df_libsize[['library_size']]
    return(library_sizes)
}

save_mean_expression_boxplot = function(df_nonzero, cnf_out2, uncorrected, corrected, font_size = 8) {
    df_nonzero_tmm = df_nonzero
    for (col in colnames(df_nonzero_tmm)) {
        tmm_normalization_factor = cnf_out2[[2]][col, 'norm.factors']
        df_nonzero_tmm[, col] = df_nonzero_tmm[, col] / tmm_normalization_factor # manually apply TMM normalization factors
    }
    mean_before = apply(df_nonzero, 2, function(x) { mean(x, na.rm = TRUE) })
    mean_after = apply(df_nonzero_tmm, 2, function(x) { mean(x, na.rm = TRUE) })
    var_before = round(var(mean_before), digits = 1)
    var_after = round(var(mean_after), digits = 1)
    txt = 'Across-species among-sample variance of mean single-copy gene raw counts before and after TMM normalization:'
    cat(txt, var_before, 'and', var_after, '\n')

    mean_sra_before = c()
    mean_sra_after = c()
    for (sp in names(corrected)) {
        mean_sra_before = c(mean_sra_before, apply(uncorrected[[sp]], 2, function(x) { mean(x, na.rm = TRUE) }))
        mean_sra_after = c(mean_sra_before, apply(corrected[[sp]], 2, function(x) { mean(x, na.rm = TRUE) }))
    }
    var_sra_before = round(var(mean_sra_before), digits = 1)
    var_sra_after = round(var(mean_sra_after), digits = 1)
    txt = 'Across-species among-sample variance of all-gene mean raw counts before and after TMM normalization:'
    cat(txt, var_sra_before, 'and', var_sra_after, '\n')

    values = c(mean_before, mean_after)
    labels = c(rep('Raw\ncounts', length(mean_before)), rep('TMM-\ncorrected\ncounts', length(mean_after)))
    df = data.frame(labels = labels, values = values)

    sra_values = c(mean_sra_before, mean_sra_after)
    sra_labels = c(rep('Raw\ncounts', length(mean_sra_before)), rep('TMM-\ncorrected\ncounts', length(mean_sra_after)))
    sra_df = data.frame(labels = sra_labels, values = sra_values)

    ps = list()
    ps[[1]] = ggplot(df, aes(x = labels, y = values)) +
        geom_boxplot(outlier.shape = NA) +
        ylab('Mean count of single-copy genes')
    ps[[2]] = ggplot(sra_df, aes(x = labels, y = values)) +
        geom_boxplot(outlier.shape = NA) +
        ylab('Mean count of all genes')

    for (i in 1:length(ps)) {
        ps[[i]] = ps[[i]] +
            ylim(0, 1000) +
            theme_bw() +
            theme(
                axis.text = element_text(size = font_size, color = 'black'),
                axis.title = element_text(size = font_size, color = 'black'),
                legend.text = element_text(size = font_size, color = 'black'),
                legend.title = element_text(size = font_size, color = 'black'),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.grid.minor.x = element_blank(),
                rect = element_rect(fill = "transparent"),
                plot.margin = unit(rep(0.1, 4), "cm"),
            )
    }

    p = ps[[1]] + ps[[2]]
    filename = file.path(dir_cstmm, 'cstmm_mean_expression_boxplot.pdf')
    ggsave(file = filename, p, height = 3.6, width = 3.6)
}

save_corrected_output_files = function(uncorrected) {
    corrected = list()
    for (sp in names(uncorrected)) {
        dat = uncorrected[[sp]]
        df_nf_sp = cnf_out2[[2]][startsWith(rownames(cnf_out2[[2]]), sp),]
        for (i in 1:length(df_nf_sp[, 1])) {
            SRR = as.character(row.names(df_nf_sp[i,]))
            tmm_normalization_factor = as.double(df_nf_sp[i, "norm.factors"])
            dat[, SRR] = dat[, SRR] / tmm_normalization_factor # manually apply TMM normalization factors
        }
        dat_out = cbind(target_id = rownames(dat), dat)
        rownames(dat_out) = NULL
        colnames(dat_out) = sub(paste0(sp, '_'), '', colnames(dat_out))
        dir_cstmm_sp = file.path(dir_cstmm, sp)
        if (!file.exists(dir_cstmm_sp)) {
            dir.create(dir_cstmm_sp)
        }
        create_eff_length_symlink(dir_count, dir_cstmm, sp)
        file_path = file.path(dir_cstmm_sp, paste0(sp, "_cstmm_counts.tsv"))
        write.table(dat_out, file_path, sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
        corrected[[sp]] = dat_out
        rownames(corrected[[sp]]) = corrected[[sp]][['target_id']]
        corrected[[sp]][, 'target_id'] = NULL
    }
    return(corrected)
}

print_excluded_run_summary = function(df_metadata) {
    exclusion_reasons = sort(unique(df_metadata[['exclusion']]))
    for (reason in exclusion_reasons) {
        num_run = sum(df_metadata[['exclusion']] == reason)
        if (reason == 'no') {
            txt_final = paste0('Number of retained runs (exclusion = "no"): ', num_run, '\n')
        } else {
            txt = paste0('Number of excluded runs (exclusion = "', reason, '"): ', num_run, '\n')
            cat(txt)
        }
    }
    cat(txt_final)
}


if (mode_tmm == 'single_species') {
    uncorrected = get_uncorrected(dir_count = dir_count, file_genecount = NA)
    stopifnot(length(names(uncorrected)) == 1)
    sp = names(uncorrected)[[1]]
    df_sog = uncorrected[[sp]]
    df_nonzero = get_df_nonzero(df_sog)
} else if (mode_tmm == 'multi_species') {
    uncorrected = get_uncorrected(dir_count = dir_count, file_genecount = file_genecount)
    df_sog = get_df_exp_single_copy_ortholog(file_genecount, file_orthogroup_table, dir_count, uncorrected)
    df_nonzero = get_df_nonzero(df_sog)
}

libraray_sizes = get_library_sizes(df_nonzero, uncorrected)
cnf_in = edgeR::DGEList(counts = df_nonzero, lib.size = libraray_sizes)
cat('Round 1: Performing TMM normalization to determine the appropriate baseline.\n')
cnf_out1 = edgeR::calcNormFactors(cnf_in, method = 'TMM', refColumn = NULL)
x = cnf_out1[[2]][['norm.factors']]
cat('Round 1: Median TMM normalization factor =', median(x), '\n')
median_value = sort(x)[ceiling(length(x) / 2)]
median_index = (1:length(x))[x == median_value]

cat('Round 2: Performing TMM normalization for output.\n')
cnf_out2 = edgeR::calcNormFactors(cnf_in, method = 'TMM', refColumn = median_index)
cat('Round 2: Median TMM normalization factor =', median(cnf_out2[[2]][['norm.factors']]), '\n')

path_metadata = file.path(dir_count, 'metadata.tsv')
df_metadata = read.table(path_metadata, header = TRUE, sep = '\t', check.names = FALSE, quote = '', comment.char = '')
df_metadata = append_tmm_stats_to_metadata(df_metadata, cnf_out2)
df_metadata = df_metadata[, !startsWith(colnames(df_metadata), 'Unnamed')]
out_path = file.path(dir_cstmm, 'metadata.tsv')
write.table(df_metadata, out_path, row.names = FALSE, sep = '\t', quote = FALSE)
print_excluded_run_summary(df_metadata)
plot_norm_factor_histogram(df_metadata = df_metadata)
plot_norm_factor_scatter(df_metadata = df_metadata)
corrected = save_corrected_output_files(uncorrected)
save_mean_expression_boxplot(df_nonzero, cnf_out2, uncorrected, corrected, font_size = 8)

cat(sprintf('Number of SRA samples for exclusion potting: %s\n', formatC(nrow(df_metadata), format = 'd', big.mark = ',')))
out_path = file.path(dir_cstmm, 'cstmm_exclusion.pdf')
save_exclusion_plot(df = df_metadata, out_path = out_path, font_size = 8)

cat('cstmm.r completed!\n')


================================================
FILE: amalgkit/curate.py
================================================
import datetime
import os
import re
import shutil
import subprocess
import sys
import warnings

from amalgkit.util import *


def get_sample_group(args, metadata):
    if args.sample_group is None:
        sample_group = metadata.df.loc[:, 'sample_group'].dropna().unique()
    else:
        sample_group = re.findall(r"[\w]+", args.sample_group)
    if (len(sample_group)==0):
        txt = 'The "sample_group" column in --metadata ({}) is not filled. Please manually edit the file.\n'
        txt += '`amalgkit curate` recognizes samples with the same string in this columns to belong the same group.\n'
        txt += 'Exiting.\n'
        sys.stderr.write(txt.format(args.metadata))
        sys.exit(1)
    print('Tissues to be included: {}'.format(', '.join(sample_group)))
    sample_group = '|'.join(sample_group)
    return sample_group

def run_curate_r_script(args, metadata, sp, input_dir):
    dist_method = args.dist_method
    mr_cut = args.mapping_rate
    correlation_threshold = args.correlation_threshold
    intermediate = args.plot_intermediate
    sample_group = get_sample_group(args, metadata)
    dir_amalgkit_script = os.path.dirname(os.path.realpath(__file__))
    r_script_path = os.path.join(dir_amalgkit_script, 'curate.r')
    r_util_path = os.path.join(dir_amalgkit_script, 'util.r')
    path_curate_input_metadata = os.path.join(input_dir, 'metadata.tsv')
    len_file = os.path.join(os.path.abspath(input_dir), sp, sp + '_eff_length.tsv')
    if 'cstmm' in input_dir:
        count_file = os.path.join(os.path.abspath(input_dir), sp, sp + '_cstmm_counts.tsv')
    else:
        count_file = os.path.join(os.path.abspath(input_dir), sp, sp + '_est_counts.tsv')
    if os.path.exists(count_file) and os.path.exists(len_file):
        print("Both counts and effective length files found: {}".format(sp), flush=True)
    else:
        if not os.path.exists(count_file):
            sys.stderr.write('Expected but undetected PATH of the count file: {}\n'.format(count_file))
        if not os.path.exists(len_file):
            sys.stderr.write('Expected but undetected PATH of the effective length file: {}\n'.format(len_file))
        sys.stderr.write('Skipping {}\n'.format(sp))
        print('Skipping {}'.format(sp), flush=True)
        return 1
    print("Starting Rscript to obtain curated {} values.".format(args.norm), flush=True)
    curate_r_exit_code = subprocess.call([
            'Rscript',
            r_script_path,
            count_file,
            path_curate_input_metadata,
            os.path.realpath(args.out_dir),
            len_file,
            dist_method,
            str(mr_cut),
            '0',
            str(intermediate),
            sample_group,
            args.sample_group_color,
            str(args.norm),
            str(args.one_outlier_per_iter),
            str(correlation_threshold),
            str(args.batch_effect_alg),
            str(args.clip_negative),
            str(args.maintain_zero),
            os.path.realpath(r_util_path),
            str(args.skip_curation)
         ])
    return curate_r_exit_code

def curate_main(args):
    check_rscript()
    if args.input_dir=='inferred':
        dir_merge = os.path.realpath(os.path.join(args.out_dir, 'merge'))
        dir_cstmm = os.path.realpath(os.path.join(args.out_dir, 'cstmm'))
        if os.path.exists(dir_cstmm):
            print('Subdirectory for amalgkit cstmm will be used as input: {}'.format(dir_cstmm))
            metadata = load_metadata(args, dir_subcommand='cstmm')
            input_dir = dir_cstmm
        else:
            print('Subdirectory for amalgkit merge will be used as input: {}'.format(dir_merge))
            metadata = load_metadata(args, dir_subcommand='merge')
            input_dir = dir_merge
    else:
        print('Input_directory: {}'.format(args.input_dir))
        metadata = load_metadata(args, dir_subcommand=os.path.basename(args.input_dir))
        input_dir = args.input_dir
    if ('tpm' in args.norm) & ('cstmm' in input_dir):
            txt = ("TPM and TMM are incompatible. "
                   "If input data are CSTMM-normalized, "
                   "please switch --norm to any of the 'fpkm' normalization methods instead.")
            sys.stderr.write(txt)
    is_selected = (metadata.df['exclusion']=='no')
    spp = metadata.df.loc[is_selected, 'scientific_name'].drop_duplicates().values
    curate_dir = os.path.join(args.out_dir, 'curate')
    if not os.path.exists(curate_dir):
        os.mkdir(curate_dir)
    print('Number of species in the selected metadata table ("exclusion"=="no"): {}'.format(len(spp)), flush=True)
    for sp in spp:
        sp = sp.replace(" ", "_")
        file_curate_completion_flag = os.path.join(curate_dir, sp, 'curate_completion_flag.txt')
        if os.path.exists(file_curate_completion_flag):
            if args.redo:
                print('Output file detected. Will be overwritten: {}'.format(sp), flush=True)
                shutil.rmtree(os.path.join(curate_dir, sp))
            else:
                print('Skipping. Output file detected: {}'.format(sp), flush=True)
                continue
        print('Starting: {}'.format(sp), flush=True)
        exit_status = run_curate_r_script(args, metadata, sp, input_dir)
        if exit_status == 0:
            with open(file_curate_completion_flag, 'w') as f:
                f.write('amalgkit curate completed at {}\n'.format(datetime.datetime.now()))



================================================
FILE: amalgkit/curate.r
================================================
#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(pcaMethods, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(colorspace, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(MASS, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(NMF, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(dendextend, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(amap, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(pvclust, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(Rtsne, quietly = TRUE)))

debug_mode = ifelse(length(commandArgs(trailingOnly = TRUE)) == 1, "debug", "batch")
log_prefix = "transcriptome_curation.r:"
cat(log_prefix, "mode =", debug_mode, "\n")
if (debug_mode == "debug") {
    out_dir = '/Users/kf/Dropbox/data/evolutionary_transcriptomics/20230527_amalgkit/amalgkit_out'
    metadata_path = '/Users/kf/Dropbox/data/evolutionary_transcriptomics/20230527_amalgkit/amalgkit_out/cstmm/metadata.tsv'
    est_counts_path = '/Users/kf/Dropbox/data/evolutionary_transcriptomics/20230527_amalgkit/amalgkit_out/cstmm/Actinidia_chinensis/Actinidia_chinensis_cstmm_counts.tsv'
    eff_length_path = '/Users/kf/Dropbox/data/evolutionary_transcriptomics/20230527_amalgkit/amalgkit_out/cstmm/Actinidia_chinensis/Actinidia_chinensis_eff_length.tsv'
    dist_method = "pearson"
    mapping_rate_cutoff = .20
    min_dif = 0
    plot_intermediate = as.logical(0)
    selected_sample_groups = c("root", "flower", "leaf")
    sample_group_colors = 'DEFAULT'
    transform_method = "log2p1-fpkm"
    one_outlier_per_iteration = as.logical(0)
    correlation_threshold = 0.3
    batch_effect_alg = 'sva'
    dist_method = "pearson"
    clip_negative = as.logical(1)
    maintain_zero = as.logical(1)
    r_util_path = '/Users/kf/Dropbox/repos/amalgkit/amalgkit/util.r'
    skip_curation_flag = FALSE
    setwd(file.path(out_dir, 'curate'))
} else if (debug_mode == "batch") {
    args = commandArgs(trailingOnly = TRUE)
    est_counts_path = args[1]
    metadata_path = args[2]
    out_dir = args[3]
    eff_length_path = args[4]
    dist_method = args[5]
    mapping_rate_cutoff = as.numeric(args[6])
    min_dif = as.numeric(args[7])
    plot_intermediate = as.logical(as.integer(args[8]))
    selected_sample_groups = strsplit(args[9], "\\|")[[1]]
    sample_group_colors = strsplit(args[10], ",")[[1]]
    transform_method = args[11]
    one_outlier_per_iteration = as.integer(args[12])
    correlation_threshold = as.numeric(args[13])
    batch_effect_alg = args[14]
    clip_negative = as.logical(as.integer(args[15]))
    maintain_zero = as.logical(as.integer(args[16]))
    r_util_path = file.path(args[17])
    skip_curation_flag = as.logical(as.integer(args[18]))
}
cat('est_counts_path:', est_counts_path, "\n")
cat('metadata_path:', metadata_path, "\n")
cat('out_dir:', out_dir, "\n")
cat('eff_length_path:', eff_length_path, "\n")
cat('dist_method:', dist_method, "\n")
cat('mapping_rate_cutoff:', mapping_rate_cutoff, "\n")
cat('min_dif:', min_dif, "\n")
cat('plot_intermediate:', plot_intermediate, "\n")
cat('selected_sample_groups:', selected_sample_groups, "\n")
cat('selected_sample_group_colors:', sample_group_colors, "\n")
cat('transform_method:', transform_method, "\n")
cat('one_outlier_per_iteration:', one_outlier_per_iteration, "\n")
cat('correlation_threshold:', correlation_threshold, "\n")
cat('batch_effect_alg:', batch_effect_alg, "\n")
cat('clip_negative:', clip_negative, "\n")
cat('maintain_zero:', maintain_zero, "\n")
cat('r_util_path:', r_util_path, "\n")
cat('skip_curation_flag:', skip_curation_flag, "\n")

source(r_util_path)

if (batch_effect_alg == "ruvseq") {
    suppressWarnings(suppressPackageStartupMessages(library(RUVSeq, quietly = TRUE)))
} else {
    suppressWarnings(suppressPackageStartupMessages(library(sva, quietly = TRUE)))
}

tc_metadata_intersect = function(tc, sra) {
    sra_run = sra[['run']]
    tc = tc[, colnames(tc) %in% sra_run, drop = FALSE]
    sra = sra[sra[['run']] %in% colnames(tc),]
    return(list(tc = tc, sra = sra))
}

remove_nonexpressed_gene = function(tc) {
    gene_sum = apply(tc, 1, sum)
    tc_ex = tc[gene_sum > 0,]
    tc_ne = tc[gene_sum == 0,]
    return(list(tc_ex = tc_ex, tc_ne = tc_ne))
}

add_color_to_curate_metadata = function(sra, selected_sample_groups, sample_group_colors) {
    sra = sra[, (!colnames(sra) %in% c("bp_color", "sp_color", "sample_group_color"))]
    if ('bioproject' %in% colnames(sra)) {
        bioproject = as.character(sra[['bioproject']])
        bioproject_u = sort(unique(bioproject))
    } else {
        bioproject_u = rep('PLACEHOLDER', nrow(sra))
    }
    scientific_name = as.character(sra[['scientific_name']])
    scientific_name_u = sort(unique(scientific_name))
    sample_group = as.character(sra[['sample_group']])
    sample_group_u = selected_sample_groups
    if (length(sample_group_colors) == 1 && sample_group_colors == "DEFAULT") {
        if (length(selected_sample_groups) <= 8) {
            sample_group_color = brewer.pal(8, "Dark2")
            sample_group_color = sample_group_color[1:length(selected_sample_groups)]
            bp_color = rainbow_hcl(length(bioproject_u), c = 50)
            sp_color = rainbow_hcl(length(scientific_name_u), c = 100)
        } else if (length(selected_sample_groups) <= 12) {
            sample_group_color = brewer.pal(length(selected_sample_groups), "Paired")
            bp_color = rainbow_hcl(length(bioproject_u), c = 50)
            sp_color = rainbow_hcl(length(scientific_name_u), c = 100)
        } else {
            sample_group_color = rainbow_hcl(length(selected_sample_groups), c = 100)
            bp_color = rainbow_hcl(length(bioproject_u), c = 50)
            sp_color = rainbow_hcl(length(scientific_name_u), c = 150)
        }
    } else {
        if (length(sample_group_colors) != length(selected_sample_groups)) {
            stop("Length of sample_group_colors must match length of selected_sample_groups")
        }
        sample_group_color = sample_group_colors
        bp_color = rainbow_hcl(length(bioproject_u), c = 50)
        sp_color = rainbow_hcl(length(scientific_name_u), c = 100)
    }
    df_sample_group = data.frame(sample_group = sample_group_u, sample_group_color = sample_group_color[1:length(sample_group_u)], stringsAsFactors = FALSE)
    df_bp = data.frame(bioproject = bioproject_u, bp_color = bp_color[1:length(bioproject_u)], stringsAsFactors = FALSE)
    df_sp = data.frame(scientific_name = scientific_name_u, sp_color = sp_color[1:length(scientific_name_u)], stringsAsFactors = FALSE)
    sra = merge(sra, df_bp, sort = FALSE, all.y = FALSE)
    sra = merge(sra, df_sp, sort = FALSE, all.y = FALSE)
    sra = merge(sra, df_sample_group, sort = FALSE, all.y = FALSE)
    return(sra)
}


sort_tc_and_metadata = function(tc, sra, sort_columns = c("sample_group", "scientific_name", "bioproject")) {
    for (column in rev(sort_columns)) {
        if (!column %in% colnames(sra)) {
            next
        }
        sra = sra[order(sra[[column]]),]
    }
    sra_intersection = sra[(sra[['run']] %in% colnames(tc)), 'run']
    tc = tc[, sra_intersection, drop = FALSE]
    return(list(tc = tc, sra = sra))
}

sort_averaged_tc = function(tc) {
    split_colnames = strsplit(colnames(tc), "_")
    genus_names = c()
    specific_names = c()
    sample_group_names = c()
    for (i in 1:length(split_colnames)) {
        genus_names = c(genus_names, split_colnames[[i]][1])
        specific_names = c(specific_names, split_colnames[[i]][2])
        sample_group_names = c(sample_group_names, split_colnames[[i]][3])
    }
    colname_order = order(sample_group_names, genus_names, specific_names)
    tc = tc[, colname_order, drop = FALSE]
    return(tc)
}

cleanY = function(y, mod, svs) {
    X = cbind(mod, svs)
    Hat = solve(t(X) %*% X) %*% t(X)
    beta = (Hat %*% t(y))
    P = ncol(mod)
    return(y - t(as.matrix(X[, -c(1:P)]) %*% beta[-c(1:P),]))
}

sample_group_mean = function(tc, sra, selected_sample_groups = NA, balance_bp = FALSE) {
    # if the data are SVA-corrected, balance_bp would not be necessary because project-specfic effects
    # were removed in SVA already.
    if (all(is.na(selected_sample_groups))) {
        sp_sample_groups = unique(sra[['sample_group']])
    } else {
        sp_sample_groups = selected_sample_groups[selected_sample_groups %in% unique(sra[['sample_group']])]
    }
    tc_ave = data.frame(matrix(rep(NA, length(sp_sample_groups) * nrow(tc)), nrow = nrow(tc)))
    colnames(tc_ave) = sp_sample_groups
    rownames(tc_ave) = rownames(tc)
    for (sample_group in sp_sample_groups) {
        exclusion_sample_group = sra[(sra[['sample_group']] == sample_group), 'exclusion']
        run_sample_group = sra[(sra[['sample_group']] == sample_group), 'run']
        run_sample_group = run_sample_group[run_sample_group %in% colnames(tc)]
        if (all(exclusion_sample_group != "no")) {
            warning_message = paste0('All samples of sample_group ', sample_group, ' are marked for exclusion. This sample_group will be omitted from further analysis.')
            selected_sample_groups = selected_sample_groups[!selected_sample_groups == sample_group]
            tc_ave = tc_ave[, !names(tc_ave) %in% c(sample_group)]
            warning(warning_message)
            next
        }
        if (length(run_sample_group) == 1) {
            exp_sample_group = tc[, run_sample_group]
        } else {
            if (balance_bp) {
                is_run = (sra[['run']] %in% colnames(tc))
                is_sample_group = (sra[['sample_group']] == sample_group)
                is_no_exclusion = (sra[['exclusion']] == "no")
                bps = unique(sra[is_run & is_sample_group & is_no_exclusion, "bioproject"])
                df_tmp = data.frame(matrix(rep(NA, nrow(tc) * length(bps)), nrow = nrow(tc), ncol = length(bps)))
                colnames(df_tmp) = bps
                for (bp in bps) {
                    is_bp = (sra[['bioproject']] == bp)
                    is_sample_group = (sra[['sample_group']] == sample_group)
                    is_no_exclusion = (sra[['exclusion']] == "no")
                    sra_ids = sra[is_bp & is_sample_group & is_no_exclusion, "run"]
                    tc_bp = tc[, sra_ids]
                    if (class(tc_bp) == "numeric") {
                        df_tmp[bp] = tc_bp
                    } else {
                        df_tmp[bp] = rowMeans(tc_bp)
                    }
                }
                exp_sample_group = rowMeans(df_tmp)
            } else {
                exp_sample_group = rowMeans(tc[, run_sample_group])
            }
        }
        tc_ave[, sample_group] = exp_sample_group
    }
    return(list(tc_ave = tc_ave, selected_sample_groups = selected_sample_groups))
}

row_tau = function(row) {
    is_nonzero = row > 0
    if (sum(is_nonzero) > 0) {
        exp_order = order(row[is_nonzero], decreasing = TRUE)
        sample_group_ordered = colnames(tc_sample_group)[is_nonzero][exp_order]
        highest = sample_group_ordered[1]
        order = paste(sample_group_ordered, collapse = "|")
        return(c(highest, order))
    } else {
        return(c(NA, NA))  # Return NA values if no nonzero elements
    }
}

sample_group2tau = function(tc_sample_group, rich.annotation = TRUE, transform_method) {
    if (rich.annotation) {
        cols = c("tau", "highest", "order")
    } else {
        cols = c("tau")
    }
    df_tau = data.frame(matrix(rep(NA, length(cols) * nrow(tc_sample_group)), nrow = nrow(tc_sample_group)))
    colnames(df_tau) = cols
    rownames(df_tau) = rownames(tc_sample_group)
    if (grepl('logn-', transform_method)) {
        tc_sample_group = exp(tc_sample_group)
    } else if (grepl('log2-', transform_method)) {
        tc_sample_group = 2**tc_sample_group
    } else if (grepl('lognp1-', transform_method)) {
        tc_sample_group = exp(tc_sample_group) - 1
    } else if (grepl('log2p1-', transform_method)) {
        tc_sample_group = 2**tc_sample_group - 1
    }
    tc_sample_group[tc_sample_group < 0] = 0
    xmax = apply(tc_sample_group, 1, function(x) { max(x, na.rm = TRUE) })
    df_tau[, 'tau'] = apply((1 - (tc_sample_group / xmax)) / (ncol(tc_sample_group) - 1), 1, sum)
    if (rich.annotation) {
        tc_sample_group[is.na(tc_sample_group)] = 0
        results = t(apply(tc_sample_group, 1, row_tau))
        df_tau[, c("highest", "order")] = results
    }
    return(df_tau)
}

check_mapping_rate = function(tc, sra, mapping_rate_cutoff) {
    if ('mapping_rate' %in% colnames(sra)) {
        cat(paste0('Mapping rate cutoff: ', mapping_rate_cutoff * 100, '%\n'))
        is_mapping_good = (sra[['mapping_rate']] > mapping_rate_cutoff * 100)
        is_mapping_good[is.na(is_mapping_good)] = TRUE
        if (any(!is_mapping_good)) {
            cat("Removed due to low mapping rate:\n")
            df_tmp = sra[!is_mapping_good,]
            for (i in rownames(df_tmp)) {
                sra_id = df_tmp[i, 'run']
                mapping_rate = df_tmp[i, 'mapping_rate']
                cat(paste0(sra_id, ': mapping rate = ', mapping_rate, '%\n'))
            }
            tc = tc[, colnames(tc) %in% sra[is_mapping_good, "run"], drop = FALSE]
        } else {
            cat("No entry removed due to low mapping rate.\n")
        }
        sra[!is_mapping_good, "exclusion"] = "low_mapping_rate"
    } else {
        cat('Mapping rate cutoff will not be applied.\n')
    }
    return(list(tc = tc, sra = sra))
}

check_within_sample_group_correlation = function(tc, sra, dist_method, min_dif, selected_sample_groups, one_out_per_iter = TRUE, correlation_threshold) {
    if (length(selected_sample_groups) == 1) {
        cat('Only one sample_group category is available. Outlier removal will be skipped.\n')
        return(list(tc = tc, sra = sra))
    }
    out = tc_metadata_intersect(tc, sra)
    tc = out[["tc"]]
    sra2 = out[["sra"]]
    sra2[, 'num_other_run_same_bp_sample_group'] = 0
    selected_sample_groups = selected_sample_groups[selected_sample_groups %in% unique(sra2[['sample_group']])]
    exclude_runs = c()
    for (sra_run in colnames(tc)) {
        is_sra = (sra2[['run']] == sra_run)
        my_sample_group = sra2[is_sra, "sample_group"]
        my_bioproject = sra2[is_sra, "bioproject"]
        is_not_my_bp = (sra2[['bioproject']] != my_bioproject)
        is_my_sample_group = (sra2[['sample_group']] == my_sample_group)
        run_other_bp = sra2[(is_not_my_bp | !is_my_sample_group), "run"]
        run_other_bp = run_other_bp[run_other_bp %in% colnames(tc)]
        tc_other_bp = tc[, run_other_bp]
        num_other_run_same_bp_sample_group = length(unique(sra2[(is_not_my_bp & is_my_sample_group), "bioproject"]))
        sra2[is_sra, "num_other_run_same_bp_sample_group"] = num_other_run_same_bp_sample_group
        sra2_other_bp <- sra2[sra2[['run']] %in% run_other_bp,]

        # If one sample_group is completely sourced from the same bioproject, we can't remove the whole bioproject for tc_ave_other_bp

        num_other_bp_same_sample_group = sum(sra2_other_bp[['sample_group']] == my_sample_group, na.rm = TRUE)
        if (num_other_bp_same_sample_group == 0) {
            tc_ave_other_bp = sample_group_mean(tc, sra2, selected_sample_groups)[['tc_ave']]
        } else {
            tc_ave_other_bp = sample_group_mean(tc_other_bp, sra2, selected_sample_groups)[['tc_ave']]
        }
        tc_ave = sample_group_mean(tc, sra2, selected_sample_groups)[['tc_ave']]
        coef = c()
        coef_other_bp = c()
        for (sample_group in selected_sample_groups) {
            tmp_coef = cor(tc[, sra_run], tc_ave[, sample_group], method = dist_method)
            if ((num_other_bp_same_sample_group == 0) & (sample_group == my_sample_group)) {
                tmp_coef_other_bp = cor(tc[, sra_run], tc_ave[, sample_group], method = dist_method)
            } else {
                tmp_coef_other_bp = cor(tc[, sra_run], tc_ave_other_bp[, sample_group], method = dist_method)
            }
            if (sample_group == my_sample_group) {
                tmp_coef = tmp_coef - min_dif
                tmp_coef_other_bp = tmp_coef_other_bp - min_dif
            }
            coef = c(coef, tmp_coef)
            coef_other_bp = c(coef_other_bp, tmp_coef_other_bp)
        }
        names(coef) = selected_sample_groups
        names(coef_other_bp) = selected_sample_groups
        if (max(coef, na.rm = TRUE) != coef[my_sample_group]) {
            cat('Registered as a candidate for exclusion. Better correlation to other sample group(s):', sra_run, '\n')
            exclude_runs = c(exclude_runs, sra_run)
        }
        if (coef_other_bp[my_sample_group] < correlation_threshold) {
            cat('Registered as a candidate for exclusion. Low within-sample-group correlation:', sra_run, '\n')
            exclude_runs = c(exclude_runs, sra_run)
        }
    }
    if (length(exclude_runs)) {
        if (one_out_per_iter == TRUE) {
            cat("Excluding only one outlier per BioProject or same sample_group. \n")
            exclude_run_bps_and_sample_group = sra2[(sra2$run %in% exclude_runs), c("bioproject", "run", "sample_group")]
            first_bp_hit = exclude_run_bps_and_sample_group[match(unique(exclude_run_bps_and_sample_group$bioproject), exclude_run_bps_and_sample_group$bioproject),]
            first_same_sample_group_hit = exclude_run_bps_and_sample_group[match(unique(exclude_run_bps_and_sample_group$sample_group), exclude_run_bps_and_sample_group$sample_group),]
            # if a first_same_sample_group_hit is part of the same bioproject as the other removal candidates, ommit the same sample_group candidates
            if (any(first_same_sample_group_hit$bioproject %in% first_bp_hit$bioproject)) {
                exclude_runs_tmp = c(first_bp_hit$run, first_same_sample_group_hit[!first_same_sample_group_hit$bioproject %in% first_bp_hit$bioproject]$run)
            }else {
                exclude_runs_tmp = c(first_bp_hit$run, first_same_sample_group_hit$run)
            }
            exclude_runs = unique(exclude_runs_tmp)
            # TODO This boolean vector should be all TRUE by definition. ???: exclude_run_bps_and_sample_group$run %in% exclude_runs
            exclude_bps = exclude_run_bps_and_sample_group[exclude_run_bps_and_sample_group$run %in% exclude_runs, "bioproject"]
        } else {
            exclude_run_bps = sra2[(sra2[['run']] %in% exclude_runs), c("bioproject", "run", "num_other_run_same_bp_sample_group")]
            exclude_bp_counts = data.frame(table(exclude_run_bps[['bioproject']]))
            exclude_run_bps = merge(exclude_run_bps, exclude_bp_counts, by.x = "bioproject", by.y = "Var1")
            exclude_run_bps = exclude_run_bps[order(exclude_run_bps[['num_other_run_same_bp_sample_group']], exclude_run_bps[['Freq']]),]
            rownames(exclude_run_bps) = 1:nrow(exclude_run_bps)
            min_other_run_same_bp_sample_group = exclude_run_bps[1, "num_other_run_same_bp_sample_group"]
            semimin_bp_count = exclude_run_bps[1, "Freq"]
            cat("Minimum number of other BioProjects within sample_group:", min_other_run_same_bp_sample_group, "\n")
            cat("Semi-minimum count of exclusion-candidate BioProjects:", semimin_bp_count, "\n")
            conditions = (exclude_run_bps[['Freq']] == semimin_bp_count)
            conditions = conditions & (exclude_run_bps[['num_other_run_same_bp_sample_group']] == min_other_run_same_bp_sample_group)
            exclude_bps = unique(exclude_run_bps[conditions, "bioproject"])
            exclude_runs = exclude_run_bps[(exclude_run_bps[['bioproject']] %in% exclude_bps), "run"]
        }
    }
    if (length(exclude_runs)) {
        cat('Partially removed BioProjects due to low within-sample_group correlation:', paste(exclude_bps, collapse = ' '), '\n')
        cat('Removed Runs due to low within-sample_group correlation:', paste(exclude_runs, collapse = ' '), '\n')
    }
    tc = tc[, !colnames(tc) %in% exclude_runs, drop = FALSE]
    sra[(sra[['run']] %in% exclude_runs), "exclusion"] = "low_within_sample_group_correlation"
    return(list(tc = tc, sra = sra))
}

batch_effect_subtraction = function(tc, sra, batch_effect_alg, transform_method, clip_negative) {
    if (ncol(tc) == 1) {
        cat('Only 1 sample is available. Skipping batch effect subtraction.\n')
        return(list(tc = tc, sva = NULL))
    }
    out = tc_metadata_intersect(tc, sra)
    tc = out[["tc"]]
    sra = out[["sra"]]
    out = sort_tc_and_metadata(tc, sra)
    tc = out[["tc"]]
    sra = out[["sra"]]
    out = remove_nonexpressed_gene(tc)
    tc = out[["tc_ex"]]
    tc_ne = out[["tc_ne"]]
    if (batch_effect_alg == "sva") {
        mod = try(model.matrix(~sample_group, data = sra))
        if ("try-error" %in% class(mod)) {
            return(list(tc = tc, sva = NULL))
        }
        mod0 = model.matrix(~1, data = sra)
        sva1 = try(sva(dat = as.matrix(tc), mod = mod, mod0 = mod0, B = 10))
        if ((class(sva1) != "try-error")) {
            cat("SVA correction was correctly performed.\n")
            tc = cleanY(y = tc, mod = mod, svs = sva1[['sv']])
            tc = rbind(tc, tc_ne)
        } else if ((class(sva1) == "try-error")) {
            cat("SVA correction failed.")
            # tc = rbind(tc, tc_ne) # Original tc shouldn't be returned as it is confusing. We may need a logic to retrieve the result of the previous round as final output.
        }
    } else if (batch_effect_alg == "combatseq") {
        bp_freq = as.data.frame(table(sra[, "bioproject"]))
        bp_freq_gt1 = bp_freq[bp_freq[, "Freq"] > 1, "Var1"]
        bp_freq_eq1 = bp_freq[bp_freq[, "Freq"] == 1, "Var1"]
        run_bp_freq_gt1 = sra[sra[, "bioproject"] %in% bp_freq_gt1, "run"]
        run_bp_freq_eq1 = sra[sra[, "bioproject"] %in% bp_freq_eq1, "run"]
        tc_combat = tc[, colnames(tc) %in% run_bp_freq_gt1]
        tcc_cn = colnames(tc_combat)
        batch = sra[(sra[, "bioproject"] %in% bp_freq_gt1), "bioproject"]
        group = sra[(sra[, "run"] %in% run_bp_freq_gt1), "sample_group"]
        tc_combat = try(ComBat_seq(as.matrix(tc_combat), batch = batch, group = group))
        if (class(tc_combat)[1] != "try-error") {
            cat("These runs are being removed, due to the bioproject only having 1 sample: \n")
            print(run_bp_freq_eq1)
            cat("Combatseq correction was correctly performed.\n")
            tc_combat = as.data.frame(tc_combat)
            colnames(tc_combat) = tcc_cn
            tc = rbind(tc_combat, tc_ne[, colnames(tc_combat)])
            sva1 = ''
        } else {
            cat("Combatseq correction failed. Trying again without group parameter. \n")
            tc_combat = tc[, colnames(tc) %in% run_bp_freq_gt1]
            tc_combat = try(ComBat_seq(as.matrix(tc_combat), batch = batch))
            if (class(tc_combat)[1] != "try-error") {
                cat("These runs are being removed, due to the bioproject only having 1 sample: \n")
                print(run_bp_freq_eq1)
                cat("Combatseq correction was correctly performed.\n")
                tc_combat = as.data.frame(tc_combat)
                colnames(tc_combat) = tcc_cn
                tc = rbind(tc_combat, tc_ne[, colnames(tc_combat)])
                sva1 = ''
            }
            else {
                stop("Combatseq correction failed.")
            }
        }
    } else if (batch_effect_alg == "ruvseq") {
        x = as.factor(sra$sample_group)
        design = try(model.matrix(~sample_group, data = sra))
        if ("try-error" %in% class(design)) {
            return(list(tc = tc, sva = NULL))
        }
        y = DGEList(counts = as.matrix(tc + 1), group = x)
        y = calcNormFactors(y, method = "upperquartile")
        y = estimateGLMCommonDisp(y, design)
        y = estimateGLMTagwiseDisp(y, design)
        fit = glmFit(y, design)
        res = residuals(fit, type = "deviance")
        seqUQ = betweenLaneNormalization(as.matrix(tc + 1), which = "upper", round = TRUE, offset = FALSE)
        controls = rep(TRUE, dim(as.matrix(tc + 1))[1])
        batch_ruv_res = try(RUVr(seqUQ, controls, k = 1, res)[[2]])
        if (class(batch_ruv_res)[1] != "try-error") {
            cat("RUVseq correction was correctly performed.\n")
            tc = rbind(batch_ruv_res, tc_ne)
            sva1 = ''
        } else {
            stop("RUVseq correction failed.")
            # tc = rbind(tc, tc_ne) # Original tc shouldn't be returned as it is confusing. We may need a logic to retrieve the result of the previous round as final output.
        }
    } else if (batch_effect_alg == "no") {
        cat("No batch effect correction was performed.\n")
        tc = rbind(tc, tc_ne)
        sva1 = ''
    } else {
        stop("Invalid batch effect correction algorithm.")
    }
    if (clip_negative) {
        if (endsWith(sub('-.*', '', transform_method), 'p1')) {
            is_negative = (tc < 0)
            txt = 'Number of negative values clipped to zero: %s/%s (%s x %s matrix)\n'
            val1 = format(sum(is_negative), big.mark = ',', scientific = FALSE)
            val2 = format(prod(dim(tc)), big.mark = ',', scientific = FALSE)
            val3 = format(dim(tc)[1], big.mark = ',', scientific = FALSE)
            val4 = format(dim(tc)[2], big.mark = ',', scientific = FALSE)
            cat(sprintf(txt, val1, val2, val3, val4))
            if (sum(is_negative) > 0) {
                tc[is_negative] = 0
            }
        } else {
            cat('`--clip_negative yes` is only applicable to `--norm log*p1-*`.\n')
        }
    }
    return(list(tc = tc, sva = sva1))
}

map_color = function(redundant_variables, c) {
    uniq_var = unique(redundant_variables)
    uniq_col = rainbow_hcl(length(uniq_var), c = c)
    df_unique = data.frame(var = uniq_var, col = uniq_col, stringsAsFactors = FALSE)
    df_redundant = data.frame(var = redundant_variables, order = seq(1, length(redundant_variables)),
                              stringsAsFactors = FALSE)
    df_redundant = merge(df_redundant, df_unique, by = "var", all.x = TRUE, stringsAsFactors = FALSE)
    df_redundant = df_redundant[order(df_redundant[['order']]),]
    return(df_redundant[['col']])
}

draw_heatmap = function(sra, tc_dist_matrix, legend = TRUE, fontsize = 7) {
    bp_fac = factor(sub(";.*", "", sra[, c("bioproject")]))
    sample_group_fac = factor(sra[, c("sample_group")])
    ann_label = data.frame(bioproject = bp_fac, sample_group = sample_group_fac)
    bp_col_uniq = unique(sra[order(sra[['bioproject']]), 'bp_color'])
    sample_group_col_uniq = unique(sra[order(sra[['sample_group']]), 'sample_group_color'])
    ann_color = list(bioproject = bp_col_uniq, sample_group = sample_group_col_uniq)
    breaks = c(0, seq(0.3, 1, 0.01))
    colnames(tc_dist_matrix) = sra[sra$run %in% colnames(tc_dist_matrix), 'run']
    head_half = substr(colnames(tc_dist_matrix), 1, 4)
    tail_half = substr(colnames(tc_dist_matrix), length(colnames(tc_dist_matrix)) - 3, length(colnames(tc_dist_matrix)))
    short_names = paste0(head_half, '..', tail_half)
    rownames(tc_dist_matrix) = ifelse(nchar(rownames(tc_dist_matrix)) <= 10, rownames(tc_dist_matrix), short_names)
    colnames(tc_dist_matrix) = ifelse(nchar(colnames(tc_dist_matrix)) <= 10, colnames(tc_dist_matrix), short_names)
    aheatmap(tc_dist_matrix, color = "-RdYlBu2:71", Rowv = NA, Colv = NA, revC = TRUE, legend = TRUE,
             breaks = breaks, annCol = ann_label, annRow = ann_label, annColors = ann_color, annLegend = legend,
             fontsize = fontsize)
}

color_children2parent = function(node) {
    if (length(node) == 2) {
        child1_color = attributes(node[[1]])[['edgePar']][["col"]]
        child2_color = attributes(node[[2]])[['edgePar']][["col"]]
        if ((!is.null(child1_color)) & (!is.null(child2_color))) {
            if (child1_color == child2_color) {
                attributes(node)[['edgePar']][["col"]] = child1_color
            }
        }
    }
    return(node)
}

draw_dendrogram = function(sra, tc_dist_dist, fontsize = 7) {
    dend = as.dendrogram(hclust(tc_dist_dist))
    dend_colors = sra[order.dendrogram(dend), 'sample_group_color']
    labels_colors(dend) = dend_colors
    dend_labels = sra[order.dendrogram(dend), 'run']
    dend = color_branches(dend, labels = dend_labels, col = dend_colors)
    dend = set(dend, "branches_lwd", 1)
    for (i in 1:ncol(tc)) {
        dend = dendrapply(dend, color_children2parent)
    }
    cex.xlab = min(fontsize, max(0.2, 0.5 / log10(nrow(sra)), na.rm = TRUE), na.rm = TRUE)
    par(cex = cex.xlab)
    plot(dend, las = 1, axes = FALSE)
    par(cex = 1)
    axis(side = 2, line = 0, las = 1)
    mtext("Distance", side = 2, line = 8.5, outer = FALSE)
    n = nrow(sra)
    symbols(
        1:n,
        rep(0, n),
        circles = rep(1, n),
        add = TRUE,
        inches = 0.02,
        xpd = TRUE,
        lwd = 1,
        bg = sra[order.dendrogram(dend), 'sample_group_color'],
        fg = sra[order.dendrogram(dend), 'bp_color']
    )
}

draw_dendrogram_ggplot = function(sra, tc_dist_dist, fontsize = 7) {

    cg_col = unique(sra[, c('sample_group', 'sample_group_color')])
    bp_col = unique(sra[, c('bioproject', 'bp_color')])
    colnames(cg_col) = c('Group', 'Color')
    colnames(bp_col) = c('Group', 'Color')
    sra_colors = rbind(cg_col, bp_col)

    group_colors <- data.table(Group = sra_colors$Group, Color = sra_colors$Color, key = "Group")
    group_colors <- transpose(group_colors, make.names = "Group")
    colnames(tc_dist_dist) <- sra$run
    hc <- hclust(tc_dist_dist)           # heirarchal clustering
    dendr <- dendro_data(hc, type = "rectangle") # convert for ggplot


    clust.df <- data.frame(label = sra$run, sample_group = factor(sra[sra$run %in% dendr$labels$label, 'sample_group']), bioproject = factor(sra[sra$run %in% dendr$labels$label, 'bioproject']))
    dendr[["labels"]] <- merge(dendr[["labels"]], clust.df, by = "label")
    ggplot() +
        geom_segment(data = dendr$segments, aes(x = x, y = y, xend = xend, yend = yend), size = .8, show.legend = FALSE) +
        geom_segment(data = merge(dendr$segments[dendr$segments$yend == 0,], dendr$labels[, c('label', 'sample_group', 'x')], by = 'x'), aes(x = x, y = y, xend = xend, yend = yend, color = sample_group), size = .8, show.legend = FALSE) +
        geom_text(data = dendr$labels, aes(x, y - .008, label = label, hjust = 0, angle = 270, color = sample_group), size = 3, show.legend = FALSE) +
        scale_y_continuous(expand = c(.2, .1)) +
        geom_point(data = dendr$labels, aes(x, y, color = bioproject), size = 3, show.legend = FALSE) +
        geom_point(data = dendr$labels, aes(x, y, color = sample_group), size = 2, show.legend = FALSE) +
        theme(axis.line.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              panel.background = element_rect(fill = "white"),
              panel.grid = element_blank()
        ) +
        scale_color_manual(values = group_colors)
}

draw_dendrogram_pvclust = function(sra, tc, nboot, pvclust_file, fontsize = 7) {

    dist_fun = function(x) {
        Dist(t(x), method = "pearson")
    }

    sp = sub(" ", "_", sra[['scientific_name']][1])
    if (file.exists(pvclust_file)) {
        if (file.info(pvclust_file)[['size']]) {
            cat("pvclust intermediate file found.\n")
            load(pvclust_file)
        }
    } else {
        cat("no pvclust intermediate file found. Start bootstrapping.\n")
        result = pvclust(tc, method.dist = dist_fun, method.hclust = "average", nboot = nboot, parallel = FALSE)  # UPGMA
        save(result, file = pvclust_file)
    }
    dend = as.dendrogram(result)
    dend_colors = sra[order.dendrogram(dend), 'sample_group_color']
    labels_colors(dend) = dend_colors
    dend_labels = sra[order.dendrogram(dend), 'run']
    dend = color_branches(dend, labels = dend_labels, col = dend_colors)
    dend = set(dend, "branches_lwd", 2)
    for (i in 1:ncol(tc)) {
        dend = dendrapply(dend, color_children2parent)
    }
    cex.xlab = min(0.2 + 1 / log10(fontsize), na.rm = TRUE)
    par(cex = cex.xlab)
    plot(dend, las = 1, ylab = "Distance", cex.axis = 1 / cex.xlab, cex.lab = 1 / cex.xlab)
    par(cex = 1)
    n = nrow(sra)
    symbols(
        1:n,
        rep(0, n),
        circles = rep(1, n),
        add = TRUE,
        inches = 0.04,
        xpd = TRUE,
        lwd = 2,
        bg = sra[order.dendrogram(dend), 'sample_group_color'],
        fg = sra[order.dendrogram(dend), 'bp_color']
    )
    text(result, print.num = FALSE, cex = 1, col.pv = "black")
}

draw_pca = function(sra, tc_dist_matrix, fontsize = 7) {
    pca = prcomp(tc_dist_matrix)
    xlabel = paste0("PC 1 (", round(summary(pca)[['importance']][2, 1] * 100, digits = 1), "%)")
    ylabel = paste0("PC 2 (", round(summary(pca)[['importance']][2, 2] * 100, digits = 1), "%)")
    plot(
        pca[['x']][, 1],
        pca[['x']][, 2],
        pch = 21,
        cex = 2,
        lwd = 1,
        bg = sra[['sample_group_color']],
        col = sra[['bp_color']],
        xlab = xlabel,
        ylab = ylabel,
        las = 1
    )
    # plot(pca$x[,1], pca$x[,2], pch=21, cex=2, lwd=2, bg=sra$sample_group_color, col=sra$bp_color, main=title,
    # xlab=xlabel, ylab=ylabel, las=1)
}

draw_mds = function(sra, tc_dist_dist, fontsize = 7) {
    try_out = tryCatch({
        isoMDS(tc_dist_dist, k = 2, maxit = 100)
    }, error = function(a) {
        return("MDS failed.")
    })
    if (mode(try_out) == "character") {
        cat("MDS failed.\n")
        plot(c(0, 1), c(0, 1), ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
    } else {
        mds <- try_out
        plot(
            mds[['points']][, 1],
            mds[['points']][, 2],
            pch = 21,
            cex = 2,
            lwd = 1,
            bg = sra[['sample_group_color']],
            col = sra[['bp_color']],
            xlab = "MDS dimension 1",
            ylab = "MDS dimension 2",
            las = 1
        )
    }
}

draw_tsne = function(sra, tc, fontsize = 7) {
    perplexity = min(30, floor(nrow(sra) / 4), na.rm = TRUE)
    try_out = tryCatch({
        Rtsne(
            as.matrix(t(tc)),
            theta = 0,
            check_duplicates = FALSE,
            verbose = FALSE,
            perplexity = perplexity,
            dims = 2
        )
    }, error = function(a) {
        return("t-SNE calculation failed.")
    })
    if (mode(try_out) == "character") {
        cat("t-SNE failed.\n")
        plot(c(0, 1), c(0, 1), ann = FALSE, bty = "n", type = "n", xaxt = "n", yaxt = "n")
        return(NULL)
    }
    out_tsne = try_out
    try_out = tryCatch({
        plot(
            out_tsne[['Y']][, 1],
            out_tsne[['Y']][, 2],
            pch = 21,
            cex = 2,
            lwd = 1,
            bg = sra[['sample_group_color']],
            col = sra[['bp_color']],
            xlab = "t-SNE dimension 1",
            ylab = "t-SNE dimension 2",
            las = 1
        )
    }, error = function(a) {
        return("t-SNE plot failed.")
    })
    if (mode(try_out) == "character") {
        cat("t-SNE failed.\n")
        plot(c(0, 1), c(0, 1), ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
    }
}

draw_sva_summary = function(sva_out, tc, sra, fontsize) {
    if ((is.null(sva_out)) | (class(sva_out) == "try-error")) {
        plot(c(0, 1), c(0, 1), ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
        df = NA
    } else {
        out = tc_metadata_intersect(tc, sra)
        tc = out[["tc"]]
        sra = out[["sra"]]
        out = sort_tc_and_metadata(tc, sra)
        tc = out[["tc"]]
        sra = out[["sra"]]
        sra[['log10_total_spots']] = log10(sra[['total_spots']])
        sra[['log10_total_bases']] = log10(sra[['total_bases']])
        cols = c("sample_group", "bioproject", "lib_layout", "lib_selection", "instrument", "mapping_rate", 'log10_total_spots', 'log10_total_bases')
        label_cols = c("Group", "BioProject", "Library layout", "Library selection", "Instrument", "Mapping rate", 'Log10 total reads', 'Log10 total bases')
        if ('tmm_normalization_factor' %in% colnames(sra)) {
            cols = c(cols, 'tmm_normalization_factor')
            label_cols = c(label_cols, 'TMM normalization factor')
        }
        num_sv = sva_out[['n.sv']]
        if (num_sv == 0) {
            cat('No surrogate variables found.\n')
            plot(c(0, 1), c(0, 1), ann = FALSE, bty = "n", type = "n", xaxt = "n", yaxt = "n")
            return(data.frame())
        }
        df = data.frame(matrix(NA, num_sv, length(cols)))
        colnames(df) = cols
        rownames(df) = paste0("SV", 1:nrow(df))
        for (i in 1:length(cols)) {
            for (j in 1:num_sv) {
                if (length(unique(sra[, cols[i]])) == 1) {
                    df[j, i] = NA
                } else {
                    lm_summary = summary(lm(sva_out[['sv']][, j] ~ sra[, cols[i]]))
                    df[j, i] = lm_summary[['adj.r.squared']]
                }
            }
        }
        colnames(df) = label_cols
        breaks = seq(0, 1, 0.02)
        colors = colorRampPalette(c("blue", "yellow", "red"))(length(breaks))
        df2 = t(df)
        df2[df2 < 0] = 0
        aheatmap(df2, color = colors, Rowv = NA, Colv = NA, revC = TRUE, breaks = breaks, fontsize = fontsize)
    }
    return(df)
}

draw_boxplot = function(sra, tc_dist_matrix, fontsize = 7) {
    is_same_bp = outer(sra[['bioproject']], sra[['bioproject']], function(x, y) {
        x == y
    })
    is_same_sample_group = outer(sra[['sample_group']], sra[['sample_group']], function(x, y) {
        x == y
    })
    plot(c(0.5, 4.5), c(0, 1), type = "n", xlab = "", ylab = "Pearson's correlation\ncoefficient", las = 1,
         xaxt = "n")
    boxplot(tc_dist_matrix[(!is_same_bp) & (!is_same_sample_group)], at = 1, add = TRUE, col = "gray", yaxt = "n")
    boxplot(tc_dist_matrix[(is_same_bp) & (!is_same_sample_group)], at = 2, add = TRUE, col = "gray", yaxt = "n")
    boxplot(tc_dist_matrix[(!is_same_bp) & (is_same_sample_group)], at = 3, add = TRUE, col = "gray", yaxt = "n")
    boxplot(tc_dist_matrix[(is_same_bp) & (is_same_sample_group)], at = 4, add = TRUE, col = "gray", yaxt = "n")
    labels = c("bw\nbw", "bw\nwi", "wi\nbw", "wi\nwi")
    axis(side = 1, at = c(1, 2, 3, 4), labels = labels, padj = 0.5)
    axis(side = 1, at = 0.35, labels = "Sample Group\nBioProject", padj = 0.5, hadj = 1, tick = FALSE)

    # Add mean PCC
    means <- c(
        mean(tc_dist_matrix[(!is_same_bp) & (!is_same_sample_group)], na.rm = TRUE),
        mean(tc_dist_matrix[(is_same_bp) & (!is_same_sample_group)], na.rm = TRUE),
        mean(tc_dist_matrix[(!is_same_bp) & (is_same_sample_group)], na.rm = TRUE),
        mean(tc_dist_matrix[(is_same_bp) & (is_same_sample_group)], na.rm = TRUE)
    )

    points(1:4, means, col = "red", pch = 16)
    lines(c(1, 2), means[1:2], col = "red")
    lines(c(3, 4), means[3:4], col = "red")

    text(x = 1.5, y = max(means[1:2]) + 0.05, labels = round(abs(means[1] - means[2]), 2), col = "red", cex = 0.8)
    text(x = 3.5, y = max(means[3:4]) + 0.05, labels = round(abs(means[3] - means[4]), 2), col = "red", cex = 0.8)

    labels = c("bw\nbw", "bw\nwi", "wi\nbw", "wi\nwi")
    axis(side = 1, at = c(1, 2, 3, 4), labels = labels, padj = 0.5)
    axis(side = 1, at = 0.35, labels = "Sample Group\nBioProject", padj = 0.5, hadj = 1, tick = FALSE)

    # Add legend in the bottom left corner
    legend("bottomleft", legend = c("mean PCC", expression(Delta ~ "mean PCC")),
           pch = c(16, NA), col = c("red", "red"),
           lty = c(NA, 1), bty = "n", cex = 0.8,
           text.width = max(strwidth(c("mean PCC", expression(Delta ~ "mean PCC")))))  # Align text

}

save_correlation = function(tc, sra, dist_method, round) {
    out = tc_metadata_intersect(tc, sra)
    tc = out[["tc"]]
    sra = out[["sra"]]
    is_same_bp = outer(sra[['bioproject']], sra[['bioproject']], function(x, y) {
        x == y
    })
    is_same_sample_group = outer(sra[['sample_group']], sra[['sample_group']], function(x, y) {
        x == y
    })


    tc_dist_matrix = cor(tc, method = dist_method)
    tc_dist_matrix[is.na(tc_dist_matrix)] = 0

    # bw = between; wi = within; order: Bioproject -> Sample Group; tc_dist_bwwi: between BPs, within CGs
    tc_dist_bwbw = tc_dist_matrix[(!is_same_bp) & (!is_same_sample_group)]
    bwbw_mea = mean(tc_dist_bwbw)
    bwbw_med = median(tc_dist_bwbw)
    bwbw_var = var(tc_dist_bwbw)

    tc_dist_wibw = tc_dist_matrix[(is_same_bp) & (!is_same_sample_group)]
    wibw_mea = mean(tc_dist_wibw)
    wibw_med = median(tc_dist_wibw)
    wibw_var = var(tc_dist_wibw)

    tc_dist_bwwi = tc_dist_matrix[(!is_same_bp) & (is_same_sample_group)]
    bwwi_mea = mean(tc_dist_bwwi)
    bwwi_med = median(tc_dist_bwwi)
    bwwi_var = var(tc_dist_bwwi)

    tc_dist_wiwi = tc_dist_matrix[(is_same_bp) & (is_same_sample_group)]
    wiwi_mea = mean(tc_dist_wiwi)
    wiwi_med = median(tc_dist_wiwi)
    wiwi_var = var(tc_dist_wiwi)

    tc_dist_stats = c(bwbw_mea, bwbw_med, bwbw_var, wibw_mea, wibw_med, wibw_var, bwwi_mea, bwwi_med, bwwi_var, wiwi_mea, wiwi_med, wiwi_var)

    # Check if dataframe exists in environment, if not create it
    if (!exists("correlation_statistics", envir = .GlobalEnv)) {
        correlation_statistics <- data.frame(matrix(tc_dist_stats, ncol = 12, dimnames = list(NULL, c("bwbw_mean", "bwbw_median", "bwbw_variance", "wibw_mean", "wibw_median", "wibw_variance", "bwwi_mean", "bwwi_median", "bwwi_variance", "wiwi_mean", "wiwi_median", "wiwi_variance"))))
        rownames(correlation_statistics) <- paste0("round_", round)
        assign("correlation_statistics", correlation_statistics, envir = .GlobalEnv)
    } else {
        correlation_statistics <- get("correlation_statistics", envir = .GlobalEnv)
        new_row <- data.frame(matrix(tc_dist_stats, ncol = 12, dimnames = list(NULL, c("bwbw_mean", "bwbw_median", "bwbw_variance", "wibw_mean", "wibw_median", "wibw_variance", "bwwi_mean", "bwwi_median", "bwwi_variance", "wiwi_mean", "wiwi_median", "wiwi_variance"))))
        rownames(new_row) <- paste0("round_", round)
        correlation_statistics <- rbind(correlation_statistics, new_row)
        assign("correlation_statistics", correlation_statistics, envir = .GlobalEnv)
    }

    return(correlation_statistics)


}


draw_tau_histogram = function(tc, sra, selected_sample_groups, fontsize = 7, transform_method) {
    df_tau = sample_group2tau(sample_group_mean(tc, sra, selected_sample_groups)[['tc_ave']], rich.annotation = FALSE, transform_method)
    hist_out = hist(df_tau[['tau']], breaks = seq(0, 1, 0.05), las = 1, xlab = "Tau (expression specificity)",
                    ylab = "Gene count", main = "", col = "gray")
    num_noexp = sum(is.na(df_tau[['tau']]))
    num_all = nrow(df_tau)
    text_noexp = paste0("Excluded due to\nno expression:\n", num_noexp, "/", num_all, " genes")
    text(0, max(hist_out[['counts']], na.rm = TRUE) * 0.85, text_noexp, pos = 4)
}

draw_exp_level_histogram = function(tc, sra, selected_sample_groups, fontsize = 7, transform_method) {
    tc_sample_group = sample_group_mean(tc, sra, selected_sample_groups)[['tc_ave']]
    xmax = apply(tc_sample_group, 1, max)
    xmax[xmax < 0] = 0
    xmax[xmax > 15] = 15
    breaks = seq(0, 15, 1)
    hist_out = hist(xmax, breaks = breaks, las = 1, xlab = paste0("Max expression (", transform_method, ")"), ylab = "Gene count",
                    main = "", col = "gray")
}

draw_legend = function(sra, new = TRUE, pos = "center", fontsize = 7, nlabel.in.col) {
    if (new) {
        plot.new()
    }
    sample_group_unique = unique(sra[['sample_group']])
    bp_unique = unique(sub(";.*", "", sra[['bioproject']]))
    sample_group_color_unique = unique(sra[['sample_group_color']])
    bp_color_unique = unique(sra[['bp_color']])
    ncol = ceiling((length(sample_group_unique) +
        length(bp_unique) +
        2) / nlabel.in.col)
    legend_text = c("Sample Group", as.character(sample_group_unique), "", "BioProject", as.character(bp_unique))
    legend_color = c(rgb(1, 1, 1, 0), rep(rgb(1, 1, 1, 0), length(sample_group_color_unique)), rgb(1, 1, 1,
                                                                                                   0), rgb(1, 1, 1, 0), bp_color_unique)
    legend_bg = c(rgb(1, 1, 1, 0), sample_group_color_unique, rgb(1, 1, 1, 0), rgb(1, 1, 1, 0), rep(rgb(1,
                                                                                                        1, 1, 0), length(bp_color_unique)))
    legend_font = c(2, rep(1, length(sample_group_color_unique)), 1, 2, rep(1, length(bp_color_unique)))
    legend(pos, legend = legend_text, pch = 21, lwd = 1, lty = 0, col = legend_color, pt.bg = legend_bg,
           text.font = legend_font, ncol = ncol, bty = "n")
}

save_plot = function(tc, sra, sva_out, dist_method, file, selected_sample_groups, sample_group_colors, fontsize = 7, transform_method, batch_effect_alg) {
    if (ncol(tc) == 1) {
        cat('Only 1 sample is available. Skipping the plot.\n')
        return()
    }
    out = tc_metadata_intersect(tc, sra)
    tc = out[["tc"]]
    sra = out[["sra"]]
    sra = add_color_to_curate_metadata(sra, selected_sample_groups, sample_group_colors)
    out = sort_tc_and_metadata(tc, sra)
    tc = out[["tc"]]
    sra = out[["sra"]]
    pdf(file.path(dir_pdf, paste0(file, ".pdf")), height = 8, width = 7.2, fonts = "Helvetica", pointsize = fontsize)
    layout_matrix = matrix(c(2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1,
                             2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1,
                             1, 1, 1, 1, 1, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 3, 3,
                             3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 9, 9, 7, 7, 7, 8, 8, 8, 9, 9, 9,
                             9, 9, 9, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 9, 9, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 9, 9, 10, 10, 10,
                             10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), 14, 12,
                           byrow = TRUE)
    layout(layout_matrix)
    colnames(tc) <- sra$run
    ##
    tc_dist_matrix = cor(tc, method = dist_method)
    tc_dist_matrix[is.na(tc_dist_matrix)] = 0
    tc_dist_dist = Dist(t(tc), method = dist_method) + 1e-09
    tc_dist_dist[is.na(tc_dist_dist)] = 1
    par(mar = c(6, 6, 1, 0))
    draw_dendrogram(sra, tc_dist_dist, fontsize)
    par(mar = c(0, 0, 0, 0))
    draw_heatmap(sra, tc_dist_matrix, legend = FALSE)
    # draw_dendrogram(sra, tc, nboot=1000, cex.xlab=0.6, pvclust_file=paste0(file, '.pvclust.RData'))
    par(mar = c(4, 4, 0.1, 1))
    draw_pca(sra, tc_dist_matrix, fontsize)
    par(mar = c(4, 4, 0.1, 1))
    draw_mds(sra, tc_dist_dist, fontsize)
    par(mar = c(4, 4, 0.1, 1))
    draw_tsne(sra, tc, fontsize)
    par(mar = c(4, 5, 0.1, 1))
    draw_boxplot(sra, tc_dist_matrix, fontsize)
    par(mar = c(4, 4, 1, 1))
    draw_exp_level_histogram(tc, sra, selected_sample_groups, fontsize, transform_method)
    par(mar = c(4, 4, 1, 1))
    draw_tau_histogram(tc, sra, selected_sample_groups, fontsize, transform_method)
    par(mar = rep(0.1, 4))
    if (batch_effect_alg == 'sva') {
        df_r2 = draw_sva_summary(sva_out, tc, sra, fontsize)
        if (!all(is.na(df_r2))) {
            write.table(df_r2, file.path(dir_tsv, paste0(file, ".r2.tsv")), sep = "\t", row.names = FALSE)
        }
    }
    par(mar = rep(0.1, 4))
    draw_legend(sra, new = TRUE, pos = "center", fontsize = fontsize, nlabel.in.col = 8)
    graphics.off()
}

transform_raw_to_fpkm = function(counts, effective_lengths, sra) {
    if ('tmm_library_size' %in% colnames(sra)) {
        cat('FPKM transformation with the original library sizes from amalgkit cstmm output.\n')
        cat('If --input_dir is specified with amalgkit cstmm output, resultant values will be TMM-FPKM.\n')
        cat('If --input_dir is specified with amalgkit merge output, resultant values will be non-TMM-corrected FPKM.\n')
        tmp = sra
        rownames(tmp) = tmp[['run']]
        library_sizes = tmp[colnames(counts), 'tmm_library_size']
    } else {
        cat('FPKM transformation with the library sizes in the input files.\n')
        cat('Irrespective of --input_dir, resultant values will be non-TMM-corrected FPKM.\n')
        library_sizes = colSums(counts)
    }
    res = counts / effective_lengths / library_sizes * 1e+09 # 1e+09 = kb * million reads
    return(as.data.frame(res))
}

transform_raw_to_tpm = function(counts, effective_lengths) {
    x <- counts / effective_lengths
    res = t(t(x) * 1e+06 / colSums(x))
    return(res)
}

apply_transformation_logic = function(tc, tc_eff_length, transform_method, batch_effect_alg,
                                      step = c('before_batch', 'before_batch_plot', 'after_batch'), sra) {
    if (batch_effect_alg == 'no') {
        if (step == 'before_batch') {
            bool_fpkm_tpm = TRUE
            bool_log = TRUE
        } else if (step == 'before_batch_plot') {
            bool_fpkm_tpm = FALSE
            bool_log = FALSE
        } else if (step == 'after_batch') {
            bool_fpkm_tpm = FALSE
            bool_log = FALSE
        }
    } else if (batch_effect_alg == 'sva') {
        if (step == 'before_batch') {
            bool_fpkm_tpm = TRUE
            bool_log = TRUE
        } else if (step == 'before_batch_plot') {
            bool_fpkm_tpm = FALSE
            bool_log = FALSE
        } else if (step == 'after_batch') {
            bool_fpkm_tpm = FALSE
            bool_log = FALSE
        }
    } else if (batch_effect_alg == 'ruvseq') {
        if (step == 'before_batch') {
            bool_fpkm_tpm = FALSE
            bool_log = FALSE
        } else if (step == 'before_batch_plot') {
            bool_fpkm_tpm = TRUE
            bool_log = TRUE
        } else if (step == 'after_batch') {
            bool_fpkm_tpm = TRUE
            bool_log = TRUE
        }
    } else if (batch_effect_alg == 'combatseq') {
        if (step == 'before_batch') {
            bool_fpkm_tpm = FALSE
            bool_log = FALSE
        } else if (step == 'before_batch_plot') {
            bool_fpkm_tpm = TRUE
            bool_log = TRUE
        } else if (step == 'after_batch') {
            bool_fpkm_tpm = TRUE
            bool_log = TRUE
        }
    }
    if (bool_fpkm_tpm == TRUE) {
        if (grepl('fpkm', transform_method)) {
            cat('Applying FPKM transformation.\n')
            tc = transform_raw_to_fpkm(tc, tc_eff_length[, colnames(tc)], sra)
        } else if (grepl('tpm', transform_method)) {
            cat('Applying TPM transformation.\n')
            tc = transform_raw_to_tpm(tc, tc_eff_length[, colnames(tc)])
        } else {
            cat('Applying neither FPKM nor TPM transformation.\n')
        }
    }
    if (bool_log == TRUE) {
        if (grepl('logn-', transform_method)) {
            cat('Applying log_n(x) normalization.\n')
            tc = log(tc)
        } else if (grepl('log2-', transform_method)) {
            cat('Applying log_2(x) normalization.\n')
            tc = log2(tc)
        } else if (grepl('lognp1-', transform_method)) {
            cat('Applying log_n(x+1) normalization.\n')
            tc = log(tc + 1)
        } else if (grepl('log2p1-', transform_method)) {
            cat('Applying log_2(x+1) normalization.\n')
            tc = log2(tc + 1)
        } else {
            cat('Applying no log normalization.\n')
        }
    }
    return(tc)
}

standardize_metadata_all = function(sra_all) {
    for (col in c('instrument', 'bioproject')) {
        if (!col %in% colnames(sra_all)) {
            next
        }
        is_missing = (sra_all[, col] == "") | (is.na(sra_all[, col]))
        sra_all[is_missing, col] = "not_provided"
    }
    return(sra_all)
}

get_species_metadata = function(sra_all, scientific_name, selected_sample_groups) {
    is_sp = (sra_all[, 'scientific_name'] == scientific_name)
    is_sample_group = (sra_all[, 'sample_group'] %in% selected_sample_groups)
    cat('Number of SRA runs for this species:', sum(is_sp), '\n')
    cat('Number of SRA runs for selected sample groups:', sum(is_sample_group), '\n')
    sra = sra_all[(is_sp & is_sample_group),]
    conditions = (sra[['exclusion']] == "no") & (!sra[['run']] %in% colnames(tc))
    if (any(conditions)) {
        cat("Failed quantification:", sra[conditions, "run"], "\n")
        sra[conditions, "exclusion"] = "failed_quantification"
    }
    return(sra)
}

exclude_inappropriate_sample_from_tc = function(tc, sra) {
    is_not_excluded = (sra[['exclusion']] == 'no')
    cat('Number of non-excluded SRA runs (exclusion=="no"):', sum(is_not_excluded), '\n')
    tc = tc[, sra[is_not_excluded, 'run'], drop = FALSE]
    return(tc)
}

exclude_inappropriate_sample_from_eff_length = function(tc_eff_length, tc) {
    tc_eff_length = tc_eff_length[, colnames(tc), drop = FALSE]
    return(tc_eff_length)
}

########cd############END OF FUNCTION DECLARATION####################################################

################################################# START OF SVA CORRECTION ########################################################


fontsize = 7

tc = read.table(est_counts_path, sep = "\t", stringsAsFactors = FALSE, header = TRUE, quote = "", fill = FALSE, row.names = 1, check.names = FALSE)
tc_eff_length = read.table(eff_length_path, sep = "\t", stringsAsFactors = FALSE, header = TRUE, quote = "", fill = FALSE, row.names = 1, check.names = FALSE)

sra_all = read.table(metadata_path, sep = "\t", header = TRUE, quote = "", fill = TRUE, comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
sra_all = standardize_metadata_all(sra_all)

scientific_name = sra_all[(sra_all[['run']] %in% colnames(tc)), "scientific_name"][1]
dir_curate = file.path(out_dir, 'curate')
dir_pdf = file.path(dir_curate, sub(" ", "_", scientific_name), 'plots')
dir.create(dir_pdf, showWarnings = FALSE, recursive = TRUE)
dir_rdata = file.path(dir_curate, sub(" ", "_", scientific_name), 'rdata')
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir_tsv = file.path(dir_curate, sub(" ", "_", scientific_name), 'tables')
dir.create(dir_tsv, showWarnings = FALSE, recursive = TRUE)
setwd(dir_curate)
cat(log_prefix, "Working at:", getwd(), "\n")

sra = get_species_metadata(sra_all, scientific_name, selected_sample_groups)
tc = exclude_inappropriate_sample_from_tc(tc, sra)
out = sort_tc_and_metadata(tc, sra); tc = out[["tc"]]; sra = out[["sra"]]
tc_eff_length = exclude_inappropriate_sample_from_eff_length(tc_eff_length, tc)
tc = apply_transformation_logic(tc, tc_eff_length, transform_method, batch_effect_alg, step = 'before_batch', sra = sra)
tc_tmp = apply_transformation_logic(tc, tc_eff_length, transform_method, batch_effect_alg, step = 'before_batch_plot', sra = sra)
is_input_zero = data.frame(tc_tmp == 0, check.names = FALSE)
file_name = file.path(dir_tsv, paste0(sub(" ", "_", scientific_name), ".uncorrected.tc.tsv"))
write_table_with_index_name(df = tc_tmp, file_path = file_name, index_name = 'target_id')
original_sample_groups = selected_sample_groups
out = sample_group_mean(tc_tmp, sra, selected_sample_groups)
tc_sample_group_uncorrected = out[['tc_ave']]
selected_sample_groups = out[['selected_sample_groups']]
if (length(selected_sample_groups) != length(original_sample_groups)) {
    sample_group_colors = sample_group_colors[match(selected_sample_groups, original_sample_groups)]
}
file_name = file.path(dir_tsv, paste0(sub(" ", "_", scientific_name), ".uncorrected.sample_group.mean.tsv"))
write_table_with_index_name(df = tc_sample_group_uncorrected, file_path = file_name, index_name = 'target_id')

if (skip_curation_flag == TRUE) {
    cat("No curation requested, finishing early.\n")
    cat("Files created: \n")
    cat(file.path(dir_tsv, paste0(sub(" ", "_", scientific_name), ".uncorrected.tc.tsv")), "\n")
    cat(file.path(dir_tsv, paste0(sub(" ", "_", scientific_name), ".uncorrected.sample_group.mean.tsv")), "\n")
    cat("Transformation applied: ", transform_method, "\n")
    cat(log_prefix, "Completed.\n")
    quit(save = 'no', status = 0)
}

cat("Removing samples with mapping rate of 0.\n")
round = 0
sva_out = NULL
tc_batch_corrected = NULL
out = check_mapping_rate(tc, sra, 0)
tc = out[["tc"]]
sra = out[["sra"]]
tc_tmp = apply_transformation_logic(tc, tc_eff_length, transform_method, batch_effect_alg, step = 'before_batch_plot', sra = sra)
save_plot(tc_tmp, sra, NULL, dist_method, paste0(sub(" ", "_", scientific_name), ".", round, ".original"),
          selected_sample_groups, sample_group_colors, fontsize, transform_method, batch_effect_alg)
save_correlation(tc_tmp, sra, dist_method, round)
out = batch_effect_subtraction(tc, sra, batch_effect_alg, transform_method, clip_negative)
tc_batch_corrected = out[["tc"]]
if (batch_effect_alg == "combatseq") {
    tc = tc[, colnames(tc_batch_corrected)]
}
sva_out = out[["sva"]]
if (!is.null(sva_out)) {
    file_name = paste0(sub(" ", "_", scientific_name), ".", batch_effect_alg, ".", round, ".RData")
    save(sva_out, file = file.path(dir_rdata, file_name))
}
tc_batch_corrected_tmp = apply_transformation_logic(tc_batch_corrected, tc_eff_length, transform_method, batch_effect_alg, step = 'after_batch', sra = sra)
save_plot(tc_batch_corrected_tmp, sra, sva_out, dist_method, paste0(sub(" ", "_", scientific_name), ".", round, ".original", ".", batch_effect_alg),
          selected_sample_groups, sample_group_colors, fontsize, transform_method, batch_effect_alg)

cat("Removing samples with the mapping rate smaller than", mapping_rate_cutoff, "\n")
round = 1
sva_out = NULL
tc_batch_corrected = NULL
out = check_mapping_rate(tc, sra, mapping_rate_cutoff)
tc = out[["tc"]]
sra = out[["sra"]]

tc_tmp = apply_transformation_logic(tc, tc_eff_length, transform_method, batch_effect_alg, step = 'before_batch_plot', sra = sra)
save_plot(tc_tmp, sra, NULL, dist_method, paste0(sub(" ", "_", scientific_name), ".", round, ".mapping_cutoff"),
          selected_sample_groups, sample_group_colors, fontsize, transform_method, batch_effect_alg)
out = batch_effect_subtraction(tc, sra, batch_effect_alg, transform_method, clip_negative)
tc_batch_corrected = out[["tc"]]
if (batch_effect_alg == "combatseq") {
    tc = tc[, colnames(tc_batch_corrected)]
}
sva_out = out[["sva"]]
if (!is.null(sva_out)) {
    save(sva_out, file = file.path(dir_rdata, paste0(sub(" ", "_", scientific_name), ".", batch_effect_alg, ".", round, ".RData")))
}
tc_batch_corrected_tmp = apply_transformation_logic(tc_batch_corrected, tc_eff_length, transform_method, batch_effect_alg, step = 'after_batch', sra = sra)
save_plot(tc_batch_corrected_tmp, sra, sva_out, dist_method, paste0(sub(" ", "_", scientific_name), ".", round, ".mapping_cutoff", ".", batch_effect_alg),
          selected_sample_groups, sample_group_colors, fontsize, transform_method, batch_effect_alg)
save_correlation(tc_batch_corrected_tmp, sra, dist_method, round)

round = 2
end_flag = 0
while (end_flag == 0) {
    cat("Iteratively checking within-sample_group correlation, round:", round, "\n")
    tc_cwtc = NULL
    num_run_before = sum(sra[['exclusion']] == "no")
    out = check_within_sample_group_correlation(tc, sra, dist_method, min_dif, selected_sample_groups, one_outlier_per_iteration, correlation_threshold)
    tc_cwtc = out[["tc"]]
    sra = out[["sra"]]
    num_run_after = sum(sra[['exclusion']] == "no")
    if ((num_run_before == num_run_after) | (plot_intermediate)) {
        sva_out = NULL
        tc_batch_corrected = NULL
        out = batch_effect_subtraction(tc[, colnames(tc_cwtc)], sra, batch_effect_alg, transform_method, clip_negative)
        tc_batch_corrected = out[["tc"]]
        if (batch_effect_alg == "combatseq") {
            tc = tc[, colnames(tc_batch_corrected)]
        }
        sva_out = out[["sva"]]
        if (!is.null(sva_out)) {
            save(sva_out, file = file.path(dir_rdata, paste0(sub(" ", "_", scientific_name), ".", batch_effect_alg, ".", round, ".RData")))
        }
        tc_cwtc_tmp = apply_transformation_logic(tc_cwtc, tc_eff_length, transform_method, batch_effect_alg, step = 'before_batch_plot', sra = sra)
        save_plot(tc_cwtc_tmp, sra, NULL, dist_method, paste0(sub(" ", "_", scientific_name), ".", round, ".correlation_cutoff"),
                  selected_sample_groups, sample_group_colors, fontsize, transform_method, batch_effect_alg)
        tc_batch_corrected_tmp = apply_transformation_logic(tc_batch_corrected, tc_eff_length, transform_method, batch_effect_alg, step = 'after_batch', sra = sra)
        file_base = paste0(sub(" ", "_", scientific_name), ".", round, ".correlation_cutoff", ".", batch_effect_alg)
        save_plot(tc_batch_corrected_tmp, sra, sva_out, dist_method, file_base, selected_sample_groups, sample_group_colors, fontsize, transform_method, batch_effect_alg)
        save_correlation(tc_batch_corrected_tmp, sra, dist_method, round)
    }
    cat("Round:", round, ": # before =", num_run_before, ": # after =", num_run_after, "\n\n")
    if (num_run_before == num_run_after) {
        end_flag = 1
    }
    tc = tc_cwtc
    round = round + 1
}

cat("Finished checking within-sample_group correlation.\n")
if (batch_effect_alg != 'sva') {
    cat("Batch-effect removal algorithm is: ", batch_effect_alg, ". Applying transformation on final batch-removed counts.\n")
    tc_batch_corrected = tc_batch_corrected_tmp
}
if (maintain_zero) {
    cat('Any zero expression levels in the input will remain as zero-values in the output tables.\n')
    tc_batch_corrected = tc_batch_corrected[order(rownames(tc_batch_corrected)),]
    is_input_zero = is_input_zero[rownames(tc_batch_corrected),]
    stopifnot(all(rownames(is_input_zero) == rownames(tc_batch_corrected)))
    for (col in colnames(tc_batch_corrected)) {
        tc_batch_corrected[is_input_zero[[col]], col] = 0
    }
}
cat("Writing summary files for", scientific_name, "\n")
file = file.path(dir_tsv, paste0(sub(" ", "_", scientific_name), ".metadata.tsv"))
write.table(sra[, colnames(sra) != 'index'], file = file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
file = file.path(dir_tsv, paste0(sub(" ", "_", scientific_name), ".", batch_effect_alg, ".tc.tsv"))
write_table_with_index_name(df = tc_batch_corrected, file_path = file, index_name = 'target_id')
out = sample_group_mean(tc_batch_corrected, sra, selected_sample_groups)
tc_sample_group = out[['tc_ave']]
file = file.path(dir_tsv, paste0(sub(" ", "_", scientific_name), ".", batch_effect_alg, ".sample_group.mean.tsv"))
write_table_with_index_name(df = tc_sample_group, file_path = file, index_name = 'target_id')
file = file.path(dir_tsv, paste0(sub(" ", "_", scientific_name), ".", batch_effect_alg, ".correlation_statistics.tsv"))
write.table(correlation_statistics, file = file, row.names = TRUE, sep = "\t", quote = FALSE)
tc_tau = sample_group2tau(tc_sample_group, rich.annotation = TRUE, transform_method)
file = file.path(dir_tsv, paste0(sub(" ", "_", scientific_name), ".", batch_effect_alg, ".tau.tsv"))
write_table_with_index_name(df = tc_tau, file_path = file, index_name = 'target_id')
cat(log_prefix, "Completed.\n")
quit(save = 'no', status = 0)



================================================
FILE: amalgkit/getfastq.py
================================================
from Bio import Entrez
import itertools
import numpy

from amalgkit.util import *

import glob
import lxml
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from urllib.error import HTTPError

def getfastq_search_term(ncbi_id, additional_search_term=None):
    # https://www.ncbi.nlm.nih.gov/books/NBK49540/
    if additional_search_term is None:
        search_term = ncbi_id
    else:
        search_term = ncbi_id + ' AND ' + additional_search_term
    return search_term

def getfastq_getxml(search_term, retmax=1000):
    entrez_db = 'sra'
    try:
        sra_handle = Entrez.esearch(db=entrez_db, term=search_term, retmax=10000000)
    except HTTPError as e:
        print(e, '- Trying Entrez.esearch() again...')
        sra_handle = Entrez.esearch(db=entrez_db, term=search_term, retmax=10000000)
    sra_record = Entrez.read(sra_handle)
    record_ids = sra_record["IdList"]
    num_record = len(record_ids)
    print('Number of SRA records:', num_record)
    root = None
    for i in numpy.arange(numpy.ceil(num_record // retmax) + 1):
        start = int(i * retmax)
        end = int(((i + 1) * retmax) - 1) if num_record >= int(((i + 1) * retmax) - 1) else num_record
        print('processing SRA records:', start, '-', end, flush=True)
        try:
            handle = Entrez.efetch(db="sra", id=record_ids[start:end], rettype="full", retmode="xml", retmax=retmax)
        except HTTPError as e:
            print(e, '- Trying Entrez.efetch() again...')
            handle = Entrez.efetch(db="sra", id=record_ids[start:end], rettype="full", retmode="xml", retmax=retmax)
        chunk = lxml.etree.parse(handle).getroot()
        if root is None:
            root = chunk
        else:
            root.append(chunk)
    xml_string = lxml.etree.tostring(root, pretty_print=True)
    for line in str(xml_string).split('\n'):
        if '<Error>' in line:
            print(line)
            raise Exception('\<Error> found in the xml. Search term: '+search_term)
    return root

def get_range(sra_stat, offset, total_sra_bp, max_bp):
    if (total_sra_bp <= max_bp):
        start = 1
        end = sra_stat['total_spot']
    else:
        if (sra_stat['total_spot'] > (sra_stat['num_read_per_sra'] + offset)):
            start = offset
            end = offset + sra_stat['num_read_per_sra']
        elif (sra_stat['total_spot'] > sra_stat['num_read_per_sra']):
            start = sra_stat['total_spot'] - sra_stat['num_read_per_sra']
            end = sra_stat['total_spot']
        elif (sra_stat['total_spot'] <= sra_stat['num_read_per_sra']):
            start = 1
            end = sra_stat['total_spot']
    return start, end

def concat_fastq(args, metadata, output_dir, g):
    layout = get_layout(args, metadata)
    inext = '.amalgkit.fastq.gz'
    infiles = list()
    for sra_id in metadata.df.loc[:, 'run']:
        infiles.append([f for f in os.listdir(output_dir) if (f.endswith(inext)) & (f.startswith(sra_id))])
    infiles = [item for sublist in infiles for item in sublist]
    num_inext_files = len(infiles)
    if (layout == 'single') & (num_inext_files == 1):
        print('Only 1', inext, 'file was detected. No concatenation will happen.', flush=True)
        if args.id is not None:
            outfile = args.id + infiles[0]
        elif args.id_list is not None:
            outfile = os.path.basename(args.id_list) + infiles[0]
        if infiles[0] != outfile:
            print('Replacing ID in the output file name:', infiles[0], outfile)
            infile_path = os.path.join(output_dir, infiles[0])
            outfile_path = os.path.join(output_dir, outfile)
            os.rename(infile_path, outfile_path)
        return None
    elif (layout == 'paired') & (num_inext_files == 2):
        print('Only 1 pair of', inext, 'files were detected. No concatenation will happen.', flush=True)
        for infile in infiles:
            if args.id is not None:
                outfile = args.id + re.sub('.*(_[1-2])', '\g<1>', infile)
            elif args.id_list is not None:
                outfile = os.path.basename(args.id_list) + re.sub('.*(_[1-2])', '\g<1>', infile)
            if infile != outfile:
                print('Replacing ID in the output file name:', infile, outfile)
                infile_path = os.path.join(output_dir, infile)
                outfile_path = os.path.join(output_dir, outfile)
                os.rename(infile_path, outfile_path)
        return None
    else:
        print('Concatenating files with the extension:', inext)
        outext = '.amalgkit.fastq.gz'
        if layout == 'single':
            subexts = ['', ]
        elif layout == 'paired':
            subexts = ['_1', '_2', ]
        for subext in subexts:
            infiles = metadata.df['run'].replace('$', subext + inext, regex=True)
            if args.id is not None:
                outfile_path = os.path.join(output_dir, args.id + subext + outext)
            elif args.id_list is not None:
                outfile_path = os.path.join(output_dir, os.path.basename(args.id_list) + subext + outext)
            if os.path.exists(outfile_path):
                os.remove(outfile_path)
            for infile in infiles:
                infile_path = os.path.join(output_dir, infile)
                assert os.path.exists(infile_path), 'Dumped fastq not found: ' + infile_path
                print('Concatenated file:', infile_path, flush=True)
                os.system('cat "' + infile_path + '" >> "' + outfile_path + '"')
            print('')
        if args.remove_tmp:
            for i in metadata.df.index:
                sra_id = metadata.df.loc[i, 'run']
                sra_stat = get_sra_stat(sra_id, metadata, g['num_bp_per_sra'])
                ext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir)
                remove_intermediate_files(sra_stat, ext=ext, work_dir=output_dir)
        return None

def remove_sra_files(metadata, amalgkit_out_dir):
    print('Starting SRA file removal.', flush=True)
    for sra_id in metadata.df['run']:
        sra_pattern = os.path.join(os.path.realpath(amalgkit_out_dir), 'getfastq', sra_id, sra_id + '.sra*')
        path_downloaded_sras = glob.glob(sra_pattern)
        if len(path_downloaded_sras) > 0:
            for path_downloaded_sra in path_downloaded_sras:
                print('Deleting SRA file: {}'.format(path_downloaded_sra))
                os.remove(path_downloaded_sra)
        else:
            print('SRA file not found. Pattern searched: {}'.format(sra_pattern))
    print('')

def get_layout(args, metadata):
    if args.layout == 'auto':
        layouts = metadata.df['lib_layout'].unique().tolist()
        if (len(layouts) != 1):
            print('Detected multiple layouts in the metadata:', layouts)
        layout = 'paired' if 'paired' in layouts else 'single'
    else:
        layout = args.layout
    return layout

def remove_old_intermediate_files(sra_id, work_dir):
    old_files = os.listdir(work_dir)
    files = [f for f in old_files if
             (f.startswith(sra_id)) & (not f.endswith('.sra')) & (os.path.isfile(os.path.join(work_dir, f)))]
    for f in files:
        f_path = os.path.join(work_dir, f)
        print('Deleting old intermediate file:', f_path)
        os.remove(f_path)


def remove_intermediate_files(sra_stat, ext, work_dir):
    file_paths = list()
    if sra_stat['layout'] == 'single':
        file_paths.append(os.path.join(work_dir, sra_stat['sra_id'] + ext))
    elif sra_stat['layout'] == 'paired':
        for i in [1, 2]:
            file_paths.append(os.path.join(work_dir, sra_stat['sra_id'] + '_' + str(i) + ext))
    for file_path in file_paths:
        if os.path.exists(file_path):
            print('Deleting intermediate file:', file_path)
            os.remove(file_path)
        else:
            print('Tried to delete but file not found:', file_path)

def download_sra(metadata, sra_stat, args, work_dir, overwrite=False):
    path_downloaded_sra = os.path.join(work_dir, sra_stat['sra_id'] + '.sra')
    individual_sra_tmp_dir = os.path.join(work_dir, sra_stat['sra_id'] + '/')
    if os.path.exists(path_downloaded_sra):
        print('Previously-downloaded sra file was detected at: {}'.format(path_downloaded_sra))
        if (overwrite):
            print('Removing', path_downloaded_sra)
            print('New sra file will be downloaded.')
            os.remove(path_downloaded_sra)
        else:
            return None
    else:
        print('Previously-downloaded sra file was not detected. New sra file will be downloaded.')

    if (args.aws) or (args.ncbi) or (args.gcp):
        sra_sources = dict()
        sra_id = sra_stat['sra_id']
        is_sra = (metadata.df['run']==sra_stat['sra_id'])
        if args.aws:
            aws_link = metadata.df.loc[is_sra,'AWS_Link'].values[0]
            if aws_link=='':
                sys.stderr.write('AWS_Link is empty and will be skipped.\n')
            else:
                sra_sources['AWS'] = aws_link
        if args.gcp:
            gcp_link = metadata.df.loc[is_sra,'GCP_Link'].values[0]
            if gcp_link=='':
                sys.stderr.write('GCP_Link is empty and will be skipped.\n')
            else:
                sra_sources['GCP'] = gcp_link
        if args.ncbi:
            ncbi_link = metadata.df.loc[is_sra,'NCBI_Link'].values[0]
            if ncbi_link=='':
                sys.stderr.write('NCBI_Link is empty and will be skipped.\n')
            else:
                sra_sources['NCBI'] = ncbi_link
        if len(sra_sources)==0:
            print('No source URL is available. Check whether --aws, --gcp, and --ncbi are properly set.')
        is_sra_download_completed = False
        for sra_source_name in sra_sources.keys():
            print("Trying to fetch {} from {}: {}".format(sra_id, sra_source_name, sra_sources[sra_source_name]))
            if str(sra_sources[sra_source_name])=='nan':
                sys.stderr.write("Skipping. No URL for {}.\n".format(sra_source_name))
                continue
            try:
                urllib.request.urlretrieve(str(sra_sources[sra_source_name]), path_downloaded_sra)
                if os.path.exists(path_downloaded_sra):
                    is_sra_download_completed = True
                    print('SRA file was downloaded with urllib.request from {}'.format(sra_source_name), flush=True)
                    break
            except urllib.error.URLError:
                sys.stderr.write("urllib.request failed SRA download from {}.\n".format(sra_source_name))
        if not is_sra_download_completed:
            sys.stderr.write("Exhausted all sources of download.\n")
        else:
            assert os.path.exists(path_downloaded_sra), 'SRA file download failed: ' + sra_stat['sra_id']
            return

    if not os.path.exists(path_downloaded_sra):
        print('Trying to download the SRA file using prefetch.')
        if os.path.exists(individual_sra_tmp_dir):
            shutil.rmtree(individual_sra_tmp_dir)
        prefetch_command = [args.prefetch_exe, '--force', 'no', '--max-size', '100G',
                            '--output-directory', './', str(sra_stat['sra_id'])]
        print('Command:', ' '.join(prefetch_command))
        prefetch_out = subprocess.run(prefetch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print('prefetch stdout:')
        print(prefetch_out.stdout.decode('utf8'))
        print('prefetch stderr:')
        print(prefetch_out.stderr.decode('utf8'))
        if (prefetch_out.returncode):
            sys.stderr.write("prefetch did not finish safely. Trying prefetch again.\n")
            prefetch_command = [args.prefetch_exe, '--force', 'no', '--max-size', '100G',
                                sra_stat['sra_id']]
            prefetch_out = subprocess.run(prefetch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print('prefetch stdout:')
            print(prefetch_out.stdout.decode('utf8'))
            print('prefetch stderr:')
            print(prefetch_out.stderr.decode('utf8'))
            if (prefetch_out.returncode !=0):
                sys.stderr.write("Again, prefetch did not finish safely.\n")
    # Move files downloaded by prefetch. This is necessary because absolute path didn't work for prefetch --output-directory
    if os.path.exists(os.path.join('./', sra_stat['sra_id'] + '/', sra_stat['sra_id'] + '.sra')):
        subprocess.run(['mv', os.path.join('./', sra_stat['sra_id'] + '/', sra_stat['sra_id'] + '.sra'), path_downloaded_sra],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        shutil.rmtree(os.path.join('./', sra_stat['sra_id'] + '/'))
    elif os.path.exists(os.path.expanduser(os.path.join('~/ncbi/public/sra/', sra_stat['sra_id'] + '.sra'))):
        subprocess.run(
            ['mv', os.path.expanduser(os.path.join('~/ncbi/public/sra/', sra_stat['sra_id'] + '.sra')), path_downloaded_sra],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Move files downloaded by ascp
    if os.path.exists(os.path.join(individual_sra_tmp_dir, sra_stat['sra_id'] + '.sra')):
        subprocess.run(['mv', os.path.join(individual_sra_tmp_dir, sra_stat['sra_id'] + '.sra'), path_downloaded_sra],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        shutil.rmtree(individual_sra_tmp_dir)
    err_txt = 'SRA file download failed for {}. Expected PATH: {}'.format(sra_stat['sra_id'], path_downloaded_sra)
    assert os.path.exists(path_downloaded_sra), err_txt

def check_getfastq_dependency(args):
    if args.pfd:
        test_pfd = subprocess.run([args.pfd_exe, '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert (test_pfd.returncode == 0), "parallel-fastq-dump PATH cannot be found: " + args.pfd_exe
        # commented out because prefetch is often not activatable in containers and no longer strictly required for getfastq.
        #test_prefetch = subprocess.run([args.prefetch_exe, '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #assert (test_prefetch.returncode == 0), "prefetch (SRA toolkit) PATH cannot be found: " + args.prefetch_exe
    if args.fastp:
        test_fp = subprocess.run([args.fastp_exe, '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert (test_fp.returncode == 0), "fastp PATH cannot be found: " + args.fastp_exe
    check_seqkit_dependency()
    return None
def run_pfd(sra_stat, args, metadata, start, end):
    path_downloaded_sra = os.path.join(sra_stat['getfastq_sra_dir'], sra_stat['sra_id'] + '.sra')
    pfd_command = ['parallel-fastq-dump', '-t', str(args.threads), '--minReadLen', str(args.min_read_length),
                   '--qual-filter-1',
                   '--skip-technical', '--split-3', '--clip', '--gzip', '--outdir', sra_stat['getfastq_sra_dir'],
                   '--tmpdir', sra_stat['getfastq_sra_dir']]
    print('Total sampled bases:', "{:,}".format(sra_stat['spot_length'] * (end - start + 1)), 'bp')
    pfd_command = pfd_command + ['--minSpotId', str(int(start)), '--maxSpotId', str(int(end))]
    # If sra_stat['sra_id'], not path_downloaded_sra, is provided, pfd couldn't find pre-downloaded .sra files
    # and start downloading it to $HOME/ncbi/public/sra/
    pfd_command = pfd_command + ['-s', path_downloaded_sra]
    print('Command:', ' '.join(pfd_command))
    pfd_out = subprocess.run(pfd_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if args.pfd_print:
        print('parallel-fastq-dump stdout:')
        print(pfd_out.stdout.decode('utf8'))
        print('parallel-fastq-dump stderr:')
        print(pfd_out.stderr.decode('utf8'))
    if (pfd_out.returncode != 0):
        sys.stderr.write("pfd did not finish safely.\n")
        sys.exit(1)
    stdout = pfd_out.stdout.decode('utf8')
    nd = [int(line.replace('Read ', '').split(' ')[0]) for line in stdout.split('\n') if line.startswith('Read')]
    nr = [int(line.replace('Rejected ', '').split(' ')[0]) for line in stdout.split('\n') if line.startswith('Rejected')]
    nw = [int(line.replace('Written ', '').split(' ')[0]) for line in stdout.split('\n') if line.startswith('Written')]
    ind_sra = metadata.df.index[metadata.df.loc[:,'run'] == sra_stat['sra_id']].values[0]
    metadata.df.at[ind_sra,'num_dumped'] += sum(nd)
    metadata.df.at[ind_sra,'num_rejected'] += sum(nr)
    metadata.df.at[ind_sra,'num_written'] += sum(nw)
    metadata.df.at[ind_sra,'bp_dumped'] += sum(nd) * sra_stat['spot_length']
    metadata.df.at[ind_sra,'bp_rejected'] += sum(nr) * sra_stat['spot_length']
    metadata.df.at[ind_sra,'bp_written'] += sum(nw) * sra_stat['spot_length']
    sra_stat = detect_layout_from_file(sra_stat)
    remove_unpaired_files(sra_stat)
    metadata.df.at[ind_sra,'layout_amalgkit'] = sra_stat['layout']
    return metadata,sra_stat

def remove_unpaired_files(sra_stat):
    if (sra_stat['layout']=='paired'):
        # Order is important in this list. More downstream should come first.
        extensions = ['.amalgkit.fastq.gz', '.rename.fastq.gz', '.fastp.fastq.gz', '.fastq.gz']
        for ext in extensions:
            single_fastq_file = os.path.join(sra_stat['getfastq_sra_dir'], sra_stat['sra_id'] + ext)
            if os.path.exists(single_fastq_file):
                print('Removing 3rd fastq file: {}'.format(single_fastq_file), flush=True)
                os.remove(single_fastq_file)

def run_fastp(sra_stat, args, output_dir, metadata):
    inext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir)
    outext = '.fastp.fastq.gz'
    if args.threads > 16:
        print('Too many threads for fastp (--threads {}). Only 16 threads will be used.'.format(args.threads))
        fastp_thread = 16
    else:
        fastp_thread = args.threads
    fp_command = ['fastp', '--thread', str(fastp_thread), '--length_required',
                  str(args.min_read_length)] + args.fastp_option.split(' ')
    if sra_stat['layout'] == 'single':
        infile = os.path.join(output_dir, sra_stat['sra_id'])
        fp_command = fp_command + ['--in1', infile + inext, '--out1', infile + outext]
    elif sra_stat['layout'] == 'paired':
        infile1 = os.path.join(output_dir, sra_stat['sra_id'] + '_1')
        infile2 = os.path.join(output_dir, sra_stat['sra_id'] + '_2')
        fp_command = fp_command + ['--in1', infile1 + inext, '--out1', infile1 + outext, '--in2', infile2 + inext,
                                   '--out2', infile2 + outext]
    fp_command = [fc for fc in fp_command if fc != '']
    print('Command:', ' '.join(fp_command))
    try:
        fp_out = subprocess.run(fp_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        raise RuntimeError("command '{}' returned with error (code {}): {}".format(e.cmd, e.returncode, e.output))
    if args.fastp_print:
        print('fastp stdout:')
        print(fp_out.stdout.decode('utf8'))
        print('fastp stderr:')
        print(fp_out.stderr.decode('utf8'))
    if args.remove_tmp:
        remove_intermediate_files(sra_stat, ext=inext, work_dir=output_dir)
    bps = fp_out.stderr.decode('utf8').split('\n')
    num_in = list()
    num_out = list()
    bp_in = list()
    bp_out = list()
    for i in range(len(bps)):
        if (' before filtering:' in bps[i]):
            num_in.append(int(bps[i + 1].replace('total reads: ', '')))
            bp_in.append(int(bps[i + 2].replace('total bases: ', '')))
        if (' after filtering:' in bps[i]) | (' aftering filtering:' in bps[i]):
            num_out.append(int(bps[i + 1].replace('total reads: ', '')))
            bp_out.append(int(bps[i + 2].replace('total bases: ', '')))
    ind_sra = metadata.df.index[metadata.df.loc[:,'run'] == sra_stat['sra_id']].values[0]
    metadata.df.at[ind_sra,'num_fastp_in'] += sum(num_in)
    metadata.df.at[ind_sra,'num_fastp_out'] += sum(num_out)
    metadata.df.at[ind_sra,'bp_fastp_in'] += sum(bp_in)
    metadata.df.at[ind_sra,'bp_fastp_out'] += sum(bp_out)
    return metadata

def rename_reads(sra_stat, args, output_dir):
    inext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir)
    outext = '.rename.fastq.gz'
    seqkit_command_base = ['seqkit', 'replace', '--threads', str(args.threads), '--pattern', '[[:space:]].*', ]
    if sra_stat['layout'] == 'single':
        inbase = os.path.join(output_dir, sra_stat['sra_id'])
        if os.path.exists(inbase + inext):
            infile = inbase + inext
            outfile = inbase + outext
            seqkit_command = seqkit_command_base + ['--replacement', "'/1'", '--out-file', outfile, infile]
            print('Command:', ' '.join(seqkit_command))
            seqkit_out = subprocess.run(seqkit_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print('SeqKit stdout:')
            print(seqkit_out.stdout.decode('utf8'))
            print('SeqKit stderr:')
            print(seqkit_out.stderr.decode('utf8'))
            if (seqkit_out.returncode != 0):
                raise Exception("SeqKit did not finish safely.")
    elif sra_stat['layout'] == 'paired':
        inbase1 = os.path.join(output_dir, sra_stat['sra_id'] + '_1')
        inbase2 = os.path.join(output_dir, sra_stat['sra_id'] + '_2')
        if os.path.exists(inbase1 + inext):
            infile1 = inbase1 + inext
            infile2 = inbase2 + inext
            outfile1 = inbase1 + outext
            outfile2 = inbase2 + outext
            seqkit_command1 = seqkit_command_base + ['--replacement', '/1', '--out-file', outfile1, infile1]
            seqkit_command2 = seqkit_command_base + ['--replacement', '/2', '--out-file', outfile2, infile2]
            print('Command1:', ' '.join(seqkit_command1))
            print('Command2:', ' '.join(seqkit_command2))
            seqkit_out1 = subprocess.run(seqkit_command1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            seqkit_out2 = subprocess.run(seqkit_command2, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print('SeqKit stdout1:')
            print(seqkit_out1.stdout.decode('utf8'))
            print('SeqKit stderr1:')
            print(seqkit_out1.stderr.decode('utf8'))
            print('SeqKit stdout2:')
            print(seqkit_out2.stdout.decode('utf8'))
            print('SeqKit stderr2:')
            print(seqkit_out2.stderr.decode('utf8'))
            if (seqkit_out1.returncode != 0) or (seqkit_out2.returncode != 0):
                raise Exception("SeqKit did not finish safely.")
    if args.remove_tmp:
        remove_intermediate_files(sra_stat, ext=inext, work_dir=output_dir)

def rename_fastq(sra_stat, output_dir, inext, outext):
    if sra_stat['layout'] == 'single':
        inbase = os.path.join(output_dir, sra_stat['sra_id'])
        os.rename(inbase + inext, inbase + outext)
    elif sra_stat['layout'] == 'paired':
        inbase1 = os.path.join(output_dir, sra_stat['sra_id'] + '_1')
        inbase2 = os.path.join(output_dir, sra_stat['sra_id'] + '_2')
        os.rename(inbase1 + inext, inbase1 + outext)
        os.rename(inbase2 + inext, inbase2 + outext)

def calc_2nd_ranges(metadata):
    sra_target_bp = metadata.df.loc[:,'bp_until_target_size']
    rate_obtained = metadata.df.loc[:,'rate_obtained']
    spot_lengths = metadata.df.loc[:,'spot_length_amalgkit']
    total_spots = metadata.df.loc[:,'total_spots']
    sra_target_reads = numpy.zeros_like(sra_target_bp)
    for i in numpy.arange(sra_target_reads.shape[0]):
        if numpy.isnan(rate_obtained[i]):
            sra_target_reads[i] = (sra_target_bp[i]/spot_lengths[i]).astype(int)+1 # If no read was extracted in 1st.
        else:
            sra_target_reads[i] = ((sra_target_bp[i]/spot_lengths[i])/rate_obtained[i]).astype(int)+1
    start_2nds = metadata.df.loc[:,'spot_end_1st'] + 1
    end_2nds = start_2nds + sra_target_reads
    pooled_missing_bp = metadata.df.loc[:,'bp_until_target_size'].sum()
    for dummy in range(1000):
        current_total_bp = 0
        for ind in end_2nds.index:
            pooled_missing_read = (pooled_missing_bp / spot_lengths.loc[ind]).astype(int)
            if ((end_2nds.loc[ind] + pooled_missing_read) < total_spots.loc[ind]):
                pooled_missing_bp = 0
                end_2nds.loc[ind] = end_2nds.loc[ind] + pooled_missing_bp
            elif (end_2nds.loc[ind] + pooled_missing_read > total_spots.loc[ind]):
                pooled_missing_bp = (end_2nds.loc[ind] + pooled_missing_read - total_spots.loc[ind]) * \
                                    spot_lengths.loc[ind]
                end_2nds.loc[ind] = total_spots.loc[ind]
            current_total_bp += end_2nds.loc[ind] * spot_lengths.loc[ind]
        all_equal_total_spots = all([e2 == ts for e2, ts in zip(end_2nds, total_spots)])
        is_enough_read = (current_total_bp >= metadata.df.loc[:,'bp_until_target_size'].sum())
        if all_equal_total_spots:
            print('Reached total spots in all SRAs.', flush=True)
            break
        if is_enough_read:
            print('Enough read numbers were assigned for the 2nd round sequence extraction.', flush=True)
            break
    metadata.df.loc[:,'spot_start_2nd'] = start_2nds
    metadata.df.loc[:,'spot_end_2nd'] = end_2nds
    return metadata

def print_read_stats(args, metadata, g, sra_stat=None, individual=False):
    if sra_stat is None:
        df = metadata.df
        print('Target size (--max_bp): {:,} bp'.format(g['max_bp']))
    else:
        df = metadata.df.loc[(metadata.df['run']==sra_stat['sra_id']),:]
        print('Individual target size: {:,} bp'.format(g['num_bp_per_sra']))
    if args.pfd:
        print('Sum of fastq_dump dumped reads: {:,} bp'.format(df['bp_dumped'].sum()))
        print('Sum of fastq_dump rejected reads: {:,} bp'.format(df['bp_rejected'].sum()))
        print('Sum of fastq_dump written reads: {:,} bp'.format(df['bp_written'].sum()))
    if args.fastp:
        print('Sum of fastp input reads: {:,} bp'.format(df['bp_fastp_in'].sum()))
        print('Sum of fastp output reads: {:,} bp'.format(df['bp_fastp_out'].sum()))
    if individual:
        print('Individual SRA IDs:', ' '.join(df['run'].values))
        read_types = list()
        keys = list()
        if args.pfd:
            read_types = read_types + ['fastq_dump dumped reads', 'fastq_dump rejected reads',
                                       'fastq_dump written reads']
            keys = keys + ['bp_dumped', 'bp_rejected', 'bp_written']
        if args.fastp:
            read_types = read_types + ['fastp input reads', 'fastp output reads']
            keys = keys + ['bp_fastp_in', 'bp_fastp_out']
        if len(read_types) > 0:
            for rt, key in zip(read_types, keys):
                values = ['{:,}'.format(s) for s in df[key].values]
                txt = ' '.join(values)
                print('Individual {} (bp): {}'.format(rt, txt))
    print('')

def getfastq_metadata(args):
    if args.id is not None:
        print('--id is specified. Downloading SRA metadata from Entrez.')
        Entrez.email = args.entrez_email
        sra_id = args.id
        search_term = getfastq_search_term(sra_id, args.entrez_additional_search_term)
        print('Entrez search term:', search_term)
        xml_root = getfastq_getxml(search_term=search_term)
        metadata = Metadata.from_xml(xml_root=xml_root)
        print('Filtering SRA entry with --layout:', args.layout)
        layout = get_layout(args, metadata)
        metadata.df = metadata.df.loc[(metadata.df['lib_layout'] == layout), :]
        if args.sci_name is not None:
            print('Filtering SRA entry with --sci_name:', args.sci_name)
            metadata.df = metadata.df.loc[(metadata.df['scientific_name'] == args.sci_name), :]
    if args.id_list is not None:
        print('--id_list is specified. Downloading SRA metadata from Entrez.')
        Entrez.email = args.entrez_email
        sra_id_list = [line.rstrip('\n') for line in open(args.id_list) if not line.startswith('#')]
        metadata_dict = dict()
        for sra_id in sra_id_list:
            search_term = getfastq_search_term(sra_id, args.entrez_additional_search_term)
            print('Entrez search term:', search_term)
            xml_root = getfastq_getxml(search_term)
            metadata_dict_tmp = Metadata.from_xml(xml_root)
            if metadata_dict_tmp.df.shape[0]==0:
                print('No associated SRA. Skipping {}'.format(sra_id))
                continue
            metadata_dict[sra_id] = metadata_dict_tmp
            print('Filtering SRA entry with --layout:', args.layout)
            layout = get_layout(args, metadata_dict[sra_id])
            metadata_dict[sra_id].df = metadata_dict[sra_id].df.loc[(metadata_dict[sra_id].df['lib_layout'] == layout), :]
            if args.sci_name is not None:
                print('Filtering SRA entry with --sci_name:', args.sci_name)
                metadata_dict[sra_id].df = metadata_dict[sra_id].df.loc[(metadata_dict[sra_id].df['scientific_name'] == args.sci_name), :]
        if len(metadata_dict)==0:
            print('No associated SRA is found with --id_list. Exiting.')
            sys.exit(1)
        metadata = list(metadata_dict.values())[0]
        metadata.df = pandas.concat([ v.df for v in metadata_dict.values() ], ignore_index=True)
    if (args.id is None)&(args.id_list is None):
        metadata = load_metadata(args)
    metadata.df['total_bases'] = metadata.df.loc[:,'total_bases'].replace('', numpy.nan).astype(float)
    metadata.df['spot_length'] = metadata.df.loc[:, 'spot_length'].replace('', numpy.nan).astype(float)
    return metadata

def is_getfastq_output_present(sra_stat):
    sra_stat = detect_layout_from_file(sra_stat)
    prefixes = [sra_stat['sra_id'], ]
    if sra_stat['layout'] == 'single':
        sub_exts = ['', ]
    elif sra_stat['layout'] == 'paired':
        sub_exts = ['_1', '_2']
    exts = ['.amalgkit.fastq.gz', ]
    is_output_present = True
    for prefix, sub_ext, ext in itertools.product(prefixes, sub_exts, exts):
        out_path1 = os.path.join(sra_stat['getfastq_sra_dir'], prefix + sub_ext + ext)
        out_path2 = os.path.join(sra_stat['getfastq_sra_dir'], prefix + sub_ext + ext + '.safely_removed')
        is_out1 = os.path.exists(out_path1)
        is_out2 = os.path.exists(out_path2)
        if is_out1:
            print('getfastq output detected: {}'.format(out_path1))
        if is_out2:
            print('getfastq output detected: {}'.format(out_path2))
        is_output_present *= (is_out1 | is_out2)
    return is_output_present

def remove_experiment_without_run(metadata):
    num_all_run = metadata.df.shape[0]
    is_missing_run = (metadata.df.loc[:, 'run'] == '')
    num_missing_run = is_missing_run.sum()
    if (num_missing_run > 0):
        print('There are {} out of {} Experiments without Run ID. Removing.'.format(num_missing_run, num_all_run))
        metadata.df = metadata.df.loc[~is_missing_run, :]
    return metadata

def initialize_columns(metadata, g):
    time_keys = ['time_start_1st', 'time_end_1st', 'time_start_2nd', 'time_end_2nd',]
    spot_keys = ['spot_start_1st', 'spot_end_1st', 'spot_start_2nd', 'spot_end_2nd',]
    keys = (['num_dumped', 'num_rejected', 'num_written', 'num_fastp_in', 'num_fastp_out','bp_amalgkit',
            'bp_dumped', 'bp_rejected', 'bp_written', 'bp_fastp_in', 'bp_fastp_out', 'bp_discarded',
            'bp_still_available', 'bp_specified_for_extraction','rate_obtained','layout_amalgkit',]
            + time_keys + spot_keys)
    for key in keys:
        if key=='layout_amalgkit':
            metadata.df.loc[:,key] = ''
        elif key=='rate_obtained':
            metadata.df.loc[:, key] = numpy.nan
        else:
            metadata.df.loc[:,key] = 0
    metadata.df.loc[:, 'bp_until_target_size'] = g['num_bp_per_sra']
    cols = ['total_spots','total_bases','size','nominal_length','nominal_sdev','spot_length']
    for col in cols:
        if any([ dtype in str(metadata.df[col].dtype) for dtype in ['str','object'] ]):
            metadata.df[col] = metadata.df.loc[:,col].astype(str).str.replace('^$', 'nan', regex=True).astype(float)
    for col in time_keys:
        metadata.df[col] = metadata.df[col].astype(float)
    return metadata

def sequence_extraction(args, sra_stat, metadata, g, start, end):
    sra_id = sra_stat['sra_id']
    ind_sra = metadata.df.index[metadata.df.loc[:,'run'] == sra_id].values[0]
    if args.pfd:
        metadata,sra_stat = run_pfd(sra_stat, args, metadata, start, end)
        bp_discarded = metadata.df.at[ind_sra,'bp_dumped'] - metadata.df.at[ind_sra,'bp_written']
        metadata.df.at[ind_sra,'bp_discarded'] += bp_discarded
    no_read_written = (metadata.df.loc[(metadata.df.loc[:,'run']==sra_id),'num_written'].values[0]==0)
    if no_read_written:
        return metadata
    if args.fastp:
        metadata = run_fastp(sra_stat, args, sra_stat['getfastq_sra_dir'], metadata)
        bp_discarded = metadata.df.at[ind_sra,'bp_dumped'] - metadata.df.at[ind_sra,'bp_fastp_out']
        metadata.df.at[ind_sra,'bp_discarded'] += bp_discarded
    if args.read_name == 'trinity':
        rename_reads(sra_stat, args, sra_stat['getfastq_sra_dir'])
    inext = get_newest_intermediate_file_extension(sra_stat, work_dir=sra_stat['getfastq_sra_dir'])
    outext = '.amalgkit.fastq.gz'
    rename_fastq(sra_stat, sra_stat['getfastq_sra_dir'], inext, outext)
    metadata.df.at[ind_sra,'bp_still_available'] = sra_stat['spot_length'] * (sra_stat['total_spot'] - end)
    bp_specified_for_extraction = sra_stat['spot_length'] * (end - start)
    metadata.df.at[ind_sra, 'bp_specified_for_extraction'] += bp_specified_for_extraction
    if args.fastp:
        metadata.df.at[ind_sra, 'bp_amalgkit'] = metadata.df.at[ind_sra,'bp_fastp_out']
    else:
        metadata.df.at[ind_sra, 'bp_amalgkit'] = metadata.df.at[ind_sra,'bp_written']
    metadata.df.at[ind_sra, 'rate_obtained'] = metadata.df.at[ind_sra, 'bp_amalgkit'] / g['num_bp_per_sra']
    metadata.df.at[ind_sra, 'bp_until_target_size'] -= metadata.df.at[ind_sra, 'bp_amalgkit']
    return metadata

def sequence_extraction_1st_round(args, sra_stat, metadata, g):
    offset = 10000  # https://edwards.sdsu.edu/research/fastq-dump/
    ind_sra = metadata.df.index[metadata.df.loc[:,'run'] == sra_stat['sra_id']].values[0]
    metadata.df.at[ind_sra, 'time_start_1st'] = time.time()
    start, end = get_range(sra_stat, offset, g['total_sra_bp'], g['max_bp'])
    metadata.df.at[ind_sra, 'spot_length_amalgkit'] = sra_stat['spot_length']
    metadata.df.at[ind_sra,'spot_start_1st'] = start
    metadata.df.at[ind_sra,'spot_end_1st'] = end
    metadata = sequence_extraction(args, sra_stat, metadata, g, start, end)
    txt = 'Time elapsed for 1st-round sequence extraction: {}, {:,.1f} sec'
    print(txt.format(sra_stat['sra_id'], int(time.time() - g['start_time'])))
    print('\n--- getfastq 1st-round sequence generation report ---')
    print_read_stats(args, metadata, g, sra_stat, individual=False)
    txt = '{:.2f}% of reads were obtained in the 1st-round sequence generation: {:,} bp out of the individual target amount of {:,} bp'
    percent_obtained = metadata.df.at[ind_sra,'rate_obtained']*100
    bp_amalgkit = metadata.df.at[ind_sra, 'bp_amalgkit']
    print(txt.format(percent_obtained, bp_amalgkit, g['num_bp_per_sra']), flush=True)
    metadata.df.at[ind_sra, 'time_end_1st'] = time.time()
    elapsed_time = metadata.df.at[ind_sra, 'time_end_1st'] - metadata.df.at[ind_sra, 'time_start_1st']
    txt = 'Time elapsed for 1st-round sequence extraction: {}, {:,.1f} sec'
    print(txt.format(sra_stat['sra_id'], elapsed_time))
    print('')
    return metadata

def sequence_extraction_2nd_round(args, sra_stat, metadata, g):
    print('Starting the 2nd-round sequence extraction.')
    ind_sra = metadata.df.index[metadata.df.loc[:,'run'] == sra_stat['sra_id']].values[0]
    metadata.df.at[ind_sra,'time_start_2nd'] = time.time()
    ext_main = '.amalgkit.fastq.gz'
    ext_1st_tmp = '.amalgkit_1st.fastq.gz'
    print('')
    sra_id = sra_stat['sra_id']
    layout = sra_stat['layout']
    start = metadata.df.at[ind_sra,'spot_start_2nd']
    end = metadata.df.at[ind_sra,'spot_end_2nd']
    if (start >= end):
        txt = '{}: All spots have been extracted in the 1st trial. Cancelling the 2nd trial. start={:,}, end={:,}'
        print(txt.format(sra_id, start, end))
        return metadata
    no_read_in_1st = (metadata.df.loc[(metadata.df.loc[:,'run']==sra_id),'bp_written'].values[0]==0)
    if no_read_in_1st:
        print('No read was extracted in 1st round. Skipping 2nd round: {}'.format(sra_id))
        return metadata
    else:
        rename_fastq(sra_stat, sra_stat['getfastq_sra_dir'], inext=ext_main, outext=ext_1st_tmp)
    metadata = sequence_extraction(args, sra_stat, metadata, g, start, end)
    if (layout == 'single'):
        subexts = ['']
    elif (layout == 'paired'):
        subexts = ['_1', '_2']
    for subext in subexts:
        added_path = os.path.join(sra_stat['getfastq_sra_dir'], sra_id + subext + ext_1st_tmp)
        adding_path = os.path.join(sra_stat['getfastq_sra_dir'], sra_id + subext + ext_main)
        assert os.path.exists(added_path), 'Dumped fastq not found: ' + added_path
        assert os.path.exists(adding_path), 'Dumped fastq not found: ' + adding_path
        os.system('cat "' + adding_path + '" >> "' + added_path + '"')
        os.remove(adding_path)
        os.rename(added_path, adding_path)
    metadata.df.at[ind_sra, 'time_end_2nd'] = time.time()
    elapsed_time = metadata.df.at[ind_sra, 'time_end_2nd'] - metadata.df.at[ind_sra, 'time_start_2nd']
    txt = 'Time elapsed for 2nd-round sequence extraction: {}, {:,} sec'
    print(txt.format(sra_stat['sra_id'], elapsed_time))
    print('')
    return metadata

def sequence_extraction_private(i, metadata, sra_stat, args):
    for col in ['read1_path','read2_path']:
        path_from = metadata.df.at[i,col]
        path_to = os.path.join(sra_stat['getfastq_sra_dir'], os.path.basename(path_from))
        path_to = path_to.replace('.fq', '.fastq')
        if not path_to.endswith('.gz'):
            path_to = path_to+'.gz' # .gz is necessary even if the original file is not compressed.
        if os.path.exists(path_from):
            if os.path.lexists(path_to):
                os.remove(path_to)
            os.symlink(src=path_from, dst=path_to)
        else:
            sys.stderr.write('Private fastq file not found: {}\n'.format(path_from))
    if args.fastp:
        metadata = run_fastp(sra_stat, args, sra_stat['getfastq_sra_dir'], metadata)
    inext = get_newest_intermediate_file_extension(sra_stat, work_dir=sra_stat['getfastq_sra_dir'])
    if (inext=='no_extension_found')&(sra_stat['layout']=='paired'):
        raise Exception('Paired-end file names may be invalid. They should contain _1 and _2 to indicate a pair: {}'.format(sra_stat['sra_id']))
    outext = '.amalgkit.fastq.gz'
    rename_fastq(sra_stat, sra_stat['getfastq_sra_dir'], inext, outext)
    return metadata

def check_metadata_validity(metadata):
    assert metadata.df.shape[0] > 0, 'No SRA entry found. Make sure whether --id or --id_list is compatible with --sci_name and --layout.'
    is_total_bases_na = metadata.df.loc[:,'total_bases'].isnull()
    is_total_bases_na |= (metadata.df.loc[:, 'total_bases']==0)
    is_total_bases_na |= (metadata.df.loc[:, 'total_bases']=='')
    if is_total_bases_na.any():
        txt = 'Empty value(s) of total_bases were detected in {}. Filling a placeholder value 999,999,999,999\n'
        sys.stderr.write(txt.format(', '.join(metadata.df.loc[is_total_bases_na, 'run'])))
        metadata.df.loc[is_total_bases_na,'total_bases'] = 999999999999
        metadata.df['total_bases'] = metadata.df.loc[:, 'total_bases'].astype(int)
    is_total_spots_na = metadata.df.loc[:, 'total_spots'].isnull()
    is_total_spots_na |=  (metadata.df.loc[:, 'total_spots']==0)
    is_total_spots_na |=  (metadata.df.loc[:, 'total_spots']=='')
    if is_total_spots_na.any():
        new_values = metadata.df.loc[is_total_spots_na,'total_bases'] / metadata.df.loc[is_total_spots_na,'spot_length']
        if is_total_spots_na.any():
            txt = 'Empty value(s) of total_spots were detected in {}. Filling a placeholder value 999,999,999,999\n'
            sys.stderr.write(txt.format(', '.join(metadata.df.loc[is_total_spots_na,'run'])))
            new_values.loc[new_values.isnull()] = 999999999999 # https://github.com/kfuku52/amalgkit/issues/110
        new_values = new_values.astype(int)
        metadata.df.loc[is_total_spots_na, 'total_spots'] = new_values
    for i in metadata.df.index:
        txt = 'Individual SRA size of {}: {:,} bp'
        print(txt.format(metadata.df.at[i, 'run'], metadata.df.at[i, 'total_bases']))
    return metadata

def initialize_global_params(args, metadata):
    g = dict()
    g['start_time'] = time.time()
    g['max_bp'] = int(args.max_bp.replace(',', ''))
    g['num_sra'] = metadata.df.shape[0]
    g['num_bp_per_sra'] = int(g['max_bp'] / g['num_sra'])
    g['total_sra_bp'] = metadata.df.loc[:,'total_bases'].sum()
    print('Number of SRAs to be processed: {:,}'.format(g['num_sra']))
    print('Total target size (--max_bp): {:,} bp'.format(g['max_bp']))
    print('The sum of SRA sizes: {:,} bp'.format(g['total_sra_bp']))
    print('Target size per SRA: {:,} bp'.format(g['num_bp_per_sra']))
    return g

def getfastq_main(args):
    check_getfastq_dependency(args)
    metadata = getfastq_metadata(args)
    metadata = remove_experiment_without_run(metadata)
    metadata = check_metadata_validity(metadata)
    g = initialize_global_params(args, metadata)
    metadata = initialize_columns(metadata, g)
    flag_private_file = False
    flag_any_output_file_present = False
    # 1st round sequence extraction
    for i in metadata.df.index:
        print('')
        sra_id = metadata.df.at[i, 'run']
        print('Processing SRA ID: {}'.format(sra_id))
        sra_stat = get_sra_stat(sra_id, metadata, g['num_bp_per_sra'])
        sra_stat['getfastq_sra_dir'] = get_getfastq_run_dir(args, sra_id)
        if (is_getfastq_output_present(sra_stat)):
            flag_any_output_file_present =True
            if not args.redo:
                txt = 'Output file(s) detected. Skipping {}. Set "--redo yes" for reanalysis.'
                print(txt.format(sra_id), flush=True)
                continue
        remove_old_intermediate_files(sra_id=sra_id, work_dir=sra_stat['getfastq_sra_dir'])
        print('Library layout:', sra_stat['layout'])
        print('Number of reads:', "{:,}".format(sra_stat['total_spot']))
        print('Single/Paired read length:', sra_stat['spot_length'], 'bp')
        print('Total bases:', "{:,}".format(int(metadata.df.loc[i, 'total_bases'])), 'bp')
        flag_private_file = False
        if 'private_file' in metadata.df.columns:
            if metadata.df.at[i,'private_file']=='yes':
                print('Processing {} as private data. --max_bp is disabled.'.format(sra_id), flush=True)
                flag_private_file = True
                sequence_extraction_private(i, metadata, sra_stat, args)
        if not flag_private_file:
            print('Processing {} as publicly available data from SRA.'.format(sra_id), flush=True)
            download_sra(metadata, sra_stat, args, sra_stat['getfastq_sra_dir'], overwrite=False)
            metadata = sequence_extraction_1st_round(args, sra_stat, metadata, g)
    # 2nd round sequence extraction
    if (not flag_private_file) & (not flag_any_output_file_present):
        g['rate_obtained_1st'] = metadata.df.loc[:,'bp_amalgkit'].sum() / g['max_bp']
        if (g['rate_obtained_1st'] < (args.tol*0.01)):
            print('Sufficient data were obtained in the 1st-round sequence extraction. Proceeding without the 2nd round.')
        else:
            txt = 'Only {:,.2f}% ({:,}/{:,}) of the target size (--max_bp) was obtained in the 1st round. Proceeding to the 2nd round read extraction.'
            print(txt.format(g['rate_obtained_1st']*100, metadata.df.loc[:,'bp_amalgkit'].sum(), g['max_bp']), flush=True)
            metadata = calc_2nd_ranges(metadata)
            for i in metadata.df.index:
                sra_id = metadata.df.at[i, 'run']
                sra_stat = get_sra_stat(sra_id, metadata, g['num_bp_per_sra'])
                sra_stat['getfastq_sra_dir'] = get_getfastq_run_dir(args, sra_id)
                metadata = sequence_extraction_2nd_round(args, sra_stat, metadata, g)
        g['rate_obtained_2nd'] = metadata.df.loc[:, 'bp_amalgkit'].sum() / g['max_bp']
        txt = '2nd round read extraction improved % bp from {:,.2f}% to {:,.2f}%'
        print(txt.format(g['rate_obtained_1st']*100, g['rate_obtained_2nd']*100), flush=True)
    # Postprocessing
    if args.remove_sra:
        remove_sra_files(metadata=metadata, amalgkit_out_dir=args.out_dir)
    else:
        print('SRA files not removed: {}'.format(sra_stat['getfastq_sra_dir']))
    if (not flag_any_output_file_present):
        print('')
        print('\n--- getfastq final report ---')
        print_read_stats(args, metadata, g, sra_stat=None, individual=True)


================================================
FILE: amalgkit/integrate.py
================================================
import pandas as pd

import glob
import os
import platform
import re
import subprocess
import warnings

from amalgkit.sanity import check_getfastq_outputs
from amalgkit.util import *

def get_fastq_stats(args):
    print("Starting integration of fastq-file metadata...")
    if not os.path.exists(args.fastq_dir):
        raise ValueError("PATH to fastq directory does not exist: {}".format(args.fastq_dir))
    all_files = os.listdir(args.fastq_dir)
    fastq_extension_regex = r'(\.fq$|\.fastq$|\.fq.gz$|\.fastq.gz$)'
    all_fastq_files = [ f for f in all_files if re.search(fastq_extension_regex, f) ]
    if len(all_fastq_files)==0:
        txt = 'No detected fastq files (with regex "{}") in: {}'.format(fastq_extension_regex, args.fastq_dir)
        raise ValueError(txt)
    id_list = list(map(os.path.basename, all_fastq_files))
    id_list = [basename.split(os.extsep)[0] for basename in id_list] # This is necessary to correctly parse single-end fastq files
    # split off last (!) underscore in case sample name includes multiple underscores.
    # Assuming naming format: sample1_1.fastq, sample1_2.fastq
    # sample_1_1.fastq, sample_1_2.fastq is also possible
    id_list = [basename.rsplit('_', 1)[0] for basename in id_list]
    # duplicates (i.e. doublets) should indicate paired-end library
    num_fastq_files = {id: id_list.count(id) for id in id_list}
    column_names = ['scientific_name', 'sample_group', 'run', 'read1_path','read2_path', 'is_sampled',
                    'is_qualified','exclusion', 'lib_layout', 'spot_length', 'total_spots', 'total_bases', 'size', 'private_file']
    tmp_metadata = pd.DataFrame(columns = column_names)
    row = 0
    for id in num_fastq_files:
        if num_fastq_files[id] == 1:
            lib_layout = 'single'
        if num_fastq_files[id] == 2:
            lib_layout = 'paired'
        if num_fastq_files[id] > 2:
            raise ValueError("found more than 2 files (", num_fastq_files[id], ") for id", id,". Please check filenames.")
        fastq_files = [ os.path.join(args.fastq_dir, f) for f in all_fastq_files if re.search(id+'[\._]', f) ]
        print("Found {} file(s) for ID {}. Lib-layout: {}".format(num_fastq_files[id], id,lib_layout), flush=True)
        print("Getting sequence statistics.", flush=True)
        tmp_file = os.path.join(args.out_dir, id+'_seqkit_stats.tmp')
        # check for file extension. seqkit is significantly slower on compressed files, but still fast on decompressed files.
        if fastq_files[0].endswith(('.fq', '.fastq')):
            is_decompressed = True
        elif fastq_files[0].endswith(('.fq.gz', '.fastq.gz')):
            is_decompressed = False
        else:
            warnings.warn("{} is not a fastq file. Skipping.".format(fastq_files[0]))
            continue

        if args.accurate_size or is_decompressed:
            print('--accurate_size set to yes. Running accurate sequence scan.')
            seqkit_command = ['seqkit', 'stats', '-T', '-j', str(args.threads), fastq_files[0]]
            seqkit_stdout = open(tmp_file, 'w')
            subprocess.run(seqkit_command, stdout=seqkit_stdout)
            seqkit_stdout.close()
            tmp_stat_df = pandas.read_csv(tmp_file, sep='\t', header=0)
            total_spots = tmp_stat_df.at[0, 'num_seqs']
        else:
            OS = platform.system()
            if OS == 'Darwin':
                zcat_command = 'zcat < '
            elif OS == 'Linux':
                zcat_command = 'zcat '
            else:
                zcat_command = 'zcat '
                sys.stderr.write('zcat may not be supported by this OS: {}\n'.format(OS))
            seqkit_command = zcat_command + fastq_files[0] + ' | head -n 4000 | seqkit stats -T -j ' + str(
                args.threads)
            total_lines_command = 'echo $(' + zcat_command + str(fastq_files[0]) + ' | wc -l)'
            seqkit_stdout = open(tmp_file, 'w')
            subprocess.run(seqkit_command, shell=True, stdout=seqkit_stdout)
            seqkit_stdout.close()
            tmp_stat_df = pandas.read_csv(tmp_file, sep='\t', header=0)
            total_lines_bytes = subprocess.check_output(total_lines_command, shell=True)
            total_lines = int(total_lines_bytes.decode().replace('\n', ''))
            total_spots = int(total_lines / 4)
        if args.remove_tmp:
            os.remove(tmp_file)
        tmp_stat_df['id'] = pandas.Series(dtype='str')
        tmp_stat_df['file1'] = pandas.Series(dtype='str')
        tmp_stat_df['file2'] = pandas.Series(dtype='str')
        tmp_stat_df.at[0,'id'] = id
        if len(fastq_files) == 2:
            tmp_stat_df.at[0, 'file2'] = fastq_files[1]
        elif len(fastq_files) == 1:
            tmp_stat_df.at[0, 'file2'] = 'unavailable'
        else:
            raise ValueError('Too many files found for set {}'.format(fastq_files))
        tmp_metadata.at[row, 'scientific_name'] = 'Please add in format: Genus species'
        tmp_metadata.at[row,'sample_group'] = 'Please add'
        tmp_metadata.at[row,'run'] = tmp_stat_df.at[0,'id']
        tmp_metadata.at[row,'read1_path'] = os.path.abspath(fastq_files[0])
        if tmp_stat_df.at[0, 'file2'] != 'unavailable':
            tmp_metadata.at[row,'read2_path'] = os.path.abspath(tmp_stat_df.at[0,'file2'])
        else:
            tmp_metadata.at[row, 'read2_path'] = 'unavailable'
        tmp_metadata.at[row,'is_sampled'] = 'yes'
        tmp_metadata.at[row,'is_qualified'] = 'yes'
        tmp_metadata.at[row, 'exclusion'] = 'no'
        tmp_metadata.at[row,'lib_layout'] = lib_layout
        tmp_metadata.at[row,'total_spots'] = total_spots
        tmp_metadata.at[row,'size'] = os.path.getsize(fastq_files[0])
        tmp_metadata.at[row, 'private_file'] = 'yes'
        tmp_metadata.at[row,'spot_length'] = int(tmp_stat_df.at[0,'avg_len'])
        total_bases = total_spots * int(tmp_stat_df.at[0,'avg_len'])
        if lib_layout == 'paired':
            total_bases = total_bases * 2
        tmp_metadata.at[row,'total_bases'] = total_bases
        row += 1
    if not os.path.exists(os.path.join(args.out_dir, 'metadata')):
        os.makedirs(os.path.join(args.out_dir, 'metadata'))
    tmp_metadata = tmp_metadata.sort_values(by='run', axis=0, ascending=True).reset_index(drop=True)
    tmp_metadata.to_csv(os.path.join(args.out_dir, 'metadata_private_fastq.tsv'), sep='\t', index=False)
    return tmp_metadata

def integrate_main(args):
    check_seqkit_dependency()
    if args.metadata=='inferred':
        relative_path = os.path.join(args.out_dir, 'metadata', 'metadata.tsv')
        metadata_path = os.path.realpath(relative_path)
    else:
        metadata_path = os.path.realpath(args.metadata)
    if os.path.exists(metadata_path):
        print('Merging existing metadata and private fastq info: {}'.format(metadata_path))
        metadata = load_metadata(args)
        metadata.df.loc[:,'private_file'] = 'no'
        print("scanning for getfastq output")
        sra_ids = metadata.df.loc[:,'run']
        data_available, data_unavailable = check_getfastq_outputs(args, sra_ids, metadata, args.out_dir)
        metadata.df.loc[metadata.df['run'].isin(data_available), 'data_available'] = 'yes'
        metadata.df.loc[metadata.df['run'].isin(data_unavailable), 'data_available'] = 'no'
        tmp_metadata = get_fastq_stats(args)
        df = pandas.concat([metadata.df, tmp_metadata])
        df.to_csv(os.path.join(args.out_dir, 'metadata', 'metadata_updated_for_private_fastq.tsv'), sep='\t', index=False)
    else:
        print('Generating a new metadata table.')
        get_fastq_stats(args)




================================================
FILE: amalgkit/merge.py
================================================
import pandas

import os
import warnings
from amalgkit.util import *

def merge_main(args):
    quant_dir = os.path.realpath(os.path.join(args.out_dir, 'quant'))
    merge_dir = os.path.realpath(os.path.join(args.out_dir, 'merge'))
    if not os.path.exists(merge_dir):
        os.makedirs(os.path.join(merge_dir))
    metadata = load_metadata(args)
    spp = metadata.df.loc[:,'scientific_name'].dropna().unique()
    for sp in spp:
        print('processing: {}'.format(sp), flush=True)
        sp_filled = sp.replace(' ', '_')
        merge_species_dir = os.path.join(os.path.join(merge_dir, sp_filled))
        is_sp = (metadata.df.loc[:,'scientific_name']==sp)
        sra_ids = metadata.df.loc[is_sp,'run'].values
        is_sampled = (metadata.df.loc[:,'exclusion']=='no')
        sampled_sra_ids = metadata.df.loc[is_sampled,'run'].values
        if len(sra_ids)==0:
            warnings.warn('No SRA Run ID found. Skipping: {}'.format(sp))
            continue
        quant_out_paths = list()
        detected_sra_ids = list()
        for sra_id in sra_ids:
            quant_out_path = os.path.join(quant_dir, sra_id, sra_id+'_abundance.tsv')
            if os.path.exists(quant_out_path):
                if not os.path.exists(merge_species_dir):
                    os.makedirs(merge_species_dir)
                quant_out_paths.append(quant_out_path)
                detected_sra_ids.append(sra_id)
            else:
                if sra_id in sampled_sra_ids:
                    print('quant outfile not found: {}'.format(quant_out_path))
        print('{:,} quant outfiles were detected.'.format(len(quant_out_paths)))
        if len(quant_out_paths) == 0:
            continue
        for col in ['eff_length','est_counts','tpm']:
            out = pandas.read_csv(quant_out_paths[0], header=0, sep='\t').loc[:,['target_id',]]
            values = [ pandas.read_csv(p, header=0, sep='\t').loc[:,[col,]] for p in quant_out_paths ]
            for sra_id,value in zip(detected_sra_ids,values):
                value.columns = [sra_id,]
            out = pandas.concat([out,]+values, axis=1, ignore_index=False, sort=False)
            outfile_name = sp_filled+'_'+col+'.tsv'
            outfile = os.path.join(merge_species_dir, outfile_name)
            print('Writing output file:', outfile)
            out.to_csv(outfile, sep='\t', index=False)
    print('Getting mapping rate from quant output and write new metadata file into merge directory.', flush=True)
    path_metadata_merge = os.path.realpath(os.path.join(args.out_dir, 'merge', 'metadata.tsv'))
    write_updated_metadata(metadata, path_metadata_merge, args)
    r_merge_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'merge.r')
    r_util_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'util.r')
    r_command = ['Rscript', r_merge_path, merge_dir, path_metadata_merge, r_util_path]
    print('Starting R script for plot generation: {}'.format(' '.join(r_command)), flush=True)
    subprocess.check_call(r_command)



================================================
FILE: amalgkit/merge.r
================================================
#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(ggplot2, quietly = TRUE)))
mode = ifelse(length(commandArgs(trailingOnly = TRUE)) == 1, 'debug', 'batch')
if (mode == "debug") {
    dir_work = '/Users/kf/Dropbox/data/evolutionary_transcriptomics/20230527_amalgkit/amalgkit_out'
    setwd(dir_work)
    dir_merge = file.path(dir_work, 'merge')
    file_metadata = file.path(dir_merge, 'metadata.tsv')
    r_util_path = '/Users/kf/Dropbox/repos/amalgkit/amalgkit/util.r'
} else if (mode == "batch") {
    args = commandArgs(trailingOnly = TRUE)
    dir_merge = args[1]
    file_metadata = args[2]
    r_util_path = args[3]
}
source(r_util_path)
font_size = 8

df = read.table(file_metadata, sep = '\t', header = TRUE, quote = '', comment.char = '', check.names = FALSE)
is_excluded = (!df[['exclusion']] == 'no')
is_mapping_rate_available = (!is.na(df[['mapping_rate']]))
cat(sprintf('Number of non-excluded SRA samples: %s\n', formatC(sum(!is_excluded), format = 'd', big.mark = ',')))
cat(sprintf('Number of excluded SRA samples: %s\n', formatC(sum(is_excluded), format = 'd', big.mark = ',')))
cat(sprintf('Number of SRA samples with available mapping rates: %s\n', formatC(sum(is_mapping_rate_available), format = 'd', big.mark = ',')))

df2 = df[((!is_excluded) & (is_mapping_rate_available)),]
cat(sprintf('Number of SRA samples for mapping_rate potting: %s\n', formatC(nrow(df2), format = 'd', big.mark = ',')))
g = ggplot(data = df2)
g = g + geom_boxplot(aes(x = scientific_name, y = mapping_rate), outlier.size = 0.3)
g = g + ylim(0, 100)
g = g + theme_bw(base_size = font_size)
g = g + labs(x = '', y = 'Mapping rate')
g = g + theme(
    axis.text = element_text(size = font_size, color = 'black'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_text(size = font_size, color = 'black'),
    #panel.grid.major.y=element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = font_size, color = 'black'),
    rect = element_rect(fill = "transparent"),
    plot.margin = unit(rep(0.1, 4), "cm")
)
num_spp = length(unique(df[['scientific_name']]))
plot_width = max(3.6, 0.11 * num_spp)
out_path = file.path(dir_merge, 'merge_mapping_rate.pdf')
ggsave(out_path, plot = g, width = plot_width, height = 3.6, units = 'in')

df2 = df[((!is_excluded)),]
cat(sprintf('Number of SRA samples for total_spots potting: %s\n', formatC(nrow(df2), format = 'd', big.mark = ',')))
g = ggplot(data = df2)
g = g + geom_boxplot(aes(x = scientific_name, y = total_spots), outlier.size = 0.3)
g = g + scale_y_log10()
g = g + theme_bw(base_size = font_size)
g = g + labs(x = '', y = 'Total spots')
g = g + theme(
    axis.text = element_text(size = font_size, color = 'black'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_text(size = font_size, color = 'black'),
    #panel.grid.major.y=element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = font_size, color = 'black'),
    rect = element_rect(fill = "transparent"),
    plot.margin = unit(rep(0.1, 4), "cm")
)
num_spp = length(unique(df[['scientific_name']]))
plot_width = max(3.6, 0.11 * num_spp)
out_path = file.path(dir_merge, 'merge_total_spots.pdf')
ggsave(out_path, plot = g, width = plot_width, height = 3.6, units = 'in')

df2 = df[((!is_excluded)),]
cat(sprintf('Number of SRA samples for total_bases potting: %s\n', formatC(nrow(df2), format = 'd', big.mark = ',')))
g = ggplot(data = df2)
g = g + geom_boxplot(aes(x = scientific_name, y = total_bases), outlier.size = 0.3)
g = g + scale_y_log10()
g = g + theme_bw(base_size = font_size)
g = g + labs(x = '', y = 'Total bases')
g = g + theme(
    axis.text = element_text(size = font_size, color = 'black'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_text(size = font_size, color = 'black'),
    #panel.grid.major.y=element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = font_size, color = 'black'),
    rect = element_rect(fill = "transparent"),
    plot.margin = unit(rep(0.1, 4), "cm")
)
num_spp = length(unique(df[['scientific_name']]))
plot_width = max(3.6, 0.11 * num_spp)
out_path = file.path(dir_merge, 'merge_total_bases.pdf')
ggsave(out_path, plot = g, width = plot_width, height = 3.6, units = 'in')

df2 = df[((!is_excluded)),]
cat(sprintf('Number of SRA samples for lib_layout potting: %s\n', formatC(nrow(df2), format = 'd', big.mark = ',')))
data_summary = aggregate(cbind(count = lib_layout) ~ scientific_name + lib_layout, df2, length)
data_summary[['total']] = ave(data_summary[['count']], data_summary[['scientific_name']], FUN = sum)
data_summary[['proportion']] = data_summary[['count']] / data_summary[['total']]
g = ggplot(data_summary, aes(x = scientific_name, y = count, fill = lib_layout))
g = g + geom_bar(stat = "identity")
g = g + labs(x = "", y = "Count", fill = "Library layout")
g = g + theme_bw(base_size = font_size)
g = g + theme(
    axis.text = element_text(size = font_size, color = 'black'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_text(size = font_size, color = 'black'),
    #panel.grid.major.y=element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = font_size, color = 'black'),
    legend.position = 'bottom',
    rect = element_rect(fill = "transparent"),
    plot.margin = unit(rep(0.1, 4), "cm")
)
num_spp = length(unique(df[['scientific_name']]))
plot_width = max(3.6, 0.11 * num_spp)
out_path = file.path(dir_merge, 'merge_library_layout.pdf')
ggsave(out_path, plot = g, width = plot_width, height = 3.6, units = 'in')

cat(sprintf('Number of SRA samples for exclusion potting: %s\n', formatC(nrow(df), format = 'd', big.mark = ',')))
out_path = file.path(dir_merge, 'merge_exclusion.pdf')
save_exclusion_plot(df = df, out_path = out_path, font_size = font_size)

cat('merge.r completed!\n')



================================================
FILE: amalgkit/metadata.py
================================================
from Bio import Entrez
import lxml.etree
import numpy
import pandas

import datetime
import os
import re
import sys
import time
import warnings

from amalgkit.util import *

from urllib.error import HTTPError

def fetch_sra_xml(search_term, retmax=1000):
    try:
        sra_handle = Entrez.esearch(db="sra", term=search_term, retmax=10000000)
    except HTTPError as e:
        print(e, '- Trying Entrez.esearch() again...')
        sra_handle = Entrez.esearch(db="sra", term=search_term, retmax=10000000)
    sra_record = Entrez.read(sra_handle)
    record_ids = sra_record["IdList"]
    num_record = len(record_ids)
    print('Number of SRA records: {:,}'.format(num_record))
    start_time = time.time()
    print('{}: SRA XML retrieval started.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    root = None
    max_retry = 10
    for i in numpy.arange(numpy.ceil(num_record//retmax)+1):
        start = int(i*retmax)
        end = int(((i+1)*retmax)-1) if num_record >= int(((i+1)*retmax)-1) else num_record
        now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print('{}: Retrieving SRA XML: {:,}-{:,} of {:,} records'.format(now, start, end, num_record), flush=True)
        for j in range(max_retry):
            try:
                handle = Entrez.efetch(db="sra", id=record_ids[start:end], rettype="full", retmode="xml", retmax=retmax)
            except HTTPError as e:
                sleep_second = 60
                print('{} - Trying Entrez.efetch() again after {:,} seconds...'.format(e, sleep_second), flush=True)
                time.sleep(sleep_second)
                continue
            try:
                chunk = lxml.etree.parse(handle).getroot()
            except:
                print('XML may be truncated. Retrying...', flush=True)
                continue
            break
        if root is None:
            root = chunk
        else:
            root.append(chunk)
    elapsed_time = int(time.time() - start_time)
    print('{}: SRA XML retrieval ended.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    print('SRA XML retrieval time: {:,.1f} sec'.format(elapsed_time), flush=True)
    xml_string = lxml.etree.tostring(root, encoding='UTF-8', pretty_print=True)
    for line in str(xml_string).split('\n'):
        if '<Error>' in line:
            print(line)
            raise Exception(',Error. found in the xml.')
    return root

def metadata_main(args):
    Entrez.email = args.entrez_email
    metadata_dir = os.path.join(args.out_dir, 'metadata')
    metadata_outfile_path = os.path.join(metadata_dir, "metadata.tsv")
    if os.path.exists(metadata_outfile_path):
        if args.redo:
            os.remove(metadata_outfile_path)
        else:
            print('Exiting. --redo is specified and the output file already exists at: {}'.format(metadata_outfile_path))
            sys.exit(0)
    for path_dir in [args.out_dir, metadata_dir]:
        if not os.path.exists(path_dir):
            print('Creating directory: {}'.format(path_dir))
            os.mkdir(path_dir)
    search_term = args.search_string
    print('Entrez search term:', search_term)
    root = fetch_sra_xml(search_term=search_term)
    metadata = Metadata.from_xml(xml_root=root)
    metadata.df['tissue'] = metadata.df['tissue'].astype(str)
    metadata.df.loc[(metadata.df['tissue']=='nan'), 'tissue'] = ''
    metadata.df.loc[:, 'sample_group'] = metadata.df.loc[:, 'tissue'].str.lower()
    metadata.reorder(omit_misc=False)
    metadata.df.to_csv(metadata_outfile_path, sep="\t", index=False)
    if metadata.df.shape[0]==0:
        txt = 'No entry was found/survived in the metadata processing. Please reconsider the --search_term specification.\n'
        sys.stderr.write(txt)



================================================
FILE: amalgkit/quant.py
================================================
import os
import re
import subprocess
import sys

from amalgkit.util import *

def quant_output_exists(sra_id, output_dir):
    out_path = os.path.join(output_dir, sra_id + '_abundance.tsv')
    is_output_present = os.path.exists(out_path)
    if is_output_present:
        print('Output file detected: {}'.format(out_path))
        return True
    else:
        print('Output file was not detected: {}'.format(out_path))
        return False

def call_kallisto(args, in_files, metadata, sra_stat, output_dir, index):
    sra_id = sra_stat['sra_id']
    lib_layout = sra_stat['layout']
    kallisto_cmd = ''
    if lib_layout == 'single':
        print("Single end reads detected. Proceeding in single mode")
        if len(in_files) != 1:
            txt = "Library layout: {} and expected 1 input file. " \
                  "Received {} input file[s]. Please check your inputs and metadata."
            raise ValueError(txt.format(lib_layout, len(in_files)))
        nominal_length = metadata.df.loc[:, 'nominal_length'].values[0]
        if nominal_length:
            print('Nominal length in metadata is unusually small ({}). Setting it to 200.'.format(nominal_length))
            if nominal_length < 200 or numpy.isnan(nominal_length):
                nominal_length = 200
        else:
            print("Could not find nominal length in metadata. Assuming fragment length.")
            nominal_length = 200
        print("Fragment length set to: {}".format(nominal_length))
        fragment_sd = nominal_length / 10
        print("Fragment length standard deviation set to: {}".format(fragment_sd))
        kallisto_cmd = ["kallisto", "quant", "--threads", str(args.threads), "--index", index, "-o", output_dir,
                        "--single", "-l", str(nominal_length), "-s", str(fragment_sd), in_files[0]]
    elif lib_layout == 'paired':
        if len(in_files) != 2:
            txt = "Library layout: {} and expected 2 input files. " \
                  "Received {} input file[s]. Please check your inputs and metadata."
            raise ValueError(txt.format(lib_layout, len(in_files)))
        print("Paired-end reads detected. Running in paired read mode.")
        kallisto_cmd = ["kallisto", "quant", "--threads", str(args.threads), "-i", index, "-o",
                        output_dir, in_files[0], in_files[1]]

    print('Command: {}'.format(' '.join(kallisto_cmd)))
    kallisto_out = subprocess.run(kallisto_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print('kallisto quant stdout:')
    print(kallisto_out.stdout.decode('utf8'))
    print('kallisto quant stderr:')
    print(kallisto_out.stderr.decode('utf8'))
    if kallisto_out.returncode != 0:
        sys.stderr.write("kallisto did not finish safely.\n")
        if 'Zero reads pseudoaligned' in kallisto_out.stderr.decode('utf8'):
            sys.stderr.write('No reads are mapped to the reference. This sample will be removed by `amalgkit curate`.')

    # move output to results with unique name
    try:
        os.rename(os.path.join(output_dir, "run_info.json"), os.path.join(output_dir, sra_id + "_run_info.json"))
        os.rename(os.path.join(output_dir, "abundance.tsv"), os.path.join(output_dir, sra_id + "_abundance.tsv"))
        os.rename(os.path.join(output_dir, "abundance.h5"), os.path.join(output_dir, sra_id + "_abundance.h5"))
    except FileNotFoundError:
        pass

    return kallisto_out


def check_kallisto_dependency():
    try:
        subprocess.run(['kallisto', '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        raise Exception("kallisto is not installed.")


def check_layout_mismatch(sra_stat, output_dir):
    if sra_stat['layout'] == 'paired':
        fastq_files = glob.glob(os.path.join(output_dir, sra_stat['sra_id'] + '*.fastq*'))
        if len(fastq_files) == 1:
            sys.stderr.write('Single-end fastq was detected even though layout = {}\n'.format(sra_stat['layout']))
            sys.stderr.write('This sample will be treated as single-end sequencing.\n')
            sra_stat['layout'] = 'single'
    return sra_stat


def run_quant(args, metadata, sra_id, index):
    # make results directory, if not already there
    output_dir = os.path.join(args.out_dir, 'quant', sra_id)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    is_quant_output_available = quant_output_exists(sra_id, output_dir)
    if is_quant_output_available:
        if args.redo:
            print('The output will be overwritten. Set "--redo no" to not overwrite results.')
        else:
            print('Continued. The output will not be overwritten. If you want to overwrite the results, set "--redo yes".')
            return
    output_dir_getfastq = os.path.join(args.out_dir, 'getfastq', sra_id)
    sra_stat = get_sra_stat(sra_id, metadata, num_bp_per_sra=None)
    sra_stat = check_layout_mismatch(sra_stat, output_dir_getfastq)
    sra_stat['getfastq_sra_dir'] = get_getfastq_run_dir(args, sra_id)
    ext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir_getfastq)
    if ext == '.safely_removed':
        print('These files have been safe-deleted. If you wish to re-obtain the .fastq file(s), run: getfastq --id ', sra_id, ' -w ', args.out_dir)
        print('Skipping.')
        return
    if ext == 'no_extension_found':
        sys.stderr.write('getfastq output not found in: {}, layout = {}\n'.format(sra_stat['getfastq_sra_dir'], sra_stat['layout']))
        txt = 'Exiting. If you wish to obtain the .fastq file(s), run: getfastq --id {}\n'
        sys.stderr.write(txt.format(sra_stat['sra_id']))
        sys.exit(1)
    in_files = glob.glob(os.path.join(args.out_dir, 'getfastq', sra_id, sra_id + "*" + ext))
    assert in_files, '{}: Fastq file not found. Check {}'.format(sra_id, output_dir)
    print('Input fastq detected:', ', '.join(in_files))
    call_kallisto(args, in_files, metadata, sra_stat, output_dir, index)
    if (args.clean_fastq & quant_output_exists(sra_id, output_dir)):
        print('Safe-deleting getfastq files.', flush=True)
        for in_file in in_files:
            print('Output file detected. Safely removing fastq:', in_file)
            os.remove(in_file)
            placeholder = open(in_file + '.safely_removed', "w")
            placeholder.write("This fastq file was safely removed after `amalgkit quant`.")
            placeholder.close()
    else:
        print('Skipping the deletion of getfastq files.', flush=True)

def get_index(args, sci_name):
    if args.index_dir is not None:
        index_dir = args.index_dir
    else:
        index_dir = os.path.join(args.out_dir, 'index')
    if not os.path.exists(index_dir) and args.build_index:
            os.mkdir(index_dir)
    if not os.path.exists(index_dir):
        raise FileNotFoundError("Could not find index folder at: {}".format(index_dir))

    index = glob.glob(os.path.join(index_dir, sci_name + '*'))
    if len(index) > 1:
        raise ValueError(
            "Found multiple index files for species. Please make sure there is only one index file for this species.")
    elif len(index) == 0:
        if args.build_index:
            print("--build_index set. Building index for {}".format(sci_name))
            if (args.fasta_dir=='inferred'):
                path_fasta_dir = os.path.join(args.out_dir, 'fasta')
            else:
                path_fasta_dir = args.fasta_dir
            fasta_files = []
            for ext in ['*.fa','*.fasta','*.fa.gz','*.fasta.gz',]:
                path_pattern = os.path.join(path_fasta_dir, sci_name + ext)
                fasta_files.extend(glob.glob(path_pattern))
            if len(fasta_files) > 1:
                txt = "Found multiple reference fasta files for this species: {}\n"
                txt += "Please make sure there is only one index file for this species.\n{}"
                raise ValueError(txt.format(', '.join(sci_name, fasta_files)))
            elif len(fasta_files) == 0:
                txt = "Could not find reference fasta file for this species: {}\n".format(sci_name)
                txt += 'If the reference fasta file is correctly placed, the column "scientific_name" of the --metadata file may need to be edited.'
                raise FileNotFoundError(txt)
            fasta_file = fasta_files[0]
            print('Reference fasta file found: {}'.format(fasta_file), flush=True)
            index_path = os.path.join(index_dir, sci_name + '.idx')
            print('Building index: {}'.format(index_path), flush=True)
            kallisto_build_cmd = ["kallisto", "index", "-i", index_path, fasta_file]
            print('Command: {}'.format(' '.join(kallisto_build_cmd)), flush=True)
            index_out = subprocess.run(kallisto_build_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print('kallisto index stdout:')
            print(index_out.stdout.decode('utf8'))
            print('kallisto index stderr:')
            print(index_out.stderr.decode('utf8'))
            index = [index_path]
        else:
            sys.stderr.write('No index file was found in: {}\n'.format(index_dir))
            sys.stderr.write('Try --fasta_dir PATH and --build_index yes\n')
            raise FileNotFoundError("Could not find index file.")
    index = index[0]
    print("Kallisto index file found: {}".format(index), flush=True)
    return index


def quant_main(args):
    check_kallisto_dependency()
    metadata = load_metadata(args)
    for i in metadata.df.index:
        print('')
        sra_id = metadata.df.at[i, 'run']
        sci_name = metadata.df.at[i, 'scientific_name']
        print('Species: {}'.format(sci_name))
        print('SRA Run ID: {}'.format(sra_id))
        sci_name = sci_name.replace(" ", "_")
        print('Looking for index folder in ', args.out_dir)
        index = get_index(args, sci_name)
        run_quant(args, metadata, sra_id, index)



================================================
FILE: amalgkit/sanity.py
================================================
import numpy as np
from amalgkit.util import *
import glob


def list_duplicates(seq):
    seen = set()
    seen_add = seen.add
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set(x for x in seq if x in seen or seen_add(x))
    # turn the set into a list (as requested)
    return list(seen_twice)


def parse_metadata(args, metadata):
    print("Checking essential entries from metadata file.")
    species = metadata.df.loc[:, 'scientific_name']
    uni_species = np.unique(species)
    if len(uni_species):
        print(len(uni_species), " species detected:")
        print(uni_species)
    else:
        txt = "{} species detected. Please check if --metadata ({}) has a 'scientific_name' column."
        raise ValueError(txt.format(len(uni_species), args.metadata))

    sra_ids = metadata.df.loc[:, 'run']
    if len(sra_ids):
        print(len(np.unique(sra_ids)), " SRA runs detected:")
        print(np.unique(sra_ids))
        # check for duplicate runs
        if len(sra_ids) > len(np.unique(sra_ids)):
            raise ValueError("Duplicate SRA IDs detected, where IDs should be unique. Please check these entries: ",
                             list_duplicates(sra_ids))
    else:
        txt = "{} SRA runs detected. Please check if --metadata ({}) has a 'run' column."
        raise ValueError(txt.format(len(np.unique(sra_ids)), args.metadata))
    return uni_species, sra_ids


def check_getfastq_outputs(args, sra_ids, metadata, output_dir):
    print("checking for getfastq outputs: ")
    if args.getfastq_dir:
        getfastq_path = args.getfastq
    else:
        getfastq_path = os.path.join(args.out_dir, "getfastq")
    data_available = []
    data_unavailable = []
    if os.path.exists(getfastq_path):
        print("amalgkit getfastq output folder detected. Checking presence of output files.")
        for sra_id in sra_ids:
            print("\n")
            print("Looking for {}".format(sra_id))
            sra_path = os.path.join(getfastq_path, sra_id)
            if os.path.exists(sra_path):
                sra_stat = get_sra_stat(sra_id, metadata)
                try:
                    ext = get_newest_intermediate_file_extension(sra_stat, sra_path)
                except FileNotFoundError:
                    print("could not find any fastq files for ", sra_id,
                          "Please make sure amalgkit getfastq ran properly")
                    data_unavailable.append(sra_id)
                    continue

                files = glob.glob(os.path.join(sra_path, sra_id + "*" + ext))
                if ext != '.safely_removed':
                    print("Found:", files)
                data_available.append(sra_id)

            else:
                print("Could not find getfastq output for: ", sra_id, "\n")
                print("Suggested command for rerun: getfastq -e email@adress.com --id ", sra_id, " -w ", args.out_dir, "--redo yes --gcp yes --aws yes --ncbi yes")
                data_unavailable.append(sra_id)

    else:
        print("Could not find getfastq output folder ", getfastq_path, ". Have you run getfastq yet?")
        data_unavailable = metadata.df['run'].tolist()

    if data_unavailable:
        print("writing SRA IDs without getfastq output to: ", os.path.join(output_dir, "SRA_IDs_without_fastq.txt"))
        file = open(os.path.join(output_dir, "SRA_IDs_without_fastq.txt"), "w")
        for sra_id in data_unavailable:
            file.write(sra_id + "\n")
        file.close()
    else:
        txt = "The getfastq output files for all SRA IDs in --metadata ({}) were found."
        print(txt.format(args.metadata))


    return data_available, data_unavailable


def check_quant_index(args, uni_species, output_dir):
    if args.index_dir:
        index_dir_path = args.index_dir
    else:
        index_dir_path = os.path.join(args.out_dir, "index")
    index_unavailable = []
    index_available = []
    if os.path.exists(index_dir_path):
        for species in uni_species:
            sci_name = species.replace(" ", "_")
            sci_name = sci_name.replace(".", "")
            index_path = os.path.join(index_dir_path, sci_name + "*")
            print("\n")
            print("Looking for index file {} for species {}".format(index_path, species))
            index_files = glob.glob(index_path)
            if not index_files:
                print("could not find anything in", index_path)
                sci_name = species.split(" ")
                # Deprecate subspecies or variants and look again
                # I.e. if Gorilla_gorilla_gorilla.idx was not found, we look for Gorilla_gorilla.idx instead.
                if len(sci_name) > 2:
                    sci_name = sci_name[0] + "_" + sci_name[1]
                    print("Ignoring subspecies.")
                    index_path = os.path.join(index_dir_path, sci_name + "*")
                    print("Looking for {}".format(index_path))
                    index_files = glob.glob(index_path)
                    if index_files:
                        print("Found ", index_files, "!")
                        index_available.append(species)
                    else:
                        print("Could not find any index files for ", species)
                        index_unavailable.append(species)
                else:
                    print("Could not find any index files for ", species)
                    index_unavailable.append(species)

            else:
                print("Found ", index_files, "!")
                index_available.append(species)

            if len(index_files) > 1:
                print("Multiple possible index files detected for ", species, ": ", index_files,
                      ". You may have to resolve ambiguity")

        if index_unavailable:
            print("writing species without index to: ", os.path.join(output_dir, "species_without_index.txt"))
            file = open(os.path.join(output_dir, "species_without_index.txt"), "w")
            for species in index_unavailable:
                file.write(species + "\n")
            file.close()
        else:
            print("index found for all species in --metadata ({})".format(args.metadata))
    else:
        print("Could not find index directory ", index_dir_path, " . Did you provide the correct Path?")

    return index_available, index_unavailable


def check_quant_output(args, sra_ids, output_dir):
    print("checking for quant outputs: ")
    quant_path = os.path.join(args.out_dir, "quant")
    data_available = []
    data_unavailable = []

    if os.path.exists(quant_path):
        print("amalgkit quant output folder detected. Checking presence of output files.")
        for sra_id in sra_ids:
            warned = []
            print("\n")
            print("Looking for {}".format(sra_id))
            sra_path = os.path.join(quant_path, sra_id)
            if os.path.exists(sra_path):
                print("Found output folder ", sra_path, " for ", sra_id)
                print("Checking for output files.")
                abundance_file = os.path.join(sra_path, sra_id + "_abundance.tsv")
                run_info_file = os.path.join(sra_path, sra_id + "_run_info.json")

                if os.path.exists(abundance_file) and os.path.exists(run_info_file):
                    print("All quant output files present for", sra_id, "!")
                    data_available.append(sra_id)
                    continue
                elif not os.path.exists(abundance_file):
                    print(abundance_file, " is missing! Please check if quant ran correctly")
                    warned = True
                elif not os.path.exists(run_info_file):
                    print(run_info_file, " is missing! Please check if quant ran correctly")
                    warned = True

                if warned:
                    data_unavailable.append(sra_id)

            else:
                print("Could not find output folder ", sra_path, " for ", sra_id)
                data_unavailable.append(sra_id)
    else:
        print("Could not find quant output folder ", quant_path, ". Have you run quant yet?")

    if data_unavailable:
        print("writing SRA IDs without quant output to: ", os.path.join(output_dir, "SRA_IDs_without_quant.txt"))
        file = open(os.path.join(output_dir, "SRA_IDs_without_quant.txt"), "w")
        for sra_id in data_unavailable:
            file.write(sra_id + "\n")
        file.close()
    else:
        print("Quant outputs found for all SRA IDs in --metadata ({})".format(args.metadata))

    return data_available, data_unavailable


def sanity_main(args):
    metadata = load_metadata(args)
    output_dir = os.path.join(args.out_dir, 'sanity')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    uni_species, sra_ids = parse_metadata(args, metadata)
    if args.getfastq or args.all:
        check_getfastq_outputs(args, sra_ids, metadata, output_dir)
    if args.index or args.all:
        check_quant_index(args, uni_species, output_dir)
    if args.quant or args.all:
        check_quant_output(args, sra_ids, output_dir)



================================================
FILE: amalgkit/select.py
================================================
import shutil

from amalgkit.util import *

def write_select_outputs(path_metadata_original, path_metadata_table, metadata_dir, metadata):
    if os.path.exists(path_metadata_original):
        print('Original metadata copy exists and will not be renewed: {}'.format(path_metadata_original), flush=True)
    else:
        print('Original metadata copy does not exist. Creating: {}'.format(path_metadata_original), flush=True)
        shutil.copyfile(path_metadata_table, path_metadata_original)
    print('Updating metadata table at: {}'.format(path_metadata_table), flush=True)
    metadata.df.to_csv(path_metadata_table, sep='\t', index=False)
    sra_qualified_pivot = metadata.pivot(n_sp_cutoff=0, qualified_only=True, sampled_only=False)
    sra_qualified_pivot.to_csv(os.path.join(metadata_dir, 'pivot_qualified.tsv'), sep='\t')
    sra_selected_pivot = metadata.pivot(n_sp_cutoff=0, qualified_only=True, sampled_only=True)
    sra_selected_pivot.to_csv(os.path.join(metadata_dir, 'pivot_selected.tsv'), sep='\t')

def select_main(args):
    metadata_dir = os.path.join(args.out_dir, 'metadata')
    if args.config_dir=='inferred':
        dir_config = os.path.join(args.out_dir, 'config')
    else:
        dir_config = args.config_dir
    check_config_dir(dir_path=dir_config, mode='select')
    path_metadata_table = os.path.join(metadata_dir, 'metadata.tsv')
    path_metadata_original = os.path.join(metadata_dir, 'metadata_original.tsv')

    metadata = load_metadata(args)
    if args.sample_group is not None:
        txt = '{}: Extracting pre-selected sample_group entries: {}'
        print(txt.format(datetime.datetime.now(), args.sample_group), flush=True)
        selected_sample_groups = args.sample_group.split(',')
        metadata.df = metadata.df.loc[metadata.df['sample_group'].isin(selected_sample_groups),:].reset_index(drop=True)
    metadata.nspot_cutoff(args.min_nspots)
    metadata.mark_redundant_biosample(args.mark_redundant_biosamples)
    metadata.remove_specialchars()
    metadata.group_attributes(dir_config)
    metadata.mark_exclude_keywords(dir_config)
    metadata.mark_treatment_terms(dir_config)
    metadata.label_sampled_data(args.max_sample)
    metadata.reorder(omit_misc=True)
    write_select_outputs(path_metadata_original, path_metadata_table, metadata_dir, metadata)



================================================
FILE: amalgkit/util.py
================================================
import json
import numpy
import pandas
import lxml.etree
import datetime

import glob
import inspect
import os
import re
import subprocess
import sys
import warnings

def strtobool(val):
    val = val.lower()
    if val in ("y", "yes", "t", "true", "on", "1"):
        return True
    elif val in ("n", "no", "f", "false", "off", "0"):
        return False
    else:
        raise ValueError(f"invalid truth value {val!r}")

class Metadata:
    column_names = ['scientific_name', 'tissue', 'sample_group', 'genotype', 'sex', 'age',
                    'treatment', 'source_name',
                    'is_sampled', 'is_qualified', 'exclusion', 'protocol', 'bioproject', 'biosample',
                    'experiment', 'run', 'sra_primary', 'sra_sample', 'sra_study', 'study_title', 'exp_title', 'design',
                    'sample_title', 'sample_description', 'lib_name', 'lib_layout', 'lib_strategy', 'lib_source',
                    'lib_selection', 'instrument', 'total_spots', 'total_bases', 'size', 'nominal_length',
                    'nominal_sdev',
                    'spot_length', 'read_index', 'read_class', 'read_type', 'base_coord', 'lab', 'center',
                    'submitter_id',
                    'pubmed_id', 'taxid', 'published_date', 'biomaterial_provider', 'cell', 'location', 'antibody',
                    'batch',
                    'misc', 'NCBI_Link', 'AWS_Link', 'GCP_Link', ]
    id_cols = ['bioproject', 'biosample', 'experiment', 'run', 'sra_primary', 'sra_sample', 'sra_study']

    def __init__(self, column_names=column_names):
        self.config_dir = ''
        self.df = pandas.DataFrame(index=[], columns=column_names)

    def reorder(self, omit_misc=False, column_names=column_names):
        if (self.df.shape[0] == 0):
            return None
        self.df.loc[:, [col for col in column_names if col not in self.df.columns]] = ''
        if omit_misc:
            self.df = self.df.loc[:, column_names]
        else:
            misc_columns = [col for col in self.df.columns if col not in column_names]
            self.df = self.df.loc[:, column_names + misc_columns]
        self.df.loc[:, 'exclusion'] = self.df.loc[:, 'exclusion'].replace('', 'no')
        # reorder sample_group to the front
        if 'sample_group' in self.df.columns:
            cols = list(self.df)
            cols.insert(1, cols.pop(cols.index('sample_group')))
            self.df = self.df.loc[:, cols]
        self.df = self.df.reset_index(drop=True)

    def from_DataFrame(df):
        metadata = Metadata()
        metadata.df = df
        metadata.reorder(omit_misc=False)
        return metadata

    def from_xml(xml_root):
        if isinstance(xml_root, lxml.etree._Element):
            xml_root = lxml.etree.ElementTree(xml_root)
        root = xml_root
        assert isinstance(root, lxml.etree._ElementTree), "Unknown input type."
        df_list = list()
        counter = 0
        metadata = Metadata()
        for entry in root.iter(tag="EXPERIMENT_PACKAGE"):
            if counter % 1000 == 0:
                now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                print('{}: Converting {:,}th sample from XML to DataFrame'.format(now, counter), flush=True)
            items = []
            bioproject = entry.findall('.//EXTERNAL_ID[@namespace="BioProject"]')
            if not len(bioproject):
                labels = entry.findall('.//LABEL')
                for label in labels:
                    text = label.text
                    if text.startswith("PRJ"):
                        bioproject = [label]
                        break
            is_single = len(entry.findall('.//LIBRARY_LAYOUT/SINGLE'))
            is_paired = len(entry.findall('.//LIBRARY_LAYOUT/PAIRED'))
            if is_single:
                library_layout = ["single"]
            elif is_paired:
                library_layout = ["paired"]
            else:
                library_layout = [""]
            values = entry.findall('.//VALUE')
            is_protected = ["No"]
            if len(values):
                for value in values:
                    text = value.text
                    if not text is None:
                        if text.endswith("PROTECTED"):
                            is_protected = ["Yes"]
                            break
            items.append(["bioproject", bioproject])
            items.append(["scientific_name", entry.xpath('./SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME')])
            items.append(["biosample", entry.findall('.//EXTERNAL_ID[@namespace="BioSample"]')])
            items.append(["experiment", entry.xpath('./EXPERIMENT/IDENTIFIERS/PRIMARY_ID')])
            items.append(["run", entry.xpath('./RUN_SET/RUN/IDENTIFIERS/PRIMARY_ID')])
            items.append(["sra_primary", entry.xpath('./SUBMISSION/IDENTIFIERS/PRIMARY_ID')])
            items.append(["sra_sample", entry.xpath('./SAMPLE/IDENTIFIERS/PRIMARY_ID')])
            items.append(["sra_study", entry.xpath('./EXPERIMENT/STUDY_REF/IDENTIFIERS/PRIMARY_ID')])
            items.append(["published_date", entry.xpath('./RUN_SET/RUN/@published')])
            items.append(["exp_title", entry.xpath('./EXPERIMENT/TITLE')])
            items.append(["design", entry.xpath('./EXPERIMENT/DESIGN/DESIGN_DESCRIPTION')])
            items.append(["lib_name", entry.xpath('./EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_NAME')])
            items.append(["lib_strategy", entry.xpath('./EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY')])
            items.append(["lib_source", entry.xpath('./EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SOURCE')])
            items.append(["lib_selection", entry.xpath('./EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SELECTION')])
            items.append(["lib_layout", library_layout])
            items.append(["nominal_length",
                          entry.xpath('./EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED/@NOMINAL_LENGTH')])
            items.append(["nominal_sdev",
                          entry.xpath('./EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED/@NOMINAL_SDEV')])
            items.append(
                ["spot_length", entry.xpath('./EXPERIMENT/DESIGN/SPOT_DESCRIPTOR/SPOT_DECODE_SPEC/SPOT_LENGTH')])
            items.append(["read_index",
                          entry.xpath('./EXPERIMENT/DESIGN/SPOT_DESCRIPTOR/SPOT_DECODE_SPEC/READ_SPEC/READ_INDEX')])
            items.append(["read_class",
                          entry.xpath('./EXPERIMENT/DESIGN/SPOT_DESCRIPTOR/SPOT_DECODE_SPEC/READ_SPEC/READ_CLASS')])
            items.append(
                ["read_type", entry.xpath('./EXPERIMENT/DESIGN/SPOT_DESCRIPTOR/SPOT_DECODE_SPEC/READ_SPEC/READ_TYPE')])
            items.append(["base_coord",
                          entry.xpath('./EXPERIMENT/DESIGN/SPOT_DESCRIPTOR/SPOT_DECODE_SPEC/READ_SPEC/BASE_COORD')])
            items.append(["instrument", entry.xpath('./EXPERIMENT/PLATFORM/ILLUMINA/INSTRUMENT_MODEL')])
            items.append(["lab", entry.xpath('./SUBMISSION/@lab_name')])
            items.append(["center", entry.xpath('./SUBMISSION/@center_name')])
            items.append(["submitter_id", entry.xpath('./SUBMISSION/IDENTIFIERS/SUBMITTER_ID')])
            items.append(["study_title", entry.xpath('./STUDY/DESCRIPTOR/STUDY_TITLE')])
            items.append(["pubmed_id", entry.xpath('./STUDY/STUDY_LINKS/STUDY_LINK/XREF_LINK/ID')])
            items.append(["sample_title", entry.xpath('./SAMPLE/TITLE')])
            items.append(["taxid", entry.xpath('./SAMPLE/SAMPLE_NAME/TAXON_ID')])
            items.append(["sample_description", entry.xpath('./SAMPLE/DESCRIPTION')])
            items.append(["total_spots", entry.xpath('./RUN_SET/RUN/@total_spots')])
            items.append(["total_bases", entry.xpath('./RUN_SET/RUN/@total_bases')])
            items.append(["size", entry.xpath('./RUN_SET/RUN/@size')])
            items.append(["NCBI_Link", entry.xpath(
                './RUN_SET/RUN/SRAFiles/SRAFile[@supertype="Primary ETL"]/Alternatives[@org="NCBI"]/@url')])
            items.append(["AWS_Link", entry.xpath(
                './RUN_SET/RUN/SRAFiles/SRAFile[@supertype="Primary ETL"]/Alternatives[@org="AWS"]/@url')])
            items.append(["GCP_Link", entry.xpath(
                './RUN_SET/RUN/SRAFiles/SRAFile[@supertype="Primary ETL"]/Alternatives[@org="GCP"]/@url')])
            row = []
            for item in items:
                try:
                    if isinstance(item[1][0], (lxml.etree._ElementUnicodeResult, int, str)):
                        row.append(str(item[1][0]))
                    else:
                        row.append(item[1][0].text)
                except:
                    row.append("")
            colnames = []
            for item in items:
                colnames.append(item[0])
            row_df = pandas.DataFrame(row).T
            row_df.columns = colnames
            sas = entry.xpath('./SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE')
            for sa in sas:
                tag = sa.xpath('./TAG')
                if not tag[0].text == None:
                    tag = tag[0].text.lower()
                    tag = re.sub(r" \(.*", "", tag)
                    tag = re.sub(r" ", "_", tag)
                    if not tag in row_df.columns:
                        value = sa.xpath('./VALUE')
                        if len(value):
                            value = value[0].text
                            if tag in colnames:
                                tag = tag + "_2"
                            sa_df = pandas.DataFrame([value])
                            sa_df.columns = [tag]
                            row_df = pandas.concat([row_df, sa_df], axis=1)
            df_list.append(row_df)
            counter += 1
        if len(df_list)==0:
            return metadata
        if len(df_list) <= 1000:
            df = pandas.concat(df_list, ignore_index=True)
        else:
            chunked = [pandas.concat(df_list[i:i+1000], ignore_index=True) for i in range(0, len(df_list), 1000)]
            df = pandas.concat(chunked, ignore_index=True)
        now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print('{}: Finished converting {:,} samples'.format(now, counter), flush=True)
        metadata.df = df
        metadata.reorder(omit_misc=False)
        return metadata

    def group_attributes(self, dir_config):
        try:
            config = pandas.read_csv(os.path.join(dir_config, 'group_attribute.config'),
                                     parse_dates=False, quotechar='"', sep='\t',
                                     header=None, index_col=None, skip_blank_lines=True, comment='#')
        except:
            config = pandas.DataFrame()
        config = config.replace(numpy.nan, '')
        for i in config.index:
            aggregate_to = config.iloc[i, 0]
            aggregate_from = config.iloc[i, 1]
            if (aggregate_from in self.df.columns) & (aggregate_from != ''):
                print('{}: Aggregating column "{}" to column "{}"'.format(datetime.datetime.now(), aggregate_from, aggregate_to), flush=True)
                if not aggregate_to in self.df.columns:
                    self.df.loc[:, aggregate_to] = ''
                is_from_empty = (self.df.loc[:, aggregate_from].isnull()) | (
                            self.df.loc[:, aggregate_from].astype(str) == '')
                is_to_empty = (self.df.loc[:, aggregate_to].isnull()) | (self.df.loc[:, aggregate_to].astype(str) == '')
                new_annotations = self.df.loc[(~is_from_empty) & (is_to_empty), aggregate_from].astype(
                    str) + '[' + aggregate_from + ']'
                self.df.loc[(~is_from_empty) & (is_to_empty), aggregate_to] = new_annotations
                new_annotations = self.df.loc[(~is_from_empty) & (~is_to_empty), aggregate_to].astype(str) + "; " + \
                                  self.df.loc[(~is_from_empty) & (~is_to_empty), aggregate_from].astype(
                                      str) + '[' + aggregate_from + ']'
                self.df.loc[(~is_from_empty) & (~is_to_empty), aggregate_to] = new_annotations
                self.df = self.df.drop(labels=aggregate_from, axis=1)
        self.reorder(omit_misc=False)

    def mark_exclude_keywords(self, dir_config):
        try:
            config = pandas.read_csv(os.path.join(dir_config, 'exclude_keyword.config'),
                                     parse_dates=False, quotechar='"', sep='\t',
                                     header=None, index_col=None, skip_blank_lines=True, comment='#')
        except:
            config = pandas.DataFrame()
        config = config.replace(numpy.nan, '')
        if config.shape[0]>0:
            print('{}: Marking SRAs with bad keywords'.format(datetime.datetime.now()), flush=True)
        for i in config.index:
            cols = config.iloc[i, 0].split(',')
            reason = config.iloc[i, 1]
            exclude_keyword = config.iloc[i, 2]
            num_detected = 0
            for col in cols:
                has_bad_keyword = self.df.loc[:, col].astype(str).str.contains(exclude_keyword, regex=True, case=False).fillna(False)
                self.df.loc[has_bad_keyword, 'exclusion'] = reason
                num_detected += has_bad_keyword.sum()
            txt = '{}: Marking {:,} SRAs with keyword "{}"'
            print(txt.format(datetime.datetime.now(), num_detected, exclude_keyword), flush=True)
        self.df.loc[~((self.df.loc[:, 'antibody'].isnull()) | (
                    self.df.loc[:, 'antibody'] == '')), 'exclusion'] = 'immunoprecipitation'
        self.df.loc[~((self.df.loc[:, 'cell'].isnull()) | (self.df.loc[:, 'cell'] == '')), 'exclusion'] = 'cell_culture'

    def mark_treatment_terms(self, dir_config):
        try:
            config = pandas.read_csv(os.path.join(dir_config, 'control_term.config'),
                                     parse_dates=False, quotechar='"', sep='\t',
                                     header=None, index_col=None, skip_blank_lines=True, comment='#')
        except:
            config = pandas.DataFrame()
        config = config.replace(numpy.nan, '')
        if config.shape[0]>0:
            print('{}: Marking SRAs with non-control terms'.format(datetime.datetime.now()), flush=True)
        for i in config.index:
            cols = config.iloc[i, 0].split(',')
            control_term = config.iloc[i, 1]
            if control_term == '':
                continue
            num_control = 0
            num_treatment = 0
            for col in cols:
                is_control = self.df.loc[:, col].astype(str).str.contains(control_term, regex=True, case=False).fillna(False)
                if not any(is_control):
                    continue
                bioprojects = self.df.loc[is_control, 'bioproject'].unique()
                for bioproject in bioprojects:
                    is_bioproject = (self.df.loc[:, 'bioproject'] == bioproject)
                    self.df.loc[(is_bioproject & -is_control), 'exclusion'] = 'non_control'
                    num_control += (is_bioproject & is_control).sum()
                    num_treatment += (is_bioproject & -is_control).sum()
            txt = '{}: Applying control term "{}" to "{}": Detected control and treatment SRAs: {:,} and {:,}'
            print(txt.format(datetime.datetime.now(), control_term, ','.join(cols), num_control, num_treatment), flush=True)

    def nspot_cutoff(self, min_nspots):
        print('{}: Marking SRAs with less than {:,} reads'.format(datetime.datetime.now(), min_nspots), flush=True)
        self.df['total_spots'] = self.df.loc[:, 'total_spots'].replace('', 0)
        self.df['total_spots'] = self.df.loc[:, 'total_spots'].fillna(0).astype(int)
        self.df.loc[-(self.df.loc[:, 'total_spots'] == 0) & (
                    self.df.loc[:, 'total_spots'] < min_nspots), 'exclusion'] = 'low_nspots'

    def mark_redundant_biosample(self, exe_flag):
        if exe_flag:
            print('{}: Marking SRAs with redundant BioSample IDs'.format(datetime.datetime.now()), flush=True)
            redundant_bool = self.df.duplicated(subset=['bioproject', 'biosample'], keep='first')
            self.df.loc[redundant_bool, 'exclusion'] = 'redundant_biosample'

    def _maximize_bioproject_sampling(self, df, target_n=10):
        while len(df.loc[(df.loc[:, 'is_sampled'] == 'yes') & (df.loc[:, 'exclusion'] == 'no'), :]) < target_n:
            if len(df) <= target_n:
                df.loc[(df.loc[:, 'exclusion'] == 'no'), 'is_sampled'] = 'yes'
                break
            else:
                df_unselected = df.loc[(df.loc[:, 'is_sampled'] == 'no') & (df.loc[:, 'exclusion'] == 'no'), :]
                bioprojects = df_unselected.loc[:, 'bioproject'].unique()
                if len(bioprojects) == 0:
                    break
                remaining_n = target_n - (df.loc[:, 'is_sampled'] == 'yes').sum()
                select_n = min([len(bioprojects), remaining_n])
                selected_bioprojects = numpy.random.choice(bioprojects, size=select_n, replace=False)
                selected_index = []
                for bioproject in selected_bioprojects:
                    is_bp = (df_unselected.loc[:, 'bioproject'] == bioproject)
                    index = numpy.random.choice(df_unselected.index[is_bp], size=1, replace=False)
                    selected_index.append(int(index))
                df.loc[selected_index, 'is_sampled'] = 'yes'
        return df

    def label_sampled_data(self, max_sample=10):
        pandas.set_option('mode.chained_assignment', None)
        txt = '{}: Selecting subsets of SRA IDs for >{:,} samples per sample_group per species'
        print(txt.format(datetime.datetime.now(), max_sample), flush=True)
        is_empty = (self.df['sample_group'] == '')
        self.df.loc[is_empty,'is_qualified'] = 'no'
        self.df.loc[is_empty,'exclusion'] = 'no_tissue_label'
        self.df['bioproject'] = self.df['bioproject'].fillna('unknown').values
        self.df['is_sampled'] = 'no'
        self.df['is_qualified'] = 'no'
        self.df.loc[(self.df.loc[:, 'exclusion'] == 'no'), 'is_qualified'] = 'yes'
        df_list = list()
        species = self.df.loc[:, 'scientific_name'].unique()
        for sp in species:
            sp_table = self.df.loc[(self.df.loc[:, 'scientific_name'] == sp), :]
            sample_groups = sp_table.loc[:, 'sample_group'].unique()
            for sample_group in sample_groups:
                sp_sample_group = sp_table.loc[(sp_table.loc[:, 'sample_group'] == sample_group), :]
                if sp_sample_group.shape[0] == 0:
                    continue
                sp_sample_group = self._maximize_bioproject_sampling(df=sp_sample_group, target_n=max_sample)
                df_list.append(sp_sample_group)
        if len(df_list) <= 100:
            self.df = pandas.concat(df_list, ignore_index=True)
        else:
            chunked = [pandas.concat(df_list[i:i+100], ignore_index=True) for i in range(0, len(df_list), 100)]
            self.df = pandas.concat(chunked, ignore_index=True)
        self.reorder(omit_misc=False)
        pandas.set_option('mode.chained_assignment', 'warn')

    def remove_specialchars(self):
        for col, dtype in zip(self.df.dtypes.index, self.df.dtypes.values):
            if any([key in str(dtype) for key in ['str', 'object']]):
                self.df.loc[:, col] = self.df[col].replace(r'\r', '', regex=True)
                self.df.loc[:, col] = self.df[col].replace(r'\n', '', regex=True)
                self.df.loc[:, col] = self.df[col].replace(r'\'', '', regex=True)
                self.df.loc[:, col] = self.df[col].replace(r'\"', '', regex=True)
                self.df.loc[:, col] = self.df[col].replace(r'\|', '', regex=True)

    def pivot(self, n_sp_cutoff=0, qualified_only=True, sampled_only=False):
        df = self.df
        if qualified_only:
            df = df.loc[(df.loc[:, 'is_qualified'] == 'yes'), :]
        if sampled_only:
            df = df.loc[(df.loc[:, 'is_sampled'] == 'yes'), :]
        df_reduced = df.loc[:, ['scientific_name', 'biosample', 'sample_group']]
        pivot = df_reduced.pivot_table(columns='sample_group', index='scientific_name', aggfunc='count')
        pivot.columns = pivot.columns.get_level_values(1)
        column_sort = pivot.count(axis='index').sort_values(ascending=False).index
        index_sort = pivot.count(axis='columns').sort_values(ascending=False).index
        pivot = pivot.loc[index_sort, column_sort]
        pivot_reduced = pivot.loc[:, pivot.count(axis='index') >= n_sp_cutoff]
        column_sort = pivot_reduced.count(axis='index').sort_values(ascending=False).index
        index_sort = pivot_reduced.count(axis='columns').sort_values(ascending=False).index
        pivot_reduced = pivot_reduced.loc[index_sort, column_sort]
        return pivot_reduced

def read_config_file(file_name, dir_path):
    try:
        df = pandas.read_csv(os.path.join(dir_path, file_name),
                             parse_dates=False, quotechar='"', sep='\t',
                             header=None, index_col=None, skip_blank_lines=True, comment='#')
    except:
        df = pandas.DataFrame([])
    if df.shape[1]==1:
        df = df.iloc[:,0]
    return df

def load_metadata(args, dir_subcommand='metadata'):
    if args.metadata=='inferred':
        relative_path = os.path.join(args.out_dir, dir_subcommand, 'metadata.tsv')
        real_path = os.path.realpath(relative_path)
    else:
        real_path = os.path.realpath(args.metadata)
    print('{}: Loading metadata from: {}'.format(datetime.datetime.now(), real_path), flush=True)
    df = pandas.read_csv(real_path, sep='\t', header=0, low_memory=False)
    metadata = Metadata.from_DataFrame(df)
    if 'batch' not in dir(args):
        return metadata
    if args.batch is None:
        return metadata
    # --batch must be handled species-wise in curate.py
    # so we need to find out where the call came from
    frm = inspect.stack()[1]
    mod = inspect.getmodule(frm[0])
    if mod.__name__ == 'amalgkit.curate':
        print('Entering --batch mode for amalgkit curate. processing 1 species', flush=True)
        txt = 'This is {:,}th job. In total, {:,} jobs will be necessary for this metadata table.'
        spp = metadata.df.loc[:, 'scientific_name'].drop_duplicates().sort_values().values
        print(txt.format(args.batch, len(spp)), flush=True)
        sp = spp[args.batch - 1]
        print('Processing species: {}'.format(sp), flush=True)
        is_sp = (metadata.df['scientific_name'] == sp)
        metadata.df = metadata.df.loc[is_sp,:].reset_index(drop=True)
        return metadata
    else:
        print('--batch is specified. Processing one SRA per job.', flush=True)
        is_sampled = numpy.array([strtobool(yn) for yn in df.loc[:, 'is_sampled']], dtype=bool)
        txt = 'This is {:,}th job. In total, {:,} jobs will be necessary for this metadata table. {:,} '
        txt += 'SRAs were excluded from the table (is_sampled==no).'
        print(txt.format(args.batch, sum(is_sampled), len(numpy.where(is_sampled == False)[0])), flush=True)
        if args.batch>sum(is_sampled):
            sys.stderr.write('--batch {} is too large. Exiting.\n'.format(args.batch))
            sys.exit(0)
        if is_sampled.sum()==0:
            print('No sample is "sampled". Please check the "is_sampled" column in the metadata. Exiting.')
            sys.exit(1)
        metadata.df = metadata.df.loc[is_sampled,:]
        metadata.df = metadata.df.reset_index()
        metadata.df = metadata.df.loc[[args.batch-1,],:]
        return metadata

def get_sra_stat(sra_id, metadata, num_bp_per_sra=None):
    sra_stat = dict()
    sra_stat['sra_id'] = sra_id
    is_sra = (metadata.df.loc[:,'run']==sra_id)
    assert is_sra.sum()==1, 'There are multiple metadata rows with the same SRA ID: '+sra_id
    sra_stat['layout'] = metadata.df.loc[is_sra,'lib_layout'].values[0]
    sra_stat['total_spot'] = int(metadata.df.loc[is_sra,'total_spots'].values[0])
    original_spot_len = metadata.df.loc[is_sra,'spot_length'].values[0]
    if (numpy.isnan(original_spot_len) | (original_spot_len==0)):
        inferred_spot_len = int(metadata.df.loc[is_sra,'total_bases'].values[0]) / int(sra_stat['total_spot'])
        sra_stat['spot_length'] = int(inferred_spot_len)
        txt = 'spot_length cannot be obtained directly from metadata. Using total_bases/total_spots instead: {:,}'
        print(txt.format(sra_stat['spot_length']))
    else:
        sra_stat['spot_length'] = int(original_spot_len)
    if num_bp_per_sra is not None:
        sra_stat['num_read_per_sra'] = int(num_bp_per_sra/sra_stat['spot_length'])
    return sra_stat

def get_newest_intermediate_file_extension(sra_stat, work_dir):
    ext_out = 'no_extension_found'
    # Order is important in this list. More downstream should come first.
    extensions = ['.amalgkit.fastq.gz','.rename.fastq.gz','.fastp.fastq.gz','.fastq.gz']
    sra_stat = detect_layout_from_file(sra_stat)
    if sra_stat['layout']=='single':
        subext = ''
    elif sra_stat['layout']=='paired':
        subext = '_1'
    files = os.listdir(work_dir)
    for ext in extensions:
        if any([ f==sra_stat['sra_id']+subext+ext for f in files ]):
            ext_out = ext
            break
    if ext_out == 'no_extension_found':
        safe_delete_files = glob.glob(os.path.join(work_dir, sra_stat['sra_id']+"*.safely_removed"))
        if len(safe_delete_files):
            txt = 'getfastq safely_removed flag was detected. `amalgkit quant` has been completed in this sample: {}\n'
            sys.stdout.write(txt.format(work_dir))
            for safe_delete_file in safe_delete_files:
                sys.stdout.write('{}\n'.format(safe_delete_file))
            return '.safely_removed'
    return ext_out

def is_there_unpaired_file(sra_stat, extensions):
    is_unpaired_file = False
    for ext in extensions:
        single_fastq_file = os.path.join(sra_stat['getfastq_sra_dir'], sra_stat['sra_id'] + ext)
        if os.path.exists(single_fastq_file):
            is_unpaired_file = True
            break
    return is_unpaired_file

def detect_layout_from_file(sra_stat):
    # Order is important in this list. More downstream should come first.
    extensions = ['.amalgkit.fastq.gz.safely_removed','.amalgkit.fastq.gz','.rename.fastq.gz','.fastp.fastq.gz','.fastq.gz']
    is_paired_end = False
    for ext in extensions:
        paired_fastq_files = [
            os.path.join(sra_stat['getfastq_sra_dir'], sra_stat['sra_id'] + '_1'+ext),
            os.path.join(sra_stat['getfastq_sra_dir'], sra_stat['sra_id'] + '_2'+ext),
        ]
        if all([os.path.exists(f) for f in paired_fastq_files]):
            is_paired_end = True
            break
    is_unpaired_file = is_there_unpaired_file(sra_stat, extensions)
    if (not is_paired_end) & is_unpaired_file:
        is_single_end = True
    else:
        is_single_end = False
    if is_single_end & (sra_stat['layout'] == 'paired'):
        txt = 'Single-end fastq was generated even though layout in the metadata = {}. '
        txt += 'This sample will be treated as single-end reads: {}\n'
        txt = txt.format(sra_stat['layout'], sra_stat['sra_id'])
        sys.stderr.write(txt)
        sra_stat['layout'] = 'single'
    if is_paired_end & (sra_stat['layout'] == 'single'):
        txt = 'Paired-end fastq was generated even though layout in the metadata = {}. '
        txt += 'This sample will be treated as paired-end reads: {}\n'
        txt = txt.format(sra_stat['layout'], sra_stat['sra_id'])
        sys.stderr.write(txt)
        sra_stat['layout'] = 'paired'
    return sra_stat

def write_updated_metadata(metadata, outpath, args):
    if os.path.exists(outpath):
        print('Updated metadata file was detected. Will be overwritten: {}'.format(outpath), flush=True)
    quant_dir = os.path.join(args.out_dir, 'quant')
    metadata = get_mapping_rate(metadata, quant_dir)
    print('Writing curate metadata containing mapping rate: {}'.format(outpath))
    metadata.df.to_csv(outpath, sep='\t', index=False)

def get_mapping_rate(metadata, quant_dir):
    if os.path.exists(quant_dir):
        print('quant directory found: {}'.format(quant_dir))
        metadata.df.loc[:, 'mapping_rate'] = numpy.nan
        sra_ids = metadata.df.loc[:, 'run'].values
        sra_dirs = [d for d in os.listdir(quant_dir) if d in sra_ids]
        print('Number of quant sub-directories that matched to metadata: {:,}'.format(len(sra_dirs)))
        for sra_id in sra_dirs:
            run_info_path = os.path.join(quant_dir, sra_id, sra_id + '_run_info.json')
            if not os.path.exists(run_info_path):
                sys.stderr.write('run_info.json not found. Skipping {}.\n'.format(sra_id))
                continue
            is_sra = (metadata.df.loc[:, 'run'] == sra_id)
            with open(run_info_path) as f:
                run_info = json.load(f)
            metadata.df.loc[is_sra, 'mapping_rate'] = run_info['p_pseudoaligned']
    else:
        txt = 'quant directory not found. Mapping rate cutoff will not be applied: {}\n'
        sys.stderr.write(txt.format(quant_dir))
    return metadata

def check_rscript():
    try:
        subprocess.run(['Rscript', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as e:
        print(e)
        print("R (Rscript) is not installed. Exiting.")
        sys.exit(1)

def orthogroup2genecount(file_orthogroup, file_genecount, spp):
    df = pandas.read_csv(file_orthogroup, sep='\t', header=0, low_memory=False)
    orthogroup_df = pandas.DataFrame({'orthogroup_id': df['busco_id'].to_numpy()})
    is_spp = df.columns.isin(spp)
    df = df.loc[:,is_spp]
    df[df.isnull()] = ''
    df[df=='-'] = ''
    gc = pandas.DataFrame(0, index=df.index, columns=df.columns)
    no_comma = (df != '') & (~df.apply(lambda x: x.str.contains(',')))
    gc[no_comma] = 1
    has_comma = df.apply(lambda x: x.str.contains(','))
    gc[has_comma] = df[has_comma].apply(lambda x: x.str.count(',') + 1)
    gc = pandas.concat([orthogroup_df, gc], axis=1)
    col_order = ['orthogroup_id'] + [col for col in gc.columns if col != 'orthogroup_id']
    gc = gc[col_order]
    gc.to_csv(file_genecount, index=False, sep='\t')

def check_ortholog_parameter_compatibility(args):
    if (args.orthogroup_table is None)&(args.dir_busco is None):
        raise Exception('One of --orthogroup_table and --dir_busco should be specified.')
    if (args.orthogroup_table is not None)&(args.dir_busco is not None):
        raise Exception('Only one of --orthogroup_table and --dir_busco should be specified.')

def generate_multisp_busco_table(dir_busco, outfile):
    print('Generating multi-species BUSCO table.', flush=True)
    col_names = ['busco_id', 'status', 'sequence', 'score', 'length', 'orthodb_url', 'description']
    species_infiles = [f for f in os.listdir(path=dir_busco) if f.endswith('.tsv')]
    species_infiles = sorted(species_infiles)
    print('BUSCO full tables for {} species were detected at: {}'.format(len(species_infiles), dir_busco), flush=True)
    for species_infile in species_infiles:
        path_to_table = os.path.join(dir_busco, species_infile)
        if not os.path.exists(path_to_table):
            warnings.warn('full_table.tsv does not exist. Skipping: '.format(species_infile))
            continue
        tmp_table = pandas.read_table(path_to_table, sep='\t', header=None, comment='#', names=col_names)
        tmp_table.loc[:, 'sequence'] = tmp_table.loc[:, 'sequence'].str.replace(r':[-\.0-9]*$', '', regex=True)
        for col in ['sequence', 'orthodb_url', 'description']:
            tmp_table[col] = tmp_table[col].fillna('').astype(str)
            tmp_table.loc[(tmp_table[col]==''), col] = '-'
        if species_infile == species_infiles[0]:
            merged_table = tmp_table.loc[:, ['busco_id', 'orthodb_url', 'description']]
            merged_table = merged_table.drop_duplicates(keep='first', inplace=False, ignore_index=True)
        else:
            is_mt_missing = (merged_table.loc[:, 'orthodb_url'] == '-')
            if is_mt_missing.sum() > 0:
                tmp_table2 = tmp_table.loc[:, ['busco_id', 'orthodb_url', 'description']]
                tmp_table2 = tmp_table2.drop_duplicates(keep='first', inplace=False, ignore_index=True)
                merged_table.loc[is_mt_missing, 'orthodb_url'] = tmp_table2.loc[is_mt_missing, 'orthodb_url']
                merged_table.loc[is_mt_missing, 'description'] = tmp_table2.loc[is_mt_missing, 'description']
        tmp_table = tmp_table.loc[:, ['busco_id', 'sequence']].groupby(['busco_id'])['sequence'].apply(
            lambda x: ','.join(x))
        tmp_table = tmp_table.reset_index()
        species_colname = species_infile
        species_colname = re.sub(r'_', 'PLACEHOLDER', species_colname)
        species_colname = re.sub(r'[-\._].*', '',  species_colname)
        species_colname = re.sub(r'PLACEHOLDER', '_', species_colname)
        tmp_table = tmp_table.rename(columns={'sequence': species_colname})
        merged_table = merged_table.merge(tmp_table, on='busco_id', how='outer')
    merged_table.to_csv(outfile, sep='\t', index=None, doublequote=False)

def check_config_dir(dir_path, mode):
    files = os.listdir(dir_path)
    if mode=='select':
        asserted_files = [
            'group_attribute.config',
            'exclude_keyword.config',
            'control_term.config',
        ]
    missing_count = 0
    for af in asserted_files:
        if af in files:
            print('Config file found: {}'.format(af))
        else:
            sys.stderr.write('Config file not found: {}\n'.format(af))
            missing_count += 1
    if (missing_count>0):
        txt = 'Please refer to the AMALGKIT Wiki for more info: https://github.com/kfuku52/amalgkit/wiki/amalgkit-metadata\n'
        sys.stderr.write(txt)

def get_getfastq_run_dir(args, sra_id):
    amalgkit_out_dir = os.path.realpath(args.out_dir)
    run_output_dir = os.path.join(amalgkit_out_dir, 'getfastq', sra_id)
    if not os.path.exists(run_output_dir):
        os.makedirs(run_output_dir)
    return run_output_dir

def check_seqkit_dependency():
    try:
        subprocess.run(['seqkit','-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("SeqKit dependency satisfied. Moving on.")
    except FileNotFoundError:
        raise Exception("SeqKit not found. Please make sure SeqKit is installed properly.")



================================================
FILE: amalgkit/util.r
================================================
get_singlecopy_bool_index = function(df_gc, spp_filled, percent_singlecopy_threshold = 50) {

    is_ge_singlecopy_threshold = function(x, num_sp, percent_singlecopy_threshold) {
        num_singlecopy_species = sum(x == 1)
        percent_singlecopy_species = num_singlecopy_species / num_sp * 100
        is_ge_singlecopy = percent_singlecopy_species >= percent_singlecopy_threshold
        return(is_ge_singlecopy)
    }

    num_sp = length(spp_filled)
    is_singlecopy = apply(df_gc[, spp_filled], 1, function(x) { is_ge_singlecopy_threshold(x, num_sp, percent_singlecopy_threshold) })
    num_sc = sum(is_singlecopy)
    txt = 'Number of single-copy orthogroups (>=%s percent species) detected for the %s species: %s\n'
    cat(sprintf(txt, percent_singlecopy_threshold, formatC(length(spp_filled), big.mark = ','), formatC(num_sc, big.mark = ',')))
    return(is_singlecopy)
}

impute_expression = function(dat, num_pc = 4) {
    is_all_na_row = apply(dat, 1, function(x) { all(is.na(x)) })
    tmp = dat[!is_all_na_row,]
    txt = 'Number of removed rows with all NA values in the expression matrix: %s\n'
    cat(sprintf(txt, formatC(sum(is_all_na_row), big.mark = ',')))
    num_na = sum(is.na(tmp))
    num_sp = ncol(tmp)
    num_gene = nrow(tmp)
    num_all = num_sp * num_gene
    txt = 'Imputing %s missing values in a total of %s observations (%s genes x %s samples).\n'
    cat(sprintf(txt, formatC(num_na, big.mark = ','), formatC(num_all, big.mark = ','), formatC(num_gene, big.mark = ','), formatC(num_sp, big.mark = ',')))
    pc = pcaMethods::pca(tmp, nPcs = num_pc, method = 'ppca')
    imputed_dat = pcaMethods::completeObs(pc)
    num_negative = sum(imputed_dat < 0)
    txt = 'Number of negative values clipped to zero in the imputed expression matrix: %s\n'
    cat(sprintf(txt, formatC(num_negative, big.mark = ',')))
    imputed_dat[imputed_dat < 0] = 0
    return(imputed_dat)
}

write_table_with_index_name = function(df, file_path, index_name = 'target_id', sort = TRUE) {
    df_index = data.frame(placeholder_name = rownames(df), stringsAsFactors = FALSE)
    colnames(df_index) = index_name
    df = cbind(df_index, df)
    if (sort) {
        df = df[order(df[[index_name]]),]
    }
    write.table(df, file = file_path, sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
}

save_exclusion_plot = function(df, out_path, font_size) {
    data_summary = aggregate(cbind(count = exclusion) ~ scientific_name + exclusion, df, length)
    data_summary[['total']] = ave(data_summary[['count']], data_summary[['scientific_name']], FUN = sum)
    data_summary[['proportion']] = data_summary[['count']] / data_summary[['total']]
    g = ggplot(data_summary, aes(x = scientific_name, y = count, fill = exclusion))
    g = g + geom_bar(stat = "identity")
    g = g + labs(x = "", y = "Count", fill = "exclusion")
    g = g + theme_bw(base_size = font_size)
    g = g + theme(
        axis.text = element_text(size = font_size, color = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = font_size, color = 'black'),
        #panel.grid.major.y=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = font_size, color = 'black'),
        legend.position = 'bottom',
        rect = element_rect(fill = "transparent"),
        plot.margin = unit(rep(0.1, 4), "cm")
    )
    num_spp = length(unique(df[['scientific_name']]))
    plot_width = max(3.6, 0.11 * num_spp)
    ggsave(out_path, plot = g, width = plot_width, height = 3.6, units = 'in')
}


================================================
FILE: amalgkit/config_dir/__init__.py
================================================



================================================
FILE: amalgkit/config_dir/base/__init__.py
================================================



================================================
FILE: amalgkit/config_dir/base/control_term.config
================================================
# case-insensitive
# regular expressions allowed
# (attribute)[TAB](control term)

"treatment"	"mock"
"exp_title"	"0 days post infection"



================================================
FILE: amalgkit/config_dir/base/exclude_keyword.config
================================================
# case-insensitive
# regular expressions allowed
# (comma-delimited attributes to be searched)[TAB](arbitrary reason of exclusion)[TAB](bad keyword)

"exp_title,study_title,design,sample_title,sample_description,lib_name,experiment,treatment,protocol,age"	"single_cell"	"single cell"
"treatment,protocol"	"small_RNA"	"miRNA"



================================================
FILE: amalgkit/config_dir/base/group_attribute.config
================================================
# case-insensitive
# regular expressions allowed
# (aggregated to)[TAB](aggregated from)

"age"	"arrayexpress-developmentalstage"
"treatment"	"infection"


================================================
FILE: amalgkit/config_dir/plantae/__init__.py
================================================



================================================
FILE: amalgkit/config_dir/plantae/control_term.config
================================================
# case-insensitive
# regular expressions allowed
# (attribute)[TAB](control term)

"treatment"	"mock"
"treatment"	"control"
"treatment"	"cntl"
"treatment"	"con\d"
"treatment"	"con_"
"treatment"	"cont\b"
"treatment"	"cont_"
"treatment"	"kontrol"
"treatment"	"contral"
"treatment"	"contron"
"treatment"	"\bCK"
"treatment"	"none"
"treatment"	"normal"
"treatment"	"natural"
"treatment"	"standard"
"treatment"	"optimal"
"treatment"	"sufficient"
"treatment"	"no(?:\s+\w+)?\s+treat"
"treatment"	"non-treat"
"treatment"	"non- treat"
"treatment"	"nontreat"
"treatment"	"untreat"
"treatment"	"without"
"treatment"	"^no$"
"treatment"	"Ad Libitum"
"treatment"	"in the absence of"
"treatment"	"Thermoneutral"
"treatment"	"21%[oxygen_concentration]"
"treatment"	"DMSO"
"treatment"	"\b0h"
"treatment"	"\b0 h"
"treatment"	"_0h"
"treatment"	"\b0day"
"treatment"	"\b0 day"
"treatment"	"\b0m"
"treatment"	"\b0 m"
"treatment"	"\b0 dpi"
"treatment"	"^0 ppm"
"treatment"	"wild"
"treatment"	"uninfected"
"treatment"	"un-infected"
"treatment"	"non infected"
"treatment"	"non-inifected"
"treatment"	"no inoculation"
"treatment"	"un-inoculated"
"treatment"	"no disease"
"treatment"	"non-diseased"
"treatment"	"unwound"
"treatment"	"un-grazing"
"treatment"	"healthy"
"treatment"	"stress-free"
"treatment"	"nostress"
"treatment"	"non-stress"
"treatment"	"unstressed"
"treatment"	"non(?:\s+\w+)?\s+stress"
"treatment"	"room temperature"
"treatment"	"ambient"
"treatment"	"water\w*\b(?!\s*deficit)"
"treatment"	"H2O\w*\b(?!\s*deficit)"

"genotype"	"wild"
"genotype"	"N2"
"genotype"	"OregonR"
"genotype"	"Oregon R"
"genotype"	"Canton S"
"genotype"	"Canton-S"
"genotype"	"y w"
"genotype"	"yw"
"genotype"	"wt"
"genotype"	"\+\/\+"
"genotype"	"columbia"
"genotype"	"col0"
"genotype"	"col-0"

"sample_title"	"wild"
"sample_title"	"\bWT"
"sample_title"	"control"
"sample_title"	"CON_"
"sample_title"	"CK_"
"sample_title"	"CK-"
"sample_title"	"CK\d"
"sample_title"	"\bCK\b"
"sample_title"	"normal"
"sample_title"	"non-treat"
"sample_title"	"non-inoculated"
"sample_title"	"ungrazed"
"sample_title"	"healthy"
"sample_title"	"water"
"sample_title"	"H2O"
"sample_title"	"\b0dpi"

"sample_description"	"control"
"sample_description"	"mock"
"sample_description"	"no treatment"
"sample_description"	"non-infected"
"sample_description"	"before"
"sample_description"	"wild"

"lib_name"	"control"
"lib_name"	"CK_"
"lib_name"	"CK-"
"lib_name"	"CK\d"
"lib_name"	"\bCK\b"
"lib_name"	"\bcon-"
"lib_name"	"ctrl"
"lib_name"	"mock"
"lib_name"	"wild"
"lib_name"	"\bWT"
"lib_name"	"(?<!\d)0h"

"exp_title"	"control"
"exp_title"	"CON_"
"exp_title"	"CK_"
"exp_title"	"CK-"
"exp_title"	"CK\d"
"exp_title"	"\bCK\b"
"exp_title"	"untreated"
"exp_title"	"mock"
"exp_title"	"without"
"exp_title"	"water"
"exp_title"	"H2O"
"exp_title"	"TimePoint1_"
"exp_title"	"0 days post infection"
"exp_title"	"\b0h"
"exp_title"	"\b0 h"
"exp_title"	"_0h"

"design"	"control"
"design"	"CK_"
"design"	"CK-"
"design"	"CK\d"
"design"	"\bCK\b"
"design"	"\bcon-"
"design"	"mock"
"design"	"normal"
"design"	"common"
"design"	"untreat"
"design"	"non-treat"
"design"	"no(?:\s+\w+)?\s+treat"
"design"	"no processing"
"design"	"without"
"design"	"water"
"design"	"H2O"
"design"	"before"
"design"	"\b0h"
"design"	"\b0 h"
"design"	"_0h"
"design"	"\b0d"
"design"	"\b0 d"
"design"	"\b0 m"
"design"	"healthy"
"design"	"ungrazed"
"design"	"no grazing"
"design"	"no-grazing"
"design"	"not inoculated"
"design"	"non-inoculated"
"design"	"uninfected"
"design"	"no incubation"
"design"	"not primed"
"design"	"not induced"
"design"	"wt"
"design"	"wild"




================================================
FILE: amalgkit/config_dir/plantae/exclude_keyword.config
================================================
# case-insensitive
# regular expressions allowed
# (comma-delimited attributes to be searched)[TAB](arbitrary reason of exclusion)[TAB](bad keyword)

"exp_title,study_title,design,sample_title,sample_description,lib_name,experiment,treatment,protocol,age"	"single_cell"	"single cell"

"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"immunoprecipitation"	"RipSeq"
"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"immunoprecipitation"	"chrom_RNAseq"
"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"RNAi"	"RNAi"
"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"small_RNA"	"shRNA RNA-seq"
"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"small_RNA"	"piRNA"
"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"small_RNA"	"smRNA"
"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"CAGE"	"CAGE"

"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol,treatment"	"treatment"	"^exposed to"

"treatment,protocol"	"small_RNA"	"miRNA"
"treatment,protocol"	"small_RNA"	"sRNA"
"treatment,protocol"	"treatment"	"uranium"

"age"	"embryonic"	"embryo"
"age"	"embryonic"	"somite"
"age"	"embryonic"	"fertiliz.*"
"age"	"embryonic"	"oocyte"

"genotype"	"transgenic"	"foxn3"



================================================
FILE: amalgkit/config_dir/plantae/group_attribute.config
================================================
# case-insensitive
# regular expressions allowed
# (aggregated to)[TAB](aggregated from)

"age"	"arrayexpress-developmentalstage"
"age"	"agedays"
"age"	"age_classification"
"age"	"age_description"
"age"	"age_in_years"
"age"	"agecat"
"age"	"dev_stage"
"age"	"ages"
"age"	"development_satge"
"age"	"development_stage"
"age"	"developmental_stage"
"age"	"developmentalstage"
"age"	"differentiation_stages"
"age"	"embryonic_day"
"age"	"embryonic_stage"
"age"	"female_age"
"age"	"sample_comment"
"age"	"stage"
"age"	"age_at_collection"
"age"	"age_classification"
"age"	"age_description"
"age"	"age_in_years"
"age"	"age_of_collection"
"age"	"age_of_fly_in_days_post_eclosion"
"age"	"agecat"
"age"	"agedays"
"age"	"age{category}"
"age"	"animal_age_at_collection"
"age"	"day"
"age"	"arrayexpress-timepoint"
"age"	"avg_age_after_second_6_months_of_open_access"
"age"	"collection_timing_in_reference_to_e_lineage"
"age"	"day_relative_to_weaning"
"age"	"dev-stage"
"age"	"developmental_time"
"age"	"duration"
"age"	"embryo_stage"
"age"	"female_birth_date"
"age"	"gestational_age"
"age"	"gestational_age_at_sample_collection"
"age"	"hotta_stage"
"age"	"male_age"
"age"	"male_birth_date"
"age"	"male_collection_date"
"age"	"male_death_date"
"age"	"timepoint"
"age"	"trimester"
"age"	"birth_date"
"age"	"collection_date"
"age"	"death_date"
"age"	"female_collection_date"
"age"	"female_death_date"
"age"	"geo_loc_name"
"age"	"sample_date"
"age"	"sampling_date"
"age"	"specimen_collection_date"
"age"	"time_-_blastomeres"
"age"	"animal_length"
"age"	"date"
"age"	"developemental_stage"
"age"	"colection_date"
"age"	"time"

"antibody"	"chip_antibody"
"antibody"	"chip_or_ip_antibody"
"antibody"	"clip-antibody"
"antibody"	"clip_antibody"
"antibody"	"rip_antibody"
"antibody"	"rna_binding_protein"

"batch"	"adapter_barcode"
"batch"	"animal"
"batch"	"animal_id"
"batch"	"animalid"
"batch"	"barcode"
"batch"	"biological_replicate"
"batch"	"biorep"
"batch"	"brain_code"
"batch"	"clutch"
"batch"	"cohort"
"batch"	"collected_by"
"batch"	"cow"
"batch"	"custom_field"
"batch"	"donator_monkey_number"
"batch"	"donor_id"
"batch"	"horse_number"
"batch"	"iclip_barcode"
"batch"	"id"
"batch"	"illumina_rna-seq_bar_codes"
"batch"	"index"
"batch"	"individual"
"batch"	"individual_id"
"batch"	"individuals"
"batch"	"internal_id"
"batch"	"labanimalnumber"
"batch"	"labexpid"
"batch"	"mahpic_non_human_primate_individual_id"
"batch"	"monkey_id"
"batch"	"non_human_primate_individual_id"
"batch"	"number"
"batch"	"oocyte/embryo_number"
"batch"	"pig_id"
"batch"	"replicate"
"batch"	"rna_processing_batch"
"batch"	"run_id"
"batch"	"sample"
"batch"	"sample-type"
"batch"	"sample_coding"
"batch"	"sample_descriptor"
"batch"	"sample_id"
"batch"	"sample_identifier"
"batch"	"sample_number"
"batch"	"series"
"batch"	"specimen_voucher"
"batch"	"study"
"batch"	"study_group"
"batch"	"study_phase"
"batch"	"subject"
"batch"	"subject_id"
"batch"	"subjectid"
"batch"	"technical_replicate"
"batch"	"tissue_abbreviation"
"batch"	"tissue_code"
"batch"	"umc_id"
"batch"	"unique_id"
"batch"	"uniqueid"
"batch"	"wur_marker_type"
"batch"	"animal_number"
"batch"	"ebi_equivalent_biosample"
"batch"	"ercc"
"batch"	"flow_cell"
"batch"	"flowcell"
"batch"	"lane"
"batch"	"lane_num"
"batch"	"num_replicons"
"batch"	"cell_barcode"
"batch"	"submission_identifier"
"batch"	"library_id"
"batch"	"library_index_sequence_used_to_demultiplex"
"batch"	"replicate_/_experiment"
"batch"	"sequencing_pool"
"batch"	"plate_col_id"
"batch"	"plate_row_id"
"batch"	"replicate_no."
"batch"	"run_no."
"batch"	"sire"
"batch"	"unique_identifier"

"biomaterial_provider"	"biosourceprovider"
"biomaterial_provider"	"female_biomaterial_provider"
"biomaterial_provider"	"male_biomaterial_provider"
"biomaterial_provider"	"library_preparation_location"
"biomaterial_provider"	"material_provider"
"biomaterial_provider"	"availability"

"bioproject"	"bioproject_id"
"bioproject"	"bioprojectid"

"cell"	"cell_type"
"cell"	"cell_line"
"cell"	"cell_class"
"cell"	"cell_description"
"cell"	"cell_subtype"
"cell"	"arrayexpress-celltype"
"cell"	"culture_collection"
"cell"	"culture_conditions"
"cell"	"cell-line"
"cell"	"cell-type"
"cell"	"cell_or_tissue_type"
"cell"	"eosinophil"
"cell"	"feature"
"cell"	"germ_layer"
"cell"	"cells_derived_from"
"cell"	"number_cells"
"cell"	"cell_organism"
"cell"	"cell_culture_protocol"
"cell"	"cell_population"
"cell"	"cell_typing"

"protocol"	"extract_protocol"
"protocol"	"female_sample_collection_protocol"
"protocol"	"protocol"
"protocol"	"sample_extraction_method"
"protocol"	"experimental_protocol"
"protocol"	"chemical_treatment_of_isolated_rna"
"protocol"	"amplification"
"protocol"	"assay"
"protocol"	"experiment_target"
"protocol"	"extraction_protocol"
"protocol"	"fixation"
"protocol"	"fraction"
"protocol"	"library_selection"
"protocol"	"library_source"
"protocol"	"library_strategy"
"protocol"	"library_type"
"protocol"	"libraryprotocol"
"protocol"	"lysis_buffer_ion_concentrations"
"protocol"	"lysis_strategy"
"protocol"	"male_sample_collection_protocol"
"protocol"	"meoh-fixed"
"protocol"	"molecule_subtype"
"protocol"	"prep_type"
"protocol"	"preperation_kit"
"protocol"	"purification_protocol"
"protocol"	"rin"
"protocol"	"rna_fraction"
"protocol"	"rna_rin_values_after_globin_depletion"
"protocol"	"rna_rin_values_before_globin_depletion"
"protocol"	"rna_subtype"
"protocol"	"sample_material"
"protocol"	"sample_storage"
"protocol"	"sample_storage_processing"
"protocol"	"small_rna_classes"
"protocol"	"specimen_collection_protocol"
"protocol"	"specimen_with_known_storage_state"
"protocol"	"tissue_state"
"protocol"	"store_cond"
"protocol"	"datatype_description"
"protocol"	"datatype"
"protocol"	"monosome_enrichment_strategy"
"protocol"	"rna_concentration_after_globin_depletion"
"protocol"	"rna_concentration_before_globin_depletion"
"protocol"	"rna_input"
"protocol"	"sequencing_type"
"protocol"	"date_run"
"protocol"	"instrument_model"
"protocol"	"ngs_platform"
"protocol"	"paired_end_seq?"
"protocol"	"paired_or_single-end"
"protocol"	"pe_read_length"
"protocol"	"platform"
"protocol"	"quality"
"protocol"	"route"
"protocol"	"read_len_orig"
"protocol"	"readtype"
"protocol"	"readtype_description"
"protocol"	"sequencer"
"protocol"	"sequencing_method"
"protocol"	"extraction"
"protocol"	"indrops_version"
"protocol"	"minimum_counts_per_cell_threshold_used_to_remove_background_barcodes"
"protocol"	"minimum_reads_per_cell_during_initial_processing"
"protocol"	"spike-in"

"nominal_length"	"mean_insert_size"

"genotype"	"arrayexpress-genotype"
"genotype"	"arrayexpress-phenotype"
"genotype"	"arrayexpress-rnai"
"genotype"	"background_strain"
"genotype"	"breed"
"genotype"	"cultivar"
"genotype"	"chicken_line"
"genotype"	"ecotype"
"genotype"	"full_genotype"
"genotype"	"genetic_background"
"genotype"	"genetic_line"
"genotype"	"genetic_modification"
"genotype"	"genotype/variaion"
"genotype"	"genotype/variation"
"genotype"	"germline_knock-down"
"genotype"	"germline_knockdown_and_other_transgenes"
"genotype"	"isolate"
"genotype"	"line"
"genotype"	"phenotype"
"genotype"	"ploidy"
"genotype"	"population"
"genotype"	"snp_allelic_state"
"genotype"	"strain"
"genotype"	"strain/background"
"genotype"	"strain_origin"
"genotype"	"strains"
"genotype"	"sub_species"
"genotype"	"subspecies"
"genotype"	"variety"
"genotype"	"breed_name"
"genotype"	"arrayexpress-strainorline"
"genotype"	"breed/line"
"genotype"	"breed_type"
"genotype"	"breeding_history"
"genotype"	"breeds"
"genotype"	"chicken_breed"
"genotype"	"clonality"
"genotype"	"marker"
"genotype"	"dstim_knockdown_status"
"genotype"	"environmental_history"
"genotype"	"fly_line"
"genotype"	"origin"
"genotype"	"phenotype/variation"
"genotype"	"pig_breed"
"genotype"	"pig_type"
"genotype"	"propagation"
"genotype"	"reporter_construct"
"genotype"	"rhesus_monkey_origin"
"genotype"	"rna_interference"
"genotype"	"strain_background"
"genotype"	"strain_description"
"genotype"	"strain_prrsv"
"genotype"	"strainorline"
"genotype"	"transfecting_vector"
"genotype"	"transgenic_strain"
"genotype"	"virus_group"
"genotype"	"arrayexpress-species"
"genotype"	"parental_strain"
"genotype"	"tax_id"
"genotype"	"transgene"
"genotype"	"genotype/phenotype"
"genotype"	"cross"
"genotype"	"dad"
"genotype"	"driver"
"genotype"	"qtl_genotype/haplotype"
"genotype"	"wild type genotype"

"location"	"birth_location"
"location"	"country"
"location"	"feeding_barn"
"location"	"female_birth_location"
"location"	"geographic_location"
"location"	"lat_lon"
"location"	"library_preparation_location_latitude"
"location"	"library_preparation_location_latitude_units"
"location"	"library_preparation_location_longitude"
"location"	"library_preparation_location_longitude_units"
"location"	"local_origin"
"location"	"male_birth_location"
"location"	"sampling_site"
"location"	"seq_center"
"location"	"sequencing_location"
"location"	"sequencing_location_latitude"
"location"	"station"
"location"	"sequencing_location_latitude_units"
"location"	"sequencing_location_longitude"
"location"	"sequencing_location_longitude_units"
"location"	"source"
"location"	"biome"
"location"	"env_biome"
"location"	"env_feature"
"location"	"env_material"
"location"	"environment"
"location"	"environmental_package"
"location"	"sampling_site"

"misc"	"ena-checklist"
"misc"	"ena-first-public"
"misc"	"ena-last-update"
"misc"	"estimated_size"
"misc"	"investigation_type"
"misc"	"lab_description"
"misc"	"library_layout"
"misc"	"mapalgorithm"
"misc"	"mapalgorithm_description"
"misc"	"new_attribute"
"misc"	"notes"
"misc"	"organism"
"misc"	"project"
"misc"	"project_name"
"misc"	"species"
"misc"	"sra_experiment_accession"
"misc"	"biosamplemodel"
"misc"	"sample_year"
"misc"	"number_of_pieces"
"misc"	"specimen_picture_url"
"misc"	"submission_description"
"misc"	"submission_title"
"misc"	"unknown"
"misc"	"body_weight"
"misc"	"file_no."
"misc"	"num_parts_in_pool"
"misc"	"processing_date"
"misc"	"Broker_name"
"misc"	"ENA_checklist"
"misc"	"INSDC_center_name"
"misc"	"INSDC_first_public"
"misc"	"INSDC_last_update"
"misc"	"INSDC_status"
"misc"	"SRA_accession"
"misc"	"Alias"



"sex"	"arrayexpress-sex"
"sex"	"gender"
"sex"	"host_sex"
"sex"	"cell_sex"
"sex"	"gender_type"
"sex"	"obsolete_sex"
"sex"	"both_sexes_pooled"

"tissue"	"arrayexpress-organismpart"
"tissue"	"bio_material"
"tissue"	"body_site"
"tissue"	"biopsy_site"
"tissue"	"brain_region"
"tissue"	"description"
"tissue"	"experiment_set"
"tissue"	"explant"
"tissue"	"isolation-source"
"tissue"	"isolation_source"
"tissue"	"label"
"tissue"	"mixedtissues"
"tissue"	"muscle_type"
"tissue"	"organ/tissue"
"tissue"	"organ_part"
"tissue"	"organism_part"
"tissue"	"organismpart"
"tissue"	"oviduct"
"tissue"	"region"
"tissue"	"sample_name"
"tissue"	"sample_origin"
"tissue"	"sampling_position"
"tissue"	"source_name"
"tissue"	"tag"
"tissue"	"tissue-type"
"tissue"	"tissue/cell"
"tissue"	"tissue/cell_type"
"tissue"	"tissue_location"
"tissue"	"tissue_or_dev_stage"
"tissue"	"tissue_origin"
"tissue"	"tissue_source"
"tissue"	"tissue_type"
"tissue"	"title"
"tissue"	"arrayexpress-organismpart"
"tissue"	"organism_part"
"tissue"	"embryo_region"
"tissue"	"organsim_part"
"tissue"	"tissue_/_cells"
"tissue"	"tissue_lib"
"tissue"	"tissue_region"

"treatment"	"agent"
"treatment"	"arrayexpress-diseasestate"
"treatment"	"arrayexpress-growthcondition"
"treatment"	"arrayexpress-immunoprecipitate"
"treatment"	"behavior"
"treatment"	"biopsy_day"
"treatment"	"challenged_with"
"treatment"	"chemical_treatment_of_isolated_rna"
"treatment"	"diet"
"treatment"	"diet_fed"
"treatment"	"dietary_group"
"treatment"	"disease"
"treatment"	"disease_status"
"treatment"	"domestication"
"treatment"	"drip"
"treatment"	"energy_balance"
"treatment"	"fed_status"
"treatment"	"feeding"
"treatment"	"feeding_type"
"treatment"	"fertility_group"
"treatment"	"food"
"treatment"	"group"
"treatment"	"health_state"
"treatment"	"infect"
"treatment"	"infected_with"
"treatment"	"infection"
"treatment"	"infection_status"
"treatment"	"inoculation"
"treatment"	"maternal_diet"
"treatment"	"meat_quality"
"treatment"	"mptp_treatment"
"treatment"	"oxidation"
"treatment"	"pathogen"
"treatment"	"percent_marinade_uptake_at_24h_post-slaughter"
"treatment"	"phenotype_sample_type"
"treatment"	"restriction_model"
"treatment"	"altitude"
"treatment"	"dpi"
"treatment"	"salmonella_shedding_status"
"treatment"	"status"
"treatment"	"stimulation"
"treatment"	"stimulus"
"treatment"	"stress"
"treatment"	"survival"
"treatment"	"transplant_type"
"treatment"	"treated_with"
"treatment"	"treatment_group"
"treatment"	"pregnancy_outcome"
"treatment"	"isol_growth_condt"
"treatment"	"agent"
"treatment"	"aliquot"
"treatment"	"ammonia_concentration"
"treatment"	"ammonia_exposure"
"treatment"	"average_altitude"
"treatment"	"blastocyst_rate"
"treatment"	"breeding_method"
"treatment"	"bull_ntm"
"treatment"	"calcium_intake"
"treatment"	"cell_status"
"treatment"	"cellular_compartment"
"treatment"	"chemical"
"treatment"	"collection_time"
"treatment"	"comment"
"treatment"	"concentration"
"treatment"	"condition"
"treatment"	"culture_condition"
"treatment"	"culture-collection"
"treatment"	"culture_type"
"treatment"	"custom_name"
"treatment"	"day_of_infection"
"treatment"	"day_post_infection"
"treatment"	"days_at_29_c"
"treatment"	"days_post_infection"
"treatment"	"differentiation_time"
"treatment"	"disease_stage"
"treatment"	"drinking_category"
"treatment"	"drug"
"treatment"	"drug_for_synchronization"
"treatment"	"egg_production"
"treatment"	"environmental_stress"
"treatment"	"enzyme_treatment"
"treatment"	"exposure_time"
"treatment"	"extracellular_component"
"treatment"	"fasted_status"
"treatment"	"fed_for"
"treatment"	"feed_efficiency"
"treatment"	"feeding_period"
"treatment"	"feeding_treatment"
"treatment"	"growth_condition"
"treatment"	"growth_protocol"
"treatment"	"harmonized_disease_score"
"treatment"	"harvest_time"
"treatment"	"health-state"
"treatment"	"health_status"
"treatment"	"health_status_at_collection"
"treatment"	"heat_hours"
"treatment"	"hematocrit"
"treatment"	"hemoglobin_measurement"
"treatment"	"host"
"treatment"	"host_tissue_sampled"
"treatment"	"hours_at_restrictive_temperature"
"treatment"	"hp-prrsv_infection"
"treatment"	"infection_route"
"treatment"	"infectious_agent"
"treatment"	"inflammatory_stimulus_treatment"
"treatment"	"injury"
"treatment"	"kinetics"
"treatment"	"lens_type"
"treatment"	"lesion_score"
"treatment"	"lesion_number"
"treatment"	"library_name"
"treatment"	"litter"
"treatment"	"loin_weight_group"
"treatment"	"lps_exposure"
"treatment"	"material"
"treatment"	"maternal_treatment"
"treatment"	"mating_status"
"treatment"	"mock"
"treatment"	"mode"
"treatment"	"morpholino"
"treatment"	"muscle_name"
"treatment"	"neutralization_group"
"treatment"	"normalized?"
"treatment"	"number_born_alive"
"treatment"	"oxygen_concentration"
"treatment"	"parity_number"
"treatment"	"passage"
"treatment"	"passage_number"
"treatment"	"passagers"
"treatment"	"passages"
"treatment"	"phase"
"treatment"	"pooled"
"treatment"	"pregnancy_status"
"treatment"	"pregnancy_time"
"treatment"	"purity"
"treatment"	"resorted"
"treatment"	"rna_treatment"
"treatment"	"sample_processing"
"treatment"	"sampleid"
"treatment"	"sampling_time_point"
"treatment"	"sampling_to_preparation_interval"
"treatment"	"selection_criteria"
"treatment"	"serologic_response_status"
"treatment"	"shear_force"
"treatment"	"skin_color"
"treatment"	"social_rank"
"treatment"	"sorted"
"treatment"	"specific_host"
"treatment"	"stim_time"
"treatment"	"stud_book_number"
"treatment"	"superovulated"
"treatment"	"temperature"
"treatment"	"thermal_treatment"
"treatment"	"time"
"treatment"	"time-point"
"treatment"	"time_point"
"treatment"	"ega-disease"
"treatment"	"time_post_infection"
"treatment"	"timepoint_postchallenge"
"treatment"	"training"
"treatment"	"training_timepoint"
"treatment"	"treatment_description"
"treatment"	"tumor"
"treatment"	"tumor_grading"
"treatment"	"vaccine"
"treatment"	"vaccine_administration"
"treatment"	"vaccine_group"
"treatment"	"viremia_level"
"treatment"	"virus"
"treatment"	"virus_dose"
"treatment"	"d"
"treatment"	"days_post-infection"
"treatment"	"ega-subject"
"treatment"	"ega-treatment"
"treatment"	"ercc_pool"
"treatment"	"host_taxid"
"treatment"	"lymphocyte"
"treatment"	"monocyte"
"treatment"	"morphology"
"treatment"	"neutrophil"
"treatment"	"platelet_count"
"treatment"	"rfi_value"
"treatment"	"sample_type"
"treatment"	"source_of_female"
"treatment"	"sperm_dna_fragmentation_index"
"treatment"	"srek_ligation"
"treatment"	"time/hpf"
"treatment"	"total_number_born"
"treatment"	"treatments"
"treatment"	"type"
"treatment"	"weight"
"treatment"	"arrayexpress-compound"
"treatment"	"arrayexpress-dose"
"treatment"	"assay_type"
"treatment"	"backfat_thickness_of_live(mm)"
"treatment"	"basophil"
"treatment"	"photorecptors_remaining"
"treatment"	"productivity"
"treatment"	"hours/days_after_injury"
"treatment"	"disease_state"
"treatment"	"mehg_exposure"
"treatment"	"operation"
"treatment"	"post-injury_time_point"
"treatment"	"embryo_source"
"treatment"	"earlobe_color"
"treatment"	"reproductive_status"
"treatment"	"specimen_size"
"treatment"	"specimen_volume"
"treatment"	"specimen_weight"
"treatment"	"liver_phenotype"
"treatment"	"clinical_information"
"treatment"	"oestrous_peroid"
"treatment"	"secondary_description"
"treatment"	"finishing_diet"
"treatment"	"gestation_duration"
"treatment"	"growth_rate"
"treatment"	"heart_weight"
"treatment"	"information"
"treatment"	"light_treatment"
"treatment"	"maternal_food"
"treatment"	"number_of_passages"
"treatment"	"paternal_experience"
"treatment"	"personal_experience"
"treatment"	"physiological_conditions"
"treatment"	"rfi_group"
"treatment"	"sample_role"
"treatment"	"slaughter_weight"
"treatment"	"target"
"treatment"	"target_source"
"treatment"	"temperature_regimen"
"treatment"	"compound"



================================================
FILE: amalgkit/config_dir/test/__init__.py
================================================



================================================
FILE: amalgkit/config_dir/test/control_term.config
================================================
# case-insensitive
# regular expressions allowed
# (attribute)[TAB](control term)

"treatment"	"mock"
"treatment"	"control"
"treatment"	"cntl"
"treatment"	"con\d"
"treatment"	"con_"
"treatment"	"cont\b"
"treatment"	"cont_"
"treatment"	"kontrol"
"treatment"	"contral"
"treatment"	"contron"
"treatment"	"\bCK"
"treatment"	"none"
"treatment"	"normal"
"treatment"	"natural"
"treatment"	"standard"
"treatment"	"optimal"
"treatment"	"sufficient"
"treatment"	"no(?:\s+\w+)?\s+treat"
"treatment"	"non-treat"
"treatment"	"non- treat"
"treatment"	"nontreat"
"treatment"	"untreat"
"treatment"	"without"
"treatment"	"^no$"
"treatment"	"Ad Libitum"
"treatment"	"in the absence of"
"treatment"	"Thermoneutral"
"treatment"	"21%[oxygen_concentration]"
"treatment"	"DMSO"
"treatment"	"\b0h"
"treatment"	"\b0 h"
"treatment"	"_0h"
"treatment"	"\b0day"
"treatment"	"\b0 day"
"treatment"	"\b0m"
"treatment"	"\b0 m"
"treatment"	"\b0 dpi"
"treatment"	"^0 ppm"
"treatment"	"wild"
"treatment"	"uninfected"
"treatment"	"un-infected"
"treatment"	"non infected"
"treatment"	"non-inifected"
"treatment"	"no inoculation"
"treatment"	"un-inoculated"
"treatment"	"no disease"
"treatment"	"non-diseased"
"treatment"	"unwound"
"treatment"	"un-grazing"
"treatment"	"healthy"
"treatment"	"stress-free"
"treatment"	"nostress"
"treatment"	"non-stress"
"treatment"	"unstressed"
"treatment"	"non(?:\s+\w+)?\s+stress"
"treatment"	"room temperature"
"treatment"	"ambient"
"treatment"	"water\w*\b(?!\s*deficit)"
"treatment"	"H2O\w*\b(?!\s*deficit)"

"genotype"	"wild"
"genotype"	"N2"
"genotype"	"OregonR"
"genotype"	"Oregon R"
"genotype"	"Canton S"
"genotype"	"Canton-S"
"genotype"	"y w"
"genotype"	"yw"
"genotype"	"wt"
"genotype"	"\+\/\+"

"sample_title"	"wild"
"sample_title"	"\bWT"
"sample_title"	"control"
"sample_title"	"CON_"
"sample_title"	"CK_"
"sample_title"	"CK-"
"sample_title"	"CK\d"
"sample_title"	"\bCK\b"
"sample_title"	"normal"
"sample_title"	"non-treat"
"sample_title"	"non-inoculated"
"sample_title"	"ungrazed"
"sample_title"	"healthy"
"sample_title"	"water"
"sample_title"	"H2O"
"sample_title"	"\b0dpi"

"sample_description"	"control"
"sample_description"	"mock"
"sample_description"	"no treatment"
"sample_description"	"non-infected"
"sample_description"	"before"
"sample_description"	"wild"

"lib_name"	"control"
"lib_name"	"CK_"
"lib_name"	"CK-"
"lib_name"	"CK\d"
"lib_name"	"\bCK\b"
"lib_name"	"\bcon-"
"lib_name"	"ctrl"
"lib_name"	"mock"
"lib_name"	"wild"
"lib_name"	"\bWT"
"lib_name"	"(?<!\d)0h"

"exp_title"	"control"
"exp_title"	"CON_"
"exp_title"	"CK_"
"exp_title"	"CK-"
"exp_title"	"CK\d"
"exp_title"	"\bCK\b"
"exp_title"	"untreated"
"exp_title"	"mock"
"exp_title"	"without"
"exp_title"	"water"
"exp_title"	"H2O"
"exp_title"	"TimePoint1_"
"exp_title"	"0 days post infection"
"exp_title"	"\b0h"
"exp_title"	"\b0 h"
"exp_title"	"_0h"

"design"	"control"
"design"	"CK_"
"design"	"CK-"
"design"	"CK\d"
"design"	"\bCK\b"
"design"	"\bcon-"
"design"	"mock"
"design"	"normal"
"design"	"common"
"design"	"untreat"
"design"	"non-treat"
"design"	"no(?:\s+\w+)?\s+treat"
"design"	"no processing"
"design"	"without"
"design"	"water"
"design"	"H2O"
"design"	"before"
"design"	"\b0h"
"design"	"\b0 h"
"design"	"_0h"
"design"	"\b0d"
"design"	"\b0 d"
"design"	"\b0 m"
"design"	"healthy"
"design"	"ungrazed"
"design"	"no grazing"
"design"	"no-grazing"
"design"	"not inoculated"
"design"	"non-inoculated"
"design"	"uninfected"
"design"	"no incubation"
"design"	"not primed"
"design"	"not induced"
"design"	"wt"
"design"	"wild"



================================================
FILE: amalgkit/config_dir/test/exclude_keyword.config
================================================
# case-insensitive
# regular expressions allowed
# (comma-delimited attributes to be searched)[TAB](arbitrary reason of exclusion)[TAB](bad keyword)

"exp_title,study_title,design,sample_title,sample_description,lib_name,experiment,treatment,protocol,age"	"single_cell"	"single cell"

"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"immunoprecipitation"	"RipSeq"
"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"immunoprecipitation"	"chrom_RNAseq"
"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"RNAi"	"RNAi"
"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"small_RNA"	"shRNA RNA-seq"
"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"small_RNA"	"piRNA"
"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"small_RNA"	"smRNA"
"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"CAGE"	"CAGE"

"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol,treatment"	"treatment"	"^exposed to"

"treatment,protocol"	"small_RNA"	"miRNA"
"treatment,protocol"	"small_RNA"	"sRNA"
"treatment,protocol"	"treatment"	"uranium"

"age"	"embryonic"	"embryo"
"age"	"embryonic"	"somite"
"age"	"embryonic"	"fertiliz.*"
"age"	"embryonic"	"oocyte"

"genotype"	"transgenic"	"foxn3"



================================================
FILE: amalgkit/config_dir/test/group_attribute.config
================================================
# case-insensitive
# regular expressions allowed
# (aggregated to)[TAB](aggregated from)

"age"	"arrayexpress-developmentalstage"
"age"	"agedays"
"age"	"age_classification"
"age"	"age_description"
"age"	"age_in_years"
"age"	"agecat"
"age"	"dev_stage"
"age"	"ages"
"age"	"development_satge"
"age"	"development_stage"
"age"	"developmental_stage"
"age"	"developmentalstage"
"age"	"differentiation_stages"
"age"	"embryonic_day"
"age"	"embryonic_stage"
"age"	"female_age"
"age"	"sample_comment"
"age"	"stage"
"age"	"age_at_collection"
"age"	"age_classification"
"age"	"age_description"
"age"	"age_in_years"
"age"	"age_of_collection"
"age"	"age_of_fly_in_days_post_eclosion"
"age"	"agecat"
"age"	"agedays"
"age"	"age{category}"
"age"	"animal_age_at_collection"
"age"	"day"
"age"	"arrayexpress-timepoint"
"age"	"avg_age_after_second_6_months_of_open_access"
"age"	"collection_timing_in_reference_to_e_lineage"
"age"	"day_relative_to_weaning"
"age"	"dev-stage"
"age"	"developmental_time"
"age"	"duration"
"age"	"embryo_stage"
"age"	"female_birth_date"
"age"	"gestational_age"
"age"	"gestational_age_at_sample_collection"
"age"	"hotta_stage"
"age"	"male_age"
"age"	"male_birth_date"
"age"	"male_collection_date"
"age"	"male_death_date"
"age"	"timepoint"
"age"	"trimester"
"age"	"birth_date"
"age"	"collection_date"
"age"	"death_date"
"age"	"female_collection_date"
"age"	"female_death_date"
"age"	"geo_loc_name"
"age"	"sample_date"
"age"	"sampling_date"
"age"	"specimen_collection_date"
"age"	"time_-_blastomeres"
"age"	"animal_length"
"age"	"date"
"age"	"developemental_stage"
"age"	"colection_date"

"antibody"	"chip_antibody"
"antibody"	"chip_or_ip_antibody"
"antibody"	"clip-antibody"
"antibody"	"clip_antibody"
"antibody"	"rip_antibody"
"antibody"	"rna_binding_protein"

"batch"	"adapter_barcode"
"batch"	"animal"
"batch"	"animal_id"
"batch"	"animalid"
"batch"	"barcode"
"batch"	"biological_replicate"
"batch"	"biorep"
"batch"	"brain_code"
"batch"	"clutch"
"batch"	"cohort"
"batch"	"collected_by"
"batch"	"cow"
"batch"	"custom_field"
"batch"	"donator_monkey_number"
"batch"	"donor_id"
"batch"	"horse_number"
"batch"	"iclip_barcode"
"batch"	"id"
"batch"	"illumina_rna-seq_bar_codes"
"batch"	"index"
"batch"	"individual"
"batch"	"individual_id"
"batch"	"individuals"
"batch"	"internal_id"
"batch"	"labanimalnumber"
"batch"	"labexpid"
"batch"	"mahpic_non_human_primate_individual_id"
"batch"	"monkey_id"
"batch"	"non_human_primate_individual_id"
"batch"	"number"
"batch"	"oocyte/embryo_number"
"batch"	"pig_id"
"batch"	"replicate"
"batch"	"rna_processing_batch"
"batch"	"run_id"
"batch"	"sample"
"batch"	"sample-type"
"batch"	"sample_coding"
"batch"	"sample_descriptor"
"batch"	"sample_id"
"batch"	"sample_identifier"
"batch"	"sample_number"
"batch"	"series"
"batch"	"specimen_voucher"
"batch"	"study"
"batch"	"study_group"
"batch"	"study_phase"
"batch"	"subject"
"batch"	"subject_id"
"batch"	"subjectid"
"batch"	"technical_replicate"
"batch"	"tissue_abbreviation"
"batch"	"tissue_code"
"batch"	"umc_id"
"batch"	"unique_id"
"batch"	"uniqueid"
"batch"	"wur_marker_type"
"batch"	"animal_number"
"batch"	"ebi_equivalent_biosample"
"batch"	"ercc"
"batch"	"flow_cell"
"batch"	"flowcell"
"batch"	"lane"
"batch"	"lane_num"
"batch"	"num_replicons"
"batch"	"cell_barcode"
"batch"	"submission_identifier"
"batch"	"library_id"
"batch"	"library_index_sequence_used_to_demultiplex"
"batch"	"replicate_/_experiment"
"batch"	"sequencing_pool"
"batch"	"plate_col_id"
"batch"	"plate_row_id"
"batch"	"replicate_no."
"batch"	"run_no."
"batch"	"sire"
"batch"	"unique_identifier"

"biomaterial_provider"	"biosourceprovider"
"biomaterial_provider"	"female_biomaterial_provider"
"biomaterial_provider"	"male_biomaterial_provider"
"biomaterial_provider"	"library_preparation_location"
"biomaterial_provider"	"material_provider"
"biomaterial_provider"	"availability"

"bioproject"	"bioproject_id"
"bioproject"	"bioprojectid"

"cell"	"cell_type"
"cell"	"cell_line"
"cell"	"cell_class"
"cell"	"cell_description"
"cell"	"cell_subtype"
"cell"	"arrayexpress-celltype"
"cell"	"culture_collection"
"cell"	"culture_conditions"
"cell"	"cell-line"
"cell"	"cell-type"
"cell"	"cell_or_tissue_type"
"cell"	"eosinophil"
"cell"	"feature"
"cell"	"germ_layer"
"cell"	"cells_derived_from"
"cell"	"number_cells"
"cell"	"cell_organism"
"cell"	"cell_culture_protocol"
"cell"	"cell_population"
"cell"	"cell_typing"

"protocol"	"extract_protocol"
"protocol"	"female_sample_collection_protocol"
"protocol"	"protocol"
"protocol"	"sample_extraction_method"
"protocol"	"experimental_protocol"
"protocol"	"chemical_treatment_of_isolated_rna"
"protocol"	"amplification"
"protocol"	"assay"
"protocol"	"experiment_target"
"protocol"	"extraction_protocol"
"protocol"	"fixation"
"protocol"	"fraction"
"protocol"	"library_selection"
"protocol"	"library_source"
"protocol"	"library_strategy"
"protocol"	"library_type"
"protocol"	"libraryprotocol"
"protocol"	"lysis_buffer_ion_concentrations"
"protocol"	"lysis_strategy"
"protocol"	"male_sample_collection_protocol"
"protocol"	"meoh-fixed"
"protocol"	"molecule_subtype"
"protocol"	"prep_type"
"protocol"	"preperation_kit"
"protocol"	"purification_protocol"
"protocol"	"rin"
"protocol"	"rna_fraction"
"protocol"	"rna_rin_values_after_globin_depletion"
"protocol"	"rna_rin_values_before_globin_depletion"
"protocol"	"rna_subtype"
"protocol"	"sample_material"
"protocol"	"sample_storage"
"protocol"	"sample_storage_processing"
"protocol"	"small_rna_classes"
"protocol"	"specimen_collection_protocol"
"protocol"	"specimen_with_known_storage_state"
"protocol"	"tissue_state"
"protocol"	"store_cond"
"protocol"	"datatype_description"
"protocol"	"datatype"
"protocol"	"monosome_enrichment_strategy"
"protocol"	"rna_concentration_after_globin_depletion"
"protocol"	"rna_concentration_before_globin_depletion"
"protocol"	"rna_input"
"protocol"	"sequencing_type"
"protocol"	"date_run"
"protocol"	"instrument_model"
"protocol"	"ngs_platform"
"protocol"	"paired_end_seq?"
"protocol"	"paired_or_single-end"
"protocol"	"pe_read_length"
"protocol"	"platform"
"protocol"	"quality"
"protocol"	"route"
"protocol"	"read_len_orig"
"protocol"	"readtype"
"protocol"	"readtype_description"
"protocol"	"sequencer"
"protocol"	"sequencing_method"
"protocol"	"extraction"
"protocol"	"indrops_version"
"protocol"	"minimum_counts_per_cell_threshold_used_to_remove_background_barcodes"
"protocol"	"minimum_reads_per_cell_during_initial_processing"
"protocol"	"spike-in"

"nominal_length"	"mean_insert_size"

"genotype"	"arrayexpress-genotype"
"genotype"	"arrayexpress-phenotype"
"genotype"	"arrayexpress-rnai"
"genotype"	"background_strain"
"genotype"	"breed"
"genotype"	"cultivar"
"genotype"	"chicken_line"
"genotype"	"ecotype"
"genotype"	"full_genotype"
"genotype"	"genetic_background"
"genotype"	"genetic_line"
"genotype"	"genetic_modification"
"genotype"	"genotype/variaion"
"genotype"	"genotype/variation"
"genotype"	"germline_knock-down"
"genotype"	"germline_knockdown_and_other_transgenes"
"genotype"	"isolate"
"genotype"	"line"
"genotype"	"phenotype"
"genotype"	"ploidy"
"genotype"	"population"
"genotype"	"snp_allelic_state"
"genotype"	"strain"
"genotype"	"strain/background"
"genotype"	"strain_origin"
"genotype"	"strains"
"genotype"	"sub_species"
"genotype"	"subspecies"
"genotype"	"variety"
"genotype"	"breed_name"
"genotype"	"arrayexpress-strainorline"
"genotype"	"breed/line"
"genotype"	"breed_type"
"genotype"	"breeding_history"
"genotype"	"breeds"
"genotype"	"chicken_breed"
"genotype"	"clonality"
"genotype"	"marker"
"genotype"	"dstim_knockdown_status"
"genotype"	"environmental_history"
"genotype"	"fly_line"
"genotype"	"origin"
"genotype"	"phenotype/variation"
"genotype"	"pig_breed"
"genotype"	"pig_type"
"genotype"	"propagation"
"genotype"	"reporter_construct"
"genotype"	"rhesus_monkey_origin"
"genotype"	"rna_interference"
"genotype"	"strain_background"
"genotype"	"strain_description"
"genotype"	"strain_prrsv"
"genotype"	"strainorline"
"genotype"	"transfecting_vector"
"genotype"	"transgenic_strain"
"genotype"	"virus_group"
"genotype"	"arrayexpress-species"
"genotype"	"parental_strain"
"genotype"	"tax_id"
"genotype"	"transgene"
"genotype"	"genotype/phenotype"
"genotype"	"cross"
"genotype"	"dad"
"genotype"	"driver"
"genotype"	"qtl_genotype/haplotype"

"location"	"birth_location"
"location"	"country"
"location"	"feeding_barn"
"location"	"female_birth_location"
"location"	"geographic_location"
"location"	"lat_lon"
"location"	"library_preparation_location_latitude"
"location"	"library_preparation_location_latitude_units"
"location"	"library_preparation_location_longitude"
"location"	"library_preparation_location_longitude_units"
"location"	"local_origin"
"location"	"male_birth_location"
"location"	"sampling_site"
"location"	"seq_center"
"location"	"sequencing_location"
"location"	"sequencing_location_latitude"
"location"	"station"
"location"	"sequencing_location_latitude_units"
"location"	"sequencing_location_longitude"
"location"	"sequencing_location_longitude_units"
"location"	"source"
"location"	"biome"
"location"	"env_biome"
"location"	"env_feature"
"location"	"env_material"
"location"	"environment"
"location"	"environmental_package"

"misc"	"ena-checklist"
"misc"	"ena-first-public"
"misc"	"ena-last-update"
"misc"	"estimated_size"
"misc"	"investigation_type"
"misc"	"lab_description"
"misc"	"library_layout"
"misc"	"mapalgorithm"
"misc"	"mapalgorithm_description"
"misc"	"new_attribute"
"misc"	"notes"
"misc"	"organism"
"misc"	"project"
"misc"	"project_name"
"misc"	"species"
"misc"	"sra_experiment_accession"
"misc"	"biosamplemodel"
"misc"	"sample_year"
"misc"	"number_of_pieces"
"misc"	"specimen_picture_url"
"misc"	"submission_description"
"misc"	"submission_title"
"misc"	"unknown"
"misc"	"body_weight"
"misc"	"file_no."
"misc"	"num_parts_in_pool"
"misc"	"processing_date"

"sex"	"arrayexpress-sex"
"sex"	"gender"
"sex"	"host_sex"
"sex"	"cell_sex"
"sex"	"gender_type"
"sex"	"obsolete_sex"
"sex"	"both_sexes_pooled"

"tissue"	"arrayexpress-organismpart"
"tissue"	"bio_material"
"tissue"	"body_site"
"tissue"	"biopsy_site"
"tissue"	"brain_region"
"tissue"	"description"
"tissue"	"experiment_set"
"tissue"	"explant"
"tissue"	"isolation-source"
"tissue"	"isolation_source"
"tissue"	"label"
"tissue"	"mixedtissues"
"tissue"	"muscle_type"
"tissue"	"organ/tissue"
"tissue"	"organ_part"
"tissue"	"organism_part"
"tissue"	"organismpart"
"tissue"	"oviduct"
"tissue"	"region"
"tissue"	"sample_name"
"tissue"	"sample_origin"
"tissue"	"sampling_position"
"tissue"	"source_name"
"tissue"	"tag"
"tissue"	"tissue-type"
"tissue"	"tissue/cell"
"tissue"	"tissue/cell_type"
"tissue"	"tissue_location"
"tissue"	"tissue_or_dev_stage"
"tissue"	"tissue_origin"
"tissue"	"tissue_source"
"tissue"	"tissue_type"
"tissue"	"title"
"tissue"	"arrayexpress-organismpart"
"tissue"	"organism_part"
"tissue"	"embryo_region"
"tissue"	"organsim_part"
"tissue"	"tissue_/_cells"
"tissue"	"tissue_lib"
"tissue"	"tissue_region"

"treatment"	"agent"
"treatment"	"arrayexpress-diseasestate"
"treatment"	"arrayexpress-growthcondition"
"treatment"	"arrayexpress-immunoprecipitate"
"treatment"	"behavior"
"treatment"	"biopsy_day"
"treatment"	"challenged_with"
"treatment"	"chemical_treatment_of_isolated_rna"
"treatment"	"diet"
"treatment"	"diet_fed"
"treatment"	"dietary_group"
"treatment"	"disease"
"treatment"	"disease_status"
"treatment"	"domestication"
"treatment"	"drip"
"treatment"	"energy_balance"
"treatment"	"fed_status"
"treatment"	"feeding"
"treatment"	"feeding_type"
"treatment"	"fertility_group"
"treatment"	"food"
"treatment"	"group"
"treatment"	"health_state"
"treatment"	"infect"
"treatment"	"infected_with"
"treatment"	"infection"
"treatment"	"infection_status"
"treatment"	"inoculation"
"treatment"	"maternal_diet"
"treatment"	"meat_quality"
"treatment"	"mptp_treatment"
"treatment"	"oxidation"
"treatment"	"pathogen"
"treatment"	"percent_marinade_uptake_at_24h_post-slaughter"
"treatment"	"phenotype_sample_type"
"treatment"	"restriction_model"
"treatment"	"altitude"
"treatment"	"dpi"
"treatment"	"salmonella_shedding_status"
"treatment"	"status"
"treatment"	"stimulation"
"treatment"	"stimulus"
"treatment"	"stress"
"treatment"	"survival"
"treatment"	"transplant_type"
"treatment"	"treated_with"
"treatment"	"treatment_group"
"treatment"	"pregnancy_outcome"
"treatment"	"isol_growth_condt"
"treatment"	"agent"
"treatment"	"aliquot"
"treatment"	"ammonia_concentration"
"treatment"	"ammonia_exposure"
"treatment"	"average_altitude"
"treatment"	"blastocyst_rate"
"treatment"	"breeding_method"
"treatment"	"bull_ntm"
"treatment"	"calcium_intake"
"treatment"	"cell_status"
"treatment"	"cellular_compartment"
"treatment"	"chemical"
"treatment"	"collection_time"
"treatment"	"comment"
"treatment"	"concentration"
"treatment"	"condition"
"treatment"	"culture_condition"
"treatment"	"culture-collection"
"treatment"	"culture_type"
"treatment"	"custom_name"
"treatment"	"day_of_infection"
"treatment"	"day_post_infection"
"treatment"	"days_at_29_c"
"treatment"	"days_post_infection"
"treatment"	"differentiation_time"
"treatment"	"disease_stage"
"treatment"	"drinking_category"
"treatment"	"drug"
"treatment"	"drug_for_synchronization"
"treatment"	"egg_production"
"treatment"	"enzyme_treatment"
"treatment"	"exposure_time"
"treatment"	"extracellular_component"
"treatment"	"fasted_status"
"treatment"	"fed_for"
"treatment"	"feed_efficiency"
"treatment"	"feeding_period"
"treatment"	"feeding_treatment"
"treatment"	"growth_condition"
"treatment"	"growth_protocol"
"treatment"	"harmonized_disease_score"
"treatment"	"harvest_time"
"treatment"	"health-state"
"treatment"	"health_status"
"treatment"	"health_status_at_collection"
"treatment"	"heat_hours"
"treatment"	"hematocrit"
"treatment"	"hemoglobin_measurement"
"treatment"	"host"
"treatment"	"host_tissue_sampled"
"treatment"	"hours_at_restrictive_temperature"
"treatment"	"hp-prrsv_infection"
"treatment"	"infection_route"
"treatment"	"infectious_agent"
"treatment"	"inflammatory_stimulus_treatment"
"treatment"	"injury"
"treatment"	"kinetics"
"treatment"	"lens_type"
"treatment"	"lesion_score"
"treatment"	"lesion_number"
"treatment"	"library_name"
"treatment"	"litter"
"treatment"	"loin_weight_group"
"treatment"	"lps_exposure"
"treatment"	"material"
"treatment"	"maternal_treatment"
"treatment"	"mating_status"
"treatment"	"mock"
"treatment"	"mode"
"treatment"	"morpholino"
"treatment"	"muscle_name"
"treatment"	"neutralization_group"
"treatment"	"normalized?"
"treatment"	"number_born_alive"
"treatment"	"oxygen_concentration"
"treatment"	"parity_number"
"treatment"	"passage"
"treatment"	"passage_number"
"treatment"	"passagers"
"treatment"	"passages"
"treatment"	"phase"
"treatment"	"pooled"
"treatment"	"pregnancy_status"
"treatment"	"pregnancy_time"
"treatment"	"purity"
"treatment"	"resorted"
"treatment"	"rna_treatment"
"treatment"	"sample_processing"
"treatment"	"sampleid"
"treatment"	"sampling_time_point"
"treatment"	"sampling_to_preparation_interval"
"treatment"	"selection_criteria"
"treatment"	"serologic_response_status"
"treatment"	"shear_force"
"treatment"	"skin_color"
"treatment"	"social_rank"
"treatment"	"sorted"
"treatment"	"specific_host"
"treatment"	"stim_time"
"treatment"	"stud_book_number"
"treatment"	"superovulated"
"treatment"	"temperature"
"treatment"	"thermal_treatment"
"treatment"	"time"
"treatment"	"time-point"
"treatment"	"time_point"
"treatment"	"ega-disease"
"treatment"	"time_post_infection"
"treatment"	"timepoint_postchallenge"
"treatment"	"training"
"treatment"	"training_timepoint"
"treatment"	"treatment_description"
"treatment"	"tumor"
"treatment"	"tumor_grading"
"treatment"	"vaccine"
"treatment"	"vaccine_administration"
"treatment"	"vaccine_group"
"treatment"	"viremia_level"
"treatment"	"virus"
"treatment"	"virus_dose"
"treatment"	"d"
"treatment"	"days_post-infection"
"treatment"	"ega-subject"
"treatment"	"ega-treatment"
"treatment"	"ercc_pool"
"treatment"	"host_taxid"
"treatment"	"lymphocyte"
"treatment"	"monocyte"
"treatment"	"morphology"
"treatment"	"neutrophil"
"treatment"	"platelet_count"
"treatment"	"rfi_value"
"treatment"	"sample_type"
"treatment"	"source_of_female"
"treatment"	"sperm_dna_fragmentation_index"
"treatment"	"srek_ligation"
"treatment"	"time/hpf"
"treatment"	"total_number_born"
"treatment"	"treatments"
"treatment"	"type"
"treatment"	"weight"
"treatment"	"arrayexpress-compound"
"treatment"	"arrayexpress-dose"
"treatment"	"assay_type"
"treatment"	"backfat_thickness_of_live(mm)"
"treatment"	"basophil"
"treatment"	"photorecptors_remaining"
"treatment"	"productivity"
"treatment"	"hours/days_after_injury"
"treatment"	"disease_state"
"treatment"	"mehg_exposure"
"treatment"	"operation"
"treatment"	"post-injury_time_point"
"treatment"	"embryo_source"
"treatment"	"earlobe_color"
"treatment"	"reproductive_status"
"treatment"	"specimen_size"
"treatment"	"specimen_volume"
"treatment"	"specimen_weight"
"treatment"	"liver_phenotype"
"treatment"	"clinical_information"
"treatment"	"oestrous_peroid"
"treatment"	"secondary_description"
"treatment"	"finishing_diet"
"treatment"	"gestation_duration"
"treatment"	"growth_rate"
"treatment"	"heart_weight"
"treatment"	"information"
"treatment"	"light_treatment"
"treatment"	"maternal_food"
"treatment"	"number_of_passages"
"treatment"	"paternal_experience"
"treatment"	"personal_experience"
"treatment"	"physiological_conditions"
"treatment"	"rfi_group"
"treatment"	"sample_role"
"treatment"	"slaughter_weight"
"treatment"	"target"
"treatment"	"target_source"
"treatment"	"temperature_regimen"




================================================
FILE: amalgkit/config_dir/vertebrate/__init__.py
================================================



================================================
FILE: amalgkit/config_dir/vertebrate/control_term.config
================================================
# case-insensitive
# regular expressions allowed
# (attribute)[TAB](control term)

"treatment"	"mock"
"treatment"	"control"
"treatment"	"cntl"
"treatment"	"con\d"
"treatment"	"con_"
"treatment"	"cont\b"
"treatment"	"cont_"
"treatment"	"kontrol"
"treatment"	"contral"
"treatment"	"contron"
"treatment"	"\bCK"
"treatment"	"none"
"treatment"	"normal"
"treatment"	"natural"
"treatment"	"standard"
"treatment"	"optimal"
"treatment"	"sufficient"
"treatment"	"no(?:\s+\w+)?\s+treat"
"treatment"	"non-treat"
"treatment"	"non- treat"
"treatment"	"nontreat"
"treatment"	"untreat"
"treatment"	"without"
"treatment"	"^no$"
"treatment"	"Ad Libitum"
"treatment"	"in the absence of"
"treatment"	"Thermoneutral"
"treatment"	"21%[oxygen_concentration]"
"treatment"	"DMSO"
"treatment"	"\b0h"
"treatment"	"\b0 h"
"treatment"	"_0h"
"treatment"	"\b0day"
"treatment"	"\b0 day"
"treatment"	"\b0m"
"treatment"	"\b0 m"
"treatment"	"\b0 dpi"
"treatment"	"^0 ppm"
"treatment"	"wild"
"treatment"	"uninfected"
"treatment"	"un-infected"
"treatment"	"non infected"
"treatment"	"non-inifected"
"treatment"	"no inoculation"
"treatment"	"un-inoculated"
"treatment"	"no disease"
"treatment"	"non-diseased"
"treatment"	"unwound"
"treatment"	"un-grazing"
"treatment"	"healthy"
"treatment"	"stress-free"
"treatment"	"nostress"
"treatment"	"non-stress"
"treatment"	"unstressed"
"treatment"	"non(?:\s+\w+)?\s+stress"
"treatment"	"room temperature"
"treatment"	"ambient"
"treatment"	"water\w*\b(?!\s*deficit)"
"treatment"	"H2O\w*\b(?!\s*deficit)"

"genotype"	"wild"
"genotype"	"N2"
"genotype"	"OregonR"
"genotype"	"Oregon R"
"genotype"	"Canton S"
"genotype"	"Canton-S"
"genotype"	"y w"
"genotype"	"yw"
"genotype"	"wt"
"genotype"	"\+\/\+"

"sample_title"	"wild"
"sample_title"	"\bWT"
"sample_title"	"control"
"sample_title"	"CON_"
"sample_title"	"CK_"
"sample_title"	"CK-"
"sample_title"	"CK\d"
"sample_title"	"\bCK\b"
"sample_title"	"normal"
"sample_title"	"non-treat"
"sample_title"	"non-inoculated"
"sample_title"	"ungrazed"
"sample_title"	"healthy"
"sample_title"	"water"
"sample_title"	"H2O"
"sample_title"	"\b0dpi"

"sample_description"	"control"
"sample_description"	"mock"
"sample_description"	"no treatment"
"sample_description"	"non-infected"
"sample_description"	"before"
"sample_description"	"wild"

"lib_name"	"control"
"lib_name"	"CK_"
"lib_name"	"CK-"
"lib_name"	"CK\d"
"lib_name"	"\bCK\b"
"lib_name"	"\bcon-"
"lib_name"	"ctrl"
"lib_name"	"mock"
"lib_name"	"wild"
"lib_name"	"\bWT"
"lib_name"	"(?<!\d)0h"

"exp_title"	"control"
"exp_title"	"CON_"
"exp_title"	"CK_"
"exp_title"	"CK-"
"exp_title"	"CK\d"
"exp_title"	"\bCK\b"
"exp_title"	"untreated"
"exp_title"	"mock"
"exp_title"	"without"
"exp_title"	"water"
"exp_title"	"H2O"
"exp_title"	"TimePoint1_"
"exp_title"	"0 days post infection"
"exp_title"	"\b0h"
"exp_title"	"\b0 h"
"exp_title"	"_0h"

"design"	"control"
"design"	"CK_"
"design"	"CK-"
"design"	"CK\d"
"design"	"\bCK\b"
"design"	"\bcon-"
"design"	"mock"
"design"	"normal"
"design"	"common"
"design"	"untreat"
"design"	"non-treat"
"design"	"no(?:\s+\w+)?\s+treat"
"design"	"no processing"
"design"	"without"
"design"	"water"
"design"	"H2O"
"design"	"before"
"design"	"\b0h"
"design"	"\b0 h"
"design"	"_0h"
"design"	"\b0d"
"design"	"\b0 d"
"design"	"\b0 m"
"design"	"healthy"
"design"	"ungrazed"
"design"	"no grazing"
"design"	"no-grazing"
"design"	"not inoculated"
"design"	"non-inoculated"
"design"	"uninfected"
"design"	"no incubation"
"design"	"not primed"
"design"	"not induced"
"design"	"wt"
"design"	"wild"



================================================
FILE: amalgkit/config_dir/vertebrate/exclude_keyword.config
================================================
# case-insensitive
# regular expressions allowed
# (comma-delimited attributes to be searched)[TAB](arbitrary reason of exclusion)[TAB](bad keyword)

"exp_title,study_title,design,sample_title,sample_description,lib_name,experiment,treatment,protocol,age"	"single_cell"	"single cell"

"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"immunoprecipitation"	"RipSeq"
"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"immunoprecipitation"	"chrom_RNAseq"
"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"RNAi"	"RNAi"
"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"small_RNA"	"shRNA RNA-seq"
"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"small_RNA"	"piRNA"
"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"small_RNA"	"smRNA"
"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol"	"CAGE"	"CAGE"

"exp_title,study_title,design,sample_title,sample_description,lib_name,protocol,treatment"	"treatment"	"^exposed to"

"treatment,protocol"	"small_RNA"	"miRNA"
"treatment,protocol"	"small_RNA"	"sRNA"
"treatment,protocol"	"treatment"	"uranium"

"age"	"embryonic"	"embryo"
"age"	"embryonic"	"somite"
"age"	"embryonic"	"fertiliz.*"
"age"	"embryonic"	"oocyte"

"genotype"	"transgenic"	"foxn3"



================================================
FILE: amalgkit/config_dir/vertebrate/group_attribute.config
================================================
# case-insensitive
# regular expressions allowed
# (aggregated to)[TAB](aggregated from)

"age"	"arrayexpress-developmentalstage"
"age"	"agedays"
"age"	"age_classification"
"age"	"age_description"
"age"	"age_in_years"
"age"	"agecat"
"age"	"dev_stage"
"age"	"ages"
"age"	"development_satge"
"age"	"development_stage"
"age"	"developmental_stage"
"age"	"developmentalstage"
"age"	"differentiation_stages"
"age"	"embryonic_day"
"age"	"embryonic_stage"
"age"	"female_age"
"age"	"sample_comment"
"age"	"stage"
"age"	"age_at_collection"
"age"	"age_classification"
"age"	"age_description"
"age"	"age_in_years"
"age"	"age_of_collection"
"age"	"age_of_fly_in_days_post_eclosion"
"age"	"agecat"
"age"	"agedays"
"age"	"age{category}"
"age"	"animal_age_at_collection"
"age"	"day"
"age"	"arrayexpress-timepoint"
"age"	"avg_age_after_second_6_months_of_open_access"
"age"	"collection_timing_in_reference_to_e_lineage"
"age"	"day_relative_to_weaning"
"age"	"dev-stage"
"age"	"developmental_time"
"age"	"duration"
"age"	"embryo_stage"
"age"	"female_birth_date"
"age"	"gestational_age"
"age"	"gestational_age_at_sample_collection"
"age"	"hotta_stage"
"age"	"male_age"
"age"	"male_birth_date"
"age"	"male_collection_date"
"age"	"male_death_date"
"age"	"timepoint"
"age"	"trimester"
"age"	"birth_date"
"age"	"collection_date"
"age"	"death_date"
"age"	"female_collection_date"
"age"	"female_death_date"
"age"	"geo_loc_name"
"age"	"sample_date"
"age"	"sampling_date"
"age"	"specimen_collection_date"
"age"	"time_-_blastomeres"
"age"	"animal_length"
"age"	"date"
"age"	"developemental_stage"
"age"	"colection_date"

"antibody"	"chip_antibody"
"antibody"	"chip_or_ip_antibody"
"antibody"	"clip-antibody"
"antibody"	"clip_antibody"
"antibody"	"rip_antibody"
"antibody"	"rna_binding_protein"

"batch"	"adapter_barcode"
"batch"	"animal"
"batch"	"animal_id"
"batch"	"animalid"
"batch"	"barcode"
"batch"	"biological_replicate"
"batch"	"biorep"
"batch"	"brain_code"
"batch"	"clutch"
"batch"	"cohort"
"batch"	"collected_by"
"batch"	"cow"
"batch"	"custom_field"
"batch"	"donator_monkey_number"
"batch"	"donor_id"
"batch"	"horse_number"
"batch"	"iclip_barcode"
"batch"	"id"
"batch"	"illumina_rna-seq_bar_codes"
"batch"	"index"
"batch"	"individual"
"batch"	"individual_id"
"batch"	"individuals"
"batch"	"internal_id"
"batch"	"labanimalnumber"
"batch"	"labexpid"
"batch"	"mahpic_non_human_primate_individual_id"
"batch"	"monkey_id"
"batch"	"non_human_primate_individual_id"
"batch"	"number"
"batch"	"oocyte/embryo_number"
"batch"	"pig_id"
"batch"	"replicate"
"batch"	"rna_processing_batch"
"batch"	"run_id"
"batch"	"sample"
"batch"	"sample-type"
"batch"	"sample_coding"
"batch"	"sample_descriptor"
"batch"	"sample_id"
"batch"	"sample_identifier"
"batch"	"sample_number"
"batch"	"series"
"batch"	"specimen_voucher"
"batch"	"study"
"batch"	"study_group"
"batch"	"study_phase"
"batch"	"subject"
"batch"	"subject_id"
"batch"	"subjectid"
"batch"	"technical_replicate"
"batch"	"tissue_abbreviation"
"batch"	"tissue_code"
"batch"	"umc_id"
"batch"	"unique_id"
"batch"	"uniqueid"
"batch"	"wur_marker_type"
"batch"	"animal_number"
"batch"	"ebi_equivalent_biosample"
"batch"	"ercc"
"batch"	"flow_cell"
"batch"	"flowcell"
"batch"	"lane"
"batch"	"lane_num"
"batch"	"num_replicons"
"batch"	"cell_barcode"
"batch"	"submission_identifier"
"batch"	"library_id"
"batch"	"library_index_sequence_used_to_demultiplex"
"batch"	"replicate_/_experiment"
"batch"	"sequencing_pool"
"batch"	"plate_col_id"
"batch"	"plate_row_id"
"batch"	"replicate_no."
"batch"	"run_no."
"batch"	"sire"
"batch"	"unique_identifier"

"biomaterial_provider"	"biosourceprovider"
"biomaterial_provider"	"female_biomaterial_provider"
"biomaterial_provider"	"male_biomaterial_provider"
"biomaterial_provider"	"library_preparation_location"
"biomaterial_provider"	"material_provider"
"biomaterial_provider"	"availability"

"bioproject"	"bioproject_id"
"bioproject"	"bioprojectid"

"cell"	"cell_type"
"cell"	"cell_line"
"cell"	"cell_class"
"cell"	"cell_description"
"cell"	"cell_subtype"
"cell"	"arrayexpress-celltype"
"cell"	"culture_collection"
"cell"	"culture_conditions"
"cell"	"cell-line"
"cell"	"cell-type"
"cell"	"cell_or_tissue_type"
"cell"	"eosinophil"
"cell"	"feature"
"cell"	"germ_layer"
"cell"	"cells_derived_from"
"cell"	"number_cells"
"cell"	"cell_organism"
"cell"	"cell_culture_protocol"
"cell"	"cell_population"
"cell"	"cell_typing"

"protocol"	"extract_protocol"
"protocol"	"female_sample_collection_protocol"
"protocol"	"protocol"
"protocol"	"sample_extraction_method"
"protocol"	"experimental_protocol"
"protocol"	"chemical_treatment_of_isolated_rna"
"protocol"	"amplification"
"protocol"	"assay"
"protocol"	"experiment_target"
"protocol"	"extraction_protocol"
"protocol"	"fixation"
"protocol"	"fraction"
"protocol"	"library_selection"
"protocol"	"library_source"
"protocol"	"library_strategy"
"protocol"	"library_type"
"protocol"	"libraryprotocol"
"protocol"	"lysis_buffer_ion_concentrations"
"protocol"	"lysis_strategy"
"protocol"	"male_sample_collection_protocol"
"protocol"	"meoh-fixed"
"protocol"	"molecule_subtype"
"protocol"	"prep_type"
"protocol"	"preperation_kit"
"protocol"	"purification_protocol"
"protocol"	"rin"
"protocol"	"rna_fraction"
"protocol"	"rna_rin_values_after_globin_depletion"
"protocol"	"rna_rin_values_before_globin_depletion"
"protocol"	"rna_subtype"
"protocol"	"sample_material"
"protocol"	"sample_storage"
"protocol"	"sample_storage_processing"
"protocol"	"small_rna_classes"
"protocol"	"specimen_collection_protocol"
"protocol"	"specimen_with_known_storage_state"
"protocol"	"tissue_state"
"protocol"	"store_cond"
"protocol"	"datatype_description"
"protocol"	"datatype"
"protocol"	"monosome_enrichment_strategy"
"protocol"	"rna_concentration_after_globin_depletion"
"protocol"	"rna_concentration_before_globin_depletion"
"protocol"	"rna_input"
"protocol"	"sequencing_type"
"protocol"	"date_run"
"protocol"	"instrument_model"
"protocol"	"ngs_platform"
"protocol"	"paired_end_seq?"
"protocol"	"paired_or_single-end"
"protocol"	"pe_read_length"
"protocol"	"platform"
"protocol"	"quality"
"protocol"	"route"
"protocol"	"read_len_orig"
"protocol"	"readtype"
"protocol"	"readtype_description"
"protocol"	"sequencer"
"protocol"	"sequencing_method"
"protocol"	"extraction"
"protocol"	"indrops_version"
"protocol"	"minimum_counts_per_cell_threshold_used_to_remove_background_barcodes"
"protocol"	"minimum_reads_per_cell_during_initial_processing"
"protocol"	"spike-in"

"nominal_length"	"mean_insert_size"

"genotype"	"arrayexpress-genotype"
"genotype"	"arrayexpress-phenotype"
"genotype"	"arrayexpress-rnai"
"genotype"	"background_strain"
"genotype"	"breed"
"genotype"	"cultivar"
"genotype"	"chicken_line"
"genotype"	"ecotype"
"genotype"	"full_genotype"
"genotype"	"genetic_background"
"genotype"	"genetic_line"
"genotype"	"genetic_modification"
"genotype"	"genotype/variaion"
"genotype"	"genotype/variation"
"genotype"	"germline_knock-down"
"genotype"	"germline_knockdown_and_other_transgenes"
"genotype"	"isolate"
"genotype"	"line"
"genotype"	"phenotype"
"genotype"	"ploidy"
"genotype"	"population"
"genotype"	"snp_allelic_state"
"genotype"	"strain"
"genotype"	"strain/background"
"genotype"	"strain_origin"
"genotype"	"strains"
"genotype"	"sub_species"
"genotype"	"subspecies"
"genotype"	"variety"
"genotype"	"breed_name"
"genotype"	"arrayexpress-strainorline"
"genotype"	"breed/line"
"genotype"	"breed_type"
"genotype"	"breeding_history"
"genotype"	"breeds"
"genotype"	"chicken_breed"
"genotype"	"clonality"
"genotype"	"marker"
"genotype"	"dstim_knockdown_status"
"genotype"	"environmental_history"
"genotype"	"fly_line"
"genotype"	"origin"
"genotype"	"phenotype/variation"
"genotype"	"pig_breed"
"genotype"	"pig_type"
"genotype"	"propagation"
"genotype"	"reporter_construct"
"genotype"	"rhesus_monkey_origin"
"genotype"	"rna_interference"
"genotype"	"strain_background"
"genotype"	"strain_description"
"genotype"	"strain_prrsv"
"genotype"	"strainorline"
"genotype"	"transfecting_vector"
"genotype"	"transgenic_strain"
"genotype"	"virus_group"
"genotype"	"arrayexpress-species"
"genotype"	"parental_strain"
"genotype"	"tax_id"
"genotype"	"transgene"
"genotype"	"genotype/phenotype"
"genotype"	"cross"
"genotype"	"dad"
"genotype"	"driver"
"genotype"	"qtl_genotype/haplotype"

"location"	"birth_location"
"location"	"country"
"location"	"feeding_barn"
"location"	"female_birth_location"
"location"	"geographic_location"
"location"	"lat_lon"
"location"	"library_preparation_location_latitude"
"location"	"library_preparation_location_latitude_units"
"location"	"library_preparation_location_longitude"
"location"	"library_preparation_location_longitude_units"
"location"	"local_origin"
"location"	"male_birth_location"
"location"	"sampling_site"
"location"	"seq_center"
"location"	"sequencing_location"
"location"	"sequencing_location_latitude"
"location"	"station"
"location"	"sequencing_location_latitude_units"
"location"	"sequencing_location_longitude"
"location"	"sequencing_location_longitude_units"
"location"	"source"
"location"	"biome"
"location"	"env_biome"
"location"	"env_feature"
"location"	"env_material"
"location"	"environment"
"location"	"environmental_package"

"misc"	"ena-checklist"
"misc"	"ena-first-public"
"misc"	"ena-last-update"
"misc"	"estimated_size"
"misc"	"investigation_type"
"misc"	"lab_description"
"misc"	"library_layout"
"misc"	"mapalgorithm"
"misc"	"mapalgorithm_description"
"misc"	"new_attribute"
"misc"	"notes"
"misc"	"organism"
"misc"	"project"
"misc"	"project_name"
"misc"	"species"
"misc"	"sra_experiment_accession"
"misc"	"biosamplemodel"
"misc"	"sample_year"
"misc"	"number_of_pieces"
"misc"	"specimen_picture_url"
"misc"	"submission_description"
"misc"	"submission_title"
"misc"	"unknown"
"misc"	"body_weight"
"misc"	"file_no."
"misc"	"num_parts_in_pool"
"misc"	"processing_date"

"sex"	"arrayexpress-sex"
"sex"	"gender"
"sex"	"host_sex"
"sex"	"cell_sex"
"sex"	"gender_type"
"sex"	"obsolete_sex"
"sex"	"both_sexes_pooled"

"tissue"	"arrayexpress-organismpart"
"tissue"	"bio_material"
"tissue"	"body_site"
"tissue"	"biopsy_site"
"tissue"	"brain_region"
"tissue"	"description"
"tissue"	"experiment_set"
"tissue"	"explant"
"tissue"	"isolation-source"
"tissue"	"isolation_source"
"tissue"	"label"
"tissue"	"mixedtissues"
"tissue"	"muscle_type"
"tissue"	"organ/tissue"
"tissue"	"organ_part"
"tissue"	"organism_part"
"tissue"	"organismpart"
"tissue"	"oviduct"
"tissue"	"region"
"tissue"	"sample_name"
"tissue"	"sample_origin"
"tissue"	"sampling_position"
"tissue"	"source_name"
"tissue"	"tag"
"tissue"	"tissue-type"
"tissue"	"tissue/cell"
"tissue"	"tissue/cell_type"
"tissue"	"tissue_location"
"tissue"	"tissue_or_dev_stage"
"tissue"	"tissue_origin"
"tissue"	"tissue_source"
"tissue"	"tissue_type"
"tissue"	"title"
"tissue"	"arrayexpress-organismpart"
"tissue"	"organism_part"
"tissue"	"embryo_region"
"tissue"	"organsim_part"
"tissue"	"tissue_/_cells"
"tissue"	"tissue_lib"
"tissue"	"tissue_region"

"treatment"	"agent"
"treatment"	"arrayexpress-diseasestate"
"treatment"	"arrayexpress-growthcondition"
"treatment"	"arrayexpress-immunoprecipitate"
"treatment"	"behavior"
"treatment"	"biopsy_day"
"treatment"	"challenged_with"
"treatment"	"chemical_treatment_of_isolated_rna"
"treatment"	"diet"
"treatment"	"diet_fed"
"treatment"	"dietary_group"
"treatment"	"disease"
"treatment"	"disease_status"
"treatment"	"domestication"
"treatment"	"drip"
"treatment"	"energy_balance"
"treatment"	"fed_status"
"treatment"	"feeding"
"treatment"	"feeding_type"
"treatment"	"fertility_group"
"treatment"	"food"
"treatment"	"group"
"treatment"	"health_state"
"treatment"	"infect"
"treatment"	"infected_with"
"treatment"	"infection"
"treatment"	"infection_status"
"treatment"	"inoculation"
"treatment"	"maternal_diet"
"treatment"	"meat_quality"
"treatment"	"mptp_treatment"
"treatment"	"oxidation"
"treatment"	"pathogen"
"treatment"	"percent_marinade_uptake_at_24h_post-slaughter"
"treatment"	"phenotype_sample_type"
"treatment"	"restriction_model"
"treatment"	"altitude"
"treatment"	"dpi"
"treatment"	"salmonella_shedding_status"
"treatment"	"status"
"treatment"	"stimulation"
"treatment"	"stimulus"
"treatment"	"stress"
"treatment"	"survival"
"treatment"	"transplant_type"
"treatment"	"treated_with"
"treatment"	"treatment_group"
"treatment"	"pregnancy_outcome"
"treatment"	"isol_growth_condt"
"treatment"	"agent"
"treatment"	"aliquot"
"treatment"	"ammonia_concentration"
"treatment"	"ammonia_exposure"
"treatment"	"average_altitude"
"treatment"	"blastocyst_rate"
"treatment"	"breeding_method"
"treatment"	"bull_ntm"
"treatment"	"calcium_intake"
"treatment"	"cell_status"
"treatment"	"cellular_compartment"
"treatment"	"chemical"
"treatment"	"collection_time"
"treatment"	"comment"
"treatment"	"concentration"
"treatment"	"condition"
"treatment"	"culture_condition"
"treatment"	"culture-collection"
"treatment"	"culture_type"
"treatment"	"custom_name"
"treatment"	"day_of_infection"
"treatment"	"day_post_infection"
"treatment"	"days_at_29_c"
"treatment"	"days_post_infection"
"treatment"	"differentiation_time"
"treatment"	"disease_stage"
"treatment"	"drinking_category"
"treatment"	"drug"
"treatment"	"drug_for_synchronization"
"treatment"	"egg_production"
"treatment"	"enzyme_treatment"
"treatment"	"exposure_time"
"treatment"	"extracellular_component"
"treatment"	"fasted_status"
"treatment"	"fed_for"
"treatment"	"feed_efficiency"
"treatment"	"feeding_period"
"treatment"	"feeding_treatment"
"treatment"	"growth_condition"
"treatment"	"growth_protocol"
"treatment"	"harmonized_disease_score"
"treatment"	"harvest_time"
"treatment"	"health-state"
"treatment"	"health_status"
"treatment"	"health_status_at_collection"
"treatment"	"heat_hours"
"treatment"	"hematocrit"
"treatment"	"hemoglobin_measurement"
"treatment"	"host"
"treatment"	"host_tissue_sampled"
"treatment"	"hours_at_restrictive_temperature"
"treatment"	"hp-prrsv_infection"
"treatment"	"infection_route"
"treatment"	"infectious_agent"
"treatment"	"inflammatory_stimulus_treatment"
"treatment"	"injury"
"treatment"	"kinetics"
"treatment"	"lens_type"
"treatment"	"lesion_score"
"treatment"	"lesion_number"
"treatment"	"library_name"
"treatment"	"litter"
"treatment"	"loin_weight_group"
"treatment"	"lps_exposure"
"treatment"	"material"
"treatment"	"maternal_treatment"
"treatment"	"mating_status"
"treatment"	"mock"
"treatment"	"mode"
"treatment"	"morpholino"
"treatment"	"muscle_name"
"treatment"	"neutralization_group"
"treatment"	"normalized?"
"treatment"	"number_born_alive"
"treatment"	"oxygen_concentration"
"treatment"	"parity_number"
"treatment"	"passage"
"treatment"	"passage_number"
"treatment"	"passagers"
"treatment"	"passages"
"treatment"	"phase"
"treatment"	"pooled"
"treatment"	"pregnancy_status"
"treatment"	"pregnancy_time"
"treatment"	"purity"
"treatment"	"resorted"
"treatment"	"rna_treatment"
"treatment"	"sample_processing"
"treatment"	"sampleid"
"treatment"	"sampling_time_point"
"treatment"	"sampling_to_preparation_interval"
"treatment"	"selection_criteria"
"treatment"	"serologic_response_status"
"treatment"	"shear_force"
"treatment"	"skin_color"
"treatment"	"social_rank"
"treatment"	"sorted"
"treatment"	"specific_host"
"treatment"	"stim_time"
"treatment"	"stud_book_number"
"treatment"	"superovulated"
"treatment"	"temperature"
"treatment"	"thermal_treatment"
"treatment"	"time"
"treatment"	"time-point"
"treatment"	"time_point"
"treatment"	"ega-disease"
"treatment"	"time_post_infection"
"treatment"	"timepoint_postchallenge"
"treatment"	"training"
"treatment"	"training_timepoint"
"treatment"	"treatment_description"
"treatment"	"tumor"
"treatment"	"tumor_grading"
"treatment"	"vaccine"
"treatment"	"vaccine_administration"
"treatment"	"vaccine_group"
"treatment"	"viremia_level"
"treatment"	"virus"
"treatment"	"virus_dose"
"treatment"	"d"
"treatment"	"days_post-infection"
"treatment"	"ega-subject"
"treatment"	"ega-treatment"
"treatment"	"ercc_pool"
"treatment"	"host_taxid"
"treatment"	"lymphocyte"
"treatment"	"monocyte"
"treatment"	"morphology"
"treatment"	"neutrophil"
"treatment"	"platelet_count"
"treatment"	"rfi_value"
"treatment"	"sample_type"
"treatment"	"source_of_female"
"treatment"	"sperm_dna_fragmentation_index"
"treatment"	"srek_ligation"
"treatment"	"time/hpf"
"treatment"	"total_number_born"
"treatment"	"treatments"
"treatment"	"type"
"treatment"	"weight"
"treatment"	"arrayexpress-compound"
"treatment"	"arrayexpress-dose"
"treatment"	"assay_type"
"treatment"	"backfat_thickness_of_live(mm)"
"treatment"	"basophil"
"treatment"	"photorecptors_remaining"
"treatment"	"productivity"
"treatment"	"hours/days_after_injury"
"treatment"	"disease_state"
"treatment"	"mehg_exposure"
"treatment"	"operation"
"treatment"	"post-injury_time_point"
"treatment"	"embryo_source"
"treatment"	"earlobe_color"
"treatment"	"reproductive_status"
"treatment"	"specimen_size"
"treatment"	"specimen_volume"
"treatment"	"specimen_weight"
"treatment"	"liver_phenotype"
"treatment"	"clinical_information"
"treatment"	"oestrous_peroid"
"treatment"	"secondary_description"
"treatment"	"finishing_diet"
"treatment"	"gestation_duration"
"treatment"	"growth_rate"
"treatment"	"heart_weight"
"treatment"	"information"
"treatment"	"light_treatment"
"treatment"	"maternal_food"
"treatment"	"number_of_passages"
"treatment"	"paternal_experience"
"treatment"	"personal_experience"
"treatment"	"physiological_conditions"
"treatment"	"rfi_group"
"treatment"	"sample_role"
"treatment"	"slaughter_weight"
"treatment"	"target"
"treatment"	"target_source"
"treatment"	"temperature_regimen"





================================================
FILE: util/batch_curate.sh
================================================
#!/bin/bash
#SBATCH --ntasks-per-node=2
#SBATCH --time=1:00:00

# Run the analysis
#wait 1

Rscript ${SCRIPTPATH} ${QUANT} ${META} ${WORK} ${LEN} ${DIST} '0' ${CUT} ${INTER} ${TISSUES}




================================================
FILE: .github/workflows/tag_from_version.yml
================================================
name: Create GitHub Tag from __version__

on:
  push:
    branches:
      - main
  workflow_dispatch:

jobs:
  create_tag:
    runs-on: ubuntu-latest

    steps:
      - name: Checkout code
        uses: actions/checkout@v3

      - name: Detect version using AST
        id: get_version
        run: |
          echo "🔍 Searching for __init__.py that contains __version__..."

          # Search for the correct __init__.py
          INIT_PATH=$(grep -rl '__version__' . --include='__init__.py' | head -n1)

          if [ -z "$INIT_PATH" ]; then
            echo "❌ __version__ not found in any __init__.py"
            exit 1
          fi

          echo "✅ Found: $INIT_PATH"

          # Create Python script line by line
          echo 'import ast'                         >  get_version.py
          echo 'import sys'                        >> get_version.py
          echo 'with open(sys.argv[1]) as f:'      >> get_version.py
          echo '    tree = ast.parse(f.read())'    >> get_version.py
          echo 'for node in ast.walk(tree):'       >> get_version.py
          echo '    if isinstance(node, ast.Assign):' >> get_version.py
          echo '        if getattr(node.targets[0], "id", None) == "__version__":' >> get_version.py
          echo '            print(node.value.s)'   >> get_version.py
          echo '            break'                 >> get_version.py

          # Extract version
          VERSION=$(python3 get_version.py "$INIT_PATH")
          PACKAGE=$(basename $(dirname "$INIT_PATH"))

          if [ -z "$VERSION" ]; then
            echo "❌ Failed to extract __version__ from $INIT_PATH"
            exit 1
          fi

          echo "📦 PACKAGE: $PACKAGE"
          echo "📌 VERSION: $VERSION"

          echo "PACKAGE=$PACKAGE" >> $GITHUB_ENV
          echo "VERSION=$VERSION" >> $GITHUB_ENV

      - name: Create tag if not exists
        run: |
          git fetch --tags
          if git rev-parse "v$VERSION" >/dev/null 2>&1; then
            echo "✅ Tag v$VERSION already exists. Skipping."
          else
            echo "🏷 Creating new tag v$VERSION"
            git config user.name "github-actions"
            git config user.email "actions@github.com"
            git tag "v$VERSION"
            git push origin "v$VERSION"
          fi

