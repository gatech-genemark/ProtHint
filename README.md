# ProtHint

Tomas Bruna, Alexandre Lomsadze, Mark Borodovsky

Georgia Institute of Technology, Atlanta, Georgia, USA

Reference: [GeneMark-EP+: eukaryotic gene prediction with self-training in the space of genes and proteins](https://academic.oup.com/nargab/article/2/2/lqaa026/5836691)

# Overview

ProtHint is a pipeline for predicting and scoring hints (in the form of
introns, start and stop codons) in the genome of interest by mapping and
spliced aligning predicted genes to a database of reference protein sequences.

# Contents

* [Installation](#installation)
    * [Perl dependencies](#perl-dependencies)
    * [Python dependencies](#python-dependencies)
    * [GeneMark-ES](#genemark-es)
    * [DIAMOND](#diamond)
    * [Spaln](#spaln)
        * [Spaln boundary scorer](#spaln-boundary-scorer)
    * [ProSplign](#prosplign-optional-installation)
        * [ProSplign intron scorer](#prosplign-intron-scorer)
* [Usage](#usage)
    * [Input](#input)
    * [Protein Database Preparation](#protein-database-preparation)
    * [Running ProtHint](#running-prothint)
    * [Output](#output)
* [About](#about)


# Installation

To install, copy the content of this distribution to desired location. To verify
the installation, run ProtHint with the sample input located in the [example](example)
folder.

Running ProtHint requires a Linux system with Bash. The following dependencies
need to be satisfied.

### Perl dependencies

Perl 5.10 or higher is required.

The following non-Core Perl modules are required:

* `MCE::Mutex`
* `threads`
* `YAML`
* `Math::Utils`

Core module `Thread::Queue` needs to be updated to a version `3.11` or higher.

These modules are available at CPAN and can be installed/updated with

    cpan MCE::Mutex threads YAML Thread::Queue Math::Utils

### Python dependencies

Python 3.3 or higher is required. No libraries outside of the Python Standard
Library are required.

### GeneMark-ES

There are two ways of using GeneMark-ES in ProtHint:

1.  Run ProtHint with `--geneMarkGtf genemark.gtf` option which specifies the
    path to a file with GeneMark-ES predictions. If this option is used,
    GeneMark-ES does not need to be installed as a part of ProtHint.

2.  Install GeneMark-ES, ProtHint will run it automatically as a part of the pipeline.

    Download and extract the contents of the GeneMark-ES suite (versions 4.30 and
    up) into the `ProtHint/dependencies/GeneMarkES` folder. GeneMark-ES suite is
    available at http://exon.gatech.edu/GeneMark/license_download.cgi
 
    To verify that GeneMark-ES is installed correctly, run the following command:
    `ProtHint/dependencies/GeneMarkES/check_install.bash`.


### DIAMOND

DIAMOND local sequence aligner (available at
https://github.com/bbuchfink/diamond) is included in this distribution
package.

In case the included version is not working, install DIAMOND from
https://github.com/bbuchfink/diamond and replace the `diamond` binary in
ProtHint/dependencies folder.

### Spaln

Spaln, space-efficient spliced alignment program (available at
https://github.com/ogotoh/spaln),  is included in this distribution package.

In case the included version is not working, install Spaln from
https://github.com/ogotoh/spaln and replace the `spaln` binary in
ProtHint/dependencies folder.

#### Spaln boundary scorer

Binary for parsing and scoring hints from Spaln's alignment output is included
in this distribution package.

In case the included binary is not working, compile it from source at https://github.com/gatech-genemark/spaln-boundary-scorer and replace the `spaln_boundary_scorer` binary
in ProtHint/dependencies folder.

### ProSplign (Optional Installation)

ProSplign program from NCBI – with the standard NCBI license
(http://www.ncbi.nlm.nih.gov/sutils/static/prosplign/prosplign.html) is
included in this distribution package.

This tool is only used when ProtHint is run with `--ProSplign` option. In
case the included version is not working, download the executable from
https://www.ncbi.nlm.nih.gov/sutils/static/prosplign/prosplign.html and
replace the `prosplign` binary in ProtHint/dependencies folder.

#### ProSplign intron scorer

Binary for parsing and scoring introns from ProSplign's alignment output is included
in this distribution package.

In case the included binary is not working, compile it from source at https://github.com/gatech-genemark/prosplign-intron-scorer and replace the `prosplign_intron_scorer` binary
in ProtHint/dependencies folder.

# Usage

## Input

ProtHint inputs consist of:

* Genomic sequence from the target species in multi-FASTA format
* Reference protein sequences in multi-FASTA format

The tool is applicable to complete as well as draft genome assemblies. Every
sequence in each multi-FASTA input needs to have a unique ID (first word of
a FASTA header is used for ID). Examples of valid FASTA headers:

    >contig10
    ID: contig10
    > seq3  genome Z
    ID: seq3
    >IV contig 25
    ID: IV


## Protein Database Preparation

We recommend to use a relevant portion of OrthoDB protein database as the
source of reference protein sequences.

For example, if your genome of interest is an insect, download arthropoda
proteins:

    wget https://v100.orthodb.org/download/odb10_arthropoda_fasta.tar.gz
    tar xvf odb10_arthropoda_fasta.tar.gz

and concatenate proteins from all species into a single file:

    cat arthropoda/Rawdata/* > proteins.fasta

For other genomes of interest, you can select the most specific OrthoDB
section from the list below and repeat the procedure desribed above.

* **Fungi**: https://v100.orthodb.org/download/odb10_fungi_fasta.tar.gz
* **Metazoa**: https://v100.orthodb.org/download/odb10_metazoa_fasta.tar.gz
    * **Arthropoda**: https://v100.orthodb.org/download/odb10_arthropoda_fasta.tar.gz
    * **Vertebrata**: https://v100.orthodb.org/download/odb10_vertebrata_fasta.tar.gz
* **Protozoa**: https://v100.orthodb.org/download/odb10_protozoa_fasta.tar.gz
* **Viridiplantae**: https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz

## Running ProtHint

To run ProtHint, use the following command:

    prothint.py genome.fasta proteins.fasta

See the [example](example) folder for a sample input and output.


To display a list of all available options, use:

    prothint.py --help

Frequently used options are:

    --workdir WORKDIR   Folder for results and temporary files. If not
                        specified, current directory is used
    --geneMarkGtf GENEMARKGTF
                        File with GeneMark-ES predictions in gtf format. If
                        this file is provided, GeneMark-ES run is skipped.
    --diamondPairs DIAMONDPAIRS
                        File with "seed gene-protein" hits generated by
                        DIAMOND. If this file is provided, DIAMOND search for
                        protein hits is skipped.


## Output

ProtHint generates two main outputs:

* `prothint.gff` Gff file with all reported hints (introns, starts and stops)
* `evidence.gff` High confidence subset of `prothint.gff` which is, for instance,
                 suitable for the GeneMark-EP Plus mode. This set is generated
                 using default thresholds in `ProtHint/bin/print_high_confidence.py`
                 script. If you wish to use different filtering criteria, re-run
                 `print_high_confidence.py` script with custom thresholds.

An output which is ready to be used in [BRAKER](https://github.com/Gaius-Augustus/BRAKER) and [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus) is also generated:

* `prothint_augustus.gff`

# About

ProtHint is developed by Tomas Bruna and Alexandre Lomsadze at [Dr. Mark
Borodovsky's Bioinformatics Lab](http://exon.gatech.edu/GeneMark/), Georgia
Institute of Technology, Atlanta, USA.
