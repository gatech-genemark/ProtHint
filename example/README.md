# ProtHint Usage Example

## Description

This folder contains an example input which serves to:

  * Demonstrate the usage of ProtHint
  * To test whether the program is correctly configured

The example input genome comprises of the first 3 Mbp in _Drosophila
melanogaster's_ 2L chromosome.

## Running the pipeline

To run ProtHint, use the following command

    ../bin/prothint.py input/genome.fasta input/proteins.fasta --geneMarkGtf input/genemark.gtf --workdir test

If everything is configured correctly, the results in the `test` folder should
match the contents of the `output` folder. Note that the order of hints in gff
files might differ.

The expected runtime for this example on a system with 8 CPUs (each 2.40GHz)
and 8GB of RAM is about **1.5 minutes**.
