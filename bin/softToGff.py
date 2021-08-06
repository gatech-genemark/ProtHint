#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2021, Georgia Institute of Technology, USA
#
# Get masking coordinates from a soft-masked genome. The output is sorted by
# contig names.
# ==============================================================


import argparse
import sys
import tempfile
import os


def getContigCoordinates(genome):
    """Extract repeat coordinates for each contig

    Args:
        genome (file): Soft-masked genome

    Returns:
        folder: File handle of the temporary folder with repeat coordinates
                for each contig.
    """
    contig = ""
    coordinate = 1
    repeatState = False
    genomeFile = open(genome, "r")
    tempDir = tempfile.TemporaryDirectory(dir=".")
    outputFile = None
    for line in genomeFile:
        if len(line) == 1:
            # Skip empty lines
            continue

        if line[0] == ">":
            if repeatState:
                # Close the repeat from the previous contig
                outputFile.write("\t".join([str(coordinate - 1), ".",
                                 ".", ".", "."]) + "\n")
            if contig != "":
                outputFile.close()

            contig = line.split()[0][1:]
            coordinate = 1
            repeatState = False
            outputFile = open(tempDir.name + "/" + contig, "w")
        else:
            if contig == "":
                sys.exit("error: the fasta file contains sequence outside of" +
                         " contigs.")

            for nt in line.strip():
                if nt.islower():
                    if not repeatState:
                        # Open a repeat
                        outputFile.write("\t".join([contig, "soft_masking",
                                         "repeat", str(coordinate)]) + "\t")
                        repeatState = True
                else:
                    if repeatState:
                        # Close the repeat
                        outputFile.write("\t".join([str(coordinate - 1), ".",
                                         ".", ".", "."]) + "\n")
                        repeatState = False
                coordinate += 1

    if repeatState:
        # Close the repeat at the file end
        outputFile.write("\t".join([str(coordinate - 1), ".",
                         ".", ".", "."]) + "\n")
        if contig != "":
            outputFile.close()

    return tempDir


def printSortedCoordinates(tempDir):
    """Print repeat coordinates, sorted by contig names
    """
    contigs = os.listdir(tempDir.name)
    contigs.sort()
    for contig in contigs:
        with open(tempDir.name + "/" + contig) as contigOut:
            for line in contigOut:
                print(line, end="")


def main():
    args = parseCmd()
    tempDir = getContigCoordinates(args.genome)
    printSortedCoordinates(tempDir)


def parseCmd():

    parser = argparse.ArgumentParser(description='Get masking coordinates \
        from a soft-masked genome. The output is sorted by contig names.')

    parser.add_argument('genome', metavar='genome.fasta.masked', type=str,
                        help='Masked genome')

    return parser.parse_args()


if __name__ == '__main__':
    main()
