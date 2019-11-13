#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2019, Georgia Institute of Technology, USA
#
# Compute and print the number of CDS segments overlapping each start. CDS
# regions which start before a start codon starting coordinate and end after a
# start codon ending coordinate are considered to be overlapping. Both input
# files (start codons and CDS coordinates in gff format) need to be sorted by
# chromosome, start and end. For example like this: sort -k1,1 -k4,4n -k5,5n
# starts.gff > starts_sorted.gff. The file with CDS coordinates can be collapsed
# (with combine_gff_records.gff script) for faster execution of this script.
# ==============================================================


import csv
import argparse


class CDS:

    def __init__(self, start, end, coverage):
        self.start = start
        self.end = end
        self.coverage = coverage


def loadCDS(cdsFileName):
    codingSegments = {}
    cdsFile = open(cdsFileName)
    cdses = csv.reader(cdsFile, delimiter='\t')
    for row in cdses:
        if not row[0] in codingSegments:
            codingSegments[row[0]] = []

        coverage = 1
        if row[5] != ".":
            coverage = int(row[5])

        codingSegments[row[0]].append(CDS(int(row[3]),
                                      int(row[4]), coverage))
    cdsFile.close()
    return codingSegments


def filterStarts(startsFileName, cdsFileName):
    codingSegments = loadCDS(cdsFileName)

    startsFile = open(startsFileName)
    starts = csv.reader(startsFile, delimiter='\t')

    prevChromosome = ""
    CDSpointer = 0
    CDSNum = 0

    for start in starts:
        chrom = start[0]
        startStart = int(start[3])
        startEnd = int(start[4])

        if prevChromosome != chrom:
            CDSpointer = 0
            CDSNum = len(codingSegments[chrom])

        while (CDSpointer < CDSNum and
               codingSegments[chrom][CDSpointer].end <= startEnd):
            # Shift CDS starting point to first CDS which can overlap
            # with any of the subsequent starts
            CDSpointer += 1

        startOverlaps = 0
        i = CDSpointer
        while i < CDSNum and codingSegments[chrom][i].start < startStart:
            # While any CDS starts before the current start
            if (codingSegments[chrom][i].start < startStart
               and codingSegments[chrom][i].end > startEnd):
                startOverlaps += codingSegments[chrom][i].coverage
            i += 1

        if start[8] == ".":
            start[8] = "CDS_overlap=" + str(startOverlaps) + ";"
        else:
            start[8] += " CDS_overlap=" + str(startOverlaps) + ";"

        print("\t".join(start))

        prevChromosome = chrom

    startsFile.close()


def main():
    args = parseCmd()
    filterStarts(args.starts, args.cds)


def parseCmd():

    parser = argparse.ArgumentParser(description='Compute and print the number of CDS  \
        segments overlapping each start. CDS regions which start before a start codon \
        starting coordinate and end after a start codon ending coordinate are considered \
        to be overlapping. Both input files (start codons and CDS coordinates in gff\
        format) need to be sorted by chromosome, start and end. For example like this: \
        sort -k1,1 -k4,4n -k5,5n starts.gff > starts_sorted.gff. The file with CDS \
        coordinates can be collapsed (with combine_gff_records.gff script) for faster\
        execution of this script.')

    parser.add_argument('starts', metavar='starts.gff', type=str,
                        help='Sorted start codons in gff format.')
    parser.add_argument('cds', metavar='cds.gff', type=str,
                        help='Sorted CDS regions in gff format. If the 6th score column \
                        contains a number, this number is treated as coverage of the \
                        given CDS region.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
