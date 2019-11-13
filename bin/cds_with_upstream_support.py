#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2019, Georgia Institute of Technology, USA
#
# This script selects a subset of CDS regions which upstream coordinate is
# either a start codon in starts.gff or neighbors an intron defined in
# starts.gff.
# ==============================================================


import csv
import argparse


intronEnds = set()
starts = set()


def loadHints(hints):
    for row in csv.reader(open(hints), delimiter='\t'):
        if row[2].lower() == "intron":
            if (row[6] == "+"):
                intronEnds.add(int(row[4]))
            elif row[6] == "-":
                intronEnds.add(int(row[3]))
        elif row[2].lower() == "start_codon":
            if (row[6] == "+"):
                starts.add(int(row[3]))
            elif(row[6] == "-"):
                starts.add(int(row[4]))

    return intronEnds, starts


def filterCDS(cds):
    for row in csv.reader(open(cds), delimiter='\t'):
        if row[2].lower() == "cds":
            if (row[6] == "+"):
                cdsStart = int(row[3])
                if ((cdsStart - 1) in intronEnds) or (cdsStart in starts):
                    print("\t".join(row))
            elif row[6] == "-":
                cdsStart = int(row[4])
                if ((cdsStart + 1) in intronEnds) or (cdsStart in starts):
                    print("\t".join(row))


def main():
    args = parseCmd()
    loadHints(args.starts)
    loadHints(args.introns)
    filterCDS(args.cds)


def parseCmd():

    parser = argparse.ArgumentParser(description='This script selects a subset \
                                     of CDS regions which upstream coordinate \
                                     is either a start codon in starts.gff or \
                                     neighbors an intron defined in starts.gff.')

    parser.add_argument('cds', metavar='cds.gff', type=str,
                        help='CDS segments to be filtered')
    parser.add_argument('starts', metavar='starts.gff', type=str,
                        help='Start codons used in filtering.')
    parser.add_argument('introns', metavar='introns.gff', type=str,
                        help='Introns used in filtering.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
