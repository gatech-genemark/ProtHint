#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2019, Georgia Institute of Technology, USA
#
# Make chained hints for augustus from Spaln parser gff output
# ==============================================================


import argparse
import csv
import re
import sys


def extractFeatureGff(text, feature):
    regex = feature + '=([^;]+);'
    return re.search(regex, text).groups()[0]


def printChains(gffFile, cutoff):
    for row in csv.reader(open(gffFile), delimiter='\t'):
        protein = extractFeatureGff(row[8], "prot")
        gene = extractFeatureGff(row[8], "seed_gene_id")
        row[1] = "ProtHint"
        row[8] = "grp=" + protein + "_" + gene + ";src=C;pri=4;"

        if row[2] == "CDS":
            row[2] = "CDSpart"
            exonLength = int(row[4]) - int(row[3]) + 1
            if exonLength < 2 * cutoff + 3:
                row[3] = str(int(row[3]) + (int(exonLength / 6)) * 3)
                row[4] = str(int(row[3]) + 2)
            else:
                row[3] = str(int(row[3]) + cutoff)
                row[4] = str(int(row[4]) - cutoff)

        if row[2] == "Intron":
            row[2] = "intron"

        if row[2] == "start_codon":
            row[2] = "start"

        if row[2] == "stop_codon":
            row[2] = "stop"

        print("\t".join(row))


def main():
    args = parseCmd()
    printChains(args.input, args.exonCutoff)


def parseCmd():

    parser = argparse.ArgumentParser(description='Make chained hints for\
                                     augustus from Spaln parser gff output')

    parser.add_argument('input', metavar='spaln.gff', type=str,
                        help='Input gff file with raw hints')

    parser.add_argument('--exonCutoff', type=int, default=15,
                        help='Cut exon on each side by this much. Must be\
                        divisible by 3 to preserve phase. Default = 15.')

    args = parser.parse_args()

    if args.exonCutoff % 3 != 0:
        sys.exit("error: exonCutoff parameter must be divisible by 3")

    return args


if __name__ == '__main__':
    main()
