#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2019, Georgia Institute of Technology, USA
#
# Flag hints which were mapped from the best DIAMOND protein hit
# ==============================================================


import argparse
import csv
import re
import sys


def extractFeatureGff(text, feature):
    regex = feature + '=([^;]+);'
    return re.search(regex, text).groups()[0]


def loadTopPairs(diamondPairs):
    topPairs = set()
    allPairs = set()
    prevGene = ""
    for row in csv.reader(open(diamondPairs), delimiter='\t'):
        gene = row[0]
        allPairs.add(row[0] + "_" + row[1])
        if gene != prevGene:
            topPairs.add(row[0] + "_" + row[1])
        prevGene = gene
    return topPairs, allPairs


def flagHints(hints, topPairs, allPairs):
    for row in csv.reader(open(hints), delimiter='\t'):
        hintProt = extractFeatureGff(row[8], "prot")
        seedGene = extractFeatureGff(row[8], "seed_gene_id")
        key = seedGene + "_" + hintProt

        if key not in allPairs:
            sys.exit('error: Gene-protein pair "' + seedGene + "-" +
                     hintProt + '" present in the Spaln output was not found '
                     'in the file with DIAMOND gene-protein pairs. This issue '
                     'can be caused by the presence of special characters in '
                     'the fasta headers of input files. Please remove any '
                     'special characters and re-run ProtHint. '
                     'See https://github.com/gatech-genemark/ProtHint#input '
                     'for more details about the input format.')

        if key in topPairs:
            row[8] += " topProt=TRUE;"
        print("\t".join(row))


def main():
    args = parseCmd()
    topPairs, allPairs = loadTopPairs(args.diamondPairs)
    flagHints(args.hints, topPairs, allPairs)


def parseCmd():

    parser = argparse.ArgumentParser(description='Flag hints which were mapped \
                                     from the best DIAMOND protein hit')

    parser.add_argument('hints', metavar='hints.gff', type=str,
                        help='File with hints, source gene-protein pair \
                        is assumed to be in the 9th colum, labeled \
                        "seed_gene_id" and "prot", respectively.')

    parser.add_argument('diamondPairs', metavar='diamondPairs', type=str,
                        help='File with DIAMOND pairs. File should be sorted \
                        by query protein ID in column 1 and by target protein \
                        score (from best to worst) in column 2.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
