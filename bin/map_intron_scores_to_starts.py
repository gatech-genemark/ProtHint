#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2020, Georgia Institute of Technology, USA
#
# For each start codon, save information about scores of the intron which
# borders the initial exon in which the start is located. If the intron is not
# present in prothint.gff, "-" is returned in place of the scores.
# ==============================================================


import argparse
import csv
import re
import sys


class Intron():
    def __init__(self, IBA, IMC):
        self.IMC = IMC
        self.IBA = IBA


def getSignature(row):
    return row[0] + "_" + row[3] + "_" + row[4] + "_" + row[6]


def extractFeatureGff(text, feature):
    regex = feature + '=([^;]+);'
    return re.search(regex, text).groups()[0]


def loadIntrons(args):
    introns = {}
    for row in csv.reader(open(args.introns), delimiter='\t'):
        if len(row) != 9:
            continue

        if row[2].lower() != "intron":
            continue

        signature = getSignature(row)

        IBA = extractFeatureGff(row[8], "al_score")
        IMC = row[5]

        introns[signature] = Intron(IBA, IMC)

    return introns


def getNextIntronCoordinates(row):
    nextIntron = extractFeatureGff(row[8], "nextIntron")
    if nextIntron == '-':
        return "-", "-"
    nextIntron = nextIntron.split("-")
    return int(nextIntron[0]), int(nextIntron[1])


def checkIntronPresence(iStart, iEnd, row, introns):
    if iStart == "-":
        return "-", "-"

    nextIntronSignature = row[0] + "_"
    if row[6] == '+':
        nextIntronSignature += str(int(row[3]) + iStart) + "_" + \
                               str(int(row[3]) + iEnd) + "_"
    elif row[6] == '-':
        nextIntronSignature += str(int(row[4]) - iEnd) + "_" + \
                               str(int(row[4]) - iStart) + "_"
    else:
        sys.exit("Error: Unexpected strand")
    nextIntronSignature += row[6]

    if nextIntronSignature in introns:
        return introns[nextIntronSignature].IBA, introns[nextIntronSignature].IMC

    return "-", "-"


def assignIntronScores(starts, introns, args):
    for row in csv.reader(open(starts), delimiter='\t'):
        if len(row) != 9:
            print("\t".join(row))
            continue

        if row[2].lower() != "start_codon":
            print("\t".join(row))
            continue

        iStart, iEnd = getNextIntronCoordinates(row)

        IBA, IMC = checkIntronPresence(iStart, iEnd, row, introns)

        print("\t".join(row) + " nextIntronIBA=" + IBA +
              "; nextIntronIMC=" + IMC + ";")


def main():
    args = parseCmd()
    introns = loadIntrons(args)
    assignIntronScores(args.starts, introns, args)


def parseCmd():

    parser = argparse.ArgumentParser(description='For each start codon, save \
        information about scores of the intron which borders the initial exon\
        in which the start is located. If the intron is not present in\
        prothint.gff, "-" is returned in place of the scores.')

    parser.add_argument('starts', metavar='starts.gff', type=str,
                        help='File with start codons, coming from the "spaln\
                              boundary scorer" program.')

    parser.add_argument('introns', metavar='prothint.gff', type=str,
                        help='Scored prothint introns.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
