#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2019, Georgia Institute of Technology, USA
#
# Add information about full protein alignment from second input file to
# features (possibly combined/collapsed) in the first gff file. If the first
# input file has combined features, only one of them needs to come from a fully
# aligned protein for the flag to be set in the output.
# ==============================================================


import argparse
import csv
import re


class Feature:

    def __init__(self, row):
        self.row = row
        self.fullProtein = "FALSE"


def extractFeature(text, feature):
    regex = feature + '=([^;]+)'
    return re.search(regex, text).groups()[0]


def getSignature(row):
    return row[0] + "_" + row[2].lower() + "_" + row[3] + "_" + row[4] + "_" + row[6]


def loadInput(inputFileName):
    features = {}
    inputFile = open(inputFileName)
    csvInput = csv.reader(inputFile, delimiter='\t')
    for row in csvInput:
        signature = getSignature(row)
        if signature not in features:
            features[signature] = Feature(row)
    inputFile.close()
    return features


def printOutput(features):
    for key in features:
        feature = features[key]

        if feature.row[8] == ".":
            feature.row[8] = "fullProteinAligned=" + feature.fullProtein + ";"
        else:
            feature.row[8] += " fullProteinAligned=" + feature.fullProtein + ";"

        print("\t".join(feature.row))


def mergeGff(inputFileName, asnInputFileName):
    features = loadInput(inputFileName)

    asnInputFile = open(asnInputFileName)
    csvAsnInput = csv.reader(asnInputFile, delimiter='\t')
    for row in csvAsnInput:
        signature = getSignature(row)
        fullProtein = extractFeature(row[8], "fullProteinAligned")
        if fullProtein == "TRUE" and signature in features:
            features[signature].fullProtein = "TRUE"
    asnInputFile.close()

    printOutput(features)


def main():
    args = parseCmd()
    mergeGff(args.input, args.asn)


def parseCmd():

    parser = argparse.ArgumentParser(description='Add information about \
        full protein alignment from second input file to features (possibly \
        combined/collapsed) in the first gff file. If the first input file has \
        combined features, only one of them needs to come from a fully \
        aligned protein for the flag to be set in the output.\n')

    parser.add_argument('input', metavar='input.gff', type=str,
                        help='Input file to which the information is added')
    parser.add_argument('asn', metavar='prosplign.gff', type=str,
                        help='ProSplign output with fullProteinAligned flag for \
                        each entry')

    return parser.parse_args()


if __name__ == '__main__':
    main()
