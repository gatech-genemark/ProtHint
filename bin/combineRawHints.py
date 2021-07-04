#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2021, Georgia Institute of Technology, USA
#
# Description
# ==============================================================


import argparse
import csv
import re


class Hint(object):

    def __init__(self, row):
        self.row = row[0:8]
        al_score = extractFeatureGff(row[8], "al_score")
        if al_score:
            self.al_score = float(al_score)
        else:
            self.al_score = -1

        self.topProt = extractFeatureGff(row[8], "topProt")
        self.splice_sites = extractFeatureGff(row[8], "splice_sites")

        if row[5] == ".":
            self.count = 1
        else:
            self.count = int(row[5])

    def update(self, row):
        al_score = extractFeatureGff(row[8], "al_score")
        if al_score and float(al_score) > self.al_score:
            self.al_score = float(al_score)

        if self.topProt is None:
            self.topProt = extractFeatureGff(row[8], "topProt")

        if self.splice_sites is None:
            self.splice_sites = extractFeatureGff(row[8], "splice_sites")

        if row[5] == ".":
            self.count += 1
        else:
            self.count += int(row[5])

    def print(self):
        self.row[5] = str(self.count)
        extraFeatures = "."
        if self.al_score != -1:
            al_score = str(self.al_score)
            # This is for compatibility with an older script
            if al_score == "0.0":
                al_score = "0"
            extraFeatures = self.addFeatureToString(extraFeatures,
                                                    "al_score=" + al_score)

        if self.splice_sites:
            extraFeatures = self.addFeatureToString(extraFeatures,
                                                    "splice_sites=" +
                                                    self.splice_sites)

        if self.topProt:
            extraFeatures = self.addFeatureToString(extraFeatures,
                                                    "topProt=" +
                                                    self.topProt)

        print("\t".join(self.row) + "\t" + extraFeatures)

    def addFeatureToString(self, featureString, newFeature):
        if featureString == ".":
            return str(newFeature) + ";"
        else:
            return featureString + " " + newFeature + ";"


def extractFeatureGff(text, feature):
    regex = "(^|[ ;])" + feature + '=([^;]+);'
    match = re.search(regex, text)
    if match:
        return match.groups()[1]
    else:
        return None


def getSignature(row):
    return row[0] + "_" + row[2] + "_" + row[3] + "_" + row[4] + "_" + row[6] \
           + "_" + row[7]


def combineHints(input):
    hints = {}

    for row in csv.reader(open(input), delimiter='\t'):
        signature = getSignature(row)
        if signature not in hints:
            hints[signature] = Hint(row)
        else:
            hints[signature].update(row)

    return hints


def printHints(hints):
    for key in hints:
        hints[key].print()


def main():
    args = parseCmd()
    hints = combineHints(args.input)
    printHints(hints)


def parseCmd():

    parser = argparse.ArgumentParser(description='Description')

    parser.add_argument('input', metavar='spaln.gff', type=str,
                        help='Input raw hints.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
