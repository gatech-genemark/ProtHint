#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2019, Georgia Institute of Technology, USA
#
# This script selects NUM_PROTEINS best proteins
# supporting each intron/start/stop.
# ==============================================================


import csv
import sys
import re


class Hint:

    def __init__(self, prot, gene, score):
        self.pair = gene + "\t" + prot
        self.score = score

    def getPair(self):
        return self.pair

    def __lt__(self, other):
        return self.score > other.score


def extractFeature(text, feature):
    regex = feature + '=([^;]+)'
    return re.search(regex, text).groups()[0]


def getSignature(row):
    return row[0] + "_" + row[2] + "_" + row[3] + "_" + row[4] + "_" + row[6]


def getBest(input, k):
    pairs = set()
    features = {}
    with open(input) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')

        for row in csv_reader:
            feature = None
            if row[2].lower() == "intron":
                feature = Hint(extractFeature(row[8], "prot"), extractFeature(row[8], "seed_gene_id"),
                               extractFeature(row[8], "al_score"))
            elif row[2].lower() == "start_codon" or row[2].lower() == "stop_codon":
                feature = Hint(extractFeature(row[8], "prot"), extractFeature(row[8], "seed_gene_id"),
                               extractFeature(row[8], "score"))
            else:
                continue

            signature = getSignature(row)

            if not signature in features:
                features[signature] = [feature]
            else:
                features[signature].append(feature)

    # Print up to best k supporting proteins for each feature
    for key, featureSet in features.items():
        n = len(featureSet)
        if n > k:
            featureSet = sorted(featureSet)
            n = k
        for i in range(0, n):
            pairs.add(featureSet[i].getPair())

    for pair in pairs:
        print(pair)


def main():
    if (len(sys.argv) != 3):
        sys.stderr.write("Error: Invalid number of arguments\n")
        sys.stderr.write("Usage: " + sys.argv[0] + " in.gff NUM_PROTEINS\n\n")
        sys.stderr.write("Description: This script selects NUM_PROTEINS best proteins \n" +
                         "             supporting each intron/start/stop. \n")
        return 1
    getBest(sys.argv[1], int(sys.argv[2]))


if __name__ == '__main__':
    main()
