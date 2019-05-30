#!/usr/bin/env python
# Author: Tomas Bruna

# Convert ProMapp hints to Augustus compatible format. Assign flag src=M to
# reliable evidence exceeding the specified scoring thresholds.
#
# optional arguments:
#   --scoredIntrons SCOREDINTRONS
#                         Introns scored with alignment score (al_score in 9th
#                         gff column)
#   --al_score AL_SCORE   Alignment score threshold for scored introns to be
#                         considerd reliable evidence.
#   --intronCoverage INTRONCOVERAGE
#                         Intron coverage threshold for introns to be considered
#                         reliable evidence

import argparse
import csv
import re


def extractFeature(text, feature):
    regex = feature + '=([^;]+)'
    search = re.search(regex, text)
    if search:
        return search.groups()[0]
    else:
        return None


def getSignature(row):
    return row[0] + "_" + row[3] + "_" + row[4] + "_" + row[6]


def convertStartStops(startStops, startThreshold, stopThreshold):
    with open(startStops) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            if row[2].lower() == "start_codon":
                coverageThreshold = startThreshold
            elif row[2].lower() == "stop_codon":
                coverageThreshold = stopThreshold
            else:
                continue
            source = 'P'
            fullProtein = extractFeature(row[8], "fullProteinAligned")
            if int(row[5]) >= coverageThreshold or fullProtein == "TRUE":
                source = 'M'
            row[8] = "src=" + source + ";mult=" + row[5] + ";pri=4;"
            if fullProtein:
                row[8] += "fullProteinAligned=" + fullProtein + ";"
            print("\t".join(row))


def convertScoredIntrons(scoredIntrons, scoreThreshold, coverageThreshold):
    with open(scoredIntrons) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            if row[2].lower() != "intron":
                continue
            source = 'P'
            al_score = float(extractFeature(row[8], "al_score"))
            fullProtein = extractFeature(row[8], "fullProteinAligned")
            if al_score >= scoreThreshold and (
               int(row[5]) >= coverageThreshold or fullProtein == "TRUE"):
                source = 'M'
            row[2] = "intron"
            row[8] = "src=" + source + ";mult=" + row[5] + ";pri=4;al_score=" \
                     + str(al_score) + ";"
            if fullProtein:
                row[8] += "fullProteinAligned=" + fullProtein + ";"
            print("\t".join(row))


def main():
    args = parseCmd()
    if args.scoredIntrons:
        convertScoredIntrons(args.scoredIntrons, args.al_score, args.intronCoverage)
    if args.startStops:
        convertStartStops(args.startStops, args.startCoverage, args.stopCoverage)


def parseCmd():

    parser = argparse.ArgumentParser(description='Convert ProMapp hints to \
                                     Augustus compatible format. Assign flag \
                                     src=M to reliable evidence exceeding the \
                                     specified scoring thresholds.')

    parser.add_argument('--scoredIntrons', type=str,
                        help='Introns scored with alignment score (al_score in \
                        9th gff column) and coverage score in 6th gff column.')
    parser.add_argument('--al_score', type=float,
                        help='Alignment score threshold for scored introns to \
                        be considerd reliable evidence.')
    parser.add_argument('--intronCoverage', type=int,
                        help='Intron coverage threshold for introns to be \
                        considered reliable evidence')

    parser.add_argument('--startStops', type=str,
                        help='Start and stop codons with coverage score in the \
                        6th gff column.')
    parser.add_argument('--startCoverage', type=int,
                        help='Start codon coverage threshold for starts to be \
                        considered reliable evidence')
    parser.add_argument('--stopCoverage', type=int,
                        help='Stop codon coverage threshold for stops to be \
                        considered reliable evidence')

    args = parser.parse_args()
    if args.scoredIntrons and (args.al_score is None or args.intronCoverage is None):
        parser.error("--scoredIntrons requires --al_score and --intronCoverage.")
    if args.startStops and (args.startCoverage is None or args.stopCoverage is None):
        parser.error("--startStops requires --startCoverage and --stopCoverage.")

    return args


if __name__ == '__main__':
    main()
