#!/usr/bin/env python
# Author: Tomas Bruna
#
# Select high confidence features from ProtHint output file.
#
# positional arguments:
#   prothint.gff          ProtHint output file.
#
# optional arguments:
#   --intronCoverage INTRONCOVERAGE
#                         Intron coverage score threshold. Print all introns
#                         with coverage >= intronCoverage. Default = 4.
#   --startCoverage STARTCOVERAGE
#                         Start coverage score threshold. Print all starts with
#                         coverage >= startCoverage. Default = 4.
#   --stopCoverage STOPCOVERAGE
#                         Stop coverage score threshold. Print all stops with
#                         coverage >= stopCoverage. Default = 4.
#   --intronAlignment INTRONALIGNMENT
#                         Intron alignment score threshold. Print all introns
#                         with al_score >= intronAlignment. Default = 0.3.
#   --startOverlap STARTOVERLAP
#                         Maximum alowed CDS overlap of a start. Print all
#                         starts with CDS overlap <= startOverlap. Default = 3
#   --addFullAligned      Add hints with fullProteinAligned flag even if they do
#                         not satisfy the coverage threshold condition.


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


def intron(row, args):
    al_score = float(extractFeature(row[8], "al_score"))
    fullProtein = extractFeature(row[8], "fullProteinAligned")
    if (((int(row[5]) >= args.intronCoverage) or
        (args.addFullAligned and fullProtein == "TRUE")) and
       (al_score >= args.intronAlignment)):
        print("\t".join(row))


def stop(row, args):
    fullProtein = extractFeature(row[8], "fullProteinAligned")
    if ((int(row[5]) >= args.stopCoverage) or
       (args.addFullAligned and fullProtein == "TRUE")):
        print("\t".join(row))


def start(row, args):
    CDS_overlap = float(extractFeature(row[8], "CDS_overlap"))
    fullProtein = extractFeature(row[8], "fullProteinAligned")
    if (((int(row[5]) >= args.startCoverage) or
        (args.addFullAligned and fullProtein == "TRUE")) and
       (CDS_overlap <= args.startOverlap)):
        print("\t".join(row))


def printHighConfidence(args):
    for row in csv.reader(open(args.input), delimiter='\t'):
        if (row[2].lower() == "intron"):
            intron(row, args)
        elif (row[2].lower() == "start_codon"):
            start(row, args)
        elif (row[2].lower() == "stop_codon"):
            stop(row, args)


def main():
    args = parseCmd()
    printHighConfidence(args)


def parseCmd():

    parser = argparse.ArgumentParser(description='Select and print high confidence features\
                                     from ProtHint output file.')

    parser.add_argument('input', metavar='prothint.gff', type=str,
                        help='ProtHint output file.')

    parser.add_argument('--intronCoverage', type=int,
                        help='Intron coverage score threshold. Print all introns \
                        with coverage >= intronCoverage. Default = 4.', default=4)
    parser.add_argument('--startCoverage', type=int,
                        help='Start coverage score threshold. Print all starts \
                        with coverage >= startCoverage. Default = 4.', default=4)
    parser.add_argument('--stopCoverage', type=int,
                        help='Stop coverage score threshold. Print all stops \
                        with coverage >= stopCoverage. Default = 4.', default=4)
    parser.add_argument('--intronAlignment', type=float,
                        help='Intron alignment score threshold. Print all introns \
                        with al_score >= intronAlignment. Default = 0.3.', default=0.3)
    parser.add_argument('--startOverlap', type=int,
                        help='Maximum alowed CDS overlap of a start. Print all starts \
                        with CDS overlap <= startOverlap. Default = 3', default=3)
    parser.add_argument('--addFullAligned', action='store_true',
                        help='Add hints with fullProteinAligned flag even if they do not \
                        satisfy the coverage threshold condition.')

    return parser.parse_args()


if __name__ == '__main__':
    main()