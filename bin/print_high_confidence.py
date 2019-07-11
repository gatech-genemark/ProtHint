#!/usr/bin/env python
# ==============================================================
# Tomas Bruna
# Copyright 2019, Georgia Institute of Technology, USA
#
# Select high confidence features from ProtHint output file.
# ==============================================================


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

    if not args.addAllSpliceSites:
        spliceSites = extractFeature(row[8], "splice_sites")
        if spliceSites is not None and spliceSites.lower() != "gt_ag":
            if not args.addGCAG or spliceSites.lower() != "gc_ag":
                return

    if (((int(row[5]) >= args.intronCoverage) or
        (args.addFullAligned and fullProtein == "TRUE")) and
       (al_score >= args.intronAlignment)):
        print("\t".join(row))


def stop(row, args):
    fullProtein = extractFeature(row[8], "fullProteinAligned")
    stopScore = extractFeature(row[8], "score")
    if (stopScore is None):
        stopScore = 1
    else:
        stopScore = float(stopScore)

    if ((int(row[5]) >= args.stopCoverage and stopScore >= args.stopScore) or
       (args.addFullAligned and fullProtein == "TRUE")):
        print("\t".join(row))


def start(row, args):
    CDS_overlap = extractFeature(row[8], "CDS_overlap")
    if (CDS_overlap is None):
        CDS_overlap = 0
    else:
        CDS_overlap = float(CDS_overlap)

    fullProtein = extractFeature(row[8], "fullProteinAligned")

    startScore = extractFeature(row[8], "score")
    if (startScore is None):
        startScore = 1
    else:
        startScore = float(startScore)

    if (((int(row[5]) >= args.startCoverage and startScore >= args.startScore) or
        (args.addFullAligned and fullProtein == "TRUE")) and
       (CDS_overlap <= args.startOverlap)):
        print("\t".join(row))


def printHighConfidence(args):
    for row in csv.reader(open(args.input), delimiter='\t'):
        row[1] = "ProtHint"

        if row[5] == ".":
            row[5] = "1"

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
                        with coverage >= startCoverage. Default = 10.', default=10)
    parser.add_argument('--stopCoverage', type=int,
                        help='Stop coverage score threshold. Print all stops \
                        with coverage >= stopCoverage. Default = 4.', default=4)
    parser.add_argument('--startScore', type=float,
                        help='Start alignment score threshold. Print all starts \
                        with score >= startScore. Ignore if the field does not exist. \
                        Default = 0.01.', default=0.01)
    parser.add_argument('--stopScore', type=float,
                        help='Stop alignment score threshold. Print all stops \
                        with score >= stopScore. Ignore if the field does not exist. \
                        Default = 0.01.', default=0.01)
    parser.add_argument('--intronAlignment', type=float,
                        help='Intron alignment score threshold. Print all introns \
                        with al_score >= intronAlignment. Default = 0.3.', default=0.3)
    parser.add_argument('--startOverlap', type=int,
                        help='Maximum alowed CDS overlap of a start. Print all starts \
                        with CDS overlap <= startOverlap. Default = 0', default=0)
    parser.add_argument('--addFullAligned', action='store_true',
                        help='Add hints with fullProteinAligned flag even if they do not \
                        satisfy the coverage threshold condition.')
    parser.add_argument('--addGCAG', action='store_true',
                        help='Add introns with GC_AG splice sites. By default, \
                        only introns with canonical GT_AG splice sites are printed.')
    parser.add_argument('--addAllSpliceSites', action='store_true',
                        help='Add introns with any splice sites.  By default, \
                        only introns with canonical GT_AG splice sites are printed')

    return parser.parse_args()


if __name__ == '__main__':
    main()
