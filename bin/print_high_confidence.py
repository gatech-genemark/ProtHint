#!/usr/bin/env python3
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
    regex = ";\s*" + feature + '=([^;]+)'
    text = "; " + text
    search = re.search(regex, text)
    if search:
        return search.groups()[0]
    else:
        return None


class Filter:

    def __init__(self, args):
        self.args = args

    def decide(self, row):
        self.al_score = extractFeature(row[8], "al_score")
        if self.al_score:
            self.al_score = float(self.al_score)
        self.fullProtein = extractFeature(row[8], "fullProteinAligned")
        self.topProtein = extractFeature(row[8], "topProt")
        self.row = row
        self.coverage = int(row[5])

        self.__determineCoverageThreshod()

        if (row[2].lower() == "intron"):
            return self.__intron()
        elif (row[2].lower() == "stop_codon"):
            return self.__stop()
        elif (row[2].lower() == "start_codon"):
            return self.__start()
        else:
            return True

    def __determineCoverageThreshod(self):
        if (self.row[2].lower() == "intron"):
            self.coverageThreshold = self.args.intronCoverage
        elif (self.row[2].lower() == "stop_codon"):
            self.coverageThreshold = self.args.stopCoverage
        elif (self.row[2].lower() == "start_codon"):
            self.coverageThreshold = self.args.startCoverage

        if ((self.args.addTopProteins and self.topProtein == "TRUE") or
           (self.args.addFullAligned and self.fullProtein == "TRUE")):
            self.coverageThreshold = 1

    def __intron(self):
        if not self.args.addAllSpliceSites:
            spliceSites = extractFeature(self.row[8], "splice_sites")
            if spliceSites is not None and spliceSites.lower() != "gt_ag":
                if not self.args.addGCAG or spliceSites.lower() != "gc_ag":
                    return False

        if (self.coverage >= self.coverageThreshold and
           self.al_score >= self.args.intronAlignment):
            return True

        return False

    def __stop(self):
        coverageThreshold = self.args.stopCoverage
        if self.args.addTopProteins and self.topProtein == "TRUE":
            coverageThreshold = 1

        if (self.al_score is None):
            self.al_score = 1

        if (self.coverage >= coverageThreshold and
           self.al_score >= self.args.stopAlignment):
            return True

        return False

    def __start(self):
        if (self.al_score is None):
            self.al_score = 1
        eScore = float(extractFeature(self.row[8], "eScore"))

        self.CDS_overlap = extractFeature(self.row[8], "CDS_overlap")
        if (self.CDS_overlap is None):
            self.CDS_overlap = 0
        else:
            self.CDS_overlap = int(self.CDS_overlap)

        if eScore >= self.args.shortStartThreshold:
            return self.__regularStart()
        else:
            return self.__shortStart()

    def __regularStart(self):
        coverageThreshold = self.args.startCoverage
        if self.args.addTopProteins and self.topProtein == "TRUE":
            coverageThreshold = 1

        if (self.coverage >= coverageThreshold and
           self.al_score >= self.args.startAlignment and
           self.CDS_overlap <= self.args.startOverlap):
            return True

        return False

    def __shortStart(self):
        coverageThreshold = self.args.shortStartCoverage
        if self.args.addTopProteins and self.topProtein == "TRUE":
            coverageThreshold = 1

        intronIBA = extractFeature(self.row[8], "nextIntronIBA")
        if (intronIBA == "-"):
            intronIBA = 0
        else:
            intronIBA = float(intronIBA)

        intronIMC = extractFeature(self.row[8], "nextIntronIMC")
        if (intronIMC == "-"):
            intronIMC = 0
        else:
            intronIMC = int(intronIMC)

        if intronIBA < self.args.intronAlignment or \
           intronIMC < self.args.intronCoverage:
            return False

        if (self.coverage >= coverageThreshold and
           self.al_score >= self.args.startAlignment and
           self.CDS_overlap <= self.args.shortStartOverlap):
            return True

        return False


def printHighConfidence(args):
    filter = Filter(args)
    for row in csv.reader(open(args.input), delimiter='\t'):
        row[1] = "ProtHint"

        if row[5] == ".":
            row[5] = "1"

        if filter.decide(row):
            print("\t".join(row))


def main():
    args = parseCmd()
    printHighConfidence(args)


def parseCmd():

    parser = argparse.ArgumentParser(description='Select and print high confidence features\
                                     from ProtHint output file.')

    parser.add_argument('input', metavar='prothint.gff', type=str,
                        help='ProtHint output file.')

    parser.add_argument('--shortStartThreshold', type=int,
                        help='Exon score threshold for treating starts as starts \
                        in short initial exons. Different scoring thresholds are \
                        applied to such low-scoring starts (starts with eScore < \
                        shortStartThreshold). On top of the different scoring \
                        thresholds, introns bordering the initial exons in \
                        which the low-scoring starts are must pass the intron \
                        filtering criteria for intronAlignment and intronCoverage. \
                        Otherise, the "short" starts are not printed. \
                        Default = 25.', default=25)
    parser.add_argument('--intronCoverage', type=int,
                        help='Intron coverage score threshold. Print all introns \
                        with coverage >= intronCoverage. Default = 4.', default=4)
    parser.add_argument('--startCoverage', type=int,
                        help='Start coverage score threshold. Print all regular \
                        starts with coverage >= startCoverage. Default = 4.', default=4)
    parser.add_argument('--shortStartCoverage', type=int,
                        help='Start coverage score threshold for starts in short \
                        initial exons, i. e. starts which have eScore < \
                        shortStartThreshold. Print all "short" starts with \
                        coverage >= shortStartCoverage. Default = 4.', default=4)
    parser.add_argument('--stopCoverage', type=int,
                        help='Stop coverage score threshold. Print all stops \
                        with coverage >= stopCoverage. Default = 4.', default=4)
    parser.add_argument('--startAlignment', type=float,
                        help='Start alignment score threshold. Print all starts \
                        with al_score >= startAlignment. Ignore if the field \
                        does not exist. Default = 0.01.', default=0.01)
    parser.add_argument('--stopAlignment', type=float,
                        help='Stop alignment score threshold. Print all stops \
                        with al_score >= stopAlignment. Ignore if the field does not exist. \
                        Default = 0.01.', default=0.01)
    parser.add_argument('--intronAlignment', type=float,
                        help='Intron alignment score threshold. Print all introns \
                        with al_score >= intronAlignment. Default = 0.25.', default=0.25)
    parser.add_argument('--startOverlap', type=int,
                        help='Maximum allowed CDS overlap of a regular start \
                        (eScore >= shortStartThreshold). Print all starts with \
                        CDS overlap <= startOverlap. Default = 0', default=0)
    parser.add_argument('--shortStartOverlap', type=int,
                        help='Maximum allowed CDS overlap of a start in a short \
                        initial exon (starts with eScore < shortStartThreshold). \
                        Print all "short" starts with CDS overlap <= shortStartOverlap. \
                        Default = 0', default=0)
    parser.add_argument('--addFullAligned', action='store_true',
                        help='Add hints with fullProteinAligned flag even if they do not \
                        satisfy the coverage threshold condition.')
    parser.add_argument('--addGCAG', action='store_true',
                        help='Add introns with GC_AG splice sites. By default, \
                        only introns with canonical GT_AG splice sites are printed.')
    parser.add_argument('--addAllSpliceSites', action='store_true',
                        help='Add introns with any splice sites.  By default, \
                        only introns with canonical GT_AG splice sites are printed')
    parser.add_argument('--addTopProteins', action='store_true',
                        help='Add hints corresponding to the top protein, no matter \
                        the coverage. Other scoring thresholds still apply.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
