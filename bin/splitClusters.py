#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2021, Georgia Institute of Technology, USA
#
# Description
# ==============================================================


import argparse
import csv
import math


class Border():

    def __init__(self, coordinate, start):
        self.coordinate = coordinate
        self.start = start
        if not self.start:
            self.coordinate += 1

    def print(self):
        print(self.coordinate, self.start)

    def __lt__(self, other):
        return self.coordinate < other.coordinate


class Block():

    def __init__(self, start, end, coverage):
        self.start = start
        self.end = end
        self.coverage = coverage
        self.state = ""
        
    def print(self):
        print(self.start - 1, self.end - 1, self.coverage, self.state)


class SubCluster():

    def __init__(self, start, index):
        self.start = start
        self.index = index
        self.seeds = []

    def close(self, end):
        self.end = end

    def addSeed(self, seed):
        self.seeds.append(seed)

    def print(self):
        print(self.start, self.end, self.index)

    def printToFile(self, output):
        chrom = self.seeds[0].chrom
        strand = self.seeds[0].strand
        clusterID = self.seeds[0].clusterID
        output.write("\t".join([chrom, "DIAMOND", "CDS", str(self.start),
                     str(self.end), "1", strand, "0", "gene_id \"" +
                     clusterID + "\";" + " transcript_id \"" + clusterID +
                     "\";"]) + "\n")


class Seed():

    def __init__(self, chrom, protein, start, end, score, strand, ID):
        self.chrom = chrom
        self.protein = protein
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.clusterID = ID

    def updateSeed(self, start, end, score):
        self.start = min(self.start, start)
        self.end = max(self.end, end)
        self.score += score

    def print(self):
        print(self.start, self.end, self.score)

    def __lt__(self, other):
        return self.score < other.score


def loadSeeds(clusterFile):
    seeds = {}
    CDSBorders = []
    for row in csv.reader(open(clusterFile), delimiter='\t'):
        seedID = row[8]
        if seedID not in seeds:
            seeds[seedID] = Seed(row[0], row[1], int(row[2]), int(row[3]),
                                 float(row[6]), row[7], row[9])
        else:
            seeds[seedID].updateSeed(int(row[2]), int(row[3]), float(row[6]))

        # CDSBorders.append(Border(int(row[2]), True))
        # CDSBorders.append(Border(int(row[3]), False))

    # CDSBorders.sort()

    return seeds, CDSBorders


def makeSeedBorders(seeds):
    borders = []
    for seed in seeds:
        borders.append(Border(seeds[seed].start, True))
        borders.append(Border(seeds[seed].end, False))

    borders.sort()
    return borders


def computeCoverage(borders):
    coverage = 0
    blocks = []
    prevCoordinate = borders[0].coordinate
    maxCoverage = 0
    meanCoverage = 0
    clusterLength = borders[-1].coordinate - borders[0].coordinate

    for border in borders:
        if border.coordinate != prevCoordinate:
            meanCoverage += ((border.coordinate - prevCoordinate)
                             / clusterLength) * coverage
            if len(blocks) != 0 and blocks[-1].coverage == coverage:
                blocks[-1].end = border.coordinate
            else:
                blocks.append(Block(prevCoordinate, border.coordinate,
                                    coverage))

                if coverage > maxCoverage:
                    maxCoverage = coverage

        if border.start:
            coverage += 1
        else:
            coverage -= 1

        prevCoordinate = border.coordinate

    return blocks, maxCoverage, meanCoverage


def labelBlocks(blocks, maxCoverage, lowThreshold, highThreshold):
    state = "start"
    lowThresholdInt = int(lowThreshold * maxCoverage)
    highThresholdInt = math.ceil(highThreshold * maxCoverage)
    subClusters = []

    for block in blocks:
        if state == "start":
            # Always start the first subcluster from the beginning because a
            # bridge needs to connect two subclusters
            subClusters.append(SubCluster(block.start, len(subClusters) + 1))
            if block.coverage <= lowThresholdInt:
                state = "edge"
            else:
                state = "subCluster"
        elif state == "subCluster":
            if block.coverage <= lowThresholdInt:
                state = "bridge"
                subClusters[-1].close(block.start - 1)
        elif state == "bridge":
            if block.coverage >= highThresholdInt:
                state = "subCluster"
                subClusters.append(SubCluster(block.start,
                                   len(subClusters) + 1))
        elif state == "edge":
            if block.coverage >= highThresholdInt:
                state = "subCluster"

        block.state = state

    # Close the last subcluster.
    subClusters[-1].close(blocks[-1].end - 1)

    for block in blocks:
        block.print()

    for subcluster in subClusters:
        subcluster.print()

    return subClusters


def assignSeedsToSubsclusters(subClusters, seeds):
    MARGIN = 100

    for subcluster in subClusters:
        print(len(seeds))
        for key in list(seeds.keys()):
            seed = seeds[key]
            if seed.start < subcluster.start - MARGIN:
                # Seeds starts before this cluster, it cannot be inside
                # any subsequent subcluster
                del seeds[key]
            elif seed.start > subcluster.end:
                # Seed starts after this subcluster, deal with it later
                continue
            else:
                # Seed starts in this subcluster. Add it if it also ends in
                # it. The seed cannot be inside any subsequent subcluster.
                if seed.end <= subcluster.end + MARGIN:
                    subcluster.addSeed(seed)
                del seeds[key]


def printSubClusters(output, subClusters):

    output = open(output, "a")

    for subCluster in subClusters:
        subCluster.printToFile(output)

    output.close()


def printPairs(output, subClusters):

    output = open(output, "a")

    for subCluster in subClusters:
        for seed in subCluster.seeds:
            output.write("\t".join([seed.clusterID,
                                    seed.protein,
                                    str(round(seed.score, 2))]) + "\n")

    output.close()


def main():
    args = parseCmd()

    seeds, CDSBorders = loadSeeds(args.input)
    seedBorders = makeSeedBorders(seeds)
    blocks, maxCoverage, meanCoverage = computeCoverage(seedBorders)

    subClusters = labelBlocks(blocks, meanCoverage,
                              args.lowThreshold, args.highThreshold)

    assignSeedsToSubsclusters(subClusters, seeds)
    printSubClusters(args.seedRegions, subClusters)
    printPairs(args.alignmentPairs, subClusters)


def parseCmd():

    parser = argparse.ArgumentParser()

    parser.add_argument('input', metavar='cluster', type=str,
                        help='')

    parser.add_argument('--lowThreshold', type=float, default=0.1)
    parser.add_argument('--highThreshold', type=float, default=0.2)

    parser.add_argument('--seedRegions', type=str, required=True,
                        help='Output file for seed regions.')

    parser.add_argument('--alignmentPairs', type=str, required=True,
                        help='Output with a list of up to --maxProteinsPerSeed\
        best scoring proteins per each seed region.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
