#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2021, Georgia Institute of Technology, USA
#
# Description
# ==============================================================


import argparse
import csv
import sys
import math
import subprocess


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


def systemCall(cmd):
    if subprocess.call(["bash", "-c", cmd]) != 0:
        sys.exit('error: Program exited due to an error in command: ' + cmd)


def getSeedString(row):
    return row[0] + "-" + row[1] + "-" + row[7]


def loadBorders(cluster):
    borders = []
    for row in csv.reader(open(cluster), delimiter='\t'):
        borders.append(Border(int(row[3]), True))
        borders.append(Border(int(row[4]), False))

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

    for block in blocks:
        if state == "start":
            if block.coverage <= lowThresholdInt:
                state = "edge"
            else:
                state = "subCluster"
        elif state == "subCluster":
            if block.coverage <= lowThresholdInt:
                state = "bridge"
        elif state == "bridge":
            if block.coverage >= highThresholdInt:
                state = "subCluster"
        elif state == "edge":
            if block.coverage >= highThresholdInt:
                state = "subCluster"

        block.state = state

    for block in blocks:
        block.print()


def main():
    args = parseCmd()
    borders = loadBorders(args.input)
    blocks, maxCoverage, meanCoverage = computeCoverage(borders)
    print(maxCoverage, meanCoverage)
    labelBlocks(blocks, maxCoverage, args.lowThreshold, args.highThreshold)


def parseCmd():

    parser = argparse.ArgumentParser()

    parser.add_argument('input', metavar='cluster', type=str,
                        help='')

    parser.add_argument('--lowThreshold', type=float, default=0.1)
    parser.add_argument('--highThreshold', type=float, default=0.2)

    return parser.parse_args()


if __name__ == '__main__':
    main()
