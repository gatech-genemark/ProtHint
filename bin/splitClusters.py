#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2021, Georgia Institute of Technology, USA
#
# Description
# ==============================================================


import argparse
import csv
import tempfile
import sys
import os
import subprocess
import math
from multiprocessing import Pool


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
        
    def print(self):
        print(self.start - 1, self.end - 1, self.coverage)


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

    for border in borders:
        if border.coordinate != prevCoordinate:
            if len(blocks) != 0 and blocks[-1].coverage == coverage:
                blocks[-1].end = border.coordinate
            else:
                blocks.append(Block(prevCoordinate, border.coordinate, coverage))

        if border.start:
            coverage += 1
        else:
            coverage -= 1

        prevCoordinate = border.coordinate

    for block in blocks:
        block.print()

def main():
    args = parseCmd()
    borders = loadBorders(args.input)
    computeCoverage(borders)


def parseCmd():

    parser = argparse.ArgumentParser()

    parser.add_argument('input', metavar='cluster', type=str,
                        help='')

    return parser.parse_args()


if __name__ == '__main__':
    main()
