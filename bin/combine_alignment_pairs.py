#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2019, Georgia Institute of Technology, USA
#
# Add top k pairs from diamond seed protein-gene pairs to the Spaln pairs
# ==============================================================

import argparse
import sys


def readSpalnPairs(args):
    spalnPairs = set()
    spalnFile = open(args.spalnPairs, "r")
    for line in spalnFile:
        spalnPairs.add(line)
    spalnFile.close()
    return spalnPairs

def addDiamondPairs(args, outPairs):
    diamondPairs = set()
    diamondFile = open(args.diamondPairs, "r")

    prevGene = ""
    protCounter = 0

    for line in diamondFile:
        gene = line.split("\t")[0]
        if gene != prevGene:
            protCounter = 1
        if protCounter <= args.k:
            outPairs.add(line)

        prevGene = gene
        protCounter += 1

    diamondFile.close()


def writeOutput(args, outPairs):
    outFile = open(args.output, "w")
    for pair in outPairs:
        outFile.write(pair)
    outFile.close()


def main():
    args = parseCmd()
    outPairs = readSpalnPairs(args)
    addDiamondPairs(args, outPairs)
    writeOutput(args, outPairs)


def parseCmd():

    parser = argparse.ArgumentParser(description='Add top k pairs from diamond \
        seed protein-gene pairs to the Spaln pairs.')

    parser.add_argument('--diamondPairs', type=str, required=True,
                        help='Diamond seed gene-protein pairs.')
    parser.add_argument('--spalnPairs', type=str, required=True,
                        help='Seed gene-protein pairs coming from Spaln.')
    parser.add_argument('--output', type=str, required=True,
                        help='Output file with combined pairs.')
    parser.add_argument('--k', type=int, default=5,
                        help='Output this many best pairs from diamond pairs. Default = 5.')
    return parser.parse_args()


if __name__ == '__main__':
    main()
