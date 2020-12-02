#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2020, Georgia Institute of Technology, USA
#
# Compare files with old and new gene seeds. Output gene seeds which are unique
# in the new seed file and hints from previous iteration of ProtHint which
# correspond to seeds which are identical in both files (with gene ids updated
# to match the ids in the new seed file.
# ==============================================================


import argparse
import csv
import re
import tempfile
import subprocess
import os


def extractFeatureGtf(text, feature):
    regex = feature + ' "([^"]+)"'
    return re.search(regex, text).groups()[0]


def extractFeatureGff(text, feature):
    regex = feature + '=([^;]+);'
    return re.search(regex, text).groups()[0]


def readGff(gffFile):
    for row in csv.reader(open(gffFile), delimiter='\t'):
        pass


class Gene:

    def __init__(self, firstCds, geneId):
        self.CDSs = [firstCds]
        self.geneId = geneId
        self.signature = firstCds[0] + "_" + firstCds[3] + "_" + \
            firstCds[4] + "_" + firstCds[6] + "_" + firstCds[7]

    def addCDS(self, cds):
        self.CDSs.append(cds)
        self.signature += "-" + cds[0] + "_" + cds[3] + "_" + \
            cds[4] + "_" + cds[6] + "_" + cds[7]

    def print(self, out):
        for cds in self.CDSs:
            out.write("\t".join(cds) + "\n")


def loadGenes(gtf):
    sortedGtf = sortGenes(gtf)
    genes = {}
    for row in csv.reader(open(sortedGtf), delimiter='\t'):
        if row[2] != "CDS":
            continue
        geneId = extractFeatureGtf(row[8], "gene_id")
        if geneId not in genes:
            genes[geneId] = Gene(row, geneId)
        else:
            genes[geneId].addCDS(row)
    os.remove(sortedGtf)

    return genes


def sortGenes(gtf):
    sortedGtf = tempfile.NamedTemporaryFile(delete=False).name
    subprocess.call(["bash", "-c", "sort -k1,1 -k4,4n -k5,5n " + gtf +
                    " > " + sortedGtf])
    return sortedGtf


def findUniqueNewGenes(newGenes, oldGenes, outFile):
    oldGeneSignatures = {}

    # Hash old genes signature to make the search fast
    for gene in oldGenes:
        oldGeneSignatures[oldGenes[gene].signature] = oldGenes[gene].geneId

    # Mapping of old to new gene ids for seed genes which are identical
    old2new = {}

    out = open(outFile, 'w')
    for gene in newGenes:
        newSignature = newGenes[gene].signature
        if newSignature not in oldGeneSignatures:
            newGenes[gene].print(out)
        else:
            old2new[oldGeneSignatures[newSignature]] = newGenes[gene].geneId
    out.close()

    return old2new


def printHintsWithIdenticalSeeds(prevSpalnGff, old2new, outFile):
    out = open(outFile, 'w')
    for row in csv.reader(open(prevSpalnGff), delimiter='\t'):
        seedGeneId = extractFeatureGff(row[8], "seed_gene_id")
        if seedGeneId in old2new:
            oldString = "seed_gene_id=" + seedGeneId + ";"
            newString = "seed_gene_id=" + old2new[seedGeneId] + ";"
            row[8] = row[8].replace(oldString, newString)
            out.write("\t".join(row) + "\n")
    out.close()


def main():
    args = parseCmd()
    oldGenes = loadGenes(args.prevGeneSeeds)
    newGenes = loadGenes(args.geneSeeds)
    old2new = findUniqueNewGenes(newGenes, oldGenes, args.uniqueNewSeedsOut)
    printHintsWithIdenticalSeeds(args.prevSpalnGff, old2new,
                                 args.identicalSpalnGffOut)


def parseCmd():

    parser = argparse.ArgumentParser(description='Compare files with old and \
                                     new gene seeds. Output gene seeds which \
                                     are unique in the new seed file and hints\
                                     from previous iteration of ProtHint\
                                     which correspond to seeds which are \
                                     identical in both files (with gene ids\
                                     updated to match the ids in the new \
                                     seed file.')

    parser.add_argument('--geneSeeds', type=str, required=True,
                        help='New seed genes. One isoform per gene is assumed')

    parser.add_argument('--prevGeneSeeds', type=str, required=True,
                        help='Seed genes from previous iteration. One isoform\
                        per gene is assumed.')

    parser.add_argument('--prevSpalnGff', type=str, required=True,
                        help='Scored hints from previous iteration')

    parser.add_argument('--uniqueNewSeedsOut', type=str, required=True,
                        help='Seed genes which are unique in the new set')

    parser.add_argument('--identicalSpalnGffOut', type=str, required=True,
                        help='Scored hints from previous iteration \
                        corresponding to seed genes which are the same in \
                        --geneSeeds and --prevGeneSeeds files. Gene ids of \
                        the hints are update to match the ids in the new \
                        --geneSeeds file')

    return parser.parse_args()


if __name__ == '__main__':
    main()
