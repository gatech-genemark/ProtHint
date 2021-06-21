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


def systemCall(cmd):
    if subprocess.call(["bash", "-c", cmd]) != 0:
        sys.exit('error: Program exited due to an ' +
                 'error in command: ' + cmd)


def preprocessInput(diamond):
    """ Reverse order of start and end coordinates for hits on the negative
    strand and sort the output by coordinates.

    Args:
        diamond: Raw DIAMOND output
    """
    flippedDiamond = tempfile.NamedTemporaryFile(mode="w", prefix="flipped",
                                                 dir=".", delete=False)

    for row in csv.reader(open(diamond), delimiter='\t'):
        strand = "+"
        if int(row[6]) > int(row[7]):
            row[6], row[7] = row[7], row[6]
            row[8], row[9] = row[9], row[8]
            strand = "-"
        row.append(strand)
        flippedDiamond.write("\t".join(row) + "\n")

    flippedDiamond.close()
    systemCall("sort -k1,1 -k7,7n -k8,8n " + flippedDiamond.name +
               " -o " + flippedDiamond.name)
    return flippedDiamond.name


def mergeOverlappingQueryRegions(preprocessedDiamond):
    """Merge overlapping (in the query sequence) hits from the same target.
    Extra care is taken to correctly process the score and target protein
    alignment positions of the merged hits.

    Args:
        diamond: Pre-processed DIAMOND output
    """
    mergedQueries = tempfile.NamedTemporaryFile(mode="w", prefix="mergedQ",
                                                dir=".", delete=False)
    targets = {}
    for row in csv.reader(open(preprocessedDiamond), delimiter='\t'):

        start = int(row[6])
        end = int(row[7])

        target = row[0] + "_" + row[1] + "_" + row[13]
        if target not in targets:
            targets[target] = row
            continue

        prevStart = int(targets[target][6])
        prevEnd = int(targets[target][7])

        if start < prevEnd:
            if end <= prevEnd:
                # The hit is fully inside the previous one; therefore,
                # completely skip this hit
                pass
            else:
                # Hits are overlapping, merge the hits into one. The score of
                # of the new hit is a sum of the two original scores in which
                # the overlapping portion is counted only once. The overlapping
                # portion is scored as an average of the two original scores.
                currLen = end - start + 1
                prevLen = prevEnd - prevStart + 1
                overlapLen = prevEnd - start + 1

                currUniquePortion = (currLen - overlapLen) / currLen
                prevUniquePortion = (prevLen - overlapLen) / prevLen

                currScore = float(row[12])
                prevScore = float(targets[target][12])

                newScore = currUniquePortion * currScore + \
                    prevUniquePortion * prevScore + \
                    1 / 2 * (1 - currUniquePortion) * currScore + \
                    1 / 2 * (1 - prevUniquePortion) * prevScore

                row[12] = str(round(newScore, 2))
                row[6] = str(prevStart)
                row[8] = targets[target][8]
                targets[target] = row
        else:
            # Print the previous hit, it cannot be overlapped anymore
            mergedQueries.write("\t".join(targets[target]) + "\n")
            targets[target] = row

    for key in targets:
        mergedQueries.write("\t".join(targets[key]) + "\n")

    mergedQueries.close()
    return mergedQueries.name


def splitTargets(mergedQueries, args):
    """Split hits to duplicated genes (or conserved regions) into
    distinct seeds. The split criteria are explained in the help of
    --maxTargetHitOverlap and --maxIntron cmd arguments.

    Args:
        mergedQueries: Processed DIAMOND output
        args: Cmd arguments
    """
    output = tempfile.NamedTemporaryFile(mode="w", prefix="splitT",
                                         dir=".", delete=False)
    seeds = {}
    prevRows = {}
    for row in csv.reader(open(mergedQueries), delimiter='\t'):
        target = row[0] + "_" + row[1] + "_" + row[13]
        if target not in seeds:
            seeds[target] = 1
        else:
            if int(row[6]) - int(prevRows[target][7]) - 1 > args.maxIntron:
                seeds[target] += 1
            else:
                if row[13] == "+":
                    overlap = int(prevRows[target][9]) - int(row[8]) + 1
                else:
                    overlap = int(row[8]) - int(prevRows[target][9]) + 1
                if overlap > 0:
                    if row[13] == "+":
                        prevLen = int(prevRows[target][9]) - \
                                  int(prevRows[target][8]) + 1
                        currLen = int(row[9]) - int(row[8]) + 1
                    else:
                        prevLen = int(prevRows[target][8]) - \
                                  int(prevRows[target][9]) + 1
                        currLen = int(row[8]) - int(row[9]) + 1

                    if overlap / prevLen > args.maxTargetHitOverlap or \
                       overlap / currLen > args.maxTargetHitOverlap:
                        seeds[target] += 1

        prevRows[target] = row
        row[1] = target + "_" + str(seeds[target])
        output.write("\t".join(row) + "\n")

    output.close()
    return output.name


def getRootCluster(clusters, i):
    """Return the root cluster of a given cluster.

    Args:
        clusters: Cluster pointers
        i: Index of the cluster of interest
    """
    while clusters[i] != i:
        i = clusters[i]
    return i


def clusterSeeds(processedDiamond):
    """Cluster overlapping seeds. Only CDS-level overlaps are considered.

    Args:
        processedDiamond: Processed DIAMOND output
    """

    clusteredCDS = tempfile.NamedTemporaryFile(mode="w", prefix="clusteredCDS",
                                               dir=".", delete=False)

    clusteredSeeds = tempfile.NamedTemporaryFile(mode="w", prefix="clustered",
                                                 dir=".", delete=False)

    systemCall("sort -k1,1 -k7,7n -k8,8n " + processedDiamond +
               " -o " + processedDiamond)

    clusterId = -1
    prevContig = ""
    currentClusterEnd = 0
    clusters = []
    seed2cluster = {}

    for row in csv.reader(open(processedDiamond), delimiter='\t'):
        contig = row[0]
        start = int(row[6])
        end = int(row[7])
        seed = row[1]

        if prevContig != contig or start > currentClusterEnd:
            clusterId += 1
            clusters.append(clusterId)
            currentClusterEnd = end
        else:
            if end > currentClusterEnd:
                currentClusterEnd = end

        if seed not in seed2cluster:
            seed2cluster[seed] = clusterId
        elif seed2cluster[seed] != clusterId:
            parentCluster = getRootCluster(clusters, clusterId)
            clusters[parentCluster] = getRootCluster(clusters,
                                                     seed2cluster[seed])

        row.append(str(clusterId))
        clusteredCDS.write("\t".join(row) + "\n")
        prevContig = contig

    clusteredCDS.close()

    for row in csv.reader(open(clusteredCDS.name), delimiter='\t'):
        row[14] = str(getRootCluster(clusters, int(row[14])))
        clusteredSeeds.write("\t".join(row) + "\n")

    os.remove(clusteredCDS.name)
    clusteredSeeds.close()
    return clusteredSeeds.name


def diamond2gff(preprocessedDiamond):
    """Convert the DIAMOND output to gff and print the result to stdout

    Args:
        diamond: Pre-processed DIAMOND output
    """
    for row in csv.reader(open(preprocessedDiamond), delimiter='\t'):
        print("\t".join([row[0], "DIAMOND", "CDS", row[6], row[7], "1",
                         row[13], ".", "gene_id=" + row[1] + ";" +
                         " transcript_id=" + row[1] + ";" +
                         " targetFrom=" + row[8] + ";" +
                         " targetTo=" + row[9] + ";" +
                         " score=" + row[12] + ";"]))


def printSeeds(preprocessedDiamond):
    """Print seeds in gff format. Seeds are created by joining hits belonging
    to the same target (after the splitting procedure). Seed score is a sum
    of scores of individual hits.

    Args:
        diamond: Pre-processed DIAMOND output
    """
    seeds = {}

    for row in csv.reader(open(preprocessedDiamond), delimiter='\t'):
        if row[1] not in seeds:
            seeds[row[1]] = [row[0], row[13], row[6], row[7], row[8], row[9],
                             row[12]]
        else:
            seeds[row[1]][3] = row[7]
            seeds[row[1]][5] = row[9]
            seeds[row[1]][6] = str(float(seeds[row[1]][6]) + float(row[12]))

    for key in seeds:
        seed = seeds[key]
        seed[6] = str(round(float(seed[6]), 2))
        print("\t".join([seed[0], "DIAMOND", "seed", seed[2], seed[3], "1",
                         seed[1], ".", "gene_id=" + key + ";" +
                         " transcript_id=" + key + ";" +
                         " targetFrom=" + seed[4] + ";" +
                         " targetTo=" + seed[5] + ";" +
                         " score=" + seed[6] + ";"]))


def main():
    args = parseCmd()
    preprocessedDiamond = preprocessInput(args.diamond)
    mergedQueries = mergeOverlappingQueryRegions(preprocessedDiamond)
    os.remove(preprocessedDiamond)
    final = splitTargets(mergedQueries, args)
    os.remove(mergedQueries)
    printSeeds(final)
    os.remove(final)


def parseCmd():

    parser = argparse.ArgumentParser(description='Create seed gene regions \
        and associated protein hits from DIAMOND BLASTx alignments. This \
        script processes hits corresponding to the same target protein by (a) \
        merging hits with overlapping query alignment (b) assigning hits \
        separated by introns into seed genes, (c) splitting hits to \
        duplicated genes (or conserved regions) into distinct seeds.')

    parser.add_argument('diamond', metavar='diamond.tab', type=str,
                        help='DIAMOND blastx output. The program assumes \
                        that DIAMOND was run with the following output \
                        option: "--outfmt 6 qseqid sseqid pident length \
                        mismatch gapopen qstart qend sstart send evalue \
                        bitscore score".')

    parser.add_argument('--maxTargetHitOverlap', type=float, default=0.75,
                        help='Maximum allowed fraction of a hit\'s target \
        alignment to be overlapped by a target alignment from a neighboring \
        hit from the same target (neighboring on the query). Hits exceeding \
        this fraction are split into distinct seeds. If a target alignment of \
        a hit starts or ends before the target alignment of a preceding hit \
        (preceding on the query), the overlap is considered to be 1.')

    parser.add_argument('--maxIntron', type=int, default=50000,
                        help='Maximum intron length. Hits further apart \
                        are split into separate seed genes.')

    args = parser.parse_args()

    if args.maxTargetHitOverlap > 1 or args.maxTargetHitOverlap < 0:
        sys.exit("Error: argument --maxTargetHitOverlap must be specified in" +
                 " range from 0 to 1.")

    return args


if __name__ == '__main__':
    main()
