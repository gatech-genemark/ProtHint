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


def getSeedString(row):
    return row[0] + "-" + row[1] + "-" + row[7]


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
        if int(row[2]) > int(row[3]):
            row[2], row[3] = row[3], row[2]
            row[4], row[5] = row[5], row[4]
            strand = "-"
        row.append(strand)
        flippedDiamond.write("\t".join(row) + "\n")

    flippedDiamond.close()
    systemCall("sort -k1,1 -k3,3n -k4,4n " + flippedDiamond.name +
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

        start = int(row[2])
        end = int(row[3])

        target = getSeedString(row)
        if target not in targets:
            targets[target] = row
            continue

        prevStart = int(targets[target][2])
        prevEnd = int(targets[target][3])

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

                currScore = float(row[6])
                prevScore = float(targets[target][6])

                newScore = currUniquePortion * currScore + \
                    prevUniquePortion * prevScore + \
                    1 / 2 * (1 - currUniquePortion) * currScore + \
                    1 / 2 * (1 - prevUniquePortion) * prevScore

                row[6] = str(round(newScore, 2))
                row[2] = str(prevStart)
                row[4] = targets[target][4]
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
        target = getSeedString(row)
        if target not in seeds:
            seeds[target] = 1
        else:
            if int(row[2]) - int(prevRows[target][3]) - 1 > args.maxIntron:
                seeds[target] += 1
            else:
                if row[7] == "+":
                    overlap = int(prevRows[target][5]) - int(row[4]) + 1
                else:
                    overlap = int(row[4]) - int(prevRows[target][5]) + 1
                if overlap > 0:
                    if row[7] == "+":
                        prevLen = int(prevRows[target][5]) - \
                                  int(prevRows[target][4]) + 1
                        currLen = int(row[5]) - int(row[4]) + 1
                    else:
                        prevLen = int(prevRows[target][4]) - \
                                  int(prevRows[target][5]) + 1
                        currLen = int(row[4]) - int(row[5]) + 1

                    if overlap / prevLen > args.maxTargetHitOverlap or \
                       overlap / currLen > args.maxTargetHitOverlap:
                        seeds[target] += 1

        prevRows[target] = row
        row.append(target + "_" + str(seeds[target]))
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

    systemCall("sort -k1,1 -k3,3n -k4,4n " + processedDiamond +
               " -o " + processedDiamond)

    clusterId = -1
    prevContig = ""
    currentClusterEnd = 0
    clusters = []
    seed2cluster = {}

    for row in csv.reader(open(processedDiamond), delimiter='\t'):
        contig = row[0]
        start = int(row[2])
        end = int(row[3])
        seed = row[8]

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
        row[9] = str(getRootCluster(clusters, int(row[9])))
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
        print("\t".join([row[0], "DIAMOND", "CDS", row[2], row[3], "1",
                         row[7], ".", "gene_id=" + row[8] + ";" +
                         " transcript_id=" + row[8] + ";" +
                         " targetFrom=" + row[4] + ";" +
                         " targetTo=" + row[5] + ";" +
                         " score=" + row[6] + ";",
                         row[9]]))


def printClusters(clusteredDiamond, seedRegions):
    """Print seed clusters in gff format. Clusters are defined by the left-
    and right-most coordinates of its members.

    Args:
        clusteredDiamond: Clustered DIAMOND output
        seedRegions: Output file
    """
    output = open(seedRegions, "w")
    seeds = {}

    for row in csv.reader(open(clusteredDiamond), delimiter='\t'):
        clusterID = row[9]
        if clusterID not in seeds:
            seeds[clusterID] = row
        elif int(row[3]) > int(seeds[clusterID][3]):
            seeds[clusterID][3] = row[3]

    for clusterID in seeds:
        seed = seeds[clusterID]
        output.write("\t".join([seed[0], "DIAMOND", "CDS", seed[2],
                     seed[3], "1", seed[7], "0", "gene_id \"" + clusterID +
                     "\";" + " transcript_id \"" + clusterID + "\";"]) + "\n")
    output.close()


def printPairs(clusteredDiamond, topN, alignmentPairs):
    """Print cluster ID and its score for each seed. Select only topN best
    scoring seeds per each cluster.

    Args:
        clusteredDiamond: Clustered DIAMOND output
        topN: How many best scoring seeds to report per cluster.
        alignmentPairs: Output file
    """
    output = open(alignmentPairs, "w")
    seeds = {}

    for row in csv.reader(open(clusteredDiamond), delimiter='\t'):
        seedID = row[8]
        if seedID not in seeds:
            seeds[seedID] = [float(row[6]), row[9], row[1]]
        else:
            seeds[seedID][0] = seeds[seedID][0] + float(row[6])
            if seeds[seedID][1] != row[9]:
                sys.exit("Error: cluster mismatch within the same seed")

    for seedID in seeds:
        seed = seeds[seedID]
        seed[0] = str(round(seed[0], 2))
        output.write("\t".join([seed[1], seed[2], seed[0]]) + "\n")

    output.close()


def main():
    args = parseCmd()
    preprocessedDiamond = preprocessInput(args.diamond)

    mergedQueries = mergeOverlappingQueryRegions(preprocessedDiamond)
    os.remove(preprocessedDiamond)

    processedDiamond = splitTargets(mergedQueries, args)
    os.remove(mergedQueries)

    clusteredDiamond = clusterSeeds(processedDiamond)
    os.remove(processedDiamond)

    printClusters(clusteredDiamond, args.seedRegions)
    printPairs(clusteredDiamond, args.maxProteinsPerSeed, args.alignmentPairs)

    os.remove(clusteredDiamond)


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

    parser.add_argument('--seedRegions', type=str, required=True,
                        help='Output file for seed regions.')

    parser.add_argument('--alignmentPairs', type=str, required=True,
                        help='Output with a list of up to --maxProteinsPerSeed\
        best scoring proteins per each seed region.')

    parser.add_argument('--maxProteinsPerSeed', type=str,
                        help='Maximum number of protein per seed region. The \
        best scoring proteins are selected.')

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
