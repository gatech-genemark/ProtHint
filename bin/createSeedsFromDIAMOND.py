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
from multiprocessing import Pool, Lock
import splitSeedCluster

seedLock = Lock()
pairsLock = Lock()


def systemCall(cmd):
    if subprocess.call(["bash", "-c", cmd]) != 0:
        sys.exit('error: Program exited due to an error in command: ' + cmd)


def getSeedString(row):
    return row[0] + "-" + row[1] + "-" + row[7]


def sortSingle(args):
    systemCall("sort " + args[1] + " " + args[0] + " > " +
               args[0] + ".sorted")


def append(inputName, outputName):
    with open(outputName, "a") as outfile:
        with open(inputName) as infile:
            outfile.write(infile.read())


def fastSort(inputFile, sortString, threads):
    """ Fast, parallel sort. Other parts of this script are not parallelized
    since they only involve a linear traversal of the input file.

    Args:
        inputFile: The file to sort
        sortString: Sorting keys
        threads: How many threads to use
    """
    sortFolder = tempfile.TemporaryDirectory(prefix="sort", dir=".")

    # Using wc -l is the fastest way
    cmd = "wc -l " + inputFile
    try:
        wcOut = subprocess.check_output(["bash", "-c", cmd]).decode()
    except subprocess.CalledProcessError:
        sys.exit('error: Program exited due to an error in command: ' + cmd)
    nLines = int(wcOut.split(" ")[0])

    linesPerFile = math.ceil(nLines / threads)
    # Again, system call is the fastest way to do this
    systemCall("split " + inputFile + " " + sortFolder.name +
               "/raw -a 5 -l " + str(linesPerFile))

    sortFiles = os.listdir(sortFolder.name)
    pool = Pool(processes=threads)
    pool.map(sortSingle, [[sortFolder.name + "/" + x,
                           sortString] for x in sortFiles])

    systemCall("sort " + sortString + " -m " + sortFolder.name + "/*.sorted" +
               " > " + inputFile)


def preprocessInput(diamond, threads):
    """ Reverse order of start and end coordinates for hits on the negative
    strand and sort the output by coordinates.

    Args:
        diamond: Raw DIAMOND output
        threads: How many threads to use in sorting
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
    fastSort(flippedDiamond.name, "-k1,1 -k3,3n -k4,4n", threads)
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

                # Make sure the new protein coordinates make sense
                if row[7] == "+":
                    row[4] = str(min(int(targets[target][4]), int(row[4])))
                    row[5] = str(max(int(targets[target][5]), int(row[5])))
                else:
                    row[4] = str(max(int(targets[target][4]), int(row[4])))
                    row[5] = str(min(int(targets[target][5]), int(row[5])))

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


def clusterSeeds(processedDiamond, threads):
    """Cluster overlapping seeds. Only CDS-level overlaps are considered.

    Args:
        processedDiamond: Processed DIAMOND output
        threads: How many threads to use in sorting
    """

    clusteredCDS = tempfile.NamedTemporaryFile(mode="w", prefix="clusteredCDS",
                                               dir=".", delete=False)

    clusteredSeeds = tempfile.NamedTemporaryFile(mode="w", prefix="clustered",
                                                 dir=".", delete=False)

    fastSort(processedDiamond, "-k8,8r -k1,1 -k3,3n -k4,4n", threads)

    clusterId = -1
    prevContig = ""
    prevStrand = ""
    currentClusterEnd = 0
    clusters = []
    seed2cluster = {}

    for row in csv.reader(open(processedDiamond), delimiter='\t'):
        contig = row[0]
        start = int(row[2])
        end = int(row[3])
        seed = row[8]
        strand = row[7]

        if prevContig != contig or start > currentClusterEnd or \
           prevStrand != strand:
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
        prevStrand = strand

    clusteredCDS.close()

    for row in csv.reader(open(clusteredCDS.name), delimiter='\t'):
        row[9] = str(getRootCluster(clusters, int(row[9])))
        clusteredSeeds.write("\t".join(row) + "\n")

    os.remove(clusteredCDS.name)
    clusteredSeeds.close()
    return clusteredSeeds.name


def diamond2gff(preprocessedDiamond, outputFile):
    """Convert the DIAMOND output to gff and print the result to a file

    Args:
        diamond: Pre-processed DIAMOND output
        outputFile: Output file
    """
    output = open(outputFile, "w")

    for row in csv.reader(open(preprocessedDiamond), delimiter='\t'):
        if len(row) == 8:
            row.append(row[1])

        output.write("\t".join([row[0], "DIAMOND", "CDS", row[2], row[3], "1",
                                row[7], ".", "gene_id=" + row[8] + ";" +
                                " transcript_id=" + row[8] + ";" +
                                " targetFrom=" + row[4] + ";" +
                                " targetTo=" + row[5] + ";" +
                                " score=" + row[6]]) + "\n")
    output.close()


def flushClusterBuffer(buffer, folder):
    for cluster in list(buffer.keys()):
        with open(folder + "/" + cluster, "a") as outfile:
            for row in buffer[cluster]:
                outfile.write("\t".join(row) + "\n")
        del buffer[cluster]


def saveClustersToFiles(clusteredDiamond, clusterFolder):
    buffer = {}
    FLUSH = 1000000
    counter = 0
    for row in csv.reader(open(clusteredDiamond), delimiter='\t'):
        cluster = row[9]
        if cluster not in buffer:
            buffer[cluster] = [row]
        else:
            buffer[cluster].append(row)
        counter += 1

        if counter == FLUSH:
            flushClusterBuffer(buffer, clusterFolder)
            counter = 0

    flushClusterBuffer(buffer, clusterFolder)


def processSingle(args):
    cluster = args[0]
    lowThreshold = args[1]
    highThreshold = args[2]
    topN = args[3]
    seedRegions = args[4]
    alignmentPairs = args[5]
    clusterFolder = args[6]

    splitSeedCluster.split(clusterFolder + "/" + cluster,
                           lowThreshold, highThreshold, topN,
                           clusterFolder + "/seeds_" + cluster,
                           clusterFolder + "/pairs_" + cluster)

    with seedLock:
        append(clusterFolder + "/seeds_" + cluster, seedRegions)
    with pairsLock:
        append(clusterFolder + "/pairs_" + cluster, alignmentPairs)

    os.remove(clusterFolder + "/" + cluster)
    os.remove(clusterFolder + "/seeds_" + cluster)
    os.remove(clusterFolder + "/pairs_" + cluster)


def processClusters(clusteredDiamond, topN, lowThreshold, highThreshold,
                    seedRegions, alignmentPairs, threads):

    clusterFolder = tempfile.TemporaryDirectory(prefix="clusters", dir=".")

    saveClustersToFiles(clusteredDiamond, clusterFolder.name)

    clusters = os.listdir(clusterFolder.name)

    if os.path.exists(seedRegions):
        os.remove(seedRegions)
    if os.path.exists(alignmentPairs):
        os.remove(alignmentPairs)

    pool = Pool(processes=threads)
    pool.map(processSingle, [[x,
                              lowThreshold,
                              highThreshold,
                              topN,
                              seedRegions,
                              alignmentPairs,
                              clusterFolder.name
                              ] for x in clusters])


def main():
    args = parseCmd()
    preprocessedDiamond = preprocessInput(args.diamond, args.threads)

    mergedQueries = mergeOverlappingQueryRegions(preprocessedDiamond)
    os.remove(preprocessedDiamond)
    if args.rawSeeds:
        diamond2gff(mergedQueries, args.rawSeeds)

    processedDiamond = splitTargets(mergedQueries, args)
    os.remove(mergedQueries)
    if args.splitSeeds:
        diamond2gff(processedDiamond, args.splitSeeds)

    clusteredDiamond = clusterSeeds(processedDiamond, args.threads)
    os.remove(processedDiamond)

    processClusters(clusteredDiamond, args.maxProteinsPerSeed,
                    args.lowThreshold, args.highThreshold,
                    args.seedRegions, args.alignmentPairs, args.threads)

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

    parser.add_argument('--maxProteinsPerSeed', type=int, default=25,
                        help='Maximum number of protein per seed region. The \
        best scoring proteins are selected. Default = 25.')

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

    parser.add_argument('--threads', type=int, default=1,
                        help='Number of threads to use.')

    parser.add_argument('--rawSeeds', type=str,
                        help='Output raw, unclustered, seeds here. This is \
        an optional output; mainly useful for debugging.')

    parser.add_argument('--splitSeeds', type=str,
                        help='Output split, unclustered, seeds here. This is \
        an optional output; mainly useful for debugging.')

    parser.add_argument('--lowThreshold', type=float, default=0.1)
    parser.add_argument('--highThreshold', type=float, default=0.2)

    args = parser.parse_args()

    if args.maxTargetHitOverlap > 1 or args.maxTargetHitOverlap < 0:
        sys.exit("Error: argument --maxTargetHitOverlap must be specified in" +
                 " range from 0 to 1.")

    return args


if __name__ == '__main__':
    main()
