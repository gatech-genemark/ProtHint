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


def splitSeeds(diamond):
    seeds = {}
    lastRows = {}
    for row in csv.reader(open(diamond), delimiter='\t'):
        protein = row[1]
        if protein not in seeds:
            seeds[protein] = 1
        else:
            if int(lastRows[protein][8]) >= int(row[8]):
                seeds[protein] += 1

        lastRows[protein] = row
        row[1] = row[1] + "_" + str(seeds[protein])
        print("\t".join(row))


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
            print("\t".join(targets[target]))
            targets[target] = row

    for key in targets:
        print("\t".join(targets[key]))


def main():
    args = parseCmd()
    preprocessedDiamond = preprocessInput(args.diamond)
    mergeOverlappingQueryRegions(preprocessedDiamond)
    os.remove(preprocessedDiamond)


def parseCmd():

    parser = argparse.ArgumentParser(description='Create seed gene regions \
        and associated protein hits from DIAMOND BLASTx alignments. This \
        script processes hits corresponding to the same target protein by (a) \
        merging hits separated by introns into seed genes, and (b) splitting \
        hits to duplicated genes (or conserved regions) into distinct seeds.')

    parser.add_argument('diamond', metavar='diamond.tab', type=str,
                        help='DIAMOND blastx output. The program assumes \
                        that DIAMOND was run with the following output \
                        option: "--outfmt 6 qseqid sseqid pident length \
                        mismatch gapopen qstart qend sstart send evalue \
                        bitscore score".')

    parser.add_argument('--maxIntron', type=int, default=50000,
                        help='Maximum intron length. Hits further apart \
                        are split into separate seed genes.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
