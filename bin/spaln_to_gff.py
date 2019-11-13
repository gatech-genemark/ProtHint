#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2019, Georgia Institute of Technology, USA
#
# Convert spaln output to gff with introns, start and stop codons. Starts are
# reported only if the alignment starts with ATG and the protein start matches
# the alignment start. Stops are reported only if the codon right after the last
# aligned exon is a stop codon and protein end matches the alignment end. Spaln
# output file is read from standard input.
# ==============================================================


import sys
import argparse


def getPair(row):
    return row[0] + "_" + row[1]


def printIntron(row, prevRow, eScore):
    if float(prevRow[10].strip()) < eScore or float(row[10].strip()) < eScore:
        return

    spliceSites = row[14].strip().split()[5]

    print(row[1] + '\tSpaln\tIntron\t' +
          str(int(prevRow[9]) + 1) + '\t' + str(int(row[8]) - 1) +
          '\t.\t+\t.\tprot=' + row[0] + '; iScore=' + row[12].strip() +
          '; Llen=' + prevRow[3].strip() + '; LScore=' + prevRow[10].strip() +
          '; Lid=' + prevRow[2].strip() + '; Rlen=' + row[3].strip() +
          '; RScore=' + row[10].strip() + '; Rid=' + row[2].strip() +
          '; Ssites=' + spliceSites + ';')


def printStart(row, eScore, geneFasta):
    protStart = int(row[6])
    exonScore = float(row[10].strip())
    genomeStart = int(row[8])

    if exonScore > eScore and protStart == 1:
        with open(geneFasta) as gene:
            gene.readline()
            seq = gene.read().replace("\n", "")
            firstCodon = seq[genomeStart - 1:genomeStart + 2]

        if firstCodon == "ATG":
            print(row[1] + '\tSpaln\tStart_codon\t' +
                  str(genomeStart) + '\t' + str(genomeStart + 2) +
                  '\t.\t+\t0\tprot=' + row[0] + '; eScore=' + row[10].strip() +
                  '; elen=' + row[3].strip() + '; eid=' + row[2].strip() + ';')


def printStop(row, eScore, geneFasta, protFasta):
    protEnd = int(row[7])
    exonScore = float(row[10].strip())
    genomeEnd = int(row[9])
    protLength = 0

    with open(protFasta) as prot:
        prot.readline()
        prot = prot.read().replace("\n", "")
        protLength = len(prot)

    if exonScore > eScore and protLength - protEnd < 5:
        with open(geneFasta) as gene:
            gene.readline()
            seq = gene.read().replace("\n", "")
            if genomeEnd + 3 > len(seq):
                return
            stopCodon = seq[genomeEnd:genomeEnd + 3]

            if stopCodon == "TAA" or stopCodon == "TAG" or stopCodon == "TGA":
                print(row[1] + '\tSpaln\tStop_codon\t' +
                      str(genomeEnd + 1) + '\t' + str(genomeEnd + 3) +
                      '\t.\t+\t0\tprot=' + row[0] + '; eScore=' + row[10].strip() +
                      '; elen=' + row[3].strip() + '; eid=' + row[2].strip() + ';')


def main():

    args = parseCmd()

    prevPair = ""
    prevRow = ""
    for line in sys.stdin:
        row = line.split("\t")
        if row[0][0] == '#' or row[0][0] == '@':
            continue

        if int(row[9]) - int(row[8]) < 0:
            continue

        pair = getPair(row)
        if (pair == prevPair):
            printIntron(row, prevRow, args.intronScore)
        else:
            printStart(row, args.startScore, args.gene)

        prevPair = pair
        prevRow = row

    if (prevRow != ""):
        printStop(prevRow, args.stopScore, args.gene, args.prot)


def parseCmd():

    parser = argparse.ArgumentParser(description='Convert spaln output to gff with introns, \
                                     start and stop codons. Starts are reported only if the \
                                     alignment starts with ATG and the protein start matches \
                                     the alignment start. Stops are reported only if the codon \
                                     right after the last aligned exon is a stop codon and protein \
                                     end matches the alignment end. Spaln output file is read \
                                     from standard input.')

    parser.add_argument('--intronScore', type=float, required=True,
                        help='Minimum score of exons downstream and upstream of \
                        intron fot intron to be reported')
    parser.add_argument('--startScore', type=float, required=True,
                        help='Minimum score of a first aligned exon for a start \
                        to be reported')
    parser.add_argument('--stopScore', type=float, required=True,
                        help='Minimum score of last aligned exon for a stop \
                        to be reported')
    parser.add_argument('--gene', type=str, required=True,
                        help='Fasta file with genomic sequence used for alignment.')
    parser.add_argument('--prot', type=str, required=True,
                        help='Fasta file with protein sequence used for alignment.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()
