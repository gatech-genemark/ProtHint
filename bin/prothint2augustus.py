#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2020, Georgia Institute of Technology, USA
#
# Convert ProtHint output to BRAKER/AUGUSTUS compatible format
# ==============================================================


import argparse
import csv
import prothint
from prothint import callDependency, systemCall
import os

def main():
    args = parseCmd()
    prothint.binDir = os.path.abspath(os.path.dirname(__file__))

    callDependency("log_reg_prothints.pl", "--prothint " + args.prothint +
                   " --out " + args.output + " > /dev/null")

    output = open(args.output, "a")

    for row in csv.reader(open(args.evidence), delimiter='\t'):
        if (row[2].lower() == "intron"):
            row[2] = "intron"
        elif (row[2].lower() == "start_codon"):
            row[2] = "start"
        elif (row[2].lower() == "stop_codon"):
            row[2] = "stop"
        row[8] = "src=M;mult=" + row[5] + ";pri=4"
        output.write("\t".join(row) + "\n")

    with open(args.chains, 'r') as f:
        for line in f:
            output.write(line)

    output.close()

def parseCmd():
    parser = argparse.ArgumentParser(description='Convert ProtHint outputs to \
                                     AUGUSTUS/BRAKER compatible format')

    parser.add_argument('prothint', metavar='prothint.gff', type=str,
                        help='ProtHint output file.')

    parser.add_argument('evidence', metavar='evidence.gff', type=str,
                        help='ProtHint evidence file.')

    parser.add_argument('chains', metavar='chains.gff', type=str,
                        help='ProtHint protein chains file.')

    parser.add_argument('output', metavar='output.gff', type=str,
                        help='Output file name.')

    return parser.parse_args()

if __name__ == '__main__':
    main()
