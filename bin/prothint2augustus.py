#!/usr/bin/env python
# Author: Tomas Bruna

# Convert ProtHint output to Augustus compatible format
#
# positional arguments:
#   prothint.gff  ProtHint output file.


import argparse
import csv


def main():
    args = parseCmd()
    for row in csv.reader(open(args.input), delimiter='\t'):
        if (row[2].lower() == "intron"):
            row[2] = "intron"
        elif (row[2].lower() == "start_codon"):
            row[2] = "start"
        elif (row[2].lower() == "stop_codon"):
            row[2] = "stop"
        row[8] = "src=P;mult=" + row[5] + ";pri=4"
        print("\t".join(row))


def parseCmd():
    parser = argparse.ArgumentParser(description='Convert ProtHint output to \
                                     Augustus compatible format')

    parser.add_argument('input', metavar='prothint.gff', type=str,
                        help='ProtHint output file.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
