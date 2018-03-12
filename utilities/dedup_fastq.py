#!/usr/bin/env python
"""
A script to find and remove duplicate fastq entries.

Assumes all reads have unique IDs. 
"""

import argparse
import sys
import gzip
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fq_file",
                        help="FastQ file to remove duplicates from. "
                             "Use '-' to read from stdin.")
    parser.add_argument("outfile",
                        help="FastQ file to write to. "
                             "If filename ends in '.gz', it will be compressed.")
    args = parser.parse_args()

    if args.fq_file == "-":
        in_handle = sys.stdin
    elif args.fq_file[-3:] == ".gz":
        in_handle = gzip.open(args.fq_file, mode="rt")
    else:
        in_handle = open(args.fq_file)

    if args.outfile[-3:] == ".gz":
        out_handle = gzip.open(args.outfile, mode="w")
    else:
        out_handle = open(args.outfile, mode="w")

    reads = set()
    duplicates = list()
    for record in SeqIO.parse(in_handle, "fastq"):
        if record.id in reads:
            duplicates.append(record.id)
        else:
            SeqIO.write(record, out_handle, "fastq")
            reads.add(record.id)

    out_handle.close()
    in_handle.close()

    dupcount = len(duplicates)
    if dupcount == 0:
        print("There were no duplicates in {}".format(args.fq_file))
    elif dupcout == 1:
        print("There was 1 duplicate in {}".format(args.fq_file))
    else:
        print("There were {} duplicates in {}".format(dupcout, args.fq_file))


if __name__ == '__main__':
    main()
