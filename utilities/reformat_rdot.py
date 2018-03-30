#!/usr/bin/env python

import argparse
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("rdot",
                        help="rdot file from lastz. "
                             "Use '-' to read from stdin.")
    parser.add_argument("outfile",
                        help="File to write csv file to.")
    args = parser.parse_args()
    if args.rdot == "-":
        in_handle = sys.stdin
    else:
        in_handle = open(args.rdot)

    if args.outfile == "-":
        out_handle = sys.stdout
    else:
        out_handle = open(args.outfile, mode="w")

    print("target", "query",
          "t_start", "t_end",
          "q_start", "q_end",
          file = out_handle,
          sep="," )
    target = ""
    query = ""
    target_start = None
    target_end = None
    query_start = None
    query_end = None
    for line in in_handle:
        fields = line.split()
        if fields[0] == "NA":
            print(target, query,
                  target_start, target_end,
                  query_start, query_end,
                  file = out_handle,
                  sep="," )
            target_start = None # new pair
        elif fields[0].isdigit():
            if target_start:
                target_end, query_end = fields
            else:
                target_start, query_start = fields
        else:
            target, query = fields
    in_handle.close()
    out_handle.close()


if __name__ == '__main__':
    main()
