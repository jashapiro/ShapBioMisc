#!/usr/bin/env python
# encoding: utf-8
"""
gzPaste.py
Acts like the default unix paste command, but for gzipped files.
Assumes that all files are the same length.

Created by Joshua Shapiro on 2009-03-04.


"""

import sys
import os
import glob
import gzip


def main():
  #get file list
  files = sys.argv[1:len(sys.argv)]
  handles = [gzip.open(f, 'rb') for f in files]
  nlines = len(handles[0].readlines())
  handles[0].seek(0) 
  gzstdout = gzip.GzipFile(fileobj = sys.stdout, compresslevel = 1)
  for i in range(nlines):
    line = "\t".join([h.readline().strip() for h in handles]) + "\n"
    gzstdout.write(line)
  gzstdout.close()
  [h.close() for h in handles]

if __name__ == '__main__':
    main()

