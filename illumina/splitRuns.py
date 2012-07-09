#!/usr/bin/env python
# encoding: utf-8
"""
splitRuns.py

Uses an info file to create input lines for Lance's barcode_splitter.py

Created by Joshua Shapiro on 2012-06-05
"""

import sys
import os
import argparse
import tempfile
import subprocess
import multiprocessing
import signal

class FlowLane(object):
  """Representation of a flowcell lane and its contents"""
  def __init__(self, flowcell, lane, 
               read_dir = None, fwd_read = None, rev_read = None, index_read = None):
    super(FlowLane, self).__init__()
    self.flowcell = flowcell
    self.lane = lane
    self.read_dir = read_dir
    self.fwd_read = fwd_read
    self.rev_read = rev_read
    self.index_read = index_read
    self.barcodes = dict()
  
  def addBarcode(self, sample, barcode):
    if barcode in self.barcodes:
      raise ValueError('barcode present twice in flowcell %s lane %s' % (self.flowcell, self.lane))
    self.barcodes[barcode] = sample
  
  def checkFiles(self, read_dir, fwd_read, rev_read, index_read):
    if read_dir:
      assert self.read_dir == read_dir
    if fwd_read:
      assert self.fwd_read == fwd_read
    if rev_read:
      assert self.rev_read == rev_read
    if index_read:
      assert self.index_read == index_read
  
  def writeIndexFile(self, filehandle):
    for barcode, sample in self.barcodes.items():
      filehandle.write('%s\t%s\n' % (sample, barcode))
    filehandle.flush()


def getArgs():
  """Get command line options."""
  parser = argparse.ArgumentParser(description = 'Splits runs by barcode reads')
  required = parser.add_argument_group('required arguments')
  required.add_argument('infofile', help = 'The file listing flowcell & lane info, file locations, and index identities')
  options = parser.add_argument_group('optional arguments')
  options.add_argument('-d', dest = "in_dir", help = 'Input directory', default = '', )
  options.add_argument('-o', dest = "out_dir", help = 'Output directory')
  options.add_argument('-t', dest = "threads", type=int, help = 'Number of worker threads', default = 1)
  args = parser.parse_args()
  return args

def readInfoFile(filename, base_dir):
  '''read index file info into a dictionary indexed by flowcell lanes'''
  lane_dict = dict()
  with open(filename, 'r') as f:
    for line in f:
      line = line.strip()
      if line[0] == "#":
        continue
      fields = line.split('\t')
      flow_lane = "_".join(fields[0:2])
      (strain, barcode, read_dir, fwd_read, rev_read, index_read) = fields[3:9]
      strain_flowlane = "_".join([strain, flow_lane])
      read_dir = os.path.join(base_dir, read_dir)
      if flow_lane in lane_dict:
        lane_dict[flow_lane].checkFiles(read_dir, fwd_read, rev_read, index_read)
        lane_dict[flow_lane].addBarcode(strain_flowlane, barcode)
      else:
        lane_dict[flow_lane] = FlowLane(fields[0], fields[1],
                                        read_dir, fwd_read, rev_read, index_read)
        lane_dict[flow_lane].addBarcode(strain_flowlane, barcode)
  return(lane_dict)

def runSplitter(lane, mismatches = 1):
  try:
    bcfile = tempfile.NamedTemporaryFile()
    lane.writeIndexFile(bcfile)
    output = "Flowcell %s, lane %s\n" % (lane.flowcell, lane.lane)
    output += subprocess.check_output(['barcode_splitter.py', '--bcfile=%s' % bcfile.name,
                                       '--mismatches=%s' % mismatches, '--idxread=2',
                                       os.path.join(lane.read_dir, lane.fwd_read),
                                       os.path.join(lane.read_dir, lane.index_read),
                                       os.path.join(lane.read_dir, lane.rev_read)
                                       ])
  finally:
    bcfile.close()
  return output
  
def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    

def main():
  args = getArgs()
  lanes = readInfoFile(args.infofile, base_dir = args.in_dir)
  if args.threads > 1:
    try:
        pool = multiprocessing.Pool(processes=args.threads, init_worker)
        outputs = pool.imap(runSplitter, lanes.values())
    except KeyboardInterrupt:
        pool.terminate()
        pool.wait()
  else:
    outputs = (runSplitter(l) for l in lanes.values())
  for output in outputs:
      print output



if __name__ == '__main__':
    sys.exit(main())
  
    
        
        
        
        
