#!/usr/bin/env python
# encoding: utf-8
"""
splitMap.py

Splits Fastq files by a leading barcode or a separate barcode, strips the barcode and an optional
trim length, then maps with bwa, converts to indexed, sorted bam, and recombines. 

Portions copied from bwa_wrapper.py

Created by Joshua Shapiro on 2011-03-31.
"""

import sys
import os
import subprocess
import shutil
import argparse
import tempfile
import itertools
import gzip


#make sure we can open as many files as we need. Fails on some systems.
#import resource
#resource.setrlimit(resource.RLIMIT_NOFILE, (1000, -1))

from Bio.SeqIO.QualityIO import FastqGeneralIterator

BASES=['A', 'C', 'G','T', 'N']

def getOptions():
  """Get command line options"""
  parser = argparse.ArgumentParser()
  parser.add_argument("-b", "--barcode", 
                    action = "store", dest = "barcode_file", default = None,
                    help = "Tab delimited barcode definition file.")
  parser.add_argument("-f","--read_fwd",
                    action = "store", dest = "read_fwd", default = None,
                    help = "Reads (forward if paired end) in fastq format")
  parser.add_argument("-g","--read_rev",
                    action = "store", dest = "read_rev", default = None,
                    help = "Reverse reads in fastq format")
  parser.add_argument("-c", "--read_bc",
                     action = "store", dest = "read_bc", default = None,
                     help = "Barcode reads in fastq format")
  parser.add_argument("-o","--outfile",
                    action = "store", dest = "out_file", default = None,
                    help = "Output file for merged bam")
  parser.add_argument("-t","--trim",
                    action = "store", dest = "trim", 
                    type = int, default = 0,
                    help = "How many bases beyond the barcode should be trimmed")
  parser.add_argument("-m","--mismatch",
                    action = "store", dest = "mismatch", 
                    type = int, default = 0,
                    help = "How many mismatches are allowed in a barcode")
  
  parser.add_argument("-u","--unmatched",
                    action = "store", dest = "unmatched", default = None,
                    help = ("File to save reads without any of the expected barcodes."
                            "Default is to discard unmmatched reads"))
  parser.add_argument("-r","--reference",
                    action = "store", dest = "ref_file", default = "",
                    help = "The reference genome to use. Must be indexed for bwa and samtools, with index files named as per standard")
  parser.add_argument("--bwa-aln-opts",
                    action = "store", dest = "bwa_aln_options", default = "",
                    help = ("A quoted string containing mapping options for bwa aln,"
                            "except for those related to files"))
  parser.add_argument("--bwa-sam-opts",
                    action = "store", dest = "bwa_sam_options", default = "",
                    help = ("A quoted string containing mapping options for bwa samse or bwa sampe, as appropriate"
                            "except for those related to files"))
  parser.add_argument("-l","--library",
                    action = "store", dest = "library", default = None,
                    help = "Sets a value for library. Default value is the basename of forward read file.")
  parser.add_argument("-z","--zipped",
                    action = "store_true", dest = "zipped", default = False,
                    help = "Input fastq files are gzipped")
  parser.add_argument("-D","--debug",
                    action = "store_true", dest = "debug", default = False,
                    help = "Debug mode... temp files are put in './temp' and not deleted")

  
  
  options = parser.parse_args()
  if not (options.barcode_file and options.read_fwd and options.out_file and options.ref_file):
    parser.print_help()
    sys.exit(1)
  #absolutize paths
  options.barcode_file = os.path.abspath(options.barcode_file)
  options.read_fwd = os.path.abspath(options.read_fwd)
  if options.read_rev:
    options.read_rev = os.path.abspath(options.read_rev)
  if options.read_bc:
    options.read_bc = os.path.abspath(options.read_bc)
  if not options.library:
    options.library = os.path.splitext(os.path.basename(options.read_fwd))[0]
  options.out_file = os.path.abspath(options.out_file)
  options.ref_file = os.path.abspath(options.ref_file)
  
  return options

class Barcode(object):
  """Barcode object that stores info about barcodes and the associated strain"""
  def __init__(self,sample_id, barcode, length = None):
    super(Barcode, self).__init__()
    self.barcode = barcode.upper()
    self.sample_id = sample_id
    if length:
      self.length = length
    else:
      self.length = len(barcode)
  
  def expand(self, length):
    """Create all possible barcodes of the specified length from a shorter barcode
    Retains info about the original barcode length for trimming, etc."""
    if length < self.length:
      raise ValueError("Barcodes can not be extended to be shorter than they already are.")
    if length == self.length:
      return([self.barcode])
    #now the expansion code
    bc_list = [self.barcode]
    for i in range(length - self.length):
      bc_list = [b+c for b in bc_list for c in BASES]
    return(bc_list)

def mismatchBarcodes(barcode_dict, mismatches):
    """Generate all possible one base mismatches"""
    expanded_dict = dict()
    for n in range(mismatches):
      for barcode, sample in barcode_dict.iteritems():
        for i in range(len(barcode)):
          bc_list = list(barcode)
          for b in BASES:
            bc_list[i] = b
            new_bc = ''.join(bc_list)
            # check if new_bc has already been generated and points to a different original
            if (new_bc in expanded_dict) and (expanded_dict[new_bc] != sample):
              if expanded_dict[new_bc] == AMBIGUOUS: 
                pass
              else:
                # check if the bc is from this expansion round (Hamming distance > n). If so, set to unmatched
                if hamming_distance(new_bc, expanded_dict[new_bcs].barcode) > n: 
                  expanded_dict[new_bc] = AMBIGUOUS
            else:
              expanded_dict[new_bc] = sample
      barcode_dict = expanded_dict.copy()
    return barcode_dict        

def hamming_distance(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))  

def parseBarcodeFile(file_path):
  """Parse the barcode file, which should just be a tabbed list.
  Returns a dict keyed by barcode sequence with barcode objects within
  and the length of the longest barcode."""
  infile = open(file_path, 'r')
  try:
    lines = infile.readlines()
  finally:
    infile.close()
  
  barcodes = dict()
  bc_info = [line.split() for line in lines if line[0] != '#']
  bc_len = max(len(b[1]) for b in bc_info)
  for sample_id, bc_seq in bc_info:
    base_bc = Barcode(sample_id, bc_seq)
    for bc in base_bc.expand(bc_len):
      #check that the barcode has not been seen before.
      if barcodes.setdefault(bc, base_bc) != base_bc:
        raise ValueError("Barcodes are not unique")
  return (barcodes, bc_len)
  

  
def runWrapper(cmd, tmp_dir, outfile):
  """wrapper to eliminate excception-handling boilerplate"""
  buffsize = 1048576
  try:
    tmp = tempfile.NamedTemporaryFile( dir=tmp_dir ).name
    tmp_stderr = open( tmp, 'wb' )
    proc = subprocess.Popen( args=cmd, shell=True, cwd=tmp_dir, stderr=tmp_stderr.fileno() )
    returncode = proc.wait()
    tmp_stderr.close()
    # get stderr, allowing for case where it's very large
    tmp_stderr = open( tmp, 'rb' )
    stderr = ''
    try:
      while True:
        stderr += tmp_stderr.read( buffsize )
        if not stderr or len( stderr ) % buffsize != 0:
          break
    except OverflowError:
      pass
    tmp_stderr.close()
    if returncode != 0:
      raise Exception, stderr
  except Exception, e:
    raise Exception, 'Error running command: ' + cmd + '\n'+ str( e ) 
  # check that there are results in the output file
  if outfile:
    output_size = os.path.getsize( outfile )
    return output_size

      

def runBWA(cmd, outfile, tmp_dir):
  """code mostly copied from galaxy: bwa_wrapper.py"""
  buffsize = 1048576
  try:
      # align
      try:
          tmp = tempfile.NamedTemporaryFile( dir=tmp_dir ).name
          tmp_stderr = open( tmp, 'wb' )
          proc = subprocess.Popen( args=cmd + " > " + outfile, shell=True, cwd=tmp_dir, stderr=tmp_stderr.fileno() )
          returncode = proc.wait()
          tmp_stderr.close()
          # get stderr, allowing for case where it's very large
          tmp_stderr = open( tmp, 'rb' )
          stderr = ''
          try:
              while True:
                  stderr += tmp_stderr.read( buffsize )
                  if not stderr or len( stderr ) % buffsize != 0:
                      break
          except OverflowError:
              pass
          tmp_stderr.close()
          if returncode != 0:
              raise Exception, stderr
      except Exception, e:
          raise Exception, 'Error generating alignments. ' + str( e ) 
      # check that there are results in the output file
      if os.path.getsize( outfile ) > 0:
        pass
        # sys.stdout.write( 'BWA run on single-end data.\n'  )
      else:
        sys.stdout.write( "No mapped reads.\n")
  except Exception, e:
      raise Exception, 'The alignment failed.\n' + str( e ) 
      
  
def samToBam(infile, bam_base, ref_file, tmp_dir):
  # Extract all alignments from the input SAM file to BAM format ( since no region is specified, all the alignments will be extracted ).
  if len( open( infile ).read() ) == 0:
      raise Exception, 'Initial SAM file empty'
  command = 'samtools view -uT %s %s | samtools sort - %s' % (ref_file, infile, bam_base)
  tmp = tempfile.NamedTemporaryFile( dir=tmp_dir ).name
  tmp_stderr = open( tmp, 'wb' )
  proc = subprocess.Popen( args=command, shell=True, cwd=tmp_dir, stderr=tmp_stderr.fileno() )
  returncode = proc.wait()
  tmp_stderr.close()
  # get stderr, allowing for case where it's very large
  tmp_stderr = open( tmp, 'rb' )
  stderr = ''
  buffsize = 1048576
  try:
      while True:
          stderr += tmp_stderr.read( buffsize )
          if not stderr or len( stderr ) % buffsize != 0:
              break
  except OverflowError:
      pass
  tmp_stderr.close()
  if returncode != 0:
      raise Exception, stderr
  


def main():
  opts = getOptions()
  if opts.debug:
    tmp_dir = os.path.abspath("temp") 
  else:  
    tmp_dir = tempfile.mkdtemp()
  barcodes, bc_len = parseBarcodeFile(opts.barcode_file)
  if opts.mismatch:
      barcodes = mismatchBarcodes(barcodes, opts.mismatch)
  #create files and filehandles for the fastq files
  fq_handles = dict()
  rev_handles = dict()
  samples = list(set([bc.sample_id for bc in barcodes.values()]))
  samples.sort()
  
  headers = ["@RG\tID:%s\tSM:%s\tLB:%s" % ("_".join([sample,opts.library]), sample, opts.library) for sample in samples]
  header_file = tempfile.NamedTemporaryFile( dir=tmp_dir ).name
  header_handle = open(header_file, 'w')
  for sample, header in itertools.izip (samples, headers):
    file_path = os.path.join(tmp_dir, sample + ".fastq")
    fq_handles[sample] = open(file_path, "w")
    header_handle.write(header + "\n") 
  header_handle.close()
  if opts.read_rev:
    for sample in samples:
      file_path = os.path.join(tmp_dir, sample + "_rev.fastq")
      rev_handles[sample] = open(file_path, "w")
  
  if opts.unmatched:
    unmatched = open(opts.unmatched)
  
  # read the fastQ file, split by barcode and write out
  i = 0
  try:
    if (not opts.read_rev) and (not opts.read_bc):
      # old form: single end, with integrated barcode.
      if opts.zipped:
        fh = gzip.open(opts.read_fwd, 'r')
      else:
        fh = open(opts.read_fwd, 'r')
      for title, seq, qual in FastqGeneralIterator(fh):
        bc_sample = barcodes.get(seq[:bc_len], None)
        if not bc_sample:
          if opts.unmatched:
            unmatched.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
          continue
        handle = fq_handles[bc_sample.sample_id]
        handle.write("@%s\n%s\n+\n%s\n" % 
                      (title, 
                      seq[bc_sample.length + opts.trim:],
                      qual[bc_sample.length + opts.trim:]))
      fh.close()
  
    elif opts.read_bc and not opts.read_rev:
      # new form: single ended, separate barcode read
      if opts.zipped:
        fh = gzip.open(opts.read_fwd, 'r')
        bh = gzip.open(opts.read_bc, 'r')
      else:
        fh = open(opts.read_fwd, 'r')
        bh = open(opts.read_bc, 'r')
      for fwd, bc in itertools.izip(FastqGeneralIterator(fh), 
                                    FastqGeneralIterator(bh)):
        ftitle, fseq, fqual = fwd
        btitle, bseq, bqual = bc
        #TODO: add check that read names line up
        bc_sample = barcodes.get(bseq[:bc_len], None)
        if not bc_sample:
          if opts.unmatched:
            unmatched.write("@%s\n%s\n+\n%s\n" % (ftitle, fseq, fqual))
          continue
        handle = fq_handles[bc_sample.sample_id]
        handle.write("@%s\n%s\n+\n%s\n" % (ftitle, fseq, fqual))
      fh.close()
      bh.close()
  
    elif opts.read_bc and opts.read_rev:
      # new form: double ended, separate barcode read
      if opts.zipped:
        fh = gzip.open(opts.read_fwd, 'r')
        rh = gzip.open(opts.read_rev, 'r')
        bh = gzip.open(opts.read_bc, 'r')
      else:
        fh = open(opts.read_fwd, 'r')
        rh = open(opts.read_rev, 'r')
        bh = open(opts.read_bc, 'r')
      for fwd, rev, bc in itertools.izip(FastqGeneralIterator(fh), 
                                         FastqGeneralIterator(rh), 
                                         FastqGeneralIterator(bh)):
        ftitle, fseq, fqual = fwd
        rtitle, rseq, rqual = rev
        btitle, bseq, bqual = bc
        #TODO: add check that read names line up
        bc_sample = barcodes.get(bseq[:bc_len], None)
        if not bc_sample:
          if opts.unmatched:
            unmatched.write("@%s\n%s\n+\n%s\n" % (ftitle, fseq, fqual))
            unmatched.write("@%s\n%s\n+\n%s\n" % (rtitle, rseq, rqual))
          continue
        handle = fq_handles[bc_sample.sample_id]
        handle.write("@%s\n%s\n+\n%s\n" % (ftitle, fseq, fqual))
        rhandle = rev_handles[bc_sample.sample_id]
        rhandle.write("@%s\n%s\n+\n%s\n" % (rtitle, rseq, rqual))
      fh.close()
      rh.close()
      bh.close()
    else:
      #must be paired end with barcode included.. not supported
      print >> stderr, "Included barcodes are not supported for paired end runs."
  finally:
    #close filehandles
    for fh in fq_handles.values():
      fh.close()
    for fh in rev_handles.values():
      fh.close()
    if opts.unmatched:
      unmatched.close()
  
  
  
  #skipping checks that the ref file is indexed for now
  ref_file_name = opts.ref_file
  #run bwa for each file:
  try:
    if opts.read_rev:
      #forward and reverse reads, need to map twice and go
      for sample, header in itertools.izip(samples, headers):
        fastq_fwd = os.path.join(tmp_dir, sample + ".fastq")
        fastq_rev = os.path.join(tmp_dir, sample + "_rev.fastq")
        temp_fwd = os.path.join(tmp_dir, sample + ".sai")
        temp_rev = os.path.join(tmp_dir, sample + "_rev.sai")
        sam = os.path.join(tmp_dir, sample + ".sam")
        sys.stdout.write( "bwa on %s\n" % sample)
        cmd = 'bwa aln %s %s %s > %s;' % (opts.bwa_aln_options, ref_file_name, fastq_fwd, temp_fwd)
        cmd += 'bwa aln %s %s %s > %s;' % (opts.bwa_aln_options, ref_file_name, fastq_rev, temp_rev)
        cmd += 'bwa sampe %s -r %s %s %s %s %s %s' % (opts.bwa_sam_options, repr(header), ref_file_name, temp_fwd, temp_rev, fastq_fwd, fastq_rev )
        #do alignments
        sys.stdout.write(cmd)
        runBWA(cmd, sam, tmp_dir)
        os.remove(fastq_fwd)
        os.remove(fastq_rev)
        os.remove(temp_fwd)
        os.remove(temp_rev)
    else:
      for sample, header in itertools.izip(samples, headers):
        fastq = os.path.join(tmp_dir, sample + ".fastq")
        sam = os.path.join(tmp_dir, sample + ".sam")
        sys.stdout.write( "bwa on %s\n" % sample)
        cmd = 'bwa aln %s %s %s | bwa samse %s -r %s %s - %s' % (
                opts.bwa_aln_options, ref_file_name, fastq, 
                opts.bwa_sam_options, repr(header), ref_file_name, fastq )
        #do alignments
        runBWA(cmd, sam, tmp_dir)
        os.remove(fastq)
      
  except:
    # clean up temp dir
    if os.path.exists( tmp_dir ) and not opts.debug:
          shutil.rmtree( tmp_dir )
    raise
  #convert sam files to bam files
  bam_files = list()
  for sample in samples:
    sam = os.path.join(tmp_dir, sample + ".sam")
    bam_base = os.path.join(tmp_dir, sample)
    try:
      samToBam(sam, bam_base, ref_file_name, tmp_dir)
      bam_files.append(bam_base + ".bam")
      os.remove(sam)
    except:
      #one failure should not kill us
      continue
    
  if len(bam_files) == 0:
    raise Exception, 'No bam files created.'
  #merge bam_files
  cmd = 'samtools merge -h %s %s %s' % ( header_file, opts.out_file,' '.join( bam_files ) )
  tmp = tempfile.NamedTemporaryFile().name
  try:
      tmp_stderr = open( tmp, 'wb' )
      proc = subprocess.Popen( args=cmd, shell=True, stderr=tmp_stderr.fileno() )
      returncode = proc.wait()
      tmp_stderr.close()
      # get stderr, allowing for case where it's very large
      tmp_stderr = open( tmp, 'rb' )
      stderr = ''
      buffsize = 1048576
      try:
          while True:
              stderr += tmp_stderr.read( buffsize )
              if not stderr or len( stderr ) % buffsize != 0:
                  break
      except OverflowError:
          pass
      tmp_stderr.close()
      if returncode != 0:
          raise Exception, stderr
      if os.path.exists( tmp ) and not opts.debug:
          os.unlink( tmp )
  except Exception, e:
      if os.path.exists( tmp ) and not opts.debug:
          os.unlink( tmp )
      raise Exception( 'Error running SAMtools merge tool\n' + str( e ) )
  if os.path.getsize( opts.out_file ) > 0:
      sys.stdout.write( '%s files merged.\n' %  len(bam_files) )
  else:
      raise Exception( 'The output file is empty, there may be an error with one of your input files.' )
  
  #final cleanup
  if os.path.exists( tmp_dir ) and not opts.debug:
    shutil.rmtree( tmp_dir )
  
if __name__ == '__main__':
  main()

