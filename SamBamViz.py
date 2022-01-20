#! /usr/bin/env python3
'''
SamBamViz: Tool for generating useful visualizations from a SAM/BAM file
'''

# imports
from datetime import datetime
from os import makedirs
from os.path import isdir, isfile
from sys import argv, stderr
import argparse
import pysam

# constants
VERSION = '0.0.1'
global LOGFILE; LOGFILE = None

# print log
def print_log(s='', end='\n'):
    tmp = "[%s] %s" % (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), s)
    if LOGFILE is not None:
        print(tmp, file=LOGFILE, end=end); LOGFILE.flush()
    print(tmp, file=stderr); stderr.flush()

# error message
def error(s=None):
    if s is None:
        print_log("ERROR")
    else:
        print_log("ERROR: %s" % s)
    exit(1)

# parse user args
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help="Input File (SAM/BAM)")
    parser.add_argument('-o', '--output', required=True, type=str, help="Output Directory")
    args = parser.parse_args()
    if isfile(args.output) or isdir(args.output):
        error("Output already exits: %s" % args.output)
    return args

# open input SAM/BAM file
def open_sam(fn):
    tmp = pysam.set_verbosity(0) # disable htslib verbosity to avoid "no index file" warning
    if fn.lower() == 'stdin':
        aln = pysam.AlignmentFile('-', 'r') # standard input --> SAM
    elif not isfile(fn):
        error("File not found: %s" % fn)
    elif fn.lower().endswith('.sam'):
        aln = pysam.AlignmentFile(fn, 'r')
    elif fn.lower().endswith('.bam'):
        aln = pysam.AlignmentFile(fn, 'rb')
    else:
        error("Invalid input alignment file extension: %s" % fn)
    pysam.set_verbosity(tmp) # re-enable htslib verbosity and finish up
    return aln

# compute stats from SAM/BAM file
def compute_stats(aln):
    # initialize counters/vars
    data = {
        'num_reads': {
            'all': 0,
            'unmapped': 0,
            'mapped': dict(), # keys = chromosomes, values = num mapped
        },
    }

    # iterate over SAM/BAM
    for read in aln:
        data['num_reads']['all'] += 1

        # unmapped reads
        if read.is_unmapped:
            data['num_reads']['unmapped'] += 1

        # mapped reads
        else:
            chrom = read.reference_name
            if chrom not in data['num_reads']['mapped']:
                data['num_reads']['mapped'][chrom] = 0
            data['num_reads']['mapped'][chrom] += 1
    return data

# main content
if __name__ == "__main__":
    # prep user input
    if len(argv) == 1:
        pass # TODO: In the future, run GUI here to fill in argv accordingly (so argparse will run fine)
    args = parse_args()
    aln = open_sam(args.input)
    makedirs(args.output, exist_ok=False)
    LOGFILE = open("%s/log.txt" % args.output, 'w')
    print_log("Executing SamBamViz v%s" % VERSION)
    print_log("Input file: %s" % args.input)
    print_log("Output directory: %s" % args.output)

    # compute stats from SAM/BAM
    print_log("Parsing input file...")
    data = compute_stats(aln)
    print_log("Finished parsing input file")
    print_log("- Total number entries: %d" % data['num_reads']['all'])
    print_log("- Number unmapped: %d" % data['num_reads']['unmapped'])
    print_log("- Number mapped: %d" % sum(data['num_reads']['mapped'].values()))
    for chrom in sorted(data['num_reads']['mapped'].keys()):
        print_log("    - %s: %d" % (chrom, data['num_reads']['mapped'][chrom]))

    # finish up
    LOGFILE.close()
