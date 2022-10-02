#! /usr/bin/env python3
'''
SamBamViz: Tool for generating useful visualizations from a SAM/BAM file
'''

# imports
from datetime import datetime
from os.path import isdir, isfile
from sys import argv, stderr
import argparse
import pysam

# constants
VERSION = '0.0.10'
NUM_READS_STATUS = 50000

# print log
def print_log(s='', end='\n'):
    tmp = "[%s] %s" % (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), s)
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
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output Base Counts (TSV)")
    parser.add_argument('-q', '--base_qual', required=False, type=float, default=0, help="Minimum base quality to include base in counts")
    parser.add_argument('-m', '--map_qual', required=False, type=float, default=0, help="Minimum mapping quality to include read in counts")
    parser.add_argument('--start_at_one', action="store_true", help="Use 1-based indexing (rather than 0-based indexing)")
    parser.add_argument('--force_bam', action="store_true", help="Force BAM Input (otherwise infer from filename)")
    parser.add_argument('--version', action="store_true", help="Show SamBamViz version")
    args = parser.parse_args()
    if args.output.lower() != 'stdout' and (isfile(args.output) or isdir(args.output)):
        error("Output already exits: %s" % args.output)
    return args

# open input SAM/BAM file
def open_sam(fn, force_bam=False):
    tmp = pysam.set_verbosity(0) # disable htslib verbosity to avoid "no index file" warning
    if fn.lower() == 'stdin':
        if force_bam:
            aln = pysam.AlignmentFile('-', 'rb')
        else:
            aln = pysam.AlignmentFile('-', 'r') # standard input default --> SAM
    elif not isfile(fn):
        error("File not found: %s" % fn)
    elif force_bam or fn.lower().endswith('.bam'):
        aln = pysam.AlignmentFile(fn, 'rb')
    elif fn.lower().endswith('.sam'):
        aln = pysam.AlignmentFile(fn, 'r')
    else:
        error("Invalid input alignment file extension: %s" % fn)
    pysam.set_verbosity(tmp) # re-enable htslib verbosity and finish up
    return aln

# compute stats from SAM/BAM file
def compute_stats(aln, verbose=False):
    # initialize counters/vars
    data = {
        'num_reads': {
            'unmapped': 0,
            'mapped': 0,
        },
        'read_length': {
            'raw': { # raw (unclipped) read length
                'unmapped': list(),
                'mapped': list(), # data['read_length']['raw']['mapped'] = `list` of raw (unclipped) lengths of reads mapping to reference genome
            },
            'clipped': { # clipped read length
                'unmapped': list(),
                'mapped': list() # data['read_length']['clipped']['mapped'] = `list` of clipped lengths of reads mapping to reference genome
            },
        },
        'insert_size': list(),
        'coverage': dict(), # data['coverage'][pos] = number of reads that covered postion `pos` of reference genome
        'nuc_count': dict(), # data['nuc_count'][pos][nuc] = number of occurrences of nucleotide `nuc` at position `pos` of reference genome
    }

    # iterate over SAM/BAM
    ref_name = None
    for read_num, read in enumerate(aln.fetch(until_eof=True)):
        # progress update
        if verbose and read_num % NUM_READS_STATUS == 0 and read_num != 0:
            print_log("Parsed %d reads..." % read_num)

        # skip low-map-quality reads
        if read.mapping_quality < args.map_qual:
            continue

        # insert size
        if read.template_length > 0:
            data['insert_size'].append(read.template_length)

        # unmapped reads
        if read.is_unmapped:
            data['num_reads']['unmapped'] += 1
            data['read_length']['raw']['unmapped'].append(read.query_length)
            data['read_length']['clipped']['unmapped'].append(read.query_alignment_length)

        # mapped reads
        else:
            # check reference name
            if ref_name is None:
                ref_name = read.reference_name
            elif ref_name != read.reference_name:
                error("BAM must only contain a single reference 'chromosome' (the viral genome)")

            # parse mapped read
            quals = read.query_qualities; start = read.query_alignment_start; end = read.query_alignment_end
            data['num_reads']['mapped'] += 1
            data['read_length']['raw']['mapped'].append(read.query_length)
            data['read_length']['clipped']['mapped'].append(read.query_alignment_length)
            for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True, with_seq=False):
                if read_pos < start or read_pos >= end or quals[read_pos] < args.base_qual:
                    continue # skip soft-clipped and low-quality bases
                if ref_pos not in data['coverage']:
                    data['coverage'][ref_pos] = 0
                    data['nuc_count'][ref_pos] = {'A':0, 'C':0, 'G':0, 'T':0, 'X':0}
                data['coverage'][ref_pos] += 1
                nuc = read.query_sequence[read_pos]
                if 'a' <= nuc <= 'z':
                    nuc = {'a':'A', 'c':'C', 'g':'G', 't':'T'}[nuc]
                elif nuc not in {'A','C','G','T'}:
                    nuc = 'X'
                data['nuc_count'][ref_pos][nuc] += 1
    return data

# write nucleotide counts at each position
def write_nuc_counts(data, outdir):
    if outdir.lower() == 'stdout':
        from sys import stdout as f
    else:
        f = open(outdir, 'w')
    f.write("Pos\tA\tC\tG\tT\tOther\tTotal\n")
    max_ref_pos = max(data['nuc_count'].keys())
    for ref_pos in range(max_ref_pos+1):
        out_pos = ref_pos + args.start_at_one
        if ref_pos in data['nuc_count']:
            c = data['nuc_count'][ref_pos]
            f.write("%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (out_pos, c['A'], c['C'], c['G'], c['T'], c['X'], sum(c.values())))
        else:
            f.write("%d\t0\t0\t0\t0\t0\t0\n" % out_pos)
    f.close()

# main content
if __name__ == "__main__":
    # print version (if relevant)
    if '--version' in argv or '-version' in argv:
        print("SamBamViz v%s" % VERSION); exit()

    # prep user input
    if len(argv) == 1:
        pass # TODO: In the future, run GUI here to fill in argv accordingly (so argparse will run fine)
    args = parse_args()
    aln = open_sam(args.input, args.force_bam)
    print_log("Executing SamBamViz v%s" % VERSION)
    if args.start_at_one:
        print_log("Using 1-based indexing")
    else:
        print_log("Using 0-based indexing")
    print_log("Input file: %s" % args.input)
    print_log("Output file: %s" % args.output)

    # compute stats from SAM/BAM
    print_log("Parsing input file...")
    data = compute_stats(aln, verbose=True)
    print_log("Finished parsing input file")

    # print number of reads
    unmapped = data['num_reads']['unmapped']
    mapped = data['num_reads']['mapped']
    print_log("- Number of reads: %d" % (unmapped + mapped))
    print_log("  - Unmapped: %d" % unmapped)
    print_log("  - Mapped: %d" % mapped)

    # print average read length
    for k, s in [('raw', "raw (unclipped)"), ('clipped', "clipped")]:
        print_log("- Average %s read length: %s" % (s, ((sum(data['read_length'][k]['unmapped']) + sum(data['read_length'][k]['mapped']))/(len(data['read_length'][k]['unmapped']) + len(data['read_length'][k]['mapped'])))))
        if len(data['read_length'][k]['unmapped']) == 0:
            print_log("  - Unmapped: N/A")
        else:
            print_log("  - Unmapped: %s" % (sum(data['read_length'][k]['unmapped'])/len(data['read_length'][k]['unmapped'])))
        print_log("  - Mapped: %s" % (sum(data['read_length'][k]['mapped'])/len(data['read_length'][k]['mapped'])))

    # print average insert size
    if len(data['insert_size']) == 0:
        print_log("- Average insert size: N/A")
    else:
        print_log("- Average insert size: %s" % (sum(data['insert_size'])/len(data['insert_size'])))

    # print average coverage
    print_log("- Average coverage: %s" % (sum(data['coverage'].values())/max(data['coverage'].keys())))

    # generate output file
    print_log("Writing nucleotide counts...")
    write_nuc_counts(data, args.output)
    print_log("Finished writing nucleotide counts")
