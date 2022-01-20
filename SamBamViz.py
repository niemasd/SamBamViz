#! /usr/bin/env python3
'''
SamBamViz: Tool for generating useful visualizations from a SAM/BAM file
'''

# imports
from datetime import datetime
from matplotlib import rcParams
from matplotlib.lines import Line2D
from os import makedirs
from os.path import isdir, isfile
from seaborn import color_palette, kdeplot, set_context, set_style
from sys import argv, stderr
import argparse
import matplotlib
import matplotlib.pyplot as plt
import pysam

# constants
VERSION = '0.0.1'
global LOGFILE; LOGFILE = None

# prep matplotlib/seaborn
matplotlib.use("Agg")
RC = {"font.size":12,"axes.titlesize":16,"axes.labelsize":14,"legend.fontsize":10,"xtick.labelsize":10,"ytick.labelsize":10}
set_context("paper", rc=RC); set_style("ticks"); rcParams['font.family'] = 'serif'

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
            'mapped': dict(), # data['num_reads']['mapped'][chrom] = number of reads that mapped to `chrom`
        },
        'read_length': {
            'raw': { # raw (unclipped) read length
                'unmapped': list(),
                'mapped': dict(), # data['read_length']['raw']['mapped'][chrom] = `list` of raw (unclipped) lengths of reads mapping to `chrom`
            },
            'clipped': { # clipped read length
                'unmapped': list(),
                'mapped': dict() # data['read_length']['clipped']['mapped'][chrom] = `list` of clipped lengths of reads mapping to `chrom`
            },
        },
        'insert_size': list(),
        'coverage': dict(), # data['coverage'][chrom][pos] = number of reads that covered postion `pos` of `chrom`
    }

    # iterate over SAM/BAM
    for read in aln:
        data['num_reads']['all'] += 1

        # insert size
        if read.template_length > 0: #and read.is_proper_pair and read.is_paired:
            data['insert_size'].append(read.template_length)

        # unmapped reads
        if read.is_unmapped:
            data['num_reads']['unmapped'] += 1
            data['read_length']['raw']['unmapped'].append(read.query_length)
            data['read_length']['clipped']['unmapped'].append(read.query_alignment_length)

        # mapped reads
        else:
            # parse chromosome name and add if new
            chrom = read.reference_name
            if chrom not in data['num_reads']['mapped']:
                data['num_reads']['mapped'][chrom] = 0
                data['read_length']['raw']['mapped'][chrom] = list()
                data['read_length']['clipped']['mapped'][chrom] = list()
                data['coverage'][chrom] = dict()
            data['num_reads']['mapped'][chrom] += 1
            data['read_length']['raw']['mapped'][chrom].append(read.query_length)
            data['read_length']['clipped']['mapped'][chrom].append(read.query_alignment_length)
            for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True, with_seq=False):
                if read_pos not in data['coverage'][chrom]:
                    data['coverage'][chrom][read_pos] = 0
                data['coverage'][chrom][read_pos] += 1
    return data

# plot read length distributions
def plot_read_length(data, outdir):
    for k, s in [('raw', "Raw (Unclipped)"), ('clipped', "Clipped")]:
        fig, ax = plt.subplots(); handles = list()
        if len(data['read_length'][k]['unmapped']) + sum(len(data['read_length'][k]['mapped'][chrom]) for chrom in data['read_length'][k]['mapped']) > max([len(data['read_length'][k]['unmapped'])] + [len(data['read_length'][k]['mapped'][chrom]) for chrom in data['read_length'][k]['mapped']]):
            handles.append(Line2D([0],[0],label="All",color='black',linestyle='-'))
            kdeplot(data['read_length'][k]['unmapped'] + [v for chrom in data['read_length'][k]['mapped'] for v in data['read_length'][k]['mapped'][chrom]], label="All", color='black', linestyle='-')
        if len(data['read_length'][k]['unmapped']) != 0:
            handles.append(Line2D([0],[0],label="Unmapped",color='red',linestyle='-'))
            kdeplot(data['read_length'][k]['unmapped'], label="Unmapped", color='red', linestyle='-')
        mapped_colors = color_palette('Greens', n_colors=len(data['read_length'][k]['mapped']))
        for i, chrom in enumerate(sorted(data['read_length'][k]['mapped'].keys())):
            if len(data['read_length'][k]['mapped'][chrom]) != 0:
                handles.append(Line2D([0],[0],label=chrom,color=mapped_colors[i],linestyle='-'))
                kdeplot(data['read_length'][k]['mapped'][chrom], label=chrom, color=mapped_colors[i], linestyle='-')
        plt.title("%s Read Length" % s)
        plt.xlabel("%s Read Length" % s)
        plt.ylabel("Kernel Density Estimate")
        plt.legend(handles=handles, frameon=True)
        fig.savefig("%s/read_length_%s.pdf" % (outdir,k), format='pdf', bbox_inches='tight'); plt.close(fig)

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

    # print number of reads
    print_log("- Number of reads: %d" % data['num_reads']['all'])
    print_log("  - Unmapped: %d" % data['num_reads']['unmapped'])
    print_log("  - Mapped: %d" % sum(data['num_reads']['mapped'].values()))
    for chrom in sorted(data['num_reads']['mapped'].keys()):
        print_log("    - %s: %d" % (chrom, data['num_reads']['mapped'][chrom]))

    # print average read length
    for k, s in [('raw', "raw (unclipped)"), ('clipped', "clipped")]:
        print_log("- Average %s read length: %s" % (s, ((sum(data['read_length'][k]['unmapped']) + sum(sum(data['read_length'][k]['mapped'][chrom]) for chrom in data['read_length'][k]['mapped']))/(len(data['read_length'][k]['unmapped']) + sum(len(data['read_length'][k]['mapped'][chrom]) for chrom in data['read_length'][k]['mapped'])))))
        if len(data['read_length'][k]['unmapped']) == 0:
            print_log("  - Unmapped: N/A")
        else:
            print_log("  - Unmapped: %s" % (sum(data['read_length'][k]['unmapped'])/len(data['read_length'][k]['unmapped'])))
        print_log("  - Mapped: %s" % (sum(sum(data['read_length'][k]['mapped'][chrom]) for chrom in data['read_length'][k]['mapped'])/sum(len(data['read_length'][k]['mapped'][chrom]) for chrom in data['read_length'][k]['mapped'])))
        for chrom in sorted(data['read_length'][k]['mapped'].keys()):
            if len(data['read_length'][k]['mapped'][chrom]) == 0:
                print_log("    - %s: N/A" % chrom)
            else:
                print_log("    - %s: %s" % (chrom, sum(data['read_length'][k]['mapped'][chrom])/len(data['read_length'][k]['mapped'][chrom])))

    # print average insert size
    if len(data['insert_size']) == 0:
        print_log("- Average insert size: N/A")
    else:
        print_log("- Average insert size: %s" % (sum(data['insert_size'])/len(data['insert_size'])))

    # print average coverage
    print_log("- Average coverage: %s" % (sum(sum(data['coverage'][chrom].values()) for chrom in data['coverage'])/sum(len(data['coverage'][chrom]) for chrom in data['coverage'])))
    for chrom in sorted(data['num_reads']['mapped'].keys()):
        print_log("  - %s: %s" % (chrom, sum(data['coverage'][chrom].values())/len(data['coverage'][chrom])))

    # generate plots
    print_log("Plotting read length distributions...")
    plot_read_length(data, args.output)
    print_log("Finished plotting read length distributions")

    # finish up
    LOGFILE.close()