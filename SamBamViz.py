#! /usr/bin/env python3
'''
SamBamViz: Tool for generating useful visualizations from a SAM/BAM file
'''

# imports
from datetime import datetime
from matplotlib import rcParams
from matplotlib.lines import Line2D
from os import cpu_count, makedirs
from os.path import isdir, isfile
from seaborn import color_palette, kdeplot, lineplot, set_context, set_style
from sys import argv, stderr
import argparse
import matplotlib
import matplotlib.pyplot as plt
import pysam

# constants
VERSION = '0.0.6'
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
    parser.add_argument('-q', '--base_qual', required=False, type=float, default=0, help="Minimum base quality to include base in counts")
    parser.add_argument('-m', '--map_qual', required=False, type=float, default=0, help="Minimum mapping quality to include read in counts")
    parser.add_argument('-t', '--threads', required=False, type=int, default=cpu_count(), help="Number of Threads (BAM decompression)")
    parser.add_argument('--ylog', action="store_true", help="Log-scale Y-Axis")
    parser.add_argument('--version', action="store_true", help="Show SamBamViz version")
    args = parser.parse_args()
    if args.threads < 1 or args.threads > cpu_count():
        error("Invalid number of threads: %d (must be between 1 and %d)" % (args.threads, cpu_count()))
    if isfile(args.output) or isdir(args.output):
        error("Output already exits: %s" % args.output)
    return args

# open input SAM/BAM file
def open_sam(fn, threads):
    tmp = pysam.set_verbosity(0) # disable htslib verbosity to avoid "no index file" warning
    if fn.lower() == 'stdin':
        aln = pysam.AlignmentFile('-', 'r', threads=threads) # standard input --> SAM
    elif not isfile(fn):
        error("File not found: %s" % fn)
    elif fn.lower().endswith('.sam'):
        aln = pysam.AlignmentFile(fn, 'r', threads=threads)
    elif fn.lower().endswith('.bam'):
        aln = pysam.AlignmentFile(fn, 'rb', threads=threads)
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
        'nuc_count': dict(), # data['nuc_count'][chrom][pos][nuc] = number of occurrences of nucleotide `nuc` at position `pos` of `chrom`
    }

    # iterate over SAM/BAM
    for read in aln.fetch(until_eof=True):
        if read.mapping_quality < args.map_qual:
            continue # skip low-map-quality reads
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
            quals = read.query_qualities
            if chrom not in data['num_reads']['mapped']:
                data['num_reads']['mapped'][chrom] = 0
                data['read_length']['raw']['mapped'][chrom] = list()
                data['read_length']['clipped']['mapped'][chrom] = list()
                data['coverage'][chrom] = dict()
                data['nuc_count'][chrom] = dict()
            data['num_reads']['mapped'][chrom] += 1
            data['read_length']['raw']['mapped'][chrom].append(read.query_length)
            data['read_length']['clipped']['mapped'][chrom].append(read.query_alignment_length)
            for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True, with_seq=False):
                if quals[read_pos] < args.base_qual:
                    continue # skip low-quality bases
                if ref_pos not in data['coverage'][chrom]:
                    data['coverage'][chrom][ref_pos] = 0
                    data['nuc_count'][chrom][ref_pos] = {'A':0, 'C':0, 'G':0, 'T':0, 'X':0}
                data['coverage'][chrom][ref_pos] += 1
                nuc = read.query_sequence[read_pos].upper()
                if nuc not in {'A','C','G','T'}:
                    nuc = 'X'
                data['nuc_count'][chrom][ref_pos][nuc] += 1
    return data

# plot read length distributions
def plot_read_length(data, outdir, ylog=False):
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
        if ylog:
            ax.set_yscale('log')
        plt.legend(handles=handles, frameon=True)
        fig.savefig("%s/read_length_%s.pdf" % (outdir,k), format='pdf', bbox_inches='tight')
        plt.close(fig)

# plot coverage over genome
def plot_coverage(data, outdir, ylog=False):
    for chrom in sorted(data['coverage'].keys()):
        if len(data['coverage'][chrom]) == 0:
            continue
        xmax = max(data['coverage'][chrom].keys())
        x = list(range(xmax+1))
        y = [0 for _ in range(len(x))]
        for pos in data['coverage'][chrom]:
            y[pos] = data['coverage'][chrom][pos]
        fig, ax = plt.subplots()
        lineplot(x=x, y=y, color='black', linestyle='-')
        plt.title("Coverage: %s" % chrom)
        plt.xlabel("Position (0-based)")
        plt.ylabel("Coverage")
        if ylog:
            ax.set_yscale('log')
        fig.savefig("%s/coverage_%s.pdf" % (outdir, chrom.replace(' ','-')), format='pdf', bbox_inches='tight')
        plt.close(fig)

# write nucleotide counts at each position
def write_nuc_counts(data, outdir):
    for chrom in sorted(data['nuc_count'].keys()):
        f = open("%s/nuc_count_%s.tsv" % (outdir, chrom.replace(' ','-')), 'w')
        f.write("Pos\tA\tC\tG\tT\tOther\tTotal\n")
        max_ref_pos = max(data['nuc_count'][chrom].keys())
        for ref_pos in range(max_ref_pos+1):
            if ref_pos in data['nuc_count'][chrom]:
                c = data['nuc_count'][chrom][ref_pos]
                f.write("%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (ref_pos, c['A'], c['C'], c['G'], c['T'], c['X'], sum(c.values())))
            else:
                f.write("%d\t0\t0\t0\t0\t0\t0\n" % ref_pos)
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
    aln = open_sam(args.input, args.threads)
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

    # generate output files
    print_log("Writing nucleotide counts...")
    write_nuc_counts(data, args.output)
    print_log("Finished writing nucleotide counts")
    print_log("Plotting read length distributions...")
    plot_read_length(data, args.output, ylog=args.ylog)
    print_log("Finished plotting read length distributions")
    print_log("Plotting coverage...")
    plot_coverage(data, args.output, ylog=args.ylog)
    print_log("Finished plotting coverage")

    # finish up
    LOGFILE.close()
