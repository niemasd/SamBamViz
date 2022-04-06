#! /usr/bin/env python3
'''
Given a nuc_count_*.tsv file output by SamBamViz, output a stacked barplot of base frequencies at each specified position.

The positions (0-indexed) should be given as comma-separated ranges.
- For example, to plot positions 2, 3, 4, 8, and 9, you could do the following: 2,3,4,8,9
- Or equivalently, with ranges: 2-4,8-9
'''
from matplotlib import rcParams
from matplotlib.lines import Line2D
from os.path import isfile
from seaborn import set_context, set_style
import argparse
import matplotlib
import matplotlib.pyplot as plt

# prep matplotlib/seaborn
matplotlib.use("Agg")
RC = {"font.size":12,"axes.titlesize":16,"axes.labelsize":14,"legend.fontsize":10,"xtick.labelsize":10,"ytick.labelsize":10}
set_context("paper", rc=RC); set_style("ticks"); rcParams['font.family'] = 'serif'

# defaults
DEFAULT_BAR_WIDTH = 0.35
DEFAULT_COLOR_A = 'red'
DEFAULT_COLOR_C = 'blue'
DEFAULT_COLOR_G = 'purple'
DEFAULT_COLOR_T = 'yellow'
DEFAULT_COLOR_X = 'black'

# parse args
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help="Input File (nuc_count_*.tsv)")
    parser.add_argument('-o', '--output', required=True, type=str, help="Output File (PDF)")
    parser.add_argument('-p', '--positions', required=True, type=str, help="Positions")
    parser.add_argument('-ca', '--color_a', required=False, type=str, default=DEFAULT_COLOR_A, help="Color: A")
    parser.add_argument('-cc', '--color_c', required=False, type=str, default=DEFAULT_COLOR_C, help="Color: C")
    parser.add_argument('-cg', '--color_g', required=False, type=str, default=DEFAULT_COLOR_G, help="Color: G")
    parser.add_argument('-ct', '--color_t', required=False, type=str, default=DEFAULT_COLOR_T, help="Color: T")
    parser.add_argument('-cx', '--color_x', required=False, type=str, default=DEFAULT_COLOR_X, help="Color: Other")
    parser.add_argument('-bw', '--bar_width', required=False, type=float, default=DEFAULT_BAR_WIDTH, help="Bar Width")
    parser.add_argument('-hl', '--hide_legend', action="store_true", help="Hide Legend")
    args = parser.parse_args()
    if not isfile(args.input):
        raise ValueError("Input file not found: %s" % args.input)
    if isfile(args.output):
        raise ValueError("Output file exists: %s" % args.output)
    if not args.output.lower().endswith('.pdf'):
        raise ValueError("Output file extension must be PDF: %s" % args.output)
    return args

# parse positions
def parse_positions(positions):
    if positions is None:
        return None
    out = list()
    try:
        for part in positions.split(','):
            if '-' in part:
                a,b = part.split('-'); out += list(range(int(a),int(b)+1))
            else:
                out.append(int(part))
    except:
        raise ValueError("Invalid positions selection: %s" % positions)
    for pos in out:
        if pos < 0:
            raise ValueError("Positions must be non-negative")
    return out

# load data
def load_data(in_fn, positions):
    # set things up
    if in_fn.lower() == 'stdin':
        from sys import stdin as in_f
    else:
        in_f = open(in_fn)
    if positions is not None:
        positions_set = set(positions)
    labels, prop_a, prop_c, prop_g, prop_t, prop_x = list(), list(), list(), list(), list(), list()

    # load data and return
    for line in in_f:
        if line.startswith('Pos\t'):
            continue # skip header line
        parts = line.split('\t')
        if positions is None or int(parts[0]) in positions_set:
            a, c, g, t, x = [int(v) for v in parts[1:6]]; tot = a + c + g + t + x
            labels.append('%s (%d)' % (parts[0].strip(), tot))
            if tot == 0:
                prop_x.append(1)
            else:
                prop_a.append(a/tot); prop_c.append(c/tot); prop_g.append(g/tot); prop_t.append(t/tot); prop_x.append(x/tot)
    return labels, prop_a, prop_c, prop_g, prop_t, prop_x

# plot bars
def plot_bars(out_fn, labels, prop_a, prop_c, prop_g, prop_t, prop_x, color_a=DEFAULT_COLOR_A, color_c=DEFAULT_COLOR_C, color_g=DEFAULT_COLOR_G, color_t=DEFAULT_COLOR_T, color_x=DEFAULT_COLOR_X, bar_width=DEFAULT_BAR_WIDTH, hide_legend=False):
    fig, ax = plt.subplots()
    ax.bar(labels, prop_a, bar_width, color=color_a, label='A')
    ax.bar(labels, prop_c, bar_width, color=color_c, label='C')
    ax.bar(labels, prop_g, bar_width, color=color_g, label='G')
    ax.bar(labels, prop_t, bar_width, color=color_t, label='T')
    ax.bar(labels, prop_x, bar_width, color=color_x, label='Other')
    plt.xticks(rotation=90)
    ax.set_xlabel("Position (Coverage)")
    ax.set_ylabel("Proportion")
    if not hide_legend:
        ax.legend(loc='upper right', bbox_to_anchor=(0.02, -0.01))
    fig.savefig(out_fn, format='pdf', bbox_inches='tight')
    plt.close(fig)

# main content
if __name__ == "__main__":
    args = parse_args()
    positions = parse_positions(args.positions)
    labels, prop_a, prop_c, prop_g, prop_t, prop_x = load_data(args.input, positions)
    plot_bars(args.output, labels, prop_a, prop_c, prop_g, prop_t, prop_x, color_a=args.color_a, color_c=args.color_c, color_g=args.color_g, color_t=args.color_t, color_x=args.color_x, bar_width=args.bar_width, hide_legend=args.hide_legend)
