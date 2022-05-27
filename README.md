# SamBamViz
SamBamViz is a tool for generating visualizations of coverage and base counts across entire viral genomes. It consists of two primary components:

* **[SamBamViz.py](SamBamViz.py):** A command-line cross-platform Python program that takes as input a BAM (sorted or unsorted!) containing reads mapped to a viral reference genome and which outputs a TSV file in which each 0-indexed position of the genome is a row, and the columns contain the counts of each of the bases
* **[Web Visualization Tool](https://niemasd.github.io/SamBamViz/):** A web tool that takes as input a SamBamViz TSV (and optionally, a reference genome FASTA) and which produces user-friendly customizeable visualizations

Example data can be found in the [`example`](example) directory.

## Installation
The command-line SamBamViz program is written in Python 3 and depends on [Pysam](https://github.com/pysam-developers/pysam). You can simply download [SamBamViz.py](SamBamViz.py) to your machine and make it executable:

```bash
pip install pysam # if pysam is not already installed
wget "https://raw.githubusercontent.com/niemasd/SamBamViz/main/SamBamViz.py"
chmod a+x SamBamViz.py
sudo mv SamBamViz.py /usr/local/bin/SamBamViz.py # optional step to install globally
```

The SamBamViz web visualization tool does not require any installation. You can access it from any modern web browser via the following URL: [https://niemasd.github.io/SamBamViz](https://niemasd.github.io/SamBamViz)

## Usage
The command-line SamBamViz program can be utilized as follows:

```
usage: SamBamViz.py [-h] [-i INPUT] [-o OUTPUT] [-q BASE_QUAL] [-m MAP_QUAL] [--force_bam] [--version]

  -h, --help                            show this help message and exit
  -i INPUT, --input INPUT               Input File (SAM/BAM) (default: stdin)
  -o OUTPUT, --output OUTPUT            Output Base Counts (TSV) (default: stdout)
  -q BASE_QUAL, --base_qual BASE_QUAL   Minimum base quality to include base in counts (default: 0)
  -m MAP_QUAL, --map_qual MAP_QUAL      Minimum mapping quality to include read in counts (default: 0)
  --force_bam                           Force BAM Input (otherwise infer from filename) (default: False)
  --version                             Show SamBamViz version (default: False)
```

SamBamViz will output a TSV file in which each row corresponds to a genome position and the columns denote "count of A", "count of C", etc. This TSV file can then be visualized using the [web visualization tool](https://niemasd.github.io/SamBamViz/).

# Notes
## Adding SamBamViz to Existing Pipelines
The command-line SamBamViz program can easily be added to any viral genomics pipeline. Because it is able to read SAM/BAM data streamed via standard input, and because it does ***not*** require the SAM/BAM to be sorted, the SAM/BAM data can be streamed to the command-line SamBamViz program directly from the read mapper to avoid unnecessary slow disk access. Here is an example simple pipeline that maps reads, and then feeds the SAM data to both SamBamViz and `samtools sort` via [`tee`](https://en.wikipedia.org/wiki/Tee_(command)):

```bash
minimap2 -t 6 -a -x sr NC_045512.2.fas sample_R1.fastq.gz sample_R2.fastq.gz | tee >(SamBamViz.py -q 20 -o sample.sambamviz.tsv) | samtools sort -o sample.sorted.bam
```

## Data Privacy
The web visualization tool runs ***client*-side** (not server-side), meaning your SamBamViz TSV data **do *not*** get sent to an external server: everything happens directly on your local machine.

## Credits
[Web visualization tool](https://niemasd.github.io/SamBamViz/) created by [Grant Cheng](https://www.linkedin.com/in/grant-cheng-52171b205/) and [Annie Liu](https://www.linkedin.com/in/anniejiaqiliu/).
