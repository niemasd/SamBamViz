# SamBamViz
SamBamViz is a tool for generating visualizations of coverage and base counts across entire viral genomes. It consists of two primary components:

* **[`SamBamViz.py`](SamBamViz.py):** A command-line cross-platform Python program that takes as input a BAM (sorted or unsorted!) containing reads mapped to a viral reference genome and which outputs a TSV file in which each 0-indexed position of the genome is a row, and the columns contain the counts of each of the bases
* **[Web Visualization Tool](https://niemasd.github.io/SamBamViz/):** A web tool that takes as input a SamBamViz TSV (and optionally, a reference genome FASTA) and which produces user-friendly customizeable visualizations

Example data can be found in the [`example`](example) directory.

# Notes
The web visualization tool runs ***client*-side** (not server-side), meaning your SamBamViz TSV data **do *not*** get sent to an external server: everything happens directly on your local machine.

Credits: [Grant Cheng](https://www.linkedin.com/in/grant-cheng-52171b205/) and [Annie Liu](https://www.linkedin.com/in/anniejiaqiliu/)
