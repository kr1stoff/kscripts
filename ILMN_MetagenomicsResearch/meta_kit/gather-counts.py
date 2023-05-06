#!/usr/bin/env python
"""
This script gathers & converts Salmon output counts into something that
edgeR can read ("counts files").

Run it in a directory above all of your Salmon output directories, and
it will create a bunch of '.counts' files that you can load into R.

See https://github.com/ngs-docs/2015-nov-adv-rna/ for background info.

C. Titus Brown, 11/2015
"""
import os, os.path
import sys
import csv

def process_quant_file(root, filename, outname):
    """
    Convert individual quant.sf files into .counts files (transcripts\tcount).
    """
    print('Loading counts from:', root, filename, file=sys.stderr)
    outfp = open(outname, 'w')
    print("transcript\tcount", file=outfp)

    d = {}
    full_file = os.path.join(root, filename)
    for line in open(full_file):
        if line.startswith('Name'):
            continue
        name, length, eff_length, tpm, count = line.strip().split('\t')
        print("%s\t%s" % (name, float(tpm)), file=outfp)


def main(start_dir='.'):
    """
    Find all the quant.sf files, convert them into properly named .counts
    files.

    Here, "proper name" means "directory.counts".
    """
    quantlist = []
    
    print('Starting in:', os.path.abspath(start_dir), file=sys.stderr)
    for root, dirs, files in os.walk(start_dir):
        for filename in files:
            if filename.endswith('quant.sf'):
                dirname = os.path.basename(root)
                outname = start_dir + os.sep + dirname + '.counts'
                # outname = dirname + '.counts'
                process_quant_file(root, filename, outname)
                quantlist.append(outname)
                
                break

    print(",\n".join([ "\"%s\"" % i for i in sorted(quantlist)]))

if __name__ == '__main__':
    if len(sys.argv) != 2: print('Just need 1 argument! [e.g. ./results/3.gene_prediction/salmon/SRR11575977]'); sys.exit(1)
    start_dir=sys.argv[1]
    main(start_dir)
