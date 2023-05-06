#!/usr/bin/env python

import os, sys


assert len(sys.argv) == 3, 'PROG <out.emapper.annotations> <unigenes_id.protein.txt>'
merge_emapper_outfile, uniqgenes_file = sys.argv[1], sys.argv[2]
# print(merge_emapper_outfile, uniqgenes_file, sep='\t')

emapper_dict = dict()
with open(merge_emapper_outfile, 'rt', encoding='utf-8') as f:
    for line in f:
        gene_id = line.strip().split('\t')[0]
        emapper_dict[gene_id] = line

with open(uniqgenes_file, 'rt', encoding='utf-8') as f:
    for line in f:
        gene_id = line.strip()
        if gene_id in emapper_dict.keys():
            print(emapper_dict[gene_id], end='')
