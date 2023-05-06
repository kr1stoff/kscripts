#!/usr/bin/env python

import os
from collections import defaultdict

##### configure #####
# sintax = "./results/reads.sintax"
# hitsuc = "./results/hits.uc"
# outdir = "./results"
# otu_lineage_reads_file = os.path.join(outdir, "otu_lineage_reads.tsv")
# genus_reads_file = os.path.join(outdir, "genus_reads.tsv")
##### configure snakemake #####
sintax = snakemake.input.sintax
hitsuc = snakemake.input.hitsuc
otu_lineage_reads_file = snakemake.output.otu_lineage_reads
genus_reads_file = snakemake.output.genus_reads

#~ otu - lineage dictionary
otu_lineage_dict = dict()
with open(sintax, "rt") as f:
    for line in f:
        llist = line.strip().split("\t")
        otu, lineage = llist[0], llist[3]
        otu_lineage_dict[otu] = lineage
#~ otu - reads dictionary
otu_reads_dict = defaultdict(int)
with open(hitsuc, "rt") as f:
    for line in f:
        otu = line.strip().split("\t")[9]
        if otu == "*":
            continue
        else:
            otu_reads_dict[otu] += 1
#~ extract genus and write otu - lineage - reads table .tsv
genus_reads_dict = defaultdict(int)
with open(otu_lineage_reads_file, "wt", newline="", encoding="utf-8") as g:
    for otu in otu_reads_dict:
        reads = otu_reads_dict[otu]
        lineage = otu_lineage_dict[otu]
        outline = "\t".join([otu, lineage, str(reads)]) + "\n"
        g.write(outline)
        #` genus - reads dictionary
        lineage_list = lineage.strip().split(",")
        genus_pool = list(filter(lambda x: x.__contains__("g:"), lineage_list))
        if genus_pool == []:
            continue
        else:
            #` remove "g:"
            genus = genus_pool[0].split(":")[1]
            genus_reads_dict[genus] += otu_reads_dict[otu]
#~ write genus - reads table .tsv
with open(genus_reads_file, "wt", newline="", encoding="utf-8") as g:
    for genus in genus_reads_dict:
        outline = "\t".join([genus, str(genus_reads_dict[genus])]) + "\n"
        g.write(outline)
