#!/usr/bin/env python

import sys, os
from collections import defaultdict
import pdb


## get arguments
def print_usage():
    print("Usage:\n\tPROG <results_dir>")
    print("\nExample:  gene_prediction_stats.py ./results")
    sys.exit(1)

def rank_number():
    rank_levels_list = ["R", "K", "P", "C", "O", "F", "G", "S"]
    rank_number_dict = defaultdict(int)
    with open(kraken_report_file, "rt") as f:
        for line in f:
            linelist = line.strip().split("\t")
            #~ linelist: [1] percent [2] number clade rooted [3] number directly 
            #~ [4] rank code [5] tax id [6] scientific name
            taxa_number, rank_code = linelist[1], linelist[3]
            if rank_code in rank_levels_list:
                rank_number_dict[rank_code] += int(taxa_number)
    rank_number_dict["K"] += rank_number_dict["P"]
    rank_number_dict["U"] = rank_number_dict["R"] - rank_number_dict["K"]
    rank_number_list = [rank_number_dict[r] for r in rank_levels_list[1:]]
    rank_number_list.insert(0, rank_number_dict["U"])
    rank_percent_list = list(map(lambda x: "{:.2%}".format(x / rank_number_dict["R"]), rank_number_list))
    return rank_percent_list, rank_number_dict

def get_seqs_num():
    with open(seqkit_stats_file, "rt") as f:
        next(f)
        linelist = f.read().strip().split('\t')
    num_seqs = linelist[3]
    return num_seqs
if __name__ == "__main__":
    args = sys.argv
    assert len(args) == 2, print_usage()
    global kraken_report_file
    kraken_report_file = args[1] + os.sep + "4.species_annotation/merge.kraken2_report"
    seqkit_stats_file = args[1] + os.sep + "3.gene_prediction/metagenemark_out/gene_catalogue.uniq_gene_merge.fa.stats"
    outfile = args[1] + os.sep + "4.species_annotation/species_annotation_stats.xls"
    for tmp_file in [kraken_report_file, seqkit_stats_file]:
        if not os.path.exists(tmp_file):
            print_usage()
    rank_percent_list, rank_number_dict = rank_number()
    gene_catalogue_number = get_seqs_num()
    column1 = ["Gene catalogue", 
                "Annotated on Kraken2", 
                "Annotated on Unclassified", 
                "Annotated on Kingdom level",
                "Annotated on Phylum level",
                "Annotated on Class level",
                "Annotated on Order level",
                "Annotated on Family level",
                "Annotated on Genus level",
                "Annotated on Species level"]
    column2 = [gene_catalogue_number, 
                "{}({:.2%})".format(rank_number_dict["R"], rank_number_dict["R"] / int(gene_catalogue_number))] +\
                rank_percent_list
    with open(outfile, "wt", encoding="utf-8", newline="") as g:
        for item in zip(column1, column2):
            outline = "\t".join(item) + "\n"
            g.write(outline)
