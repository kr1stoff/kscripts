#!/usr/bin/env python

import sys, os
import subprocess


## get arguments
def print_usage():
    print("Usage:\n\tPROG <config> <results_dir>")
    print("\nExample:  gene_prediction_stats.py ./config ./results")
    sys.exit(1)

def calc_gc_percent():
    total_length, gc_number = [0, 0]
    with open(seqtk_comp_file, "rt") as f:
        for line in f:
            linelist = line.strip().split("\t")
            total_length += int(linelist[1])
            gc_number += (int(linelist[3]) + int(linelist[4]))
    gc_percent = "{:.2%}".format(gc_number/total_length)
    return gc_percent

def get_seqs_num_len():
    with open(seqkit_stats_file, "rt") as f:
        next(f)
        linelist = f.read().strip().split('\t')
    num_seqs, sum_len, avg_len = linelist[3:5] + [linelist[6]]
    sum_len_Mb = str(round(int(sum_len) / (10 ** 6), ndigits=2))
    return num_seqs, sum_len_Mb, avg_len

def get_orfs_num():
    #~ whole
    cml = "grep -c '^>' " + orfs_fasta_file
    subp_return = subprocess.run(cml, shell=True, capture_output=True, check=True, encoding="utf-8")
    orfs_number = subp_return.stdout.strip()
    #~ average
    with open(sample_names_file, "rt") as f:
        linelists = [line for line in f]
    sample_number = len(linelists)
    orfs_number_avg = str(int(int(orfs_number) / sample_number))
    return orfs_number, orfs_number_avg

if __name__ == "__main__":
    args = sys.argv
    assert len(args) == 3, print_usage()
    global seqtk_comp_file, seqkit_stats_file, orfs_fasta_file, sample_names_file
    seqtk_comp_file = args[2] + os.sep + "3.gene_prediction/metagenemark_out/gene_catalogue.uniq_gene_merge.fa.comp"
    seqkit_stats_file = args[2] + os.sep + "3.gene_prediction/metagenemark_out/gene_catalogue.uniq_gene_merge.fa.stats"
    orfs_fasta_file = args[2] + os.sep + "3.gene_prediction/metagenemark_out/mgm.gene_L100_merge.fa"
    sample_names_file = args[1] + os.sep +  "sample_names.txt"
    for tmp_file in [seqtk_comp_file, seqkit_stats_file, orfs_fasta_file, sample_names_file]:
        if not os.path.exists(tmp_file):
            print_usage()
    #~ run
    gc_percent = calc_gc_percent()
    num_seqs, sum_len_Mb, avg_len = get_seqs_num_len()
    orfs_number, orfs_number_avg = get_orfs_num()
    column1 = ["Total ORFs", "Average ORFs", "Gene catalogue", "Total length (Mbp)", "Average length (bp)", "GC percent"]
    column2 = [orfs_number, orfs_number_avg, num_seqs, sum_len_Mb, avg_len, gc_percent]
    outfile = args[2] + os.sep + "3.gene_prediction/metagenemark_out/gene_prediction_stats.xls"
    with open(outfile, "wt", encoding="utf-8", newline="") as g:
        for item in zip(column1, column2):
            g.write("\t".join(item) + "\n")
