#!/usr/bin/env python

import os, sys
import pdb


def print_uasge():
    print("PROG <merged_table_file> <group_file>")
    print("\tmerge_table_file\t\tmerged function table file, as './results/5.function_annotation/KEGG/kegg.ko.merge_table.xls'")
    print("\tgroup_file\t\t\tgroup file, as './config/group.txt'")
    sys.exit(1)

def run():
    with open(merged_table_file, 'rt') as f:
        linelists = f.readlines()
    with open(group_file, 'rt') as f:
        next(f)
        group_dict = {line.strip().split('\t')[0]: line.strip().split('\t')[1] for line in f}
    sample_list = linelists[0].strip().split('\t')[1:]
    group_list = [group_dict[samp] for samp in sample_list]
    #~ sum number list
    sum_number_list = [0] * len(sample_list)
    # pdb.set_trace()
    for line in linelists[1:]:
        linelist = line.strip().split('\t')
        for idx in range(1,len(linelist)):
            sum_number_list[idx - 1] += int(linelist[idx])
    #~ write lefse_in text 
    outfile = merged_table_file.replace('merge_table.xls','lefse.txt')
    g = open(outfile, 'wt', encoding='utf-8', newline='')
    #~ line 1 group, line2 sample
    g.write('\t'.join(['class'] + group_list) + '\n')
    g.write('\t'.join(['Function'] + sample_list) + '\n')
    for line in linelists[1:]:
        linelist = line.strip().split('\t')
        for idx in range(0,len(sample_list)):
            linelist[idx + 1] = '{}'.format(int(linelist[idx + 1]) / sum_number_list[idx])
        g.write('\t'.join(linelist) + '\n')
    g.close()

if __name__ == '__main__':
    args = sys.argv
    if len(args) != 3 or not os.path.exists(args[1]) or not os.path.exists(args[2]):
        print_uasge()
    global merged_table_file, group_file
    merged_table_file, group_file = args[1:]
    run()
