#! /usr/bin/env python

import os
import argparse
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(
        description='Generate taxa table by taxa_table_species.tsv.'
        )
    parser.add_argument('-s', '--sample', required=True, help='sample list, ./config/sample_names.txt')
    parser.add_argument('-w', '--workpath', help='work directory, default: ./results/4.species_annotation', default='./results/4.species_annotation')
    args = parser.parse_args()
    return args

def run(sample_file, workpath):
    """
    sample_file  -->  'config/sample_names.txt'
    workpath  -->  './results/4.species_annotation'
    """
    taxa_table_dict = dict()
    contrast_dict = dict()
    with open(sample_file, 'rt') as f:
        for line in f:
            sample_name = line.strip()
            single_taxa_table = workpath + os.sep + sample_name + os.sep + 'taxa_table_species.tsv'
            taxa_table_dict[sample_name] = {}
            with open(single_taxa_table, 'rt') as f:
                next(f)
                for line in f:
                    tax_id, latin, number = line.strip().split('\t')
                    taxa_table_dict[sample_name][tax_id] = number
                    if tax_id not in contrast_dict.keys(): contrast_dict[tax_id] = latin
    #~ pandas merge samples table
    taxa_table_df = pd.DataFrame(taxa_table_dict)
    taxa_table_df.fillna(0, inplace=True)
    tax_id_df = pd.DataFrame(list(taxa_table_df.index), index=taxa_table_df.index, columns=['tax_id'])
    taxa_table_df = pd.concat([tax_id_df, taxa_table_df], axis=1)
    #~ DataFrame to csv taxa_table.txt
    taxa_table_df.to_csv(workpath + os.sep + 'taxa_table.txt', sep='\t', index=False, encoding='utf-8')
    #~ write contrast.txt 
    with open(workpath + os.sep + 'contrast.txt', 'wt', encoding='utf-8', newline='') as g:
        g.write('tax_id\tlatin_name\n')
        for tax_id in contrast_dict.keys():
            latin_name = contrast_dict[tax_id]
            line = '{}\t{}\n'.format(tax_id, latin_name)
            g.write(line)


if __name__ == "__main__":
    args = get_args()
    sample_file, workpath = args.sample, args.workpath
    run(sample_file, workpath)
