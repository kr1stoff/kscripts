#!/usr/bin/env python

import os, sys
import argparse
import pandas as pd

#~ find phylum function
sys.path.append('/nfs-test/mxf/Database/NCBI/utils')
import taxonomy_utils


def run_1(genomecov_out, kraken_out, report, outfile):
    """
    report  -->  report_clean  
    kraken_out  -->  kraken2_classified_clean.txt
    genomecov_out  -->  out.depth_gt2.coverage.tab  
    outfile -->  ./gene_annotation_by_kraken2.txt 
    """
    number_dict = dict()
    with open(genomecov_out, 'rt') as f:
        for line in f:
            genemark_id, number = line.strip().split('\t')
            number_dict[genemark_id] = str(int(float(number)))
    #~ kraken report
    report_dict = dict()
    with open(report, 'rt') as f:
        next(f)
        for line in f:
            tax_id,latin,rank = line.strip().split('\t')
            report_dict[tax_id] = [rank, latin]
    #~ taxa table list
    taxa_table_list = list()
    with open(kraken_out, 'rt') as f, open(outfile, 'wt', encoding='utf-8', newline='') as g:
        g.write('#genemark_id\ttax_id\tlatin\tphylum\trank\tnumber\n')
        next(f)
        for line in f:
            genemark_id, tax_id = line.strip().split('\t')
            rank, latin = report_dict[tax_id]
            #~ just remain 'genemark_id in number_dict'
            if genemark_id not in number_dict.keys(): continue
            number = number_dict[genemark_id]
            #~ phylum 
            phylum_name = taxonomy_utils.find_phynum(tax_id)
            line = '\t'.join([genemark_id, tax_id, latin, phylum_name, rank, number]) + '\n'
            g.write(line)
            taxa_table_list.append([tax_id, latin, rank, number])
    return taxa_table_list

def generate_taxa_table(taxa_table_list, outpath):
    taxa_table_df = pd.DataFrame(taxa_table_list, columns=['#tax_id','latin','rank','number'])
    # taxa_table_df.dtypes
    taxa_table_df = taxa_table_df.astype({'number': int})
    taxa_table_out = taxa_table_df.groupby(['#tax_id', 'latin', 'rank'], as_index=False).sum()
    taxa_table_out['number'] = list(map(lambda x: int(x), taxa_table_out['number']))
    taxa_table_out.to_csv(outpath + os.sep + 'taxa_table_whole.tsv', sep='\t', index=False, encoding='utf-8')
    spec_table_out = taxa_table_out[taxa_table_out['rank'].str.contains('S')][['#tax_id','latin','number']]
    spec_table_out.to_csv(outpath + os.sep + 'taxa_table_species.tsv', sep='\t', index=False, encoding='utf-8')

def get_args():
    parser = argparse.ArgumentParser(
        description='Annotate metagenemark genes by kraken2 output file.'
        )
    parser.add_argument('-n', '--number', required=True, help='Number of gene fragments, bowtie aliments result | bedtools genomecov | number > 2')
    parser.add_argument('-k', '--kraken', required=True, help='Kraken2 output file, just retain classified sequence, kraken2_classified_clean.txt')
    parser.add_argument('-r', '--report', required=True, help='Kraken2 output file, report_clean.')
    parser.add_argument('-o', '--outpath', help='Output directory, default: ./', default='./')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()
    genomecov_out, kraken_out, report, outpath = args.number, args.kraken, args.report, args.outpath
    gene_annotation_outfile = outpath + os.sep + 'gene_annotation_by_kraken2.txt'
    taxa_table_list = run_1(genomecov_out,kraken_out, report, gene_annotation_outfile)
    generate_taxa_table(taxa_table_list, outpath)
