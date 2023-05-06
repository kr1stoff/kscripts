#!/usr/bin/env python
## /nfs-test/mxf/Software/Kscripts/taxonomy_annotate_fullrank.py lefse.contrast.txt -o lefse.tax_assignment.txt

import argparse
from pprint import pprint
import pdb


def taxonomy_db():
    namesdmp = r"/nfs-test/mxf/Database/Microbiome/Taxonamy_trans/names.scientific_name.dmp"
    nodesdmp = r"/nfs-test/mxf/Database/Microbiome/Taxonamy_trans/nodes.clean.dmp"
    global namesdmp_dict, nodesdmp_dict, rank_dict
    namesdmp_dict = dict()
    with open(namesdmp, 'rt', encoding='utf-8') as f:
        for line in f:
            l_list = line.strip().split('\t')
            namesdmp_dict[l_list[0]] = l_list[1]
    nodesdmp_dict = dict()
    rank_dict = dict()
    with open(nodesdmp, 'rt', encoding='utf-8') as f:
        for line in f:
            l_list = line.strip().split('\t')
            ## rank dictionary
            rank_dict[l_list[0]] = l_list[2]
            if l_list[1] in ['1', '131567', '33154']:
                continue
            nodesdmp_dict[l_list[0]] = l_list[1]

def translate_whole_rank(infile, outfile):
    ## Hold 7 ranks, create 7 ranks list & insert dictionary.
    rank7_list = ['superkingdom','kingdom','phylum','class','order','family','genus','species']
    insert_enum_list = list(enumerate(rank7_list))[2:]
    insert_enum = {l[1]:l[0]-1 for l in insert_enum_list}
    def func_rep1(tax_id):
        lvn_id = nodesdmp_dict[tax_id]
        lvn_name = namesdmp_dict[lvn_id] if lvn_id in namesdmp_dict.keys() else 'undefined'
        if rank_dict[lvn_id] in rank7_list:
            rank_list.append(lvn_name)
            rank_current_list.append(rank_dict[lvn_id])
        return lvn_id
    with open(infile, 'rt', encoding='utf-8') as f, open(outfile, 'wt', encoding='utf-8', newline='') as g:
        g.write('#tax_id\trank_whole_name\n')
        for line in f:
            if line.startswith('#'):
                continue
            tax_id = line.strip().split('\t')[0]
            scientific_name = namesdmp_dict[tax_id]
            rank_current_list = list()
            rank_list = [scientific_name]
            ## if Genus
            if tax_id in nodesdmp_dict.keys():
                lv2_id = func_rep1(tax_id)
                ## if Family
                if lv2_id in nodesdmp_dict.keys():
                    lv3_id = func_rep1(lv2_id)
                    ## if Order
                    if lv3_id in nodesdmp_dict.keys():
                        lv4_id = func_rep1(lv3_id)
                        ## if Class
                        if lv4_id in nodesdmp_dict.keys():
                            lv5_id = func_rep1(lv4_id)
                            ## if Phylum 
                            if lv5_id in nodesdmp_dict.keys():
                                lv6_id = func_rep1(lv5_id)
                                ## if Kingdom 
                                if lv6_id in nodesdmp_dict.keys():
                                    lv7_id = func_rep1(lv6_id)
                                    ## group ???
                                    if lv7_id in nodesdmp_dict.keys():
                                        lv8_id = func_rep1(lv7_id)
                                        ## species group .h..
                                        if lv8_id in nodesdmp_dict.keys():
                                            lv9_id = func_rep1(lv8_id)
                                            ## species subgroup f...
                                            if lv9_id in nodesdmp_dict.keys():
                                                lvx_id = func_rep1(lv9_id)
                                                ## sub.. clade etc..
                                                if lvx_id in nodesdmp_dict.keys():
                                                    lvx1_id = func_rep1(lvx_id)
                                                    ## x2
                                                    if lvx1_id in nodesdmp_dict.keys():
                                                        lvx2_id = func_rep1(lvx1_id)
                                                        ## x3  
                                                        if lvx2_id in nodesdmp_dict.keys():
                                                            lvx3_id = func_rep1(lvx2_id)
            ## fix rank list
            rank_list.reverse()
            for r in rank7_list[2:-1]:
                if r not in rank_current_list:
                    rank_list.insert(insert_enum[r], 'undefined')
            rank_whole_name = ';'.join(rank_list)
            outline = '{}\t{}\n'.format(tax_id, rank_whole_name)
            g.write(outline)
        # pdb.set_trace()

def read_arguments():
    parser = argparse.ArgumentParser(
        description='Annotate taxonomy IDs, output "Kingdom;Phylum;Class;Order;Family;Genus;Species" if complete!'
    )
    parser.add_argument('input', help='Input file, column 1 must be TaxID.')
    parser.add_argument('-o', '--outpath', help='Output file path. [default: "./lefse.tax_assignment.txt"]', default='./lefse.tax_assignment.txt')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    taxonomy_db()
    args = read_arguments()
    infile = args.input
    outfile = args.outpath
    translate_whole_rank(infile, outfile)
    # pprint(rank_dict)