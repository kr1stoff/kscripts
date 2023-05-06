#!/usr/bin/env python

import os, sys
import time
import multiprocessing  
from typing import List, Dict, Tuple
# import pandas as pd


def print_usage():
    print('[ERROR] - ', time.asctime(time.localtime(time.time())))
    print('generate_function_table.py <sample_names.txt> <results/5.function_annotation>')
    sys.exit(1)

def getargs():
    if len(sys.argv) != 3:
        print_usage()
    sample_names, workpath = sys.argv[1], sys.argv[2]
    assert os.path.isdir(workpath), "Don't exist directory: " + workpath
    return(sample_names, workpath)

def merge_table(file:str='kegg.pathway_lv1.xls', col:int=3) -> Tuple:
    with open(sample_names, 'rt') as f:
        sample_list = [line.strip() for line in f]
    temp_func_dict = dict()
    temp_func_list = list()
    for samp in sample_list:
        temp_func_dict[samp] = {}
        kegg_lv1_fraction_file = '{}/{}/{}'.format(workpath, samp, file)
        with open(kegg_lv1_fraction_file, 'rt') as f:
            next(f)
            for line in f:
                if col == 3:
                    func, _, number = line.strip().split('\t')
                else:
                    func, number = line.strip().split('\t')
                temp_func_dict[samp][func] = number
                if func not in temp_func_list:
                    temp_func_list.append(func)
    for samp in temp_func_dict.keys():
        for func in temp_func_list:
            if func not in temp_func_dict[samp].keys():
                temp_func_dict[samp][func] = '0'
    return(temp_func_list, temp_func_dict)

def write_table(func_list:List, func_dict:Dict, outfile:str):
    outlines_list = list()
    outheader_list = ['func']
    for func in func_list:
        outline_list = [func]
        for samp in func_dict.keys():
            outheader_list.append(samp)
            outline_list.append(func_dict[samp][func])
        outline = '\t'.join(outline_list) + '\n'
        outlines_list.append(outline)
    outlines_list.insert(0, '\t'.join(outheader_list[:len(func_dict.keys())+1]) + '\n')
    with open(outfile, 'wt', encoding='utf-8', newline='') as g:
        for line in outlines_list:
            g.write(line)

def my_mkdir(path:str):
    if not os.path.exists(path):
        os.mkdir(path)

#~ design for multiprocessing Pool map
def run_sub1(args_tuple:Tuple):
    infile, outfile, col = args_tuple
    temp_func_list, temp_func_dict = merge_table(file=infile, col=col)
    write_table(temp_func_list, temp_func_dict, workpath + os.sep + outfile)

def run():
    ## mkdir output
    my_mkdir(workpath + os.sep + 'KEGG')
    my_mkdir(workpath + os.sep + 'eggNOG')
    my_mkdir(workpath + os.sep + 'CAZy')
    my_mkdir(workpath + os.sep + 'GO')
    my_mkdir(workpath + os.sep + 'CARD')
    #~ multiprosses
    all_task = [
    ('kegg.pathway_lv1.xls', 'KEGG/kegg.pathway_lv1.merge_table.xls', 3),
    ('eggnog.function_lv2.xls', 'eggNOG/eggnog.function_lv2.merge_table.xls', 3),
    ('cazy_func_lv1.xls', 'CAZy/cazy_func_lv1.merge_table.xls', 3),
    ('kegg.ko.xls', 'KEGG/kegg.ko.merge_table.xls', 2),
    ('eggnog.og.xls', 'eggNOG/eggnog.og.merge_table.xls', 2),
    ('cazy_func_lv2.xls', 'CAZy/cazy_func_lv2.merge_table.xls', 2),
    ('go.xls', 'GO/go.merge_table.xls', 2),
    ('aro.xls', 'CARD/aro.merge_table.xls', 2)
    ]
    with multiprocessing.Pool(processes=8) as pool:
        pool.map(run_sub1, all_task)

if __name__ == '__main__':
    global sample_names, workpath
    sample_names, workpath = getargs()
    run()
    