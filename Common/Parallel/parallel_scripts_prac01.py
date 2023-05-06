#!/home/ionadmin/software/miniconda3/envs/data_analysis/bin/python
# -*- coding:utf-8 -*-

import os, sys
import argparse
from subprocess import run
from multiprocessing import Pool


def some_func(cmd):
    print("Command Line: " + cmd)
    a = run(cmd, shell=True, check=True, capture_output=True, encoding='utf-8')
    print(a.stdout)

def main():
    cmd_list = ['echo NO.{}'.format(i) for i in range(16)]
    # print(cmd_list)
    pool = Pool(5)
    # pool.map(some_func, cmd_list)
    pool.map_async(some_func, cmd_list)
    print('Waiting for all subprocesses done...')
    pool.close()
    pool.join()
    print('All subprocesses done.')

main()
