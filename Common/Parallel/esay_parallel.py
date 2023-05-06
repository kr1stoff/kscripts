#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import sys, os
from multiprocessing import Pool


## 最简单的 OS command 命令执行语句
def run_cmd(cmd):
    os.system(cmd)

## 随便生成一些命令
cmd_list = list()
for fl in os.listdir('/xx/xx'):
    cm = 'ls ' + fl
    cmd_list.append(cm)

## Pool() 线程拉满，如果进程数大于 cpu_count()，有多少线程有用多少线程
## Pool(processes=n), n 自选线程数
if __name__ == '__main__':
    with Pool() as p:
        p.map(run_cmd, cmd_list)
