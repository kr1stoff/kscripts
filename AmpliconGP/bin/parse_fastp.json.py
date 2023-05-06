#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/7/26
# 统计fastp数据 详情

import logging
import os,re
import click
import json
import pandas as pd
from numpy import *
import sys



#运行日志格式设置：时间+运行提示日志
logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)
help = dict(help_option_names=['-h', '--help'])

def parse_json(id,injson,out,mode):
    with open(injson) as f:
        data = json.load(f)

        #基础过滤前后reads统计
        total_reads = int(data['summary']['before_filtering']['total_reads']) 
        clean_reads = int(data['summary']['after_filtering']['total_reads'])
        total_bases = int(data['summary']['before_filtering']['total_bases']) 
        clean_bases = int(data['summary']['after_filtering']['total_bases'])

        raw_q20 = float(data['summary']['before_filtering']['q20_rate'])
        raw_q30 = float(data['summary']['before_filtering']['q30_rate'])
        clean_q20 = float(data['summary']['after_filtering']['q20_rate'])
        clean_q30 = float(data['summary']['after_filtering']['q30_rate'])
        raw_q20_base = format(int(data['summary']['before_filtering']['q20_bases']),',')
        raw_q30_base = format(int(data['summary']['before_filtering']['q30_bases']),',')
        clean_q20_base = format(int(data['summary']['after_filtering']['q20_bases']),',')
        clean_q30_base = format(int(data['summary']['after_filtering']['q30_bases']),',')

        # 各碱基
        raw_A = mean(data['read1_before_filtering']['content_curves']['A'])
        raw_G = mean(data['read1_before_filtering']['content_curves']['G'])
        raw_C = mean(data['read1_before_filtering']['content_curves']['C'])
        raw_T = mean(data['read1_before_filtering']['content_curves']['T'])
        raw_N = mean(data['read1_before_filtering']['content_curves']['N'])
        clean_A = mean(data['read1_after_filtering']['content_curves']['A'])
        clean_G = mean(data['read1_after_filtering']['content_curves']['G'])
        clean_C = mean(data['read1_after_filtering']['content_curves']['C'])
        clean_T = mean(data['read1_after_filtering']['content_curves']['T'])
        clean_N = mean(data['read1_after_filtering']['content_curves']['N'])


        raw_mean_length = int(data['summary']['before_filtering']['read1_mean_length'])
        clean_mean_length = int(data['summary']['after_filtering']['read1_mean_length'])

        if mode == 'PE':
            raw_read2_mean_length =  int(data['summary']['before_filtering']['read2_mean_length'])
            clean__read2_mean_length = int(data['summary']['after_filtering']['read2_mean_length'])

            raw_mean_length = int((raw_mean_length + raw_read2_mean_length)/2)
            clean_mean_length = int((clean_mean_length + clean__read2_mean_length)/2)

            # 各碱基
            raw_A = (raw_A + mean(data['read2_before_filtering']['content_curves']['A']))/2
            raw_G = (raw_G + mean(data['read2_before_filtering']['content_curves']['G']))/2
            raw_C = (raw_C + mean(data['read2_before_filtering']['content_curves']['C']))/2
            raw_T = (raw_T + mean(data['read2_before_filtering']['content_curves']['T']))/2
            raw_N = (raw_N + mean(data['read2_before_filtering']['content_curves']['N']))/2
            clean_A = (clean_A + mean(data['read2_after_filtering']['content_curves']['A']))/2
            clean_G = (clean_G + mean(data['read2_after_filtering']['content_curves']['G']))/2
            clean_C = (clean_C + mean(data['read2_after_filtering']['content_curves']['C']))/2
            clean_T = (clean_T + mean(data['read2_after_filtering']['content_curves']['T']))/2
            clean_N = (clean_N + mean(data['read2_after_filtering']['content_curves']['N']))/2

        #过滤具体统计
        low_quality_reads = int(data['filtering_result']['low_quality_reads'])
        N_reads = int(data['filtering_result']['too_many_N_reads'])
        low_complexity_reads = int(data['filtering_result']['low_complexity_reads'])
        dup_rate = float(data['duplication']['rate'])/100
        too_short_reads=int(data['filtering_result']['too_short_reads'])

        try:
            adapter_reads = float(data['adapter_cutting']['adapter_trimmed_reads'])
            adapter_rate = (adapter_reads/int(total_reads)) if total_reads > 0 else 0
        except: adapter_reads,adapter_rate = 0,0


    ###################
    #写入基础统计结果
    with open(f"{out}/{id}.basic.stat.txt",'w',encoding='utf-8') as w:
        w.write(f"过滤质控\t过滤前\t过滤后\n")
        w.write(f"总reads数\t{format(total_reads,',')}\t{format(clean_reads,',')}\n")
        w.write(f"总碱基数\t{format(total_bases,',')}\t{format(clean_bases,',')}\n")
        w.write(f"序列平均长度\t{raw_mean_length}\t{clean_mean_length}\n")
        w.write(f"A碱基\t{format(int(raw_A*total_bases),',')} ({raw_A:.2%})\t{format(int(raw_A*total_bases),',')} ({clean_A:.2%})\n")
        w.write(f"G碱基\t{format(int(raw_G*total_bases),',')} ({raw_G:.2%})\t{format(int(raw_G*total_bases),',')} ({clean_G:.2%})\n")
        w.write(f"C碱基\t{format(int(raw_C*total_bases),',')} ({raw_C:.2%})\t{format(int(raw_C*total_bases),',')} ({clean_C:.2%})\n")
        w.write(f"T碱基\t{format(int(raw_T*total_bases),',')} ({raw_T:.2%})\t{format(int(raw_T*total_bases),',')} ({clean_T:.2%})\n")
        w.write(f"N碱基\t{format(int(raw_N*total_bases),',')} ({raw_N:.2%})\t{format(int(raw_N*total_bases),',')} ({clean_N:.2%})\n")
        w.write(f"Q20\t{raw_q20_base} ({raw_q20:.2%})\t{clean_q20_base} ({clean_q20:.2%})\n")
        w.write(f"Q30\t{raw_q30_base} ({raw_q30:.2%})\t{clean_q30_base} ({clean_q30:.2%})\n")
        w.write(f"低质量reads数\t{format(low_quality_reads,',')}\t-\n")
        w.write(f"含N碱基过多的reads数\t{format(N_reads,',')}\t-\n")
        w.write(f"低复杂度reads数\t{format(low_complexity_reads,',')}\t-\n")
        w.write(f"重复reads比例\t{dup_rate:.2%}\t{dup_rate:.2%}\n")
        w.write(f"adapter\t{format(adapter_reads,',')} ({adapter_rate:.2%})\t-\n")


    #写入其他指标
    with open(f"{out}/{id}.detail.stat.txt",'w',encoding='utf-8_sig') as w:
        w.write(f"样本编号\t低质量reads数\t含N碱基过多的reads数\t低复杂度reads数\t重复reads比例\tadapter\n")
        w.write(f"{id}\t{format(low_quality_reads,',')}\t{format(N_reads,',')}\t{format(low_complexity_reads,',')}\t{dup_rate:.2%}\t{format(adapter_reads,',')} ({adapter_rate:.2%})\n")
    logging.info(f"Output {out}/{id}.basic.stat.txt {id}.detail.stat.txt")



#参数
@click.command(context_settings=help)
@click.option('-id',required=True,type=click.STRING,help="sample id")
@click.option('-i', '--injson',required=True,type=click.Path(),help="fastp.json路径")
@click.option('-m','--mode',required=True,type=click.Choice(['PE','SE']),default='PE',show_default=True,help="测序模式选择,单/双端")



def main(id,injson,mode):
    """
    Parse fastp json file
    """
    indir=os.path.dirname((os.path.abspath(injson)))

    parse_json(id,injson,indir,mode)



#调用全局函数
if __name__ == "__main__":
    main()
