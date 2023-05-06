#!/usr/bin/env python

import os
import sys
import logging
import argparse
import yaml
from SAMTrimPrimer import SAMTrimPrimer
import pdb


# 设置运行日志
logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)s - %(asctime)s - %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S"
)

class FQTrimPrimer():
    def __init__(self, fastq1, fastq2, reference, outfq, primers, minlen, offset, include, parallel):
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.outfq = outfq
        self.ref = reference
        self.primers = primers
        self.minlen = minlen
        self.offset = offset
        self.include = include
        self.parallel = parallel
        self.my_config()

    def my_config(self):
        logging.info("读取输入配置文件")
        # path
        self.dir_out = os.path.dirname(self.outfq) if os.path.dirname(self.outfq) != "" else "."
        self.dir_tmp = f"{self.dir_out}/ktmp"
        os.makedirs(self.dir_tmp, exist_ok=True)
        self.dir_conf = os.path.join(os.path.dirname(__file__), "../conf")
        self.software_config()

    def software_config(self):
        with open(f"{self.dir_conf}/software.yaml", "rt") as fh:
            self.dict_input = yaml.safe_load(fh)
        self.vsearch = self.dict_input["vsearch"]
        self.bwa = self.dict_input["bwa"]

    def prepare_fq(self):
        """如果双端FQ,vsearch合并,单端软链接"""
        logging.info("前处理输入FQ文件")
        if self.fastq2:
            cml = f"{self.vsearch} --fastq_mergepairs {self.fastq1} --reverse {self.fastq2} "\
                f"--fastqout {self.dir_tmp}/prepared.fq"
        else:
            cml = f"ln -sf {os.path.abspath(self.fastq1)} {self.dir_tmp}/prepared.fq"

        logging.debug(cml)
        os.system(cml)

    def alignment(self):
        logging.info("比对")
        _threads = self.parallel if self.parallel != 1 else 8
        cml = f"{self.bwa} mem -t {_threads} -M -Y -R '@RG\\tID:trim\\tSM:trim' "\
            f"{self.ref} {self.dir_tmp}/prepared.fq > {self.dir_tmp}/raw.sam"
        logging.debug(cml)
        os.system(cml)

    def trim_sam(self):
        logging.info("剪切SAM引物生成FQ")
        trim = SAMTrimPrimer(
            insam=f"{self.dir_tmp}/raw.sam",
            outfq=self.outfq,
            primers=self.primers,
            offset=self.offset,
            include=self.include,
            parallel=self.parallel,
            minlen=self.minlen
        )
        trim.execute()

    def remove_tmp(self):
        logging.info("删除中间文件")

    def execute(self):
        self.prepare_fq()
        self.alignment()
        self.trim_sam()

def get_args():
    parser = argparse.ArgumentParser()
    parser.description = "去引物程序,原始FQ到去引物后FQ.小数据量单进程,大数据量多进程."
    parser.add_argument("-i", "--fastq1", required=True, help="单端FQ或双端FQ1.")
    parser.add_argument("-r", "--reference", required=True, help="参考基因组,与引物信息文件参考一致,建好索引.")
    parser.add_argument("-p", "--primers", required=True, 
        help="引物信息文件,五列chrom,left_start,left_end,right_start,right_end,详见fgbio TrimPrimers.")
    parser.add_argument("-o", "--outfq", default="./out.ktrim.fq", 
        help="[可选]输出FQ文件. [default: ./out.ktrim.fq]")
    parser.add_argument("-I", "--fastq2", help="[可选]双端FQ2.")
    parser.add_argument("-l", "--minlen", default=40, type=int, 
        help="[可选]允许剪切后read长度的最小阈值,低于该值丢弃. [default: 40]")
    parser.add_argument("--offset", default=1, type=int, help="[可选]允许primer与read间的偏移值. [default: 1]")
    parser.add_argument("--include", action="store_true", 
        help="[可选]默认丢弃掉无primer的read.如果加上该参数,则保留没有primer的read.")
    parser.add_argument("-w", "--parallel", type=int, default=1, help="[可选]最大并行数. [default: 1]")
    args = parser.parse_args()
    return args

if __name__ =="__main__":
    args = get_args()
    pipe = FQTrimPrimer(
        fastq1=args.fastq1,
        fastq2=args.fastq2,
        reference=args.reference,
        outfq=args.outfq,
        primers=args.primers,
        offset=args.offset,
        include=args.include,
        parallel=args.parallel,
        minlen=args.minlen
    )
    pipe.execute()
