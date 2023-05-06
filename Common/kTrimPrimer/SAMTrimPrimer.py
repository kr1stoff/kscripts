#!/usr/bin/env python

import os
import sys
import logging
import argparse
import shutil
import re
from collections import Counter
from glob import glob
from subprocess import run
from multiprocessing import Pool
from Bio.Seq import Seq
import pdb


#设置运行日志
logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)s - %(asctime)s - %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S"
)


class SAMTrimPrimer():
    """输入args, 从输入SAM到去引物后的FQ和统计文件"""
    def __init__(self, insam, outfq, primers, offset, include, parallel, minlen):
        self.in_sam = insam
        self.outfq = outfq
        self.primer_tab = primers
        self.offset = offset
        self.include = include
        self.parallel = parallel
        self.minlen = minlen
        self.my_config()
        self.get_primer_info() #引物列表字典 

    def my_config(self):
        self.dir_out = os.path.dirname(self.outfq) if os.path.dirname(self.outfq) != "" else "."
        self.dir_tmp = f"{self.dir_out}/ktmp"
        os.makedirs(self.dir_tmp, exist_ok=True)
        self.basenm = os.path.basename(self.in_sam)

    def get_primer_info(self):
        """1.获取引物信息转到array, 2.创建引物截切字典"""
        self.primers = [] #[["NC123",1,2,3,4], ["NC123",5,6,7,8]...]
        self.dict_trim = {} #{"NC123-1-2-3-4":{"reads":11, "bases":111}, ... }
        with open(self.primer_tab, "rt") as fh:
            next(fh)
            for line in fh:
                listl = line.strip().split("\t")
                self.primers.append([listl[0]] + list(map(int, listl[1:])))
                primer = "-".join(listl)
                self.dict_trim[primer] = {}
                self.dict_trim[primer]["reads"] = 0
                self.dict_trim[primer]["bases"] = 0

    def update_dict_trim(self, key, base:int, dict_trim):
        if key in dict_trim:
            dict_trim[key]["reads"] += 1
            dict_trim[key]["bases"] += base
        else:
            logging.warning(f"update_dict_trim --> {key}")

    def find_best_primer(self, rname, start, end, primers, flag:int):
        """判断最优匹配的引物"""
        sam_ref_name = rname.split(".")[0]
        dict_better_primer = {} #每条read扫描所有引物对,任何有匹配到的引物对和切割信息都记录
        for lrse in primers: #lrse: left right start end
            primer_ref_name = lrse[0].split(".")[0]
            left_start = lrse[1] - self.offset #左引物起始
            left_end = lrse[2] #左引物终止
            right_start = lrse[3] #右引物起始
            right_end = lrse[4] + self.offset #右引物终止
            #1)read左右两端都在引物上 2)read左端在引物上 3)read右端在引物上
            #[230314 逻辑优化] 单端测序如果read没有两边都识别到同一对引物, 默认引物都是从5'端开始切
            if (sam_ref_name == primer_ref_name) and (left_start <= start <= left_end) \
                and (right_start <= end <= right_end): 
                trim_left = left_end - start + 1 #[220714 NC_045512.2:860-884 是25bp]
                trim_right = end - right_start + 1
            elif (sam_ref_name == primer_ref_name) and (left_start <= start <= left_end) \
                and (flag == 0): 
                trim_left = left_end - start + 1
                trim_right = 0
            elif (sam_ref_name == primer_ref_name) and (right_start <= end <= right_end) \
                and (flag & 16 == 16): 
                trim_left = 0
                trim_right = end - right_start + 1
            else:
                trim_left = trim_right = 0
            trim_len = trim_left + trim_right
            #read匹配到primer 
            if trim_len != 0: 
                dict_better_primer[trim_len] = lrse + [[trim_left, trim_right]]
        #如果一条read有多对引物, 认为切割长度最长的是最佳引物
        if dict_better_primer != {}:
            best_primer = dict_better_primer[max(dict_better_primer.keys())]
        else:
            best_primer = None
        return best_primer #best_primer ["NC123",1,2,3,4,[trim_left, trim_right]]

    def trim_primer(self, trim_left, trim_right, seq, qual):
        """返回trim后的seq,qual"""
        seq_len = len(seq)
        slice_right = seq_len - trim_right
        new_seq = seq[trim_left:slice_right]
        new_qual = qual[trim_left:slice_right]
        return new_seq, new_qual

    def flag_seq_qual(self, flag:int, new_seq, new_qual):
        """根据flag,处理seq和qual"""
        bool_continue=False
        seq, qual = new_seq, new_qual
        #flag
        if (flag & 4 == 4) or (flag & 256 == 256) or (flag & 2048 == 2048): #没比对上,不是最优比对,补充序列
            bool_continue = True
        elif flag == 0: #正向比对占位
            pass
        elif flag & 16 == 16: #反向改
            seq = Seq(new_seq).reverse_complement() #序列反向互补
            qual = new_qual[::-1]  #质量反向
        else:   #其他未知情况
            bool_continue = True
            logging.info(f"其他未知情况 FLAG: {flag}")
        #短read过滤
        if len(seq) < self.minlen: 
            bool_continue = True 
        return bool_continue, seq, qual

    def clip_read_by_cigar(self, cigar, iseq, iqual):
        """软剪切部分直接切掉, soft clipping只有三种情况,只可能在两端"""
        patt1 = re.compile("(\d+)S\d+[A-Z]") #左面软剪切
        patt2 = re.compile("\d+[A-Z](\d+)S") #右面软剪切
        patt3 = re.compile("(\d+)S.*?(\d+)S") #两边软剪切
        count_cigar = Counter(cigar)
        oseq, oqual = iseq, iqual
        #判断S在CIGAR上的情况
        if (count_cigar["S"] == 1) and (re.findall(patt1, cigar)): #左
            pos = int(re.findall(patt1, cigar)[0])
            oseq = iseq[pos:]
            oqual = iqual[pos:]
        elif (count_cigar["S"] == 1) and (re.findall(patt2, cigar)): #右
            pos = int(re.findall(patt2, cigar)[0])
            oseq = iseq[:-pos]
            oqual = iqual[:-pos]
        elif (count_cigar["S"] == 2) and (re.findall(patt3, cigar)): #两边
            pos1, pos2 = list(map(int, re.findall(patt3, cigar)[0]))
            oseq = iseq[pos1:-pos2]
            oqual = iqual[pos1:-pos2]
        elif count_cigar["S"] == 0: #CIGAR没有S, 不处理
            pass
        else:
            logging.debug(f"位置的软剪切情况: {cigar}")
        return oseq, oqual

    def parse_sam_line(self, sline):
        """
        解析SAM记录行
        0.QNAME, 1.FLAG, 2.RNAME, 3.POS, 5. CIGAR, 9.SEQ, 10.QUAL
        """
        listl = sline.strip().split("\t")
        readid, flag, rname = listl[:3]
        cigar = listl[5]
        _seq, _qual = listl[9:11]
        seq, qual = self.clip_read_by_cigar(cigar, _seq, _qual)
        start = int(listl[3])
        end = start + len(seq)
        return readid, flag, rname, seq, qual, start, end

    def my_trim(self, io_tuple):
        """
        执行单样本引物剪切, 输入io_tuple (insam,outfq)元组
        """
        in_sam, outfq = io_tuple
        fh = open(in_sam, "rt")
        gh = open(outfq, "wt", encoding="utf-8", newline="")
        wh = open(f"{outfq}.waste.txt", "wt", encoding="utf-8", newline="")
        for line in fh:
            _dict_trim = self.dict_trim
            if line.startswith("@"): continue  #跳过 SAM header
            readid, flag, rname, seq, qual, start, end = self.parse_sam_line(line)
            if len(seq) < self.minlen: continue #[221011] 短read过滤
            best_primer = self.find_best_primer(rname, start, end, self.primers, int(flag))
            if best_primer != None:
                trim_left, trim_right = best_primer[-1]
                new_seq, new_qual = self.trim_primer(trim_left, trim_right, seq, qual)
                primer_key = "-".join(list(map(str, best_primer[:-1])))
                bool_continue, seq, qual = self.flag_seq_qual(int(flag), new_seq, new_qual)
                if bool_continue: continue #不是最佳比对都过滤掉
                self.update_dict_trim(key=primer_key, base=trim_left+trim_right, dict_trim=_dict_trim)
                gh.write(f"@{readid} 1:N:0:0\n{seq}\n+\n{qual}\n")
            elif self.include: #加'--include'参数时
                gh.write(f"@{readid} 1:N:0:0\n{seq}\n+\n{qual}\n")
            else: #短read过滤continue会丢失一些waste reads
                wh.write(f"{readid}\t{flag}\t{best_primer}\n")
            self.trim_stat(_dict_trim, outfq)
        fh.close()
        gh.close()
        wh.close()

    def trim_stat(self, dict_trim, prifix):
        #引物剪切统计
        with open(f"{prifix}.primer_stats.tsv", "wt", encoding="utf-8", newline="") as gh:
            gh.write(f"参考\t起始F\t终止F\t起始R\t终止R\t序列数\t碱基数\n")
            for key in dict_trim:
                reads = str(dict_trim[key]["reads"])
                bases = str(dict_trim[key]["bases"])
                gh.write("\t".join(key.split("-") + [reads, bases]) + "\n")

    def split_sam(self):
        """split拆分SAM"""
        cml = f"split -a 3 -d -l 100000 {self.in_sam} {self.dir_tmp}/{self.basenm}_"
        logging.debug(cml)
        os.system(cml)
        self.pieces = glob(f"{self.dir_tmp}/{self.basenm}_[0-9][0-9][0-9]")

    def trim_primer_pieces(self):
        """多进程运行SAM piece文件"""
        io_tuples = ((pie, f"{pie}.fq") for pie in self.pieces)
        with Pool(self.parallel) as ph:
            ph.map(self.my_trim, io_tuples)

    def merge_fq_stats(self):
        """合并FQ和统计文件"""
        #fq
        cml = f"cat {self.dir_tmp}/{self.basenm}_*.fq > {self.outfq}"
        run(cml, shell=True)
        #waste reads
        cml = f"cat {self.dir_tmp}/{self.basenm}_*.waste.txt > {self.outfq}.waste.txt"
        run(cml, shell=True)
        #stats
        dict_stats = {}
        file_stats = glob(f"{self.dir_tmp}/{self.basenm}_*.primer_stats.tsv")
        for fil in file_stats:
            with open(fil, "rt") as fh:
                next(fh) #header
                for line in fh:
                    listl = line.strip().split("\t")
                    pos_5col = "-".join(listl[:5])
                    reads, base = list(map(int, listl[5:]))
                    if pos_5col not in dict_stats:
                        dict_stats[pos_5col] = {}
                        dict_stats[pos_5col]["reads"] = reads
                        dict_stats[pos_5col]["bases"] = base
                    else:
                        dict_stats[pos_5col]["reads"] += reads
                        dict_stats[pos_5col]["bases"] += base
        self.trim_stat(dict_stats, self.outfq)

    def execute(self):
        """执行分为单样本和多样本并行模式"""
        if self.parallel == 1:
            self.my_trim((self.in_sam, self.outfq))
        else:
            self.split_sam()
            self.trim_primer_pieces()
            self.merge_fq_stats()

def get_args():
    parser = argparse.ArgumentParser()
    parser.description = "去引物程序, 仅支持单端SAM. 小数据量单进程, 大数据量多进程."
    parser.add_argument("-i", "--insam", required=True, help="输入SAM文件.")
    parser.add_argument("-p", "--primers", required=True, 
        help="引物信息文件, 五列chrom,left_start,left_end,right_start,right_end, 详见 fgbio TrimPrimers.")
    parser.add_argument("-o", "--outfq", default="out.ktrim.fq", 
        help="[可选]输出FQ文件. [default: out.ktrim.fq]")
    parser.add_argument("-l", "--minlen", default=40, type=int, 
        help="[可选]允许剪切后read长度的最小阈值,低于该值丢弃. [default: 40]")
    parser.add_argument("--offset", default=1, type=int, help="[可选]允许primer与read间的偏移值. [default: 1]")
    parser.add_argument("--include", action="store_true", 
        help="[可选]默认丢弃掉无primer的read. 如果加上该参数, 则保留没有primer的read.")
    parser.add_argument("-w", "--parallel", type=int, default=1, help="[可选]最大并行数. [default: 1]")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_args()
    if (not os.path.isfile(args.insam)) or (not os.path.isfile(args.primers)):
        logging.error(f"SAM或PRIMERS文件不存在: {args.insam} {args.primers}")
        sys.exit(1)
    pipe = SAMTrimPrimer(
        insam=args.insam,
        outfq=args.outfq,
        primers=args.primers,
        offset=args.offset,
        include=args.include,
        parallel=args.parallel,
        minlen=args.minlen
    )
    pipe.execute()
