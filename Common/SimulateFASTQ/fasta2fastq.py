#!/usr/bin/env python

# @Author   : mengxf
# @Version  : 1.0
# @Create   : 2022/05/06
# @Modified : 2022/05/06

# @用意     : contig/scafford/whole genome 模拟FASTQ数据
# @例子     : 模拟生成1000x以上的SE50新冠FQ数据
#	        : fasta2fastq.py -i ON421856.1.fa -s SE50 -o ON421856.1.fq --min_cov 1000 

import re
import argparse
from Bio import SeqIO


def get_args():
    """获取命令行参数"""
    parser = argparse.ArgumentParser(
        description="Contigs/Scaffords FASTA convert to SE/PE FASTQ."
    )
    parser.add_argument("-i", "--input", help="Input Contigs/Scaffords FASTA File.", required=True)
    parser.add_argument("-s", "--strategy", help="Sequencing strategy, such as [SE100, PE150, PE250...]", required=True)
    parser.add_argument("-o", "--outfq1", help="Output FASTQ, if PE, FASTQ1.", required=True)
    parser.add_argument("-O", "--outfq2", help="Output FASTQ2, if PE.")
    parser.add_argument("--min_cov", help="Min coverage depth. [default: 20]", default=20)
    parser.add_argument("--insert_len", help="Insert DNA length. If not provided, the program automatically predicts.")
    args = parser.parse_args()
    return args

class SimulateFQ():
    """模拟生成单双端FASTQ文件"""
    def __init__(self):
        self.args = get_args()
        self.infile = self.args.input
        self.outfq1 = self.args.outfq1
        self.mincov = int(self.args.min_cov)
        self.strategy = self.args.strategy
        self.assign_seq_strategy()
        self.assign_insert_len()
        self.counter = 0

    def assign_seq_strategy(self):
        """检查测序策略,并拆分"""
        pattern = re.compile("^[PS]E[1-9]\d{1,2}$")
        if re.match(pattern, self.strategy):
            self.seq_type = self.strategy[:2]
            self.seq_length = int(self.strategy[2:])
        else:
            raise Exception(f"错误的测序策略! - {self.strategy}")

    def assign_insert_len(self):
        """插入片段长度"""
        if self.args.insert_len is None:
            self.insert_len = self.seq_length*2-50
        else:
            try:
                self.insert_len = int(self.args.insert_len)
            except ValueError:
                raise ValueError(f"不是数字! - {self.insert_len}")

    def simulate_se(self):
        """模拟单端数据"""
        fake_read_cov = 2 * self.mincov
        stride = self.seq_length / (0.5 * fake_read_cov)
        with open(self.outfq1, "wt", encoding="utf-8", newline="") as gh:
            for record in SeqIO.parse(self.infile, "fasta"):
                seq = str(record.seq).upper()
                seq_rev = str(record.seq.reverse_complement()).upper()
                seq_len = len(seq)
                min_len = min(seq_len, self.seq_length)
                item = -min_len
                while item < seq_len + min_len:
                    item += stride
                    start = max(item, 0)
                    start = min(start, seq_len - min_len)
                    for oseq in (seq, seq_rev):
                        self.counter += 1
                        gh.write(f"@read{self.counter}\n{oseq[int(start):int(start+min_len)]}\n+\n{'H'*min_len}\n")

    def simulate_pe(self):
        """模拟双端数据"""
        fake_read_cov = 2 * self.mincov
        stride = self.insert_len / (0.5 * fake_read_cov)
        with open(self.outfq1, "wt", encoding="utf-8", newline="") as g1h,\
            open(self.args.outfq2, "wt", encoding="utf-8", newline="") as g2h:
            for record in SeqIO.parse(self.infile, "fasta"):
                seq = str(record.seq).upper()
                seq_rev = str(record.seq.reverse_complement()).upper()
                seq_len = len(seq)
                min_len = min(seq_len, self.insert_len)
                item = -min_len
                while item < seq_len + min_len:
                    item += stride
                    start = max(item, 0)
                    start = min(start, seq_len - min_len)
                    for oseq in (seq, seq_rev):
                        self.counter += 1
                        seq1 = oseq[int(start):int(start+self.seq_length)]
                        g1h.write(f"@read{self.counter}/1\n{seq1}\n+\n{'H'*len(seq1)}\n")
                        seq2 = oseq[int(start+min_len-self.seq_length):int(start+min_len)]
                        g2h.write(f"@read{self.counter}/2\n{seq2}\n+\n{'H'*len(seq2)}\n")

    def simulate_execute(self):
        """执行模拟"""
        if self.seq_type == "SE":
            self.simulate_se()
        elif self.seq_type == "PE":
            if self.args.outfq2 is None:
                raise Exception("模拟双端需要参数'-O/--outfq2'!") 
            else:
                self.simulate_pe()
        else:
            raise Exception(f"测序类型不是SE/PE! - {self.seq_type}")

if __name__ == "__main__":
    simu = SimulateFQ()
    simu.simulate_execute()
