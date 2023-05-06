#!/usr/bin/env python

import sys, os
import math


usage = """
PROG <VCF file> <output prefix>
    [0] VCF file           Structual Variants merged file, created by SURVIVOR. 
    [1] output prefix      Output file prefix, as sample id/name.
"""

def print_usage():
    print(usage)
    sys.exit(1)
args = sys.argv
#~ arguments number is 2, vcf file must be subsistent
if len(args) != 3 or not os.path.exists(args[1]):
    print_usage()
vcf_file, qual_depth_vaf_genot_file = args[1], args[2] + ".qual_dp_vaf_genot.txt"
#~ if output exists, delete
if os.path.exists(qual_depth_vaf_genot_file):
    os.remove(qual_depth_vaf_genot_file)

# #~ run
# with open(vcf_file, 'rt') as f, open(qual_depth_vaf_genot_file, 'wt', newline="", encoding="utf-8") as g:
#     g.write("QUAL\tDEPTH\tVAF\tGENOTYPE\n")
#     for line in f.readlines():
#         if line.startswith('#'):
#             continue
#         linelist = line.strip().split("\t")
#         qual = int(linelist[5])
#         caller1_info, caller2_info = linelist[9:11]
#         caller1_info_list = caller1_info.strip().split(":")
#         genot = caller1_info_list[0]
#         caller2_info_list = caller2_info.strip().split(":")
#         dp1 = tuple(map(int, caller1_info_list[3].split(",")))
#         dp2 = tuple(map(int, caller2_info_list[3].split(",")))
#         dp0 = [math.ceil((dp1[i] + dp2[i]) / 2) for i in range(2)]
#         outline = "{}\t{},{}\t{}\t{}\n".format(qual, dp0[0], dp0[1], round(dp0[1]/sum(dp0), ndigits=4), genot)
#         ## outline
#         # QUAL	DEPTH	VAF	    GENOTYPE
#         # 1	    6,1	    0.1429	0/0
#         # 4	    3,1	    0.25	0/1
#         # 30	    1,6	    0.8571	1/1
#         g.write(outline)
#~ cuteSV
with open(vcf_file, 'rt') as f, open(qual_depth_vaf_genot_file, 'wt', newline="", encoding="utf-8") as g:
    g.write("QUAL\tDEPTH\tVAF\tGENOTYPE\n")
    for line in f.readlines():
        if line.startswith('#'):
            continue
        linelist = line.strip().split("\t")
        qual = linelist[5]
        info = linelist[9]
        info_list = info.strip().split(":")
        genot = info_list[0]
        if info_list[1] == ".":
            continue
        dr, dv = list(map(int, info_list[1:3]))
        outline = "{}\t{},{}\t{}\t{}\n".format(qual, dr, dv, round(dv/(dr+dv), ndigits=4), genot)
        ## outline
        # QUAL	DEPTH	VAF	    GENOTYPE
        # 1	    6,1	    0.1429	0/0
        # 4	    3,1	    0.25	0/1
        # 30	1,6	    0.8571	1/1
        g.write(outline)