#!/usr/bin/env python

import sys, os
import math
from collections import Counter
import pdb


#K add arguments
usage = """
PROG <VCF file> <output path> <ample ID>
    [0] VCF file            Structual Variants merged file, created by SURVIVOR. 
    [1] output VCF file     Output vcf file, as a123/merged_3callers_filt_remake.vcf.
    [2] Sample ID           Sample name/id.
"""
def print_usage():
    print(usage)
    sys.exit(1)
args = list(sys.argv)
#k arguments number is 2, vcf file must be subsistent
if len(args) != 4 or not os.path.exists(args[1]):
    print_usage()
vcf_file, out_vcf_file, sample_id = args[1:]
#k if output exists, delete
if os.path.exists(out_vcf_file):
    os.remove(out_vcf_file)

#K parameters threshold adjustment
#k AF < upper, WT; AF > lower, HOM; Others, Het
typing_af_threshold = [0.25, 0.75]
#k DP < lower, lowDP; DP > upper; highDP
low_high_depth_threshold = [2, 30]
#k QUAL < qual, lowQUAL
qual_threshold = 1
#k AF < af and GT == "0/0", poorGT
filter_gt_af_threshold = 0.1

#K main
delete_ID_list = ["SUPP", "SUPP_VEC", "SVMETHOD", "STRANDS"]
with open(vcf_file, "rt") as f, open(out_vcf_file, "wt", newline="", encoding="utf-8") as g:
    for line in f:
        #k header
        if line.startswith("##"):
            flag = 0
            for dID in delete_ID_list:
                if dID in line:
                    flag = 1
            if flag == 0:
                g.write(line)
        elif line.startswith("#CHROM"):
            linelist = line.strip().split("\t")
            outline = "\t".join(linelist[:9] + [sample_id]) + "\n"
            g.write(outline)
        #k body
        else:
            linelist = line.strip().split("\t")
            #k 6 columns: chrom, pos, id, ref, alt, qual    
            basic_cols = linelist[:6]
            qual = float(linelist[5])
            filt, info = linelist[6:8]
            cute, snif, svim = linelist[9:12]
            #k info
            info_list = info.strip().split(";")
            info_new_list = [it for it in info_list if it.split("=")[0] not in delete_ID_list]
            info_new = ";".join(info_new_list)
            #k format
            fmt = "GT:DR:ST"
            #k sample detail
            #~ 0/1    4,3    +-
            dr_list, dv_list, st_list = list(), list(), list()
            def parse_detail(detail=cute):
                xlist = detail.strip().split(":")
                dp, st = xlist[3:5]
                st_list.append(st)
                if "." not in dp and dp != "0,0":
                    dr, dv = list(map(int, dp.strip().split(",")))
                    dr_list.append(dr)
                    dv_list.append(dv)
            #~cute
            parse_detail()
            #~ sniffles
            parse_detail(detail=snif)
            #~ svim
            parse_detail(detail=svim)
            #~ merge
            dr_mean = math.ceil(sum(dr_list)/len(dr_list))
            dv_mean = math.ceil(sum(dv_list)/len(dv_list))
            af = dv_mean / (dr_mean + dv_mean)
            #~ genotype
            if af < typing_af_threshold[0]:
                genot = "0/0"
            elif typing_af_threshold[0] <= af < typing_af_threshold[1]:
                genot = "0/1"
            else:
                genot = "1/1"
            #~ strands
            counter = dict(Counter(st_list))
            strands = sorted(counter.items(), key=lambda x:x[1], reverse=True)[0][0]
            sample_detail = "{}:{},{}:{}".format(genot, dr_mean, dv_mean, strands)
            #k filter:  lowDR, highDR, lowQUAL, poorGT
            #~ lowDR, highDR
            filt = "PASS"
            if dr_mean + dv_mean < low_high_depth_threshold[0]:
                filt = "lowDR"
            elif dr_mean + dv_mean > low_high_depth_threshold[1]:
                filt = "highDR"
            #~ lowQUAL
            if qual < qual_threshold:
                if filt == "PASS":
                    filt = "lowQUAL"
                else:
                    filt = filt + ";lowQUAL"
            #~ poorGT
            if genot == "0/0" and af < filter_gt_af_threshold:
                if filt == "PASS":
                    filt = "poorGT"
                else:
                    filt = filt + ";poorGT"
            #k make up
            #~ chrom, pos, id, ref, alt, qual, filter, info, format, sample
            outlist = basic_cols + [filt, info_new, fmt, sample_detail]
            g.write("\t".join(outlist) + "\n")
# pdb.set_trace()
