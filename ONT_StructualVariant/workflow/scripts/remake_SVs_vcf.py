#!/usr/bin/env python

import sys, os
import re


## configure
vcf_file = snakemake.input[0]
out_vcf_file = snakemake.output[0]
sample_id = snakemake.params[0]
#~ test
# vcf_file = sys.argv[1]
# out_vcf_file = sys.argv[2]
# sample_id = sys.argv[3]
#~ check VCF file exist
assert os.path.exists(vcf_file), "No exists VCF file!"
#~ if output exists, delete
if os.path.exists(out_vcf_file):
    os.remove(out_vcf_file)

## parameters threshold adjustment
#~ AF < upper, WT; AF > lower, HOM; Others, Het
typing_af_threshold = [0.25, 0.75]
#~ DP < lower, lowDP; DP > upper; highDP
low_high_depth_threshold = [2, 30]
#~ QUAL < qual, lowQUAL
qual_threshold = 1
#~ AF < af and GT == "0/0", poorGT
filter_gt_af_threshold = 0.1

## main
with open(vcf_file, "rt") as f, open(out_vcf_file, "wt", newline="", encoding="utf-8") as g:
    for line in f:
        #~ header
        if line.startswith("##"):
            g.write(line)
        elif line.startswith("#CHROM"):
            linelist = line.strip().split("\t")
            outline = "\t".join(linelist[:9] + [sample_id]) + "\n"
            g.write(outline)
        #~ body
        else:
            linelist = line.strip().split("\t")
            #`CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NULL
            #~ 6 columns: chrom, pos, id, ref, alt, qual    
            basic_cols = linelist[:6]
            qual = float(linelist[5])
            filt, info, format, detail = linelist[6:]
            #~ format
            fmt = "GT:DR:DV"
            #~ sample detail
            #` 0/1    4,3    +-
            xlist = detail.strip().split(":")
            dr, dv = list(map(int, xlist[1:3]))
            af = dv / (dr + dv)
            #` genotype
            if af < typing_af_threshold[0]:
                genot = "0/0"
            elif typing_af_threshold[0] <= af < typing_af_threshold[1]:
                genot = "0/1"
            else:
                genot = "1/1"
            #` detail
            sample_detail = "{}:{}:{}".format(genot, dr, dv)
            #~ filter:  lowDR, highDR, lowQUAL, poorGT, q5
            #` lowDR, highDR
            # filt = "PASS"
            if dr + dv < low_high_depth_threshold[0]:
                filt = "lowDR"
            elif dr + dv > low_high_depth_threshold[1]:
                filt = "highDR"
            #` lowQUAL
            if qual < qual_threshold:
                if filt == "PASS":
                    filt = "lowQUAL"
                else:
                    filt = filt + ";lowQUAL"
            #` poorGT
            if genot == "0/0" and af < filter_gt_af_threshold:
                if filt == "PASS":
                    filt = "poorGT"
                else:
                    filt = filt + ";poorGT"
            #~ make up
            #` chrom, pos, id, ref, alt, qual, filter, info, format, sample
            outlist = basic_cols + [filt, info, fmt, sample_detail]
            g.write("\t".join(outlist) + "\n")