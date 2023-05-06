#!/usr/bin/env python

import os
import time
import yaml
import argparse
from datetime import datetime
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation


def handle_qualifiers(quali:str, feat):
    """
    处理feature.qualifiers中的信息, 存在就选取, 不存在就用'-'. ','分隔qualifiers
    quali   ::  feature的限定词,也就是feature.qualifiers字典的key
    feat    ::  feature对象
    """
    if quali in feat.qualifiers: 
        qualis = feat.qualifiers[quali]
    else:
        qualis = ["-"]
    return " | ".join(qualis)

def get_taxon(feat):
    """获取record的taxon id"""
    if ("db_xref" in feat.qualifiers) and ("taxon:" in feat.qualifiers["db_xref"][0]):
        return feat.qualifiers["db_xref"][0].replace("taxon:", "")
    else:
        return "-"

def gbk2dict(genbank):
    """genbank格式转成字典"""
    recs = SeqIO.parse(genbank, format="genbank")
    outdict = {}
    dict_strand = {1:"+", -1:"-", 0:"?"}
    for rec in recs:
        ver = rec.id
        outdict[ver] = {}
        outdict[ver]["ACCESSION"] = rec.name
        outdict[ver]["VERSION"] = ver
        outdict[ver]["DATE"] = datetime.strftime(datetime.strptime(rec.annotations['date'], '%d-%b-%Y'), "%Y-%m-%d")
        outdict[ver]["ORGANISM"] = rec.annotations["organism"] 
        outdict[ver]["JOURNAL"] = " | ".join([anno.title for anno in rec.annotations["references"]]) # \s|\s 分隔journal
        outdict[ver]["FEATURE"] = []
        for feat in rec.features:
            ftype = feat.type
            if ftype == "source":
                outdict[ver]["HOST"]  = handle_qualifiers("host", feat)
                outdict[ver]["COUNTRY"] = handle_qualifiers("country", feat)
                outdict[ver]["TAXON"] = get_taxon(feat)
            elif ftype in ["gene", "CDS", "mat_peptide"]:
                outgenes = handle_qualifiers("gene", feat)
                outprods = handle_qualifiers("product", feat)
                locs = []
                for loc in feat.location.parts:
                    start = loc.start
                    end = loc.end
                    strand = dict_strand[loc.strand]
                    loc = f"{start}-{end}:{strand}"
                    locs.append(loc)
                outlocs = ",".join(locs) # 可能没有位置，一个位置，或者多个位置
                outdict[ver]["FEATURE"].append({"TYPE":ftype,"GENE":outgenes, "LOCATION":outlocs, "PRODUCT":outprods})
    return outdict

def dict2output(outdict, output, byaml:bool):
    """gbk字典转表格, yaml为可选输出"""
    # table
    gh = open(output, "wt", encoding="utf-8", newline="")
    headers = ["ACCESSION","VERSION","DATE","ORGANISM","TAXON","HOST","COUNTRY","TYPE","GENE","LOCATION","PRODUCT","JOURNAL"]
    gh.write("\t".join(headers) + "\n")
    for ver in outdict:
        for feat in outdict[ver]["FEATURE"]:
            outlist = [outdict[ver]["ACCESSION"], outdict[ver]["VERSION"], outdict[ver]["DATE"], 
                    outdict[ver]["ORGANISM"], outdict[ver]["TAXON"], outdict[ver]["HOST"], outdict[ver]["COUNTRY"],
                    feat["TYPE"], feat["GENE"], feat["LOCATION"], feat["PRODUCT"], outdict[ver]["JOURNAL"]]
            gh.write("\t".join(outlist) + "\n")
    gh.close()
    # yaml
    if byaml:
        with open(f"{os.path.splitext(output)[0]}.yaml", "wt", encoding="utf-8", newline="") as gh:
            gh.write(yaml.dump(outdict, allow_unicode=True))

def get_args():
    parser = argparse.ArgumentParser()
    parser.description = "读入genbank格式文件, 输出gene,cds,mat_peptide相关的feature信息."
    parser.add_argument("-i", "--ingbk", required=True, help="输入gbk文件.")
    parser.add_argument("-o", "--output", default="./out.tsv", help="[可选]输出表文件. [默认:./out.tsv]")
    parser.add_argument("-y", "--yaml", action="store_true", help="[可选]是否输出yaml文件.默认不输出.")
    return parser.parse_args()
    
    
if __name__ == "__main__":
    args = get_args()
    outdict = gbk2dict(args.ingbk)
    dict2output(outdict, args.output, byaml=args.yaml)
