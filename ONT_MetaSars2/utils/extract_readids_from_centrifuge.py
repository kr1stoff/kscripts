#!/usr/bin/env python

import os
import sys


def get_sars_descendants():
    # sars 694009 descendant tax_ids list
    sars_descendant_taxidlist = list()
    #~ centrifuge refseq 106 Sar-Cov-2 tax_ids
    sc2_106_taxidlist = list(range(9000001, 9000107))
    sars_descendant_taxidlist += sc2_106_taxidlist
    with open(sars_694009_descendant) as fh:
        for line in fh:
            tax_id = line.strip().split(" ")[0]
            sars_descendant_taxidlist.append(tax_id)
    return sars_descendant_taxidlist

def centrifuge2taxids(sars_descendant_taxidlist):
    # classification output & output sars_descendant_readid.txt for seqtk
    # centrifuge 
    uniq_readid_list = list()
    gh = open(sars_descendant_readid_file, "w", encoding="utf-8", newline="")
    with open(centrifuge_class_out) as fh:
        header      = next(fh)
        header_list = header.strip().split("\t")
        header_enum = list(enumerate(header_list))
        header_dict = {en[1]:en[0] for en in header_enum}
        for line in fh:
            linelist    = line.strip().split("\t")
            readid      = linelist[header_dict["readID"]]
            taxid       = linelist[header_dict["taxID"]]
            # remove readid repetition
            if readid in uniq_readid_list:
                continue
            else:
                uniq_readid_list.append(readid)
            if taxid in sars_descendant_taxidlist:
                gh.write(readid + "\n")
    gh.close()

# position 
script_dir                  = os.path.split(sys.path[0])[0]
sars_694009_descendant      = script_dir + os.sep + "libs/SARS_694009.taxonkit_list.txt"
# IO
assert len(sys.argv) == 2, \
    "ERROR - extract_readids_from_centrifuge.py - Just need 1 Argumeengs, read_classifications.tsv"
# centrifuge_class_out        = "temp/barcode01_read_classifications.tsv"
centrifuge_class_out        = sys.argv[1]
output_dir                  = os.path.dirname(centrifuge_class_out)
sars_descendant_readid_file = output_dir + os.sep + "sars_descendant_readid.txt"
# run
sars_descendant_taxidlist = get_sars_descendants()
centrifuge2taxids(sars_descendant_taxidlist)