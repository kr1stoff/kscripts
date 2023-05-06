#!/usr/bin/env python

import sys
import os
import subprocess
import yaml


# Metagenomics Sars-Cov-2 analysis environments
env_yaml = "metagenomics_sars_scripts/libs/meta_sars2_ont.yaml"
with open(env_yaml) as fh:
    env_dict = yaml.safe_load(fh)

# main
cmd = "bash {shell} {sample} {fastq} {results}".format(
    shell = "metagenomics_sars_scripts/sars_cov2_meta_ont.sh",
    sample = "barcode01",
    fastq = "sars-barcoded-samples/barcode01/fastq_runid_0000000000000000000000000000000000000000_0.fastq.gz",
    results = "results/barcode01"
    )
sres = subprocess.run(cmd, shell=True, capture_output=True, check=True, encoding='utf-8', env=env_dict)
ofile = "results/logs/barcode01.sars_cov2_meta_ont.o"
efile = "results/logs/barcode01.sars_cov2_meta_ont.e"
with open(ofile, "w", encoding="utf-8", newline="") as gh:
    gh.write(sres.stdout)
with open(efile, "w", encoding="utf-8", newline="") as gh:
    gh.write(sres.stderr)
