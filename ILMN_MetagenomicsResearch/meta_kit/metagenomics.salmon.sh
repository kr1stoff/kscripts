#!/usr/bin/env bash
#~ ${1}  -->  SRR11575977[sample_name]

if test $# -ne 1;then echo 'Just need 1 argument! [e.g. SRR11575977]';exit -1;fi

mkdir -p ${Outpath}/3.gene_prediction/salmon/${1}
salmon index -t ${Outpath}/3.gene_prediction/metagenemark_out/${1}/gene_catalogue.raw.fa --type quasi -i ${Outpath}/3.gene_prediction/salmon/${1}/salmon_index
salmon quant -i ${Outpath}/3.gene_prediction/salmon/${1}/salmon_index -l A -1 ${Outpath}/1.quality_control/cleandata/${1}_1_kneaddata_paired_1.fastq -2 ${Outpath}/1.quality_control/cleandata/${1}_1_kneaddata_paired_2.fastq --validateMappings -o ${Outpath}/3.gene_prediction/salmon/${1}/salmon_quant
${Workpath}/scripts/meta_kit/gather-counts.py ${Outpath}/3.gene_prediction/salmon/${1}
awk '$2>2' ${Outpath}/3.gene_prediction/salmon/${1}/salmon_quant.counts > ${Outpath}/3.gene_prediction/salmon/${1}/salmon_quant.depth_gt2.counts
