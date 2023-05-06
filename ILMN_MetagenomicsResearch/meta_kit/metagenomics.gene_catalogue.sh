#!/usr/bin/env bash
#~ ${1}  -->  SRR11575977[sample_name]

if test $# -ne 1;then echo 'Just need 1 argument! [e.g. SRR11575977]';exit -1;fi

seqtk subseq ${Outpath}/3.gene_prediction/metagenemark_out/${1}/gene_catalogue.raw.fa ${Outpath}/3.gene_prediction/bowtie2_out/${1}/unigenes_id.txt > ${Outpath}/3.gene_prediction/metagenemark_out/${1}/gene_catalogue.fa
seqtk seq ${Outpath}/3.gene_prediction/metagenemark_out/${1}/mgm.meta.protein.fa | sed 's/|.*//g' > ${Outpath}/3.gene_prediction/metagenemark_out/${1}/mgm.meta.protein_clean_header.fa
grep '^>' ${Outpath}/3.gene_prediction/metagenemark_out/${1}/gene_catalogue.fa | sed 's/>//g' > ${Outpath}/3.gene_prediction/metagenemark_out/${1}/gene_id.txt
seqtk subseq ${Outpath}/3.gene_prediction/metagenemark_out/${1}/mgm.meta.protein_clean_header.fa ${Outpath}/3.gene_prediction/metagenemark_out/${1}/gene_id.txt > ${Outpath}/3.gene_prediction/metagenemark_out/${1}/gene_catalogue.protein.fa
