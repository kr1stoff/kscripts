#!/usr/bin/env bash
#~ ${1}  -->  SRR11575977[sample_name]

if test $# -ne 1;then echo 'Just need 1 argument! [e.g. SRR11575977]';exit -1;fi

mkdir -p ${Outpath}/3.gene_prediction/metagenemark_out/${1}
${metagenemark}/gmhmmp -m ${metagenemark}/MetaGeneMark_v1.mod -A ${Outpath}/3.gene_prediction/metagenemark_out/${1}/mgm.meta.protein.fa -D ${Outpath}/3.gene_prediction/metagenemark_out/${1}/mgm.meta.gene.fa -o ${Outpath}/3.gene_prediction/metagenemark_out/${1}/mgm.meta.gff -g 11 -f G ${Outpath}/2.assamble/megahit_out/${1}/final.contigs.fa

## [Probably will be use!]
##~ duplicate k141_
#sed 's/|.*//g' ${Outpath}/3.gene_prediction/metagenemark_out/${1}/mgm.meta.gene.fa | seqtk seq -L 100 > ${Outpath}/3.gene_prediction/metagenemark_out/${1}/mgm.gene.L100.fa
##~ CD-HIT
#cd-hit-est -c 0.95 -G 0 -aS 0.9 -g 1 -d 0 -i ${Outpath}/3.gene_prediction/metagenemark_out/${1}/mgm.gene.L100.fa -o ${Outpath}/3.gene_prediction/metagenemark_out/${1}/gene_catalogue.raw.fa

## Merge to create a reference, gene & protein
sed 's/|.*/_'${1}'/g' ${Outpath}/3.gene_prediction/metagenemark_out/${1}/mgm.meta.gene.fa | seqtk seq -L 100 > ${Outpath}/3.gene_prediction/metagenemark_out/${1}/mgm.gene.L100.fa
sed 's/|.*/_'${1}'/g' ${Outpath}/3.gene_prediction/metagenemark_out/${1}/mgm.meta.protein.fa > ${Outpath}/3.gene_prediction/metagenemark_out/${1}/mgm.protein.newheader.fa
