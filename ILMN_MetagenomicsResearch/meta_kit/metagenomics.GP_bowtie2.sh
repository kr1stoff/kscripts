#!/usr/bin/env bash
#~ ${1}  -->  SRR11575977[sample_name]

if test $# -ne 1;then echo 'Just need 1 argument! [e.g. SRR11575977]';exit -1;fi

mkdir -p ${Outpath}/3.gene_prediction/bowtie2_out/${1}
#~ uodate merged all reference sequences
#bowtie2-build --threads ${Threads} ${Outpath}/3.gene_prediction/metagenemark_out/${1}/gene_catalogue.raw.fa ${Outpath}/3.gene_prediction/bowtie2_out/${1}/gene_catalogue
#bowtie2 -p ${Threads} --end-to-end --sensitive -I 200 -x ${Outpath}/3.gene_prediction/bowtie2_out/${1}/gene_catalogue -1 ${Outpath}/1.quality_control/cleandata/${1}_1_kneaddata_paired_1.fastq -2 ${Outpath}/1.quality_control/cleandata/${1}_1_kneaddata_paired_2.fastq | samtools view -@ ${Threads} -Sb - | samtools sort -@ ${Threads} - > ${Outpath}/3.gene_prediction/bowtie2_out/${1}/out.sort.bam
bowtie2 -p ${Threads} --end-to-end --sensitive -I 200 -x ${Outpath}/3.gene_prediction/bowtie2_out/gene_catalogue.merge -1 ${Outpath}/1.quality_control/cleandata/${1}_1_kneaddata_paired_1.fastq -2 ${Outpath}/1.quality_control/cleandata/${1}_1_kneaddata_paired_2.fastq | samtools view -@ 4 -Sb - | samtools sort -@ 4 - > ${Outpath}/3.gene_prediction/bowtie2_out/${1}/out.sort.bam
samtools index ${Outpath}/3.gene_prediction/bowtie2_out/${1}/out.sort.bam

#~ average depth of gene catalogue reference 
bedtools genomecov -ibam ${Outpath}/3.gene_prediction/bowtie2_out/${1}/out.sort.bam > ${Outpath}/3.gene_prediction/bowtie2_out/${1}/out.genomecov.txt
${MetagenomicsResearch}/meta_kit/calculate-contig-coverage.py ${Outpath}/3.gene_prediction/bowtie2_out/${1}/out.genomecov.txt  ## out.genomecov.coverage.tab
awk '$2>2' ${Outpath}/3.gene_prediction/bowtie2_out/${1}/out.genomecov.coverage.tab > ${Outpath}/3.gene_prediction/bowtie2_out/${1}/out.depth_gt2.coverage.tab
cut -f1 ${Outpath}/3.gene_prediction/bowtie2_out/${1}/out.depth_gt2.coverage.tab > ${Outpath}/3.gene_prediction/bowtie2_out/${1}/unigenes_id.txt

#~ reads count of gene catalogue reference 
samtools idxstats  -@ 4 ${Outpath}/3.gene_prediction/bowtie2_out/${1}/out.sort.bam | awk -F '\t' '$3>2' | cut -f1,3 > ${Outpath}/3.gene_prediction/bowtie2_out/${1}/out.reads_gt2.coverage.tab
cut -f1 ${Outpath}/3.gene_prediction/bowtie2_out/${1}/out.reads_gt2.coverage.tab > ${Outpath}/3.gene_prediction/bowtie2_out/${1}/unigenes_id.protein.txt
