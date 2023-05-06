#!/bin/bash

#-----------------------------------
# conda environment: biosoft_py37
# notices:
#   1.zcat multiple fastq files, and minimap using the pip characters
#-----------------------------------

############# database #############
hg38_reference=/nfs-test/mxf/Database/references/homo_sapiens/hg38.fa
hg38_minimap2_index=${hg38_reference/fa/mmi}
ucsc_dir=/nfs-test/mxf/Database/UCSC/hg38
############## sample ##############
threads=32
home_dir=/nfs-test/mxf/Project/3.Structural_variants/202105_ZD_LQQ
sample_work_dir=/nfs-test/mxf/Project/3.Structural_variants/202105_ZD_LQQ/data/20210503XF
fastq_dir=${sample_work_dir}/fastq
####################################

## create minimap2 hg38 index 
# minimap2 -ax map-ont --MD -Y $hg38_reference -d $hg38_minimap2_index
## Alignment 
zcat ${fastq_dir}/*fastq.gz | \
  minimap2 -t $threads -ax map-ont --MD -Y ${hg38_minimap2_index} - | \
  samtools sort -@ $threads -O BAM -o ${sample_work_dir}/minimap2.bam - && \
samtools index -@ $threads ${sample_work_dir}/minimap2.bam
## no bed file, create bed from bam
bedtools bamtobed -i ${sample_work_dir}/minimap2.bam | \
  bedtools merge -i - | \
  awk '$1!~/_random|chrUn_/' > ${sample_work_dir}/minimap2_bam.bed

## Call SVs by 3 tools
#~ sniffles
sniffles -t $threads --min_support 1 --min_length 50 --num_reads_report -1 \
  --min_seq_size 500 --genotype \
  -m ${sample_work_dir}/minimap2.bam \
  -v ${sample_work_dir}/sniffles.vcf
#~ cuteSV
mkdir ${sample_work_dir}/cutesv_calls
cuteSV -t $threads --min_size 50 --max_size 100000 --retain_work_dir \
  --report_readid --min_support 1 --genotype \
  --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 \
  --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
  ${sample_work_dir}/minimap2.bam $hg38_reference \
  ${sample_work_dir}/cutesv.vcf ${sample_work_dir}/cutesv_calls
#~ svim
mkdir ${sample_work_dir}/svim_calls
svim alignment --min_sv_size 50 --minimum_depth 1 \
  --read_names --sequence_alleles \
  ${sample_work_dir}/svim_calls ${sample_work_dir}/minimap2.bam $hg38_reference
## Merge SVs
ls ${sample_work_dir}/sniffles.vcf ${sample_work_dir}/cutesv.vcf \
  ${sample_work_dir}/svim_calls/variants.vcf > ${sample_work_dir}/vcf_names
#~ Minimum number of supporting caller 1
SURVIVOR merge ${sample_work_dir}/vcf_names 1000 1 1 1 0 50 ${sample_work_dir}/merged_1callers.vcf
#~ Minimum number of supporting caller 2
SURVIVOR merge ${sample_work_dir}/vcf_names 1000 2 1 1 0 50 ${sample_work_dir}/merged_2callers.vcf
#~ Minimum number of supporting caller 3
SURVIVOR merge ${sample_work_dir}/vcf_names 1000 3 1 1 0 50 ${sample_work_dir}/merged_3callers.vcf

## Filter
grep -v -E '_random|chrUn_' ${sample_work_dir}/merged.vcf | \
  #~ [filter] target region 
  bedtools intersect -header -a - -b ${sample_work_dir}/minimap2_bam.bed | \
  #~ [filter] UCSC gap region
  bedtools intersect -header -f 0.5 -v -a - -b ${ucsc_dir}/UCSC_gap.hg38.bed | \
  bedtools intersect -header -f 0.5 -v -a - -b ${ucsc_dir}/UCSC_centromeres.hg38.bed | \
  bedtools intersect -header -f 0.5 -v -a - -b ${ucsc_dir}/UCSC_telomere.hg38.bed | \
  bedtools intersect -header -f 0.5 -v -a - -b ${ucsc_dir}/UCSC_segmental_dups.hg38.bed | \
  bedtools intersect -header -f 0.5 -v -a - -b ${ucsc_dir}/UCSC_single_repeats.hg38.bed > \
  ${sample_work_dir}/merged_filt.vcf

## AnnotSV
$ANNOTSV/bin/AnnotSV -genomeBuild GRCh38 -SVinputInfo 0 -annotationMode full \
  -SVinputFile ${sample_work_dir}/merged_filt.vcf \
  -outputFile ${sample_work_dir}/merged.hg38_annotsv.tsv
