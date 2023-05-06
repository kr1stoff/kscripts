#!/usr/bin/env bash

################################# CONFIURATION ################################
export Workpath=/nfs-test/mxf/Project/2.metagenomics/202010_metapipe
export Rawdata=${Workpath}/rawdata
export Outpath=${Workpath}/results
export SampleNamesFile=${Workpath}/config/sample_names.txt
export Threads=16
export NumberOfParallel=4
let Max_Threads=$Threads*$NumberOfParallel
export Max_Threads
export metagenemark=/nfs-test/mxf/Software/MetaGeneMark/MetaGeneMark/mgm
################################### DATABASE ##################################
export kraken_db=/nfs-test/mxf/Database/kraken2/200924
export eggnog_db=/nfs-test/mxf/Database/eggNOG
export local_eggnog_datadir=/nfs-test/mxf/Database/eggNOG
##################################### RUN #####################################
source /nfs-test/mxf/Software/miniconda3/bin/activate base  ## base
conda activate microbe_analysis  ## microbe_analysis
if test -d ${Workpath}/results;then rm -r ${Workpath}/results;fi
### 1. Quality Control and Filter
mkdir -p ${Outpath}/1.quality_control/qc_raw
mkdir -p ${Outpath}/1.quality_control/cleandata
mkdir  -p ${Outpath}/1.quality_control/qc_clean
## fastqc rawdata & multiqc
fastqc -t ${Max_Threads} ${Rawdata}/*.fastq.gz -o ${Outpath}/1.quality_control/qc_raw
multiqc ${Outpath}/1.quality_control/qc_raw -o ${Outpath}/1.quality_control/qc_raw
## trim rawdate & remove host,rRNA
#~ parallel
parallel -j ${NumberOfParallel} --xapply "kneaddata -i {1} -i {2} -o ${Outpath}/1.quality_control/cleandata \
  -t ${Threads} --remove-intermediate-output \
  --trimmomatic /nfs-test/mxf/Software/miniconda3/envs/microbe_analysis/share/trimmomatic \
  --trimmomatic-options 'ILLUMINACLIP:/nfs-test/mxf/Software/miniconda3/envs/microbe_analysis/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:15 MINLEN:36' \
  --bowtie2-options '--very-sensitive --dovetail --reorder' \
  -db /nfs-test/mxf/Database/UCSC/bowtie2_index/hg19 \
  -db /nfs-test/mxf/Database/NCBI/RefSeq/homo_sapiens/hs_rRNA" 
  ::: ${Rawdata}/*1.fastq.gz ::: ${Rawdata}/*2.fastq.gz
kneaddata_read_count_table --input ${Outpath}/1.quality_control/cleandata --output ${Outpath}/1.quality_control/cleandata/kneaddata_read_count.xls \
  && cut -f1-5,16-17 ${Outpath}/1.quality_control/cleandata/kneaddata_read_count.xls | \
  sed 's/ pair//g' > ${Outpath}/1.quality_control/cleandata/kneaddata_read_count.simplified.xls
## fastqc cleandata & multiqc
fastqc -t ${Max_Threads} ${Outpath}/1.quality_control/cleandata/*_paired_[12].fastq -o ${Outpath}/1.quality_control/qc_clean
multiqc ${Outpath}/1.quality_control/qc_clean -o ${Outpath}/1.quality_control/qc_clean

### 2. Assemble
## megahit
mkdir -p ${Outpath}/2.assamble/megahit_out/
mkdir ${Outpath}/2.assamble/quast_out/
#~ echo clean fastq
parallel -j ${NumberOfParallel} --xapply "megahit -t ${Threads} -1 ${Outpath}/1.quality_control/cleandata/{1}_1_kneaddata_paired_1.fastq \
  -2 ${Outpath}/1.quality_control/cleandata/{1}_1_kneaddata_paired_2.fastq -o ${Outpath}/2.assamble/megahit_out/{1}" :::: ${SampleNamesFile}
## quast
parallel -j ${NumberOfParallel} --xapply "quast ${Outpath}/2.assamble/megahit_out/{1}/final.contigs.fa -o ${Outpath}/2.assamble/quast_out/{1}" :::: ${SampleNamesFile}

### 3. Gene Prediction
mkdir -p ${Outpath}/3.gene_prediction/metagenemark_out
mkdir -p ${Outpath}/3.gene_prediction/bowtie2_out
## MetaGeneMark
parallel -j ${NumberOfParallel} --xapply ${Workpath}/scripts/meta_kit/metagenomics.metagenemark.sh {1} :::: ${SampleNamesFile}
#~ ## Merge to create a reference
cat ${Outpath}/3.gene_prediction/metagenemark_out/*/mgm.gene.L100.fa > ${Outpath}/3.gene_prediction/metagenemark_out/mgm.gene_L100_merge.fa
cat ${Outpath}/3.gene_prediction/metagenemark_out/*/mgm.protein.newheader.fa > ${Outpath}/3.gene_prediction/metagenemark_out/mgm.protein.raw_merge.fa
#~ cd-hit-est cluster  50+ min
cd-hit-est -c 0.95 -G 0 -aS 0.9 -g 1 -d 0 -i ${Outpath}/3.gene_prediction/metagenemark_out/mgm.gene_L100_merge.fa -o ${Outpath}/3.gene_prediction/metagenemark_out/gene_catalogue.merge.fa
## GP(Gene Prediction) Bowtie2  1.reads count  2.depth count
bowtie2-build --threads ${Threads} ${Outpath}/3.gene_prediction/metagenemark_out/gene_catalogue.merge.fa ${Outpath}/3.gene_prediction/bowtie2_out/gene_catalogue.merge
parallel -j ${NumberOfParallel} --xapply ${Workpath}/scripts/meta_kit/metagenomics.GP_bowtie2.sh {1} :::: ${SampleNamesFile}
## Gene Catalogue (Unigenes) - filter (1.seqtk -L 100; 2.coverage depth > 2)
#~ 1.gene_catalogue.fa  2.gene_catalogue.protein.fa
#~ MERGED -- 1.gene_catalogue.uniq_gene_merge.fa  2.gene_catalogue.uniq_protein_merge.fa
cat ${Outpath}/3.gene_prediction/bowtie2_out/*/unigenes_id.txt | sort | uniq > ${Outpath}/3.gene_prediction/metagenemark_out/unigenes_id.gene_merge.txt
seqtk subseq ${Outpath}/3.gene_prediction/metagenemark_out/mgm.gene_L100_merge.fa ${Outpath}/3.gene_prediction/metagenemark_out/unigenes_id.gene_merge.txt > \
  ${Outpath}/3.gene_prediction/metagenemark_out/gene_catalogue.uniq_gene_merge.fa
cat ${Outpath}/3.gene_prediction/bowtie2_out/*/unigenes_id.protein.txt | sort | uniq > ${Outpath}/3.gene_prediction/metagenemark_out/unigenes_id.protein_merge.txt
seqtk subseq ${Outpath}/3.gene_prediction/metagenemark_out/mgm.protein.raw_merge.fa ${Outpath}/3.gene_prediction/metagenemark_out/unigenes_id.protein_merge.txt > \
  ${Outpath}/3.gene_prediction/metagenemark_out/gene_catalogue.uniq_protein_merge.fa
## gene prediction stats
#~ Gene catalogue number,  Total length,  Average length
seqkit stats -b -T results/3.gene_prediction/metagenemark_out/gene_catalogue.uniq_gene_merge.fa > results/3.gene_prediction/metagenemark_out/gene_catalogue.uniq_gene_merge.fa.stats
#~ Gene catalogue nucleotide composition
#~ Output format: chr, length, #A, #C, #G, #T, #2, #3, #4, #CpG, #tv, #ts, #CpG-ts
seqtk comp results/3.gene_prediction/metagenemark_out/gene_catalogue.uniq_gene_merge.fa > results/3.gene_prediction/metagenemark_out/gene_catalogue.uniq_gene_merge.fa.comp
${Workpath}/scripts/meta_kit/gene_prediction_stats.py ${Workpath}/config ${Outpath}

### 4. Species Annotation
mkdir -p ${Outpath}/4.species_annotation
#~ MERGED -- kraken2 annotate merge species
kraken2 --db ${kraken_db} --threads ${Threads} \
  --report ${Outpath}/4.species_annotation/merge.kraken2_report \
  --output ${Outpath}/4.species_annotation/merge.kraken2.txt \
  ${Outpath}/3.gene_prediction/metagenemark_out/gene_catalogue.uniq_gene_merge.fa
grep '^C' ${Outpath}/4.species_annotation/merge.kraken2.txt | \
  cut -f 2,3 > ${Outpath}/4.species_annotation/merge.kraken2_classified_clean.txt && \
  sed -i '1i#genemark_id\ttax_id' ${Outpath}/4.species_annotation/merge.kraken2_classified_clean.txt
awk 'BEGIN{FS="\t";OFS="\t"}{gsub(/^[ ]+/,"",$6);print $5,$6,$4}' ${Outpath}/4.species_annotation/merge.kraken2_report > \
  ${Outpath}/4.species_annotation/merge.kraken2_report_clean \
  && sed -i '1i#tax_id\tlatin_name\trank' ${Outpath}/4.species_annotation/merge.kraken2_report_clean
for samp in $(cat ${SampleNamesFile});do mkdir -p ${Outpath}/4.species_annotation/${samp};done
parallel -j ${NumberOfParallel} --xapply "${Workpath}/scripts/meta_kit/annotate_genemark_by_kraken2out.py \
  -n ${Outpath}/3.gene_prediction/bowtie2_out/{1}/out.reads_gt2.coverage.tab \
  -k ${Outpath}/4.species_annotation/merge.kraken2_classified_clean.txt \
  -r ${Outpath}/4.species_annotation/merge.kraken2_report_clean \
  -o ${Outpath}/4.species_annotation/{1}" :::: ${SampleNamesFile}
#~ run 16S Research Pipeline
${Workpath}/scripts/meta_kit/generate_taxa_table.py -s ${SampleNamesFile} -w ${Outpath}/4.species_annotation
/nfs-test/mxf/Software/Kscripts/run_16S_analysis_keyan_meta.sh ${Workpath}/config ${Outpath}/4.species_annotation
conda deactivate  ## microbe_analysis
## species annotation stats
${Workpath}/scripts/meta_kit/species_annotation_stats.py ${Outpath}

### 5. Function Pathway  &&  Antibiotic Resistance Ontology Annotation
${Workpath}/scripts/meta_kit/metagenomics.emapper_rgi.sh
#~ single sample
parallel -j ${NumberOfParallel} --xapply ${Workpath}/scripts/meta_kit/metagenomics.function_annotation.sh {1} :::: ${SampleNamesFile}
#~ multiple sample
${Workpath}/scripts/meta_kit/generate_function_table.py ${SampleNamesFile} ${Outpath}/5.function_annotation
/nfs-test/mxf/Software/miniconda3/envs/R363/bin/Rscript ${Workpath}/scripts/meta_kit/function_analysis.r ${Workpath}/config ${Outpath}/5.function_annotation
/nfs-test/mxf/Software/miniconda3/envs/R363/bin/Rscript ${Workpath}/scripts/meta_kit/function_analysis_card.r ${Workpath}/config ${Outpath}
#~ lefse
${Workpath}/scripts/meta_kit/lefse.meta_pipe.sh ${Outpath}/5.function_annotation ${Workpath}/config/group.txt

conda deactivate  ## base