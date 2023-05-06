#!/usr/bin/env bash

######## Options ########
printUsage() {
    echo -e 'PROGRAM <sample_name> <rawdata_path> <output_dir>'
    echo -e '\tsample_name    ::   Sample name. e.g. barcode01'
    echo -e '\trawdata_path   ::   Nanopore FASTQ path. e.g. fastq_runid_0.fastq.gz'
    echo -e '\toutput_dir     ::   Single sample results direcotory. e.g. "."\n'
    echo -e 'Example:  PROGRAM barcode01 fastq_runid_0.fastq.gz temp'
    exit 1
}
#` check arguments number
if [ $# != 3 ];then
    echo -e 'Need 3 Arguments!\n'
    printUsage
fi

# Export Conda Envirionments
export PATH="$MICROBIAL_WGS_BIN:$PATH:$CENTRIFUGE_BIN:$NANOPORE_BIN"
# SCRIPTS_PATH=/sdbb/bioinfor/mengxf/Project/4.nanopore_pipe/metagenomics_sars_scripts
SCRIPTS_PATH=$(dirname $0)

# Sample information
# SAMPLENAME="barcode01"
# INPUT="sars-barcoded-samples/barcode01/fastq_runid_0000000000000000000000000000000000000000_0.fastq.gz"
# RESULTS="results/barcode01"
SAMPLENAME=$1
INPUT=$2
RESULTS=$3

# File tree
mkdir -p ${RESULTS}/1.qc/qc_report
mkdir -p ${RESULTS}/1.qc/clean_data
mkdir -p ${RESULTS}/1.qc/remove_hosts
mkdir -p ${RESULTS}/2.map/target_sars2
mkdir -p ${RESULTS}/2.map/map_qc
mkdir -p ${RESULTS}/3.variants
mkdir -p ${RESULTS}/4.consensus
mkdir -p ${RESULTS}/5.lineages/pangolin
mkdir -p ${RESULTS}/5.lineages/nextclade

# QC & Clean Data
# nanoplot
NanoPlot --plots dot -f png \
    --fastq $INPUT \
    -o ${RESULTS}/1.qc/qc_report
# filtlong
filtlong --min_length 200 --min_mean_q 7 \
    $INPUT \
    > ${RESULTS}/1.qc/clean_data/${SAMPLENAME}.clean.fq
#~ remove hosts
minimap2 -t $THREADS -ax map-ont \
    $HG38HOST_MMI \
    ${RESULTS}/1.qc/clean_data/${SAMPLENAME}.clean.fq \
    > ${RESULTS}/1.qc/remove_hosts/${SAMPLENAME}_hg38.sam 
samtools view -@ 4 -h -f 4 \
    ${RESULTS}/1.qc/remove_hosts/${SAMPLENAME}_hg38.sam \
    > ${RESULTS}/1.qc/remove_hosts/${SAMPLENAME}_hg38_unmapped.sam
picard SamToFastq \
    --INPUT ${RESULTS}/1.qc/remove_hosts/${SAMPLENAME}_hg38_unmapped.sam \
    --FASTQ ${RESULTS}/1.qc/remove_hosts/${SAMPLENAME}_hg38_unmapped.fq

# Mapping
#~ centrifuge find reads
centrifuge -p $THREADS -q \
    --host-taxids 9606 \
    --ignore-quals \
    -S ${RESULTS}/2.map/target_sars2/${SAMPLENAME}_read_classifications.tsv \
    --report-file ${RESULTS}/2.map/target_sars2/${SAMPLENAME}_centrifuge_report.tsv \
    -x $CENTRIFUGE_DB \
    -U ${RESULTS}/1.qc/remove_hosts/${SAMPLENAME}_hg38_unmapped.fq
#~ extract sars reads. out: <sars_descendant_readid.txt>
python ${SCRIPTS_PATH}/utils/extract_readids_from_centrifuge.py \
    ${RESULTS}/2.map/target_sars2/${SAMPLENAME}_read_classifications.tsv
seqkit grep -j 4 \
    -f ${RESULTS}/2.map/target_sars2/sars_descendant_readid.txt \
    ${RESULTS}/1.qc/remove_hosts/${SAMPLENAME}_hg38_unmapped.fq \
    > ${RESULTS}/2.map/target_sars2/${SAMPLENAME}_sars.fq
#~ minimap2 to NC_045512.2
minimap2 -t $THREADS -ax map-ont \
    $SARS2_MMI \
    ${RESULTS}/2.map/target_sars2/${SAMPLENAME}_sars.fq \
    | samtools view -@ 4 -hbS - \
    | samtools sort -@ 4 -o ${RESULTS}/2.map/${SAMPLENAME}_sars2.sorted.bam -
samtools index ${RESULTS}/2.map/${SAMPLENAME}_sars2.sorted.bam
#~ low depth region
samtools depth -a ${RESULTS}/2.map/${SAMPLENAME}_sars2.sorted.bam \
    > ${RESULTS}/2.map/map_qc/${SAMPLENAME}_sars2.sorted.bam.depth
# output. 1 bam_stats.txt, 2 genome_coverage_depth.png, 3 genome_coverage_depth_ylim1000.png
Rscript ${SCRIPTS_PATH}/utils/genome_coverage.R \
    ${RESULTS}/2.map/map_qc/${SAMPLENAME}_sars2.sorted.bam.depth
#~ low coverage regions, region > 20bp
bedtools genomecov -bga \
    -ibam ${RESULTS}/2.map/${SAMPLENAME}_sars2.sorted.bam \
    | awk -v cov="$MIN_COVERAGE" '$4<cov' \
    | bedtools merge -i - \
    | awk '$3-$2>20' \
    1> ${RESULTS}/2.map/${SAMPLENAME}.lowcovmask.bed

# Variants
#~ maskfasta before consensus
bedtools maskfasta \
    -fi $SARS2_FA \
    -bed ${RESULTS}/2.map/${SAMPLENAME}.lowcovmask.bed \
    -fo ${RESULTS}/4.consensus/${SAMPLENAME}.reference_masked.fa
#~ medaka network 
medaka consensus \
    ${RESULTS}/2.map/${SAMPLENAME}_sars2.sorted.bam \
    ${RESULTS}/2.map/target_sars2/consensus_probs.hdf \
    --threads 4
#~ medaka variants
medaka variant \
    ${RESULTS}/4.consensus/${SAMPLENAME}.reference_masked.fa \
    ${RESULTS}/2.map/target_sars2/consensus_probs.hdf \
    ${RESULTS}/3.variants/${SAMPLENAME}.vcf
medaka tools annotate \
    ${RESULTS}/3.variants/${SAMPLENAME}.vcf \
    ${RESULTS}/4.consensus/${SAMPLENAME}.reference_masked.fa \
    ${RESULTS}/2.map/${SAMPLENAME}_sars2.sorted.bam \
    ${RESULTS}/3.variants/${SAMPLENAME}.anno.vcf
#~ medaka consensus
medaka stitch --threads 4 \
    ${RESULTS}/2.map/target_sars2/consensus_probs.hdf \
    ${RESULTS}/4.consensus/${SAMPLENAME}.reference_masked.fa \
    ${RESULTS}/4.consensus/${SAMPLENAME}.consensus.fasta
sed -i 's/'${RLABEL}'/'${SAMPLENAME}'/g' ${RESULTS}/4.consensus/${SAMPLENAME}.consensus.fasta
#~ snpEff annotation, database NC_045512.2
snpEff NC_045512.2 \
    ${RESULTS}/3.variants/${SAMPLENAME}.anno.vcf \
    -htmlStats ${RESULTS}/3.variants/snpEff_summary.html \
    -csvStats ${RESULTS}/3.variants/snpEff_summary.csv \
    > ${RESULTS}/3.variants/${SAMPLENAME}.snpeff.vcf

# Panglin & Nextclade
pangolin ${RESULTS}/4.consensus/${SAMPLENAME}.consensus.fasta \
    -o ${RESULTS}/5.lineages/pangolin
nextclade -j 8 \
    --in-order \
    --input-fasta ${RESULTS}/4.consensus/${SAMPLENAME}.consensus.fasta \
    --input-dataset $NEXTCLADE_SARS2 \
    --output-json ${RESULTS}/5.lineages/nextclade/${SAMPLENAME}.nextclade.json \
    --output-csv ${RESULTS}/5.lineages/nextclade/${SAMPLENAME}.nextclade.csv \
    --output-tsv ${RESULTS}/5.lineages/nextclade/${SAMPLENAME}.nextclade.tsv \
    --output-tree ${RESULTS}/5.lineages/nextclade/${SAMPLENAME}.nextclade.auspice.json \
    --output-dir ${RESULTS}/5.lineages/nextclade/ \
    --output-basename ${SAMPLENAME}
