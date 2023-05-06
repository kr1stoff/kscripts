#!/usr/bin/env bash
#~ ${1}  -->  SRR11575977[sample_name]

if test $# -ne 1;then echo 'Just need 1 argument! [e.g. SRR11575977]';exit -1;fi

mkdir -p ${Outpath}/4.species_annotation/${1}

#kraken2 --db /nfs-test/mxf/Database/kraken2/200924/ --threads ${Threads} --report ${Outpath}/4.species_annotation/${1}/report --output ${Outpath}/4.species_annotation/${1}/kraken2.txt ${Outpath}/3.gene_prediction/metagenemark_out/${1}/gene_catalogue.fa
#grep '^C' ${Outpath}/4.species_annotation/${1}/kraken2.txt | cut -f 2,3 > ${Outpath}/4.species_annotation/${1}/kraken2_classified_clean.txt && sed -i '1i#genemark_id\ttax_id' ${Outpath}/4.species_annotation/${1}/kraken2_classified_clean.txt
#~ [file]report  -->  $3: number of fragment directly to taxon;    $5: NCBI taxonomic ID number;    $6: Indented scientific name
#awk 'BEGIN{FS="\t";OFS="\t"}{gsub(/^[ ]+/,"",$6);print $5,$6,$4}' ${Outpath}/4.species_annotation/${1}/report > ${Outpath}/4.species_annotation/${1}/report_clean && sed -i '1i#tax_id\tlatin_name\trank' ${Outpath}/4.species_annotation/${1}/report_clean
${Workpath}/scripts/meta_kit/annotate_genemark_by_kraken2out.py -n ${Outpath}/3.gene_prediction/bowtie2_out/${1}/out.depth_gt2.coverage.tab -k ${Outpath}/4.species_annotation/${1}/kraken2_classified_clean.txt -r ${Outpath}/4.species_annotation/${1}/report_clean -o ${Outpath}/4.species_annotation/${1}
