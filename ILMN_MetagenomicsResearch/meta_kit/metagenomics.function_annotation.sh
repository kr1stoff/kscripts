   #!/usr/bin/env bash
#~ ${1}  -->  SRR11575977[sample_name]

if test $# -ne 1;then echo 'Just need 1 argument! [e.g. SRR11575977]';exit -1;fi

mkdir -p ${Outpath}/5.function_annotation/${1}

## split emapper annotation results by each sample
#~ eggNOG, KEGG, CAZy, GO
python3 ${Workpath}/scripts/meta_kit/merge_emapper_depth.py ${Outpath}/5.function_annotation/out.emapper.annotations ${Outpath}/3.gene_prediction/bowtie2_out/${1}/unigenes_id.protein.txt > ${Outpath}/5.function_annotation/${1}/out.emapper_simgle.annotations
#~ CARD
python3 ${Workpath}/scripts/meta_kit/merge_emapper_depth.py ${Outpath}/5.function_annotation/gene_catalogue.uniq_protein_merge.rgi.txt ${Outpath}/3.gene_prediction/bowtie2_out/${1}/unigenes_id.protein.txt > ${Outpath}/5.function_annotation/${1}/rgi_single.txt
## function annotation [GO, KEGG, eggNOG, CAZy]
/nfs-test/mxf/Software/miniconda3/envs/R363/bin/Rscript ${Workpath}/scripts/meta_kit/function_annotation.r ${Outpath} ${1}
#~ function annotation CARD
/nfs-test/mxf/Software/miniconda3/envs/R363/bin/Rscript ${Workpath}/scripts/meta_kit/function_annotation_card.r  ${Outpath} ${1}